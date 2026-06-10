#!/usr/bin/env Rscript
# 52_cosmos_mechanistic_paths.R — COSMOS Mechanistic Pathway Analysis
#
# Uses COSMOS (Causal Oriented Search of Multi-Omic Space; Dugourd et al.,
# Genome Medicine 2021) to find signaling paths connecting:
#   - Upstream: GWAS-causal genes (genetically determined perturbations)
#   - Downstream: Transcriptomically dysregulated genes (disease state)
#
# COSMOS operates on gene-level summary statistics and prior knowledge
# networks (PKN), requiring NO matched samples across modalities. It
# formulates an integer linear program (ILP) to find the most parsimonious
# signaling network connecting upstream perturbations to downstream effects
# through known signaling, regulatory, and metabolic interactions.
#
# This provides *mechanistic interpretation* of why convergent genes converge:
# the current atlas identifies WHAT converges, COSMOS explains HOW.
#
# Usage: sbatch run_multiomics_integration.sh  (or Rscript on login node)
# Compute: login node OK (~10-30 min)
# Requires: cosmosR, CARNIVAL, lpSolve, decoupleR, OmnipathR, data.table

suppressPackageStartupMessages({
  library(data.table)
  library(cosmosR)
  library(decoupleR)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ME   <- file.path(BASE, "RNA-seq/results/multi_evidence")
OUTDIR <- file.path(ME, "cosmos_mechanistic")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== COSMOS Mechanistic Pathway Analysis ===\n")
cat("Output directory:", OUTDIR, "\n\n")

# ── 1. Load data ─────────────────────────────────────────────────────────
dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"))

atlas <- fread(file.path(ME, "multi_evidence_atlas.csv"))

# Map ensembl to symbol
dream[, ensembl_clean := sub("\\..*", "", gene)]
symbol_map <- atlas[!is.na(human_symbol) & human_symbol != "",
                     .(ensembl_clean = sub("\\..*", "", ensembl_id), human_symbol)]
symbol_map <- symbol_map[!duplicated(ensembl_clean)]
dream <- merge(dream, symbol_map, by = "ensembl_clean", all.x = TRUE)
dream <- dream[!is.na(human_symbol) & !duplicated(human_symbol)]
cat("Dream results with symbols:", nrow(dream), "genes\n")

# ── 2. Define upstream inputs (GWAS-causal genes) ────────────────────────
# Include ALL COLOC sources for comprehensive causal gene set
# (FinnGen removed 2026-04-08: GWAS archived, duplicate of Whitfield 2023)
coloc_cols <- intersect(c("broadaway_coloc_pp4", "best_liver_enzyme_pp4",
                           "ukbb_alt_coloc_pp4", "bbj_alt_coloc_pp4",
                           "bbj_ast_coloc_pp4", "bbj_ggt_coloc_pp4"),
                         names(atlas))
coloc_genes <- character(0)
for (cc in coloc_cols) {
  hits <- atlas[!is.na(get(cc)) & get(cc) > 0.5, human_symbol]
  coloc_genes <- union(coloc_genes, hits)
  cat(sprintf("  %s > 0.5: %d genes\n", cc, length(hits)))
}
# MR-gene list removed 2026-04-22 — MR ditched from paper.
upstream_genes <- unique(coloc_genes)
cat("Upstream GWAS-causal genes (COLOC PP4>0.5):", length(upstream_genes), "\n")

# ── 3. Define downstream measurements (DE t-statistics) ──────────────────
sig_degs <- dream[padj < 0.1 & abs(t) > 2]
cat("Significant DEGs (padj<0.1, |t|>2):", nrow(sig_degs), "\n")
downstream_tstat <- setNames(sig_degs$t, sig_degs$human_symbol)

# ── 4. Compute TF activities as signaling readout ────────────────────────
cat("\nComputing TF activities via decoupleR...\n")

tstat_all <- setNames(dream$t, dream$human_symbol)
tstat_mat <- matrix(tstat_all, ncol = 1,
                     dimnames = list(names(tstat_all), "MASLD"))

tf_net_decoupler <- as.data.table(get_dorothea(organism = "human", levels = c("A", "B", "C")))
tf_acts <- as.data.table(run_ulm(mat = tstat_mat, net = tf_net_decoupler, .source = "source",
                    .target = "target", .mor = "mor", minsize = 5))

tf_activities <- setNames(tf_acts$score, tf_acts$source)
cat("TF activities computed:", length(tf_activities), "TFs\n")

# ── 5. Load COSMOS prior knowledge network + TF regulon ──────────────────
cat("\nLoading COSMOS meta-network...\n")
data("meta_network", package = "cosmosR")
cat("Meta-network:", nrow(meta_network), "interactions\n")

tf_regulon <- load_tf_regulon_dorothea()
tf_regulon <- as.data.table(tf_regulon)
cat("COSMOS TF regulon:", nrow(tf_regulon), "TF-target pairs,",
    length(unique(tf_regulon$tf)), "TFs\n")

# ── 6. Format COSMOS inputs ──────────────────────────────────────────────
# Upstream: GWAS genes — direction previously sourced from MR beta has been
# dropped 2026-04-22 (MR ditched from paper). All upstream inputs now use the
# conservative neutral +1 default (COLOC PP4>0.5 means the locus is causal
# but direction is unknown from COLOC alone). This assumption is documented
# in the paper's methods. Alternatives (excluding these genes or imputing
# direction from transcriptomic logFC) would either lose all COLOC inputs
# for COSMOS or introduce circularity with the transcriptomic discovery set.
upstream_input <- numeric(0)
n_no_dir <- 0
for (g in upstream_genes) {
  upstream_input[g] <- 1
  n_no_dir <- n_no_dir + 1
}
upstream_input <- upstream_input[upstream_input != 0]
cat("\nUpstream inputs:", length(upstream_input), "total\n")
cat("  COLOC-only (neutral +1):", n_no_dir, "\n")

# Filter upstream to genes present in PKN (COSMOS requirement)
pkn_genes <- unique(c(meta_network$source, meta_network$target))
upstream_in_pkn <- upstream_input[names(upstream_input) %in% pkn_genes]
upstream_not_in_pkn <- setdiff(names(upstream_input), pkn_genes)
cat("  In PKN:", length(upstream_in_pkn), "/ Not in PKN:", length(upstream_not_in_pkn), "\n")

# Limit inputs for tractability with lpSolve
# Full run: set N_UP=135, N_DOWN=300, DEPTH=8 and use SLURM
N_UP   <- as.integer(Sys.getenv("COSMOS_N_UP", "50"))
N_DOWN <- as.integer(Sys.getenv("COSMOS_N_DOWN", "100"))
DEPTH  <- as.integer(Sys.getenv("COSMOS_DEPTH", "5"))
cat("COSMOS limits: N_UP=", N_UP, ", N_DOWN=", N_DOWN, ", DEPTH=", DEPTH, "\n")

# Sort upstream by max PP4 (strongest COLOC evidence first)
up_evidence <- sapply(names(upstream_in_pkn), function(g) {
  pp4 <- atlas[human_symbol == g, best_liver_enzyme_pp4]
  if (length(pp4) == 0 || is.na(pp4)) return(0)
  pp4
})
upstream_input <- upstream_in_pkn[order(-up_evidence)][seq_len(min(N_UP, length(upstream_in_pkn)))]

# Downstream: top N_DOWN DE genes by absolute t-statistic (filtered to PKN)
top_degs <- sort(abs(downstream_tstat), decreasing = TRUE)
keep_genes <- names(top_degs)[names(top_degs) %in% pkn_genes]
if (length(keep_genes) > N_DOWN) keep_genes <- keep_genes[1:N_DOWN]
downstream_input <- downstream_tstat[keep_genes]
cat("Downstream inputs (top DE genes in PKN):", length(downstream_input), "\n")

# Remove genes that are both upstream and downstream
overlap_genes <- intersect(names(upstream_input), names(downstream_input))
if (length(overlap_genes) > 0) {
  cat("Removing", length(overlap_genes), "genes present in both up/downstream\n")
  downstream_input <- downstream_input[!names(downstream_input) %in% overlap_genes]
}

# Full diff expression data for network filtering
diff_expr_data <- tstat_all

# CARNIVAL options
carnival_opts <- default_CARNIVAL_options(solver = "lpSolve")
carnival_opts$outputFolder <- OUTDIR
carnival_opts$workdir <- OUTDIR

# ── 7. Run COSMOS (signaling-to-metabolism direction) ─────────────────────
cat("\n--- Running COSMOS signaling→metabolism ---\n")
cat("This finds signaling paths from GWAS-causal upstream perturbations\n")
cat("to transcriptomically dysregulated downstream effects.\n\n")

cosmos_success <- FALSE

tryCatch({
  cosmos_input <- preprocess_COSMOS_signaling_to_metabolism(
    meta_network = meta_network,
    tf_regulon = tf_regulon,
    signaling_data = upstream_input,
    metabolic_data = downstream_input,
    diff_expression_data = diff_expr_data,
    maximum_network_depth = DEPTH,
    remove_unexpressed_nodes = TRUE,
    filter_tf_gene_interaction_by_optimization = TRUE,
    CARNIVAL_options = carnival_opts
  )

  cat("COSMOS preprocessed successfully\n")
  cat("  Network nodes:", length(unique(c(cosmos_input$meta_network$source,
                                           cosmos_input$meta_network$target))), "\n")
  cat("  Network edges:", nrow(cosmos_input$meta_network), "\n")

  # Run COSMOS
  cat("Running COSMOS optimization (this may take several minutes)...\n")
  cosmos_result <- run_COSMOS_signaling_to_metabolism(
    data = cosmos_input,
    CARNIVAL_options = carnival_opts
  )

  cat("COSMOS optimization complete\n")

  # Format results
  cosmos_formatted <- format_COSMOS_res(cosmos_result)
  solution_edges <- as.data.table(cosmos_formatted[[1]])
  solution_nodes <- as.data.table(cosmos_formatted[[2]])

  cat("\nSolution network:\n")
  cat("  Nodes:", nrow(solution_nodes), "\n")
  cat("  Edges:", nrow(solution_edges), "\n")

  fwrite(solution_edges, file.path(OUTDIR, "cosmos_solution_edges.csv"))
  fwrite(solution_nodes, file.path(OUTDIR, "cosmos_solution_nodes.csv"))
  cosmos_success <- TRUE

  # ── 8. Analyze solution network ──────────────────────────────────────────
  if (nrow(solution_edges) > 0) {
    # Node degree in solution
    out_degree <- solution_edges[, .N, by = source]
    setnames(out_degree, c("node", "out_degree"))
    in_degree <- solution_edges[, .N, by = target]
    setnames(in_degree, c("node", "in_degree"))

    hub_dt <- merge(out_degree, in_degree, by = "node", all = TRUE)
    hub_dt[is.na(out_degree), out_degree := 0]
    hub_dt[is.na(in_degree), in_degree := 0]
    hub_dt[, total_degree := out_degree + in_degree]
    setorder(hub_dt, -total_degree)

    # Annotate with atlas evidence
    hub_dt <- merge(hub_dt,
      atlas[, .(human_symbol, bulk_logFC, bulk_padj, layers_active,
                is_conserved, dgidb_druggable)],
      by.x = "node", by.y = "human_symbol", all.x = TRUE)

    fwrite(hub_dt, file.path(OUTDIR, "cosmos_hub_nodes.csv"))

    cat("\nTop 20 signaling hubs in COSMOS solution:\n")
    print(hub_dt[1:min(20, nrow(hub_dt)),
      .(node, degree = total_degree, layers = layers_active,
        CC = is_conserved, druggable = dgidb_druggable)])

    # Bridging nodes
    bridge_nodes <- setdiff(hub_dt$node,
                             c(names(upstream_input), names(downstream_input)))
    bridge_dt <- hub_dt[node %in% bridge_nodes]
    setorder(bridge_dt, -total_degree)
    fwrite(bridge_dt, file.path(OUTDIR, "cosmos_bridge_nodes.csv"))

    cat("\nTop bridging nodes (mechanistic intermediaries):\n")
    print(bridge_dt[1:min(15, nrow(bridge_dt)),
      .(node, degree = total_degree, layers = layers_active, druggable = dgidb_druggable)])

    # Enrichment analysis
    n_sol <- nrow(hub_dt)
    n_atlas <- nrow(atlas)
    for (feat in c("dgidb_druggable", "is_conserved")) {
      n_sol_feat <- hub_dt[get(feat) == TRUE, .N]
      n_atlas_feat <- atlas[get(feat) == TRUE, .N]
      OR <- (n_sol_feat / n_sol) / (n_atlas_feat / n_atlas)
      p <- phyper(n_sol_feat - 1, n_atlas_feat, n_atlas - n_atlas_feat,
                   n_sol, lower.tail = FALSE)
      cat(sprintf("\n%s enrichment: OR=%.2f, p=%.3g (%d/%d vs %d/%d)\n",
                  feat, OR, p, n_sol_feat, n_sol, n_atlas_feat, n_atlas))
    }
  }

}, error = function(e) {
  cat("COSMOS error:", conditionMessage(e), "\n")
  cat("Proceeding with TF-mediated signaling analysis only.\n\n")
})

# ── 9. TF-mediated GWAS↔DE bridging analysis ─────────────────────────────
# This runs regardless of whether COSMOS succeeded, as a standalone analysis
cat("\n--- TF-mediated GWAS↔DE bridging analysis ---\n")
cat("Identifying TFs that regulate both GWAS-causal and DE genes.\n")
cat("Using hypergeometric test to assess whether overlap exceeds chance.\n")

# Use the same DoRothEA regulon as decoupleR (Script 51) for consistency
# tf_regulon from cosmosR has 271 TFs; tf_net_decoupler has 429 TFs
# Use tf_net_decoupler for consistency
N_universe <- nrow(atlas)
N_gwas <- length(upstream_genes)
N_deg  <- length(names(downstream_tstat))

tf_gwas_reg <- rbindlist(lapply(names(tf_activities[abs(tf_activities) > 1]), function(tf_name) {
  targets <- tf_net_decoupler[source == tf_name, target]
  targets_in_atlas <- intersect(targets, atlas$human_symbol)
  n_targets <- length(targets_in_atlas)
  if (n_targets < 5) return(NULL)  # skip TFs with <5 targets in atlas

  n_gwas <- sum(targets_in_atlas %in% upstream_genes)
  n_deg  <- sum(targets_in_atlas %in% names(downstream_tstat))

  # Hypergeometric test: is GWAS overlap greater than chance?
  # Expected by chance: n_targets * N_gwas / N_universe
  expected_gwas <- n_targets * N_gwas / N_universe
  gwas_pval <- phyper(n_gwas - 1, N_gwas, N_universe - N_gwas,
                       n_targets, lower.tail = FALSE)

  data.table(
    tf = tf_name,
    tf_activity = tf_activities[tf_name],
    n_total_targets = n_targets,
    n_gwas_targets = n_gwas,
    expected_gwas = round(expected_gwas, 1),
    gwas_enrichment = n_gwas / max(expected_gwas, 0.1),
    gwas_hypergeom_pval = gwas_pval,
    n_deg_targets = n_deg,
    gwas_targets = paste(intersect(targets_in_atlas, upstream_genes), collapse = ";"),
    bridges_gwas_to_de = n_gwas > 0 & n_deg > 5
  )
}))

if (nrow(tf_gwas_reg) > 0) {
  # Multiple testing correction on GWAS enrichment
  tf_gwas_reg[, gwas_hypergeom_padj := p.adjust(gwas_hypergeom_pval, method = "BH")]
  setorder(tf_gwas_reg, gwas_hypergeom_pval)
  fwrite(tf_gwas_reg, file.path(OUTDIR, "tf_gwas_de_bridge.csv"))

  bridging_tfs <- tf_gwas_reg[bridges_gwas_to_de == TRUE]
  sig_bridges  <- tf_gwas_reg[bridges_gwas_to_de == TRUE & gwas_hypergeom_padj < 0.1]
  cat("\nBridging TFs (any GWAS + DE overlap):", nrow(bridging_tfs), "\n")
  cat("Statistically significant bridges (GWAS enrichment padj < 0.1):", nrow(sig_bridges), "\n")

  cat("\nTop 15 bridging TFs by GWAS enrichment p-value:\n")
  print(tf_gwas_reg[bridges_gwas_to_de == TRUE][1:min(15, nrow(bridging_tfs)),
    .(tf, activity = round(tf_activity, 2), n_gwas = n_gwas_targets,
      expected = expected_gwas, enrichment = round(gwas_enrichment, 1),
      padj = signif(gwas_hypergeom_padj, 3), n_deg = n_deg_targets)])

  if (nrow(sig_bridges) > 0) {
    cat("\nStatistically significant GWAS-enriched bridging TFs:\n")
    print(sig_bridges[, .(tf, n_gwas = n_gwas_targets, expected = expected_gwas,
      enrichment = round(gwas_enrichment, 1), padj = signif(gwas_hypergeom_padj, 3),
      gwas_genes = gwas_targets)])
  } else {
    cat("\nNo bridging TFs have statistically significant GWAS target enrichment.\n")
    cat("This is expected: with", N_gwas, "GWAS genes across", N_universe,
        "atlas genes (", round(100 * N_gwas / N_universe, 1), "%),\n")
    cat("individual TF regulons are generally too small for significant enrichment.\n")
    cat("The bridging pattern is informative at the aggregate level, not per-TF.\n")
  }
} else {
  cat("No TFs with GWAS-bridging targets found\n")
}

# ── 10. GWAS→TF→DE path enumeration ──────────────────────────────────────
# For each bridging TF: list the mechanistic path GWAS gene → TF → DE genes
cat("\n--- GWAS→TF→DE mechanistic paths ---\n")

if (nrow(tf_gwas_reg) > 0 && any(tf_gwas_reg$bridges_gwas_to_de)) {
  path_dt <- rbindlist(lapply(which(tf_gwas_reg$bridges_gwas_to_de), function(i) {
    tf_name <- tf_gwas_reg$tf[i]
    targets <- tf_regulon[tf == tf_name, target]
    gwas_t  <- intersect(targets, upstream_genes)
    deg_t   <- intersect(targets, names(downstream_tstat[abs(downstream_tstat) > 3]))

    if (length(gwas_t) > 0 && length(deg_t) > 0) {
      # Top DE targets by t-statistic
      top_deg <- head(names(sort(abs(downstream_tstat[deg_t]), decreasing = TRUE)), 10)
      data.table(
        tf = tf_name,
        tf_activity = tf_activities[tf_name],
        gwas_gene = paste(gwas_t, collapse = ";"),
        n_de_targets = length(deg_t),
        top_de_targets = paste(top_deg, collapse = ";"),
        path_description = paste0(
          paste(gwas_t, collapse = "/"), " → ", tf_name, " → ",
          paste(head(top_deg, 5), collapse = "/"), "...")
      )
    }
  }))

  if (nrow(path_dt) > 0) {
    fwrite(path_dt, file.path(OUTDIR, "gwas_tf_de_paths.csv"))
    cat("Mechanistic paths found:", nrow(path_dt), "\n")
    for (i in seq_len(min(10, nrow(path_dt)))) {
      cat(sprintf("  Path %d: %s (TF activity=%.1f, %d DE targets)\n",
                  i, path_dt$path_description[i],
                  path_dt$tf_activity[i], path_dt$n_de_targets[i]))
    }
  }
}

# ── Summary ───────────────────────────────────────────────────────────────
cat("\n=== RESULTS SUMMARY ===\n")
cat("Upstream GWAS-causal genes:", length(upstream_input), "\n")
cat("Downstream DE genes:", length(downstream_input), "\n")
cat("TF activities computed:", length(tf_activities), "\n")
cat("Meta-network interactions:", nrow(meta_network), "\n")
cat("COSMOS ILP solution:", ifelse(cosmos_success, "SUCCESS", "FAILED (TF bridging analysis completed)"), "\n")

sol_files <- list.files(OUTDIR, pattern = "\\.csv$", full.names = FALSE)
if (length(sol_files) > 0) {
  cat("Output files:\n")
  for (f in sol_files) cat("  ", f, "\n")
}

cat("\nAll outputs written to:", OUTDIR, "\n")
cat("Done.\n")
