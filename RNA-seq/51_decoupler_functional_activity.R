#!/usr/bin/env Rscript
# 51_decoupler_functional_activity.R — Functional Activity Space via decoupleR
#
# Projects dream t-statistics into TF activity and pathway activity space
# using decoupleR (Badia-i-Mompel et al., Bioinformatics Advances 2022).
# This enables principled cross-modality comparison in a shared functional
# space: TF activities inferred from bulk DE can be compared against
# SCENIC+ regulon activities, GWAS-implicated genes, and spatial data.
#
# Produces:
#   1. TF activity scores (ULM — recommended single method) from DoRothEA regulons
#   2. Pathway activity scores (ULM) from PROGENy footprints
#   3. TF-source convergence: TFs whose targets overlap DE + GWAS + regulons
#   4. Functional activity atlas (TF × pathway per gene annotation)
#
# Usage: Rscript 51_decoupler_functional_activity.R
# Compute: login node OK (~5-10 min)
# Requires: decoupleR 2.12+, progeny, OmnipathR, data.table

suppressPackageStartupMessages({
  library(data.table)
  library(decoupleR)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ME   <- file.path(BASE, "RNA-seq/results/multi_evidence")
OUTDIR <- file.path(ME, "functional_activity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== decoupleR Functional Activity Space ===\n")
cat("Output directory:", OUTDIR, "\n\n")

# ── 1. Load dream t-statistics ────────────────────────────────────────────
dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"))
cat("Dream results loaded:", nrow(dream), "genes\n")

# Build named t-statistic vector (gene symbol → t-stat)
# Need gene symbols — map from ensembl via atlas
atlas <- fread(file.path(ME, "multi_evidence_atlas.csv"))
cat("Atlas loaded:", nrow(atlas), "genes\n")

# Map dream genes to symbols
dream[, ensembl_clean := sub("\\..*", "", gene)]
symbol_map <- atlas[!is.na(human_symbol) & human_symbol != "",
                     .(ensembl_clean = sub("\\..*", "", ensembl_id), human_symbol)]
symbol_map <- symbol_map[!duplicated(ensembl_clean)]

dream <- merge(dream, symbol_map, by = "ensembl_clean", all.x = TRUE)
dream_named <- dream[!is.na(human_symbol) & human_symbol != ""]
dream_named <- dream_named[!duplicated(human_symbol)]  # keep first if duplicates

# Create t-statistic matrix: genes (rows) × 1 condition (column)
# decoupleR expects features in rows, samples in columns
tstat_vec <- setNames(dream_named$t, dream_named$human_symbol)
tstat_mat <- matrix(tstat_vec, ncol = 1, dimnames = list(names(tstat_vec), "MASLD_vs_Control"))
cat("T-statistic matrix:", nrow(tstat_mat), "genes ×", ncol(tstat_mat), "condition\n\n")

# ── 2. Load prior knowledge networks ─────────────────────────────────────
cat("Loading DoRothEA TF-target network...\n")
tf_net <- as.data.table(get_dorothea(organism = "human", levels = c("A", "B", "C")))
cat("  DoRothEA:", nrow(tf_net), "interactions,", length(unique(tf_net$source)), "TFs\n")

cat("Loading PROGENy pathway footprints...\n")
pw_net <- as.data.table(get_progeny(organism = "human", top = 500))
cat("  PROGENy:", nrow(pw_net), "gene-pathway links,",
    length(unique(pw_net$source)), "pathways\n\n")

# ── 3. Run decoupleR: TF activity inference ──────────────────────────────
# Use ULM (Univariate Linear Model) — the recommended decoupleR method for
# single-sample activity inference. Avoids invalid averaging across methods
# that produce scores on incompatible scales (ULM=t-stat, WMEAN=normalized
# weighted mean, WSUM=normalized weighted sum).
cat("Running TF activity inference (ULM)...\n")

tf_ulm <- as.data.table(run_ulm(mat = tstat_mat, net = tf_net, .source = "source",
                   .target = "target", .mor = "mor", minsize = 5))

tf_consensus <- tf_ulm[, .(
  tf = source,
  condition = condition,
  score = score,
  p_value = p_value,
  statistic = statistic
)]

# Rank by absolute activity
tf_consensus[, abs_score := abs(score)]
tf_consensus[, rank := rank(-abs_score, ties.method = "first")]
tf_consensus[, direction := fifelse(score > 0, "activated", "repressed")]
tf_consensus[, padj := p.adjust(p_value, method = "BH")]

setorder(tf_consensus, rank)
fwrite(tf_consensus, file.path(OUTDIR, "tf_activity_scores.csv"))

cat("\nTop 20 TFs by activity:\n")
print(tf_consensus[1:20, .(tf, score = round(score, 3),
                             padj = signif(padj, 3), direction)])

n_sig_tf <- tf_consensus[padj < 0.05, .N]
cat("\nSignificant TFs (padj < 0.05):", n_sig_tf, "of", nrow(tf_consensus), "\n")
cat("  Activated:", tf_consensus[padj < 0.05 & direction == "activated", .N], "\n")
cat("  Repressed:", tf_consensus[padj < 0.05 & direction == "repressed", .N], "\n\n")

# ── 4. Run decoupleR: Pathway activity inference ─────────────────────────
cat("Running pathway activity inference (ULM)...\n")

pw_ulm <- as.data.table(run_ulm(mat = tstat_mat, net = pw_net, .source = "source",
                    .target = "target", .mor = "weight", minsize = 5))

pw_dt <- pw_ulm
setnames(pw_dt, "source", "pathway")
pw_dt[, direction := fifelse(score > 0, "activated", "repressed")]
pw_dt[, padj := p.adjust(p_value, method = "BH")]
setorder(pw_dt, p_value)

fwrite(pw_dt, file.path(OUTDIR, "pathway_activity_scores.csv"))

cat("\nPathway activities (all 14 PROGENy pathways):\n")
print(pw_dt[, .(pathway, score = round(score, 3), padj = signif(padj, 3), direction)])

# ── 5. Cross-reference with GWAS-implicated genes ────────────────────────
cat("\n--- TF-GWAS convergence analysis ---\n")

# Get GWAS-causal genes from atlas — COLOC sources only
# broadaway_coloc_pp4, best_liver_enzyme_pp4, ukbb_alt_coloc_pp4, bbj COLOC
# (FinnGen removed 2026-04-08: GWAS archived, duplicate of Whitfield 2023)
# (MR significance removed 2026-04-22: MR ditched from paper; TWAS+COLOC+INTACT)
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
causal_genes <- coloc_genes
cat("Total GWAS causal genes (any COLOC PP4>0.5):", length(causal_genes), "\n")

# For each significant TF: what fraction of its DoRothEA targets are GWAS-causal?
sig_tfs <- tf_consensus[padj < 0.05, tf]
tf_gwas_overlap <- rbindlist(lapply(sig_tfs, function(tf_name) {
  targets <- tf_net[source == tf_name, target]
  targets_in_atlas <- intersect(targets, atlas$human_symbol)
  n_causal <- sum(targets_in_atlas %in% causal_genes)
  n_targets <- length(targets_in_atlas)

  # Hypergeometric test
  N_universe <- nrow(atlas)
  K_causal <- length(intersect(causal_genes, atlas$human_symbol))
  pval <- phyper(n_causal - 1, K_causal, N_universe - K_causal,
                  n_targets, lower.tail = FALSE)

  data.table(
    tf = tf_name,
    n_targets = n_targets,
    n_causal_targets = n_causal,
    frac_causal = n_causal / max(n_targets, 1),
    hypergeom_pval = pval,
    tf_activity_score = tf_consensus[tf == tf_name, score],
    tf_direction = tf_consensus[tf == tf_name, direction]
  )
}))

tf_gwas_overlap[, hypergeom_padj := p.adjust(hypergeom_pval, method = "BH")]
setorder(tf_gwas_overlap, hypergeom_pval)

fwrite(tf_gwas_overlap, file.path(OUTDIR, "tf_gwas_target_overlap.csv"))

n_gwas_enriched <- tf_gwas_overlap[hypergeom_padj < 0.1, .N]
cat("\nTFs with GWAS-enriched target sets (padj < 0.1):", n_gwas_enriched, "\n")
if (n_gwas_enriched > 0) {
  cat("Top TFs with GWAS-enriched targets:\n")
  print(tf_gwas_overlap[hypergeom_padj < 0.1,
    .(tf, n_causal = n_causal_targets, n_targets, padj = signif(hypergeom_padj, 3))])
}

# ── 6. Cross-reference with SCENIC+ regulon activities ───────────────────
cat("\n--- decoupleR vs SCENIC+ comparison ---\n")

# Load SCENIC+ regulon convergence data
conv_file <- file.path(BASE, "RNA-seq/results/convergence/regulon_drug_gwas_convergence.csv")
if (file.exists(conv_file)) {
  scenic <- fread(conv_file)
  # Get unique TF-level SCENIC+ activities
  scenic_tfs <- scenic[!duplicated(tf_name),
                        .(tf_name, scenic_activity_diff = regulon_activity_diff,
                          scenic_padj = activity_padj)]

  # Merge with decoupleR TF activities
  comparison <- merge(tf_consensus[, .(tf, decoupler_score = score, decoupler_padj = padj)],
                       scenic_tfs, by.x = "tf", by.y = "tf_name", all = FALSE)

  if (nrow(comparison) > 2) {
    rho <- cor(comparison$decoupler_score, comparison$scenic_activity_diff,
               method = "spearman", use = "complete.obs")
    r <- cor(comparison$decoupler_score, comparison$scenic_activity_diff,
             method = "pearson", use = "complete.obs")
    cat("TFs in both decoupleR and SCENIC+:", nrow(comparison), "\n")
    cat("Spearman rho:", round(rho, 3), "\n")
    cat("Pearson r:", round(r, 3), "\n")

    # Direction concordance with binomial test
    both_dir <- comparison[!is.na(scenic_activity_diff) & scenic_activity_diff != 0]
    concordant <- sum(sign(both_dir$decoupler_score) == sign(both_dir$scenic_activity_diff))
    binom_p <- binom.test(concordant, nrow(both_dir), p = 0.5)$p.value
    cat("Direction concordance:", concordant, "/", nrow(both_dir),
        "(", round(100 * concordant / nrow(both_dir), 1), "%)\n")
    cat("Binomial test (H0: 50% concordance): p =", signif(binom_p, 3), "\n")
    cat("NOTE: Low rho with significant direction concordance indicates\n")
    cat("  agreement in sign but not magnitude — expected since decoupleR\n")
    cat("  scores derive from bulk t-statistics while SCENIC+ derives from\n")
    cat("  scATAC regulon activity in a different cohort.\n")

    fwrite(comparison, file.path(OUTDIR, "decoupler_vs_scenic_comparison.csv"))
  } else {
    cat("Too few overlapping TFs for comparison\n")
  }
} else {
  cat("SCENIC+ convergence file not found — skipping comparison\n")
}

# ── 7. Annotate atlas genes with TF/pathway functional context ───────────
cat("\n--- Annotating atlas genes with functional activity context ---\n")

# For each gene in atlas: which significant TFs regulate it (via DoRothEA)?
# And which pathways is it a footprint gene for (via PROGENy)?

# TF regulation annotation
sig_tf_net <- as.data.table(tf_net)[source %in% sig_tfs]
gene_tf_annot <- sig_tf_net[, .(
  n_regulating_tfs = .N,
  regulating_tfs = paste(source, collapse = ";"),
  mean_tf_mor = mean(mor)
), by = target]
setnames(gene_tf_annot, "target", "human_symbol")

# Pathway annotation from PROGENy
sig_pathways <- pw_dt[padj < 0.05, pathway]
sig_pw_net <- as.data.table(pw_net)[source %in% sig_pathways]
gene_pw_annot <- sig_pw_net[, .(
  n_pathway_memberships = .N,
  pathways = paste(unique(source), collapse = ";"),
  mean_pathway_weight = mean(abs(weight))
), by = target]
setnames(gene_pw_annot, "target", "human_symbol")

# Merge annotations with atlas
func_atlas <- merge(atlas[, .(human_symbol, ensembl_id, bulk_logFC, bulk_padj,
                               bulk_tstat, layers_active, sources_active,
                               is_conserved, dgidb_druggable, opentargets_drug)],
                     gene_tf_annot, by = "human_symbol", all.x = TRUE)
func_atlas <- merge(func_atlas, gene_pw_annot, by = "human_symbol", all.x = TRUE)

# Fill NAs
func_atlas[is.na(n_regulating_tfs), n_regulating_tfs := 0]
func_atlas[is.na(n_pathway_memberships), n_pathway_memberships := 0]

# Functional activity score: genes regulated by more sig TFs + in more sig pathways
func_atlas[, functional_score := n_regulating_tfs + 2 * n_pathway_memberships]
func_atlas[, functional_rank := rank(-functional_score, ties.method = "first")]

setorder(func_atlas, functional_rank)
fwrite(func_atlas, file.path(OUTDIR, "gene_functional_annotation.csv"))

cat("Genes with TF regulation context:", func_atlas[n_regulating_tfs > 0, .N], "\n")
cat("Genes with pathway membership:", func_atlas[n_pathway_memberships > 0, .N], "\n")

# ── 8. TF-Pathway-Source convergence table ────────────────────────────────
cat("\n--- Building TF-Pathway-Source convergence ---\n")

# For each significant TF: its pathway context + source evidence
tf_pathway_source <- rbindlist(lapply(sig_tfs, function(tf_name) {
  # Get TF's targets
  targets <- tf_net[source == tf_name, target]
  targets_in_atlas <- intersect(targets, atlas$human_symbol)

  # How many targets are DEGs?
  # Exploratory annotation threshold; primary DEGs: padj<0.05 + |logFC|>0.5 (Script 05b)
  n_deg <- atlas[human_symbol %in% targets_in_atlas & !is.na(bulk_padj) & bulk_padj < 0.1, .N]

  # How many targets are GWAS-causal?
  n_gwas <- sum(targets_in_atlas %in% causal_genes)

  # How many targets are Conserved?
  n_cc <- atlas[human_symbol %in% targets_in_atlas & is_conserved == TRUE, .N]

  # How many targets are druggable?
  n_drug <- atlas[human_symbol %in% targets_in_atlas & dgidb_druggable == TRUE, .N]

  # Is the TF itself a DEG?
  # Exploratory annotation threshold; primary DEGs: padj<0.05 + |logFC|>0.5 (Script 05b)
  tf_is_deg <- atlas[human_symbol == tf_name & !is.na(bulk_padj) & bulk_padj < 0.1, .N] > 0

  # Is the TF itself GWAS-causal?
  tf_is_gwas <- tf_name %in% causal_genes

  data.table(
    tf = tf_name,
    activity_score = tf_consensus[tf == tf_name, score],
    activity_padj = tf_consensus[tf == tf_name, padj],
    direction = tf_consensus[tf == tf_name, direction],
    n_targets_in_atlas = length(targets_in_atlas),
    n_deg_targets = n_deg,
    frac_deg = n_deg / max(length(targets_in_atlas), 1),
    n_gwas_targets = n_gwas,
    n_conserved_targets = n_cc,
    n_druggable_targets = n_drug,
    tf_is_deg = tf_is_deg,
    tf_is_gwas = tf_is_gwas,
    convergence_layers = sum(c(n_deg > 0, n_gwas > 0, n_cc > 0, tf_is_deg, tf_is_gwas))
  )
}))

setorder(tf_pathway_source, -convergence_layers, activity_padj)
fwrite(tf_pathway_source, file.path(OUTDIR, "tf_convergence_summary.csv"))

cat("\nTop TFs by convergence:\n")
print(tf_pathway_source[1:15,
  .(tf, direction, n_deg = n_deg_targets, n_gwas = n_gwas_targets,
    tf_deg = tf_is_deg, tf_gwas = tf_is_gwas, layers = convergence_layers)])

# ── 9. Key drug target TF context ────────────────────────────────────────
cat("\n--- Key drug target TF context ---\n")
key_targets <- c("THRB", "NR1H4", "PPARA", "PPARG", "GLP1R", "PNPLA3", "TM6SF2")
for (g in key_targets) {
  regs <- gene_tf_annot[human_symbol == g]
  pws  <- gene_pw_annot[human_symbol == g]
  if (nrow(regs) > 0 || nrow(pws) > 0) {
    cat(sprintf("  %s: %d regulating TFs", g,
                ifelse(nrow(regs) > 0, regs$n_regulating_tfs, 0)))
    if (nrow(regs) > 0) cat(sprintf(" (%s)", regs$regulating_tfs))
    if (nrow(pws) > 0) cat(sprintf(", %d pathway(s) (%s)", pws$n_pathway_memberships, pws$pathways))
    cat("\n")
  } else {
    cat(sprintf("  %s: no TF/pathway annotation\n", g))
  }
}

# ── Summary ───────────────────────────────────────────────────────────────
cat("\n=== RESULTS SUMMARY ===\n")
cat("TF activities:", nrow(tf_consensus), "TFs scored,", n_sig_tf, "significant (padj<0.05)\n")
cat("Pathway activities:", nrow(pw_dt), "pathways,", pw_dt[padj < 0.05, .N], "significant\n")
cat("TF-GWAS enrichment:", n_gwas_enriched, "TFs with enriched causal targets\n")
cat("Genes with functional context:", func_atlas[n_regulating_tfs > 0 | n_pathway_memberships > 0, .N], "\n")
cat("\nOutputs written to:", OUTDIR, "\n")
cat("  tf_activity_scores.csv\n")
cat("  pathway_activity_scores.csv\n")
cat("  tf_gwas_target_overlap.csv\n")
cat("  decoupler_vs_scenic_comparison.csv\n")
cat("  gene_functional_annotation.csv\n")
cat("  tf_convergence_summary.csv\n")
cat("\nDone.\n")
