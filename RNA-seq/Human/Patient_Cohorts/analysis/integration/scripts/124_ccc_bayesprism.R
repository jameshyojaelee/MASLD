#!/usr/bin/env Rscript
# 124_ccc_bayesprism.R
# ---------------------------------------------------------------------------
# Cell-Cell Communication Inference from BayesPrism Expression
#
# Computes ligand-receptor interaction scores per sample using
# BayesPrism-deconvolved cell-type-specific expression. Tests which
# L-R pairs change at each fibrosis transition.
#
# Pipeline:
#   1. Load LIANA consensus L-R database (exported from liana-py)
#   2. Load BayesPrism cell-type expression for sender/receiver types
#   3. For each sample × L-R pair: score = geomean(ligand_sender, receptor_receiver)
#      weighted by sender proportion (accounts for cell abundance)
#   4. Wilcoxon test at each fibrosis transition (F0→F1, F1→F2, F2→F3, F3→F4)
#   5. Annotate L-R pairs with signaling pathway from MSigDB
#
# Input:
#   - results/progression/cibersortx_celltype_expression/bayesprism_*.csv.gz
#   - results/progression/cibersortx_celltype_expression/bayesprism_proportions.csv
#   - results/staging_classifier/modeling_metadata.csv
#   - data/liana_consensus_lr.csv (exported from liana-py; generated if absent)
#
# Output (all to results/progression/):
#   - ccc_lr_scores.csv         (sample × L-R pair × axis interaction scores)
#   - ccc_transition_de.csv     (transition DE for each L-R × axis)
#   - ccc_summary.csv           (overview: axes, pairs, pathway enrichment)
#
# SLURM: cpu, 16 CPUs, 120G RAM, 48h
# Env:   micromamba activate rnaseq
#
# Usage:
#   sbatch --job-name=stg124_ccc \
#          --partition=cpu --cpus-per-task=16 --mem=120G --time=48:00:00 \
#          --output=logs/124_ccc_%j.out \
#          --error=logs/124_ccc_%j.err \
#          --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && \
#                  micromamba activate rnaseq && \
#                  cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \
#                  Rscript 124_ccc_bayesprism.R'"
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
  library(msigdbr)
})

set.seed(42)
cat("=== 124: Cell-Cell Communication from BayesPrism ===\n")
cat(sprintf("Started: %s\n", Sys.time()))

# ── Paths ──────────────────────────────────────────────────────────────────
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INTEG <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
BP_DIR <- file.path(INTEG, "results/progression/cibersortx_celltype_expression")
OUTDIR <- file.path(INTEG, "results/progression")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(INTEG, "scripts/logs"), showWarnings = FALSE, recursive = TRUE)

# Data directory for cached L-R resource
DATA_DIR <- file.path(INTEG, "data")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

n_cores <- min(parallel::detectCores(), 16)
cat(sprintf("  Using %d cores\n", n_cores))

# ── 1. Load or generate L-R database ─────────────────────────────────────
cat("\n=== 1. Loading Ligand-Receptor Database ===\n")

lr_file <- file.path(DATA_DIR, "liana_consensus_lr.csv")

if (file.exists(lr_file)) {
  lr_db <- fread(lr_file)
  cat(sprintf("  Loaded cached L-R database: %d pairs\n", nrow(lr_db)))
} else {
  # Export from LIANA raw interactions (which used LIANA's consensus resource)
  # Build from the existing scRNA LIANA results
  liana_raw <- file.path(BASE,
    "Analysis/SingleCell/results_gpu_v2/fig2_data/liana_raw_interactions.csv")

  if (file.exists(liana_raw)) {
    cat("  Building L-R database from existing LIANA results...\n")
    liana <- fread(liana_raw)
    lr_db <- unique(liana[, .(ligand = ligand_complex, receptor = receptor_complex)])
    cat(sprintf("  Extracted %d unique L-R pairs from LIANA results\n", nrow(lr_db)))
  } else {
    stop("No L-R database found. Run: python3 -c \"from liana.resource import ",
         "select_resource; select_resource('consensus').to_csv('",
         lr_file, "', index=False)\"")
  }

  fwrite(lr_db, lr_file)
  cat(sprintf("  Saved L-R database to %s\n", lr_file))
}

# ── 2. Load metadata ─────────────────────────────────────────────────────
cat("\n=== 2. Loading Metadata ===\n")
meta <- fread(file.path(INTEG, "results/staging_classifier/modeling_metadata.csv"))
cat(sprintf("  Metadata: %d samples\n", nrow(meta)))

# Fibrosis distribution
fib_tab <- table(meta$fibrosis_stage, useNA = "ifany")
cat("  Fibrosis stage distribution:\n")
for (nm in names(fib_tab)) {
  cat(sprintf("    F%s: %d\n", nm, fib_tab[nm]))
}

# ── 3. Load BayesPrism proportions ───────────────────────────────────────
cat("\n=== 3. Loading Cell Proportions ===\n")
proportions <- fread(file.path(BP_DIR, "bayesprism_proportions.csv"))
cat(sprintf("  Proportions: %d samples x %d cell types\n",
            nrow(proportions), ncol(proportions) - 1))

# ── 4. Define communication axes ─────────────────────────────────────────
cat("\n=== 4. Defining Communication Axes ===\n")

# Map between BayesPrism cell type names and biologically relevant axes
# BayesPrism names: Hepatocyte, Stellate, Macrophage, Endothelial, Cholangiocyte, etc.
communication_axes <- list(
  Hepatocyte_to_Stellate    = list(sender = "Hepatocyte",  receiver = "Stellate"),
  Hepatocyte_to_Macrophage  = list(sender = "Hepatocyte",  receiver = "Macrophage"),
  Stellate_to_Hepatocyte    = list(sender = "Stellate",    receiver = "Hepatocyte"),
  Macrophage_to_Hepatocyte  = list(sender = "Macrophage",  receiver = "Hepatocyte"),
  Endothelial_to_Stellate   = list(sender = "Endothelial", receiver = "Stellate"),
  Endothelial_to_Hepatocyte = list(sender = "Endothelial", receiver = "Hepatocyte"),
  Macrophage_to_Stellate    = list(sender = "Macrophage",  receiver = "Stellate"),
  Stellate_to_Macrophage    = list(sender = "Stellate",    receiver = "Macrophage")
)

cat(sprintf("  %d communication axes:\n", length(communication_axes)))
for (ax in names(communication_axes)) {
  cat(sprintf("    %s: %s → %s\n", ax,
              communication_axes[[ax]]$sender,
              communication_axes[[ax]]$receiver))
}

# ── 5. Load cell-type expression matrices ─────────────────────────────────
cat("\n=== 5. Loading BayesPrism Cell-Type Expression ===\n")

# Identify which cell types we need
needed_celltypes <- unique(c(
  sapply(communication_axes, `[[`, "sender"),
  sapply(communication_axes, `[[`, "receiver")
))
cat(sprintf("  Cell types needed: %s\n", paste(needed_celltypes, collapse = ", ")))

# Load each cell-type expression matrix
expr_matrices <- list()
gene_universe <- NULL

for (ct in needed_celltypes) {
  ct_file <- file.path(BP_DIR, paste0("bayesprism_", ct, ".csv.gz"))
  if (!file.exists(ct_file)) {
    cat(sprintf("  WARNING: %s not found, skipping\n", ct_file))
    next
  }

  ct_dt <- fread(ct_file)
  sids <- ct_dt$sample_id
  ct_mat <- as.matrix(ct_dt[, -1, with = FALSE])
  rownames(ct_mat) <- sids
  expr_matrices[[ct]] <- ct_mat

  if (is.null(gene_universe)) {
    gene_universe <- colnames(ct_mat)
  } else {
    gene_universe <- intersect(gene_universe, colnames(ct_mat))
  }

  cat(sprintf("  %s: %d samples x %d genes\n", ct, nrow(ct_mat), ncol(ct_mat)))
}

cat(sprintf("  Gene universe (intersection): %d genes\n", length(gene_universe)))

# Align all matrices to common gene set
for (ct in names(expr_matrices)) {
  expr_matrices[[ct]] <- expr_matrices[[ct]][, gene_universe, drop = FALSE]
}

# ── 6. Filter L-R database to feasible pairs ─────────────────────────────
cat("\n=== 6. Filtering L-R Pairs to BayesPrism Gene Universe ===\n")

# Helper: split complex names and check gene availability
split_complex <- function(x) {
  strsplit(as.character(x), "_")[[1]]
}

# Check which pairs have ALL subunit genes in the universe
lr_db[, ligand_genes := lapply(ligand, split_complex)]
lr_db[, receptor_genes := lapply(receptor, split_complex)]

lr_db[, lig_feasible := sapply(ligand_genes, function(gs) all(gs %in% gene_universe))]
lr_db[, rec_feasible := sapply(receptor_genes, function(gs) all(gs %in% gene_universe))]
lr_db[, feasible := lig_feasible & rec_feasible]

n_feasible <- sum(lr_db$feasible)
cat(sprintf("  Total L-R pairs: %d\n", nrow(lr_db)))
cat(sprintf("  Feasible (all genes in BayesPrism): %d (%.1f%%)\n",
            n_feasible, n_feasible / nrow(lr_db) * 100))

lr_feasible <- lr_db[feasible == TRUE, .(ligand, receptor)]
cat(sprintf("  Using %d L-R pairs for analysis\n", nrow(lr_feasible)))

if (nrow(lr_feasible) < 50) {
  warning("Fewer than 50 feasible L-R pairs. Results may lack power.")
}

# ── 7. Compute L-R interaction scores ────────────────────────────────────
cat("\n=== 7. Computing Interaction Scores ===\n")
cat("  Score = sqrt(mean_ligand_in_sender * mean_receptor_in_receiver) * sender_proportion\n")

# Align sample IDs across all matrices and metadata
common_samples <- Reduce(intersect, lapply(expr_matrices, rownames))
common_samples <- intersect(common_samples, meta$sample_id)
common_samples <- intersect(common_samples, proportions$sample_id)
cat(sprintf("  Common samples across all sources: %d\n", length(common_samples)))

meta_aligned <- meta[match(common_samples, meta$sample_id), ]

# Prepare proportions (aligned to common_samples)
prop_aligned <- as.data.frame(proportions)
rownames(prop_aligned) <- prop_aligned$sample_id
prop_aligned <- prop_aligned[common_samples, -1, drop = FALSE]

# For each axis, compute all L-R scores across all samples
compute_axis_scores <- function(axis_name, axis_spec, lr_pairs, expr_mats,
                                props, sample_ids) {
  sender <- axis_spec$sender
  receiver <- axis_spec$receiver

  sender_expr <- expr_mats[[sender]][sample_ids, , drop = FALSE]
  receiver_expr <- expr_mats[[receiver]][sample_ids, , drop = FALSE]
  sender_prop <- props[sample_ids, sender]

  n_samples <- length(sample_ids)
  n_pairs <- nrow(lr_pairs)

  # Pre-allocate score matrix: samples x LR pairs
  scores <- matrix(NA_real_, nrow = n_samples, ncol = n_pairs)

  for (p in seq_len(n_pairs)) {
    lig_genes <- split_complex(lr_pairs$ligand[p])
    rec_genes <- split_complex(lr_pairs$receptor[p])

    # Ligand expression: geometric mean across complex subunits in sender
    if (length(lig_genes) == 1) {
      lig_expr <- sender_expr[, lig_genes[1]]
    } else {
      lig_expr <- exp(rowMeans(log(pmax(sender_expr[, lig_genes, drop = FALSE], 1e-3))))
    }

    # Receptor expression: geometric mean across complex subunits in receiver
    if (length(rec_genes) == 1) {
      rec_expr <- receiver_expr[, rec_genes[1]]
    } else {
      rec_expr <- exp(rowMeans(log(pmax(receiver_expr[, rec_genes, drop = FALSE], 1e-3))))
    }

    # Interaction score = geometric mean of ligand and receptor, weighted by sender proportion
    raw_score <- sqrt(pmax(lig_expr, 0) * pmax(rec_expr, 0))
    scores[, p] <- raw_score * sender_prop
  }

  colnames(scores) <- paste0(lr_pairs$ligand, "::", lr_pairs$receptor)
  rownames(scores) <- sample_ids

  # Return as data.table in long format (minimal columns to control file size)
  dt_list <- list()
  for (p in seq_len(n_pairs)) {
    dt_list[[p]] <- data.table(
      sample_id = sample_ids,
      axis = axis_name,
      lr_pair = paste0(lr_pairs$ligand[p], "::", lr_pairs$receptor[p]),
      score = scores[, p]
    )
  }
  rbindlist(dt_list)
}

# Run across all axes
cat("  Computing scores per axis...\n")
all_scores <- list()

for (ax_name in names(communication_axes)) {
  ax <- communication_axes[[ax_name]]

  # Check that both cell types are available
  if (!(ax$sender %in% names(expr_matrices)) ||
      !(ax$receiver %in% names(expr_matrices))) {
    cat(sprintf("    %s: SKIPPED (missing expression matrix)\n", ax_name))
    next
  }

  t0 <- Sys.time()
  all_scores[[ax_name]] <- compute_axis_scores(
    ax_name, ax, lr_feasible, expr_matrices, prop_aligned, common_samples
  )
  dt <- difftime(Sys.time(), t0, units = "secs")
  n_nonzero <- sum(all_scores[[ax_name]]$score > 0, na.rm = TRUE)
  cat(sprintf("    %s: %d scores computed (%.0f%% non-zero) [%.1fs]\n",
              ax_name, nrow(all_scores[[ax_name]]),
              n_nonzero / nrow(all_scores[[ax_name]]) * 100, dt))
}

scores_dt <- rbindlist(all_scores)
cat(sprintf("  Total score records: %s\n", format(nrow(scores_dt), big.mark = ",")))

# ── 8. Save L-R scores (wide format, one row per sample) ─────────────────
cat("\n=== 8. Saving L-R Scores ===\n")

# Save in long format (more flexible for downstream analysis)
fwrite(scores_dt, file.path(OUTDIR, "ccc_lr_scores.csv"))
cat(sprintf("  Saved ccc_lr_scores.csv (%s rows)\n",
            format(nrow(scores_dt), big.mark = ",")))

# ── 9. Transition DE: Wilcoxon test at each fibrosis step ────────────────
cat("\n=== 9. Transition DE for L-R Pairs ===\n")

fib_transitions <- list(
  F0_to_F1 = c(0, 1),
  F1_to_F2 = c(1, 2),
  F2_to_F3 = c(2, 3),
  F3_to_F4 = c(3, 4)
)

# Minimum score threshold: skip LR pairs where >90% of samples are zero
MIN_NONZERO_FRAC <- 0.10

# Get unique axis-LR combinations
axes_lr <- unique(scores_dt[, .(axis, lr_pair)])
cat(sprintf("  Testing %d axis-LR combinations across %d transitions\n",
            nrow(axes_lr), length(fib_transitions)))

# Build a wide score matrix per axis for faster access
scores_wide <- list()
for (ax_name in unique(scores_dt$axis)) {
  sub <- scores_dt[axis == ax_name]
  wide <- dcast(sub, sample_id ~ lr_pair, value.var = "score")
  sids <- wide$sample_id
  mat <- as.matrix(wide[, -1, with = FALSE])
  rownames(mat) <- sids
  scores_wide[[ax_name]] <- mat
}

# Run transition DE
run_transition_de <- function(trans_name, stages) {
  s_from <- stages[1]
  s_to <- stages[2]

  idx_from <- which(meta_aligned$fibrosis_stage == s_from)
  idx_to <- which(meta_aligned$fibrosis_stage == s_to)

  if (length(idx_from) < 10 || length(idx_to) < 10) {
    cat(sprintf("    %s: SKIPPED (n=%d, %d)\n", trans_name,
                length(idx_from), length(idx_to)))
    return(NULL)
  }

  sids_from <- common_samples[idx_from]
  sids_to <- common_samples[idx_to]

  results <- list()

  for (ax_name in names(scores_wide)) {
    mat <- scores_wide[[ax_name]]
    lr_names <- colnames(mat)

    # Parallelize over L-R pairs
    de_results <- mclapply(seq_along(lr_names), function(p) {
      lr_name <- lr_names[p]
      x_from <- mat[sids_from, p]
      x_to <- mat[sids_to, p]

      # Skip if too few non-zero values
      nonzero_frac <- mean(c(x_from, x_to) > 0)
      if (nonzero_frac < MIN_NONZERO_FRAC) {
        return(NULL)
      }

      # Skip if no variance
      if (sd(c(x_from, x_to)) < 1e-10) {
        return(NULL)
      }

      wt <- tryCatch(
        wilcox.test(x_to, x_from, alternative = "two.sided"),
        error = function(e) list(p.value = 1)
      )

      mean_from <- mean(x_from, na.rm = TRUE)
      mean_to <- mean(x_to, na.rm = TRUE)
      if (mean_from > 0 && mean_to > 0) {
        lfc <- log2(mean_to / mean_from)
      } else {
        lfc <- mean_to - mean_from
      }

      data.table(
        axis = ax_name,
        lr_pair = lr_name,
        transition = trans_name,
        mean_from = mean_from,
        mean_to = mean_to,
        lfc = lfc,
        pval = wt$p.value,
        n_from = length(idx_from),
        n_to = length(idx_to),
        nonzero_frac = nonzero_frac
      )
    }, mc.cores = n_cores)

    results <- c(results, de_results[!sapply(de_results, is.null)])
  }

  if (length(results) > 0) rbindlist(results) else NULL
}

# Run all transitions
de_results_list <- list()
for (trans_name in names(fib_transitions)) {
  cat(sprintf("  Testing %s (F%d vs F%d)...\n",
              trans_name,
              fib_transitions[[trans_name]][1],
              fib_transitions[[trans_name]][2]))
  t0 <- Sys.time()
  de_results_list[[trans_name]] <- run_transition_de(
    trans_name, fib_transitions[[trans_name]]
  )
  dt <- difftime(Sys.time(), t0, units = "secs")
  if (!is.null(de_results_list[[trans_name]])) {
    cat(sprintf("    %d tests completed [%.1fs]\n",
                nrow(de_results_list[[trans_name]]), dt))
  }
}

transition_de <- rbindlist(de_results_list[!sapply(de_results_list, is.null)])

# Multiple testing correction: BH per transition × axis
transition_de[, padj := p.adjust(pval, method = "BH"),
              by = .(transition, axis)]

# Split lr_pair back into ligand/receptor columns
transition_de[, c("ligand", "receptor") := tstrsplit(lr_pair, "::", fixed = TRUE)]

# Count significant
n_sig <- transition_de[padj < 0.1 & abs(lfc) > 0.25, .N]
cat(sprintf("\n  Total transition DE results: %d\n", nrow(transition_de)))
cat(sprintf("  Significant (padj<0.1 & |LFC|>0.25): %d\n", n_sig))

# Print per-transition summary
cat("\n  Per-transition summary:\n")
trans_summary <- transition_de[, .(
  n_tested = .N,
  n_sig = sum(padj < 0.1 & abs(lfc) > 0.25, na.rm = TRUE),
  n_up = sum(padj < 0.1 & lfc > 0.25, na.rm = TRUE),
  n_down = sum(padj < 0.1 & lfc < -0.25, na.rm = TRUE),
  median_lfc = median(lfc, na.rm = TRUE)
), by = .(transition)]

for (i in seq_len(nrow(trans_summary))) {
  cat(sprintf("    %s: %d sig / %d tested (up=%d, down=%d)\n",
              trans_summary$transition[i],
              trans_summary$n_sig[i],
              trans_summary$n_tested[i],
              trans_summary$n_up[i],
              trans_summary$n_down[i]))
}

# Per-axis breakdown
cat("\n  Per-axis summary:\n")
axis_summary <- transition_de[, .(
  n_sig = sum(padj < 0.1 & abs(lfc) > 0.25, na.rm = TRUE),
  n_tested = .N
), by = .(axis)]

for (i in seq_len(nrow(axis_summary))) {
  cat(sprintf("    %-30s: %4d sig / %d tested\n",
              axis_summary$axis[i],
              axis_summary$n_sig[i],
              axis_summary$n_tested[i]))
}

# ── 10. Pathway annotation for L-R pairs ─────────────────────────────────
cat("\n=== 10. Pathway Annotation ===\n")

# Build pathway map: gene → pathway (Hallmark + KEGG + Reactome signaling)
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
kegg <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS")
reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")

# Focus on signaling-relevant pathways
signaling_keywords <- c("TGF", "WNT", "NOTCH", "TNF", "HEDGEHOG", "HIPPO",
                         "NFkB", "NFKB", "JAK_STAT", "JAK.STAT",
                         "INTERLEUKIN", "INTERFERON", "COMPLEMENT",
                         "APOPTOSIS", "INFLAMMATORY", "CHEMOKINE",
                         "PDGF", "EGF", "FGF", "VEGF", "BMP",
                         "COLLAGEN", "ECM", "INTEGRIN", "ADHESION",
                         "TOLL_LIKE", "TOLL.LIKE", "INNATE_IMMUNE",
                         "FIBROSIS", "STELLATE", "HEPATIC",
                         "FATTY_ACID", "LIPID", "CHOLESTEROL",
                         "OXIDATIVE", "FERROPTOSIS")

sig_pattern <- paste(signaling_keywords, collapse = "|")

all_pathways <- rbind(
  as.data.table(hallmark)[, .(gs_name, gene_symbol)],
  as.data.table(kegg)[grepl(sig_pattern, gs_name, ignore.case = TRUE),
       .(gs_name, gene_symbol)],
  as.data.table(reactome)[grepl(sig_pattern, gs_name, ignore.case = TRUE),
           .(gs_name, gene_symbol)]
)
all_pathways <- unique(all_pathways)
cat(sprintf("  Signaling pathways loaded: %d genes across %d gene sets\n",
            length(unique(all_pathways$gene_symbol)),
            length(unique(all_pathways$gs_name))))

# Build gene -> pathway map (take top 3 pathways per gene, prioritize Hallmark)
gene_pathway_map <- all_pathways[, .(
  pathways = paste(head(unique(gs_name), 3), collapse = ";")
), by = gene_symbol]

# Annotate each LR pair with pathway
annotate_lr_pathway <- function(ligand, receptor) {
  lig_genes <- split_complex(ligand)
  rec_genes <- split_complex(receptor)
  all_genes <- unique(c(lig_genes, rec_genes))

  pw <- gene_pathway_map[gene_symbol %in% all_genes]
  if (nrow(pw) > 0) {
    all_pw <- unique(unlist(strsplit(pw$pathways, ";")))
    # Simplify: return top pathway
    return(paste(head(all_pw, 3), collapse = ";"))
  }
  return(NA_character_)
}

# Annotate transition DE results
cat("  Annotating L-R pairs with signaling pathways...\n")
unique_pairs <- unique(transition_de[, .(ligand, receptor)])
unique_pairs[, pathway := mapply(annotate_lr_pathway, ligand, receptor)]

n_annotated <- sum(!is.na(unique_pairs$pathway))
cat(sprintf("  Annotated: %d / %d pairs (%.1f%%)\n",
            n_annotated, nrow(unique_pairs),
            n_annotated / nrow(unique_pairs) * 100))

# Merge pathway annotation back
transition_de <- merge(transition_de, unique_pairs[, .(ligand, receptor, pathway)],
                       by = c("ligand", "receptor"), all.x = TRUE)

# ── 11. Save transition DE results ───────────────────────────────────────
cat("\n=== 11. Saving Results ===\n")

# Order by significance
transition_de[, abs_lfc := abs(lfc)]
setorder(transition_de, padj, -abs_lfc)
transition_de[, abs_lfc := NULL]
fwrite(transition_de, file.path(OUTDIR, "ccc_transition_de.csv"))
cat(sprintf("  Saved ccc_transition_de.csv (%d rows)\n", nrow(transition_de)))

# ── 12. Summary table ────────────────────────────────────────────────────
cat("\n=== 12. Building Summary ===\n")

# Overall summary
summary_overall <- data.table(
  metric = c("total_lr_pairs_in_db",
             "feasible_lr_pairs",
             "communication_axes",
             "total_samples",
             "total_score_records",
             "total_transition_tests",
             "sig_transitions_padj01_lfc025",
             "pct_lr_annotated_pathway"),
  value = c(nrow(lr_db),
            nrow(lr_feasible),
            length(communication_axes),
            length(common_samples),
            nrow(scores_dt),
            nrow(transition_de),
            n_sig,
            round(n_annotated / nrow(unique_pairs) * 100, 1))
)

# Top rewired L-R pairs per transition
top_rewired <- transition_de[padj < 0.1 & abs(lfc) > 0.25][
  order(padj)][, head(.SD, 10), by = transition]

# Per-axis × per-transition DEG counts
axis_trans_summary <- transition_de[, .(
  n_sig = sum(padj < 0.1 & abs(lfc) > 0.25, na.rm = TRUE),
  n_up = sum(padj < 0.1 & lfc > 0.25, na.rm = TRUE),
  n_down = sum(padj < 0.1 & lfc < -0.25, na.rm = TRUE),
  top_pair = {
    sig <- .SD[padj < 0.1 & abs(lfc) > 0.25]
    if (nrow(sig) > 0) sig$lr_pair[which.min(sig$padj)] else NA_character_
  }
), by = .(axis, transition)]

# Pathway enrichment of rewired pairs
if (nrow(transition_de[padj < 0.1 & abs(lfc) > 0.25]) > 0) {
  sig_pairs <- transition_de[padj < 0.1 & abs(lfc) > 0.25]
  sig_pathways <- sig_pairs[!is.na(pathway)]
  if (nrow(sig_pathways) > 0) {
    pw_split <- unlist(strsplit(sig_pathways$pathway, ";"))
    pw_counts <- sort(table(pw_split), decreasing = TRUE)
    cat("\n  Top pathways among rewired L-R pairs:\n")
    for (i in seq_len(min(10, length(pw_counts)))) {
      cat(sprintf("    %s: %d\n", names(pw_counts)[i], pw_counts[i]))
    }
  }
}

# Save summary
fwrite(summary_overall, file.path(OUTDIR, "ccc_summary.csv"))
cat("  Saved ccc_summary.csv\n")

fwrite(axis_trans_summary, file.path(OUTDIR, "ccc_axis_transition_summary.csv"))
cat("  Saved ccc_axis_transition_summary.csv\n")

if (nrow(top_rewired) > 0) {
  fwrite(top_rewired, file.path(OUTDIR, "ccc_top_rewired.csv"))
  cat(sprintf("  Saved ccc_top_rewired.csv (%d rows)\n", nrow(top_rewired)))
}

# ── 13. Print biological highlights ──────────────────────────────────────
cat("\n=== Biological Highlights ===\n")

for (trans in names(fib_transitions)) {
  sig <- transition_de[transition == trans & padj < 0.1 & abs(lfc) > 0.25]
  if (nrow(sig) == 0) {
    cat(sprintf("  %s: No significant rewiring detected\n", trans))
    next
  }

  cat(sprintf("\n  %s (%d rewired L-R pairs):\n", trans, nrow(sig)))

  # Top gained (upregulated) interactions
  gained <- sig[lfc > 0][order(padj)][1:min(3, sum(lfc > 0))]
  if (nrow(gained) > 0) {
    cat("    Gained:\n")
    for (j in seq_len(nrow(gained))) {
      cat(sprintf("      %s [%s] LFC=%.2f, padj=%.2e\n",
                  gained$lr_pair[j], gained$axis[j],
                  gained$lfc[j], gained$padj[j]))
    }
  }

  # Top lost (downregulated) interactions
  lost <- sig[lfc < 0][order(padj)][1:min(3, sum(lfc < 0))]
  if (nrow(lost) > 0) {
    cat("    Lost:\n")
    for (j in seq_len(nrow(lost))) {
      cat(sprintf("      %s [%s] LFC=%.2f, padj=%.2e\n",
                  lost$lr_pair[j], lost$axis[j],
                  lost$lfc[j], lost$padj[j]))
    }
  }
}

cat(sprintf("\n=== 124: COMPLETE (%s) ===\n", Sys.time()))
