#!/usr/bin/env Rscript
# 115_transition_programs.R
# ---------------------------------------------------------------------------
# Stage-Transition Gene Programs: Identify genes that change specifically at
# EACH disease transition, not just "disease vs control".
#
# Pipeline:
#   1. Load merged DGE (1,444 samples) + staging metadata
#   2. Define adjacent transitions:
#      Fibrosis: F0→F1, F1→F2, F2→F3, F3→F4
#      NAS groups: NAS0-1→NAS2-4, NAS2-4→NAS5, NAS5→NAS6-8
#   3. For each transition, run dream contrast (stage_i+1 vs stage_i)
#   4. Compute transition specificity index (tau) — genes changing at
#      ONLY one specific transition
#   5. Pathway enrichment per transition program (fgsea, Hallmark + KEGG)
#   6. Save per-transition results + unified transition matrix
#
# Input:
#   - results/integration/merged_dge.rds (1,444 samples)
#   - results/staging_classifier/modeling_metadata.csv
#
# Output (all to results/progression/):
#   - transition_fib_dream_results.csv   (all fibrosis transition DE)
#   - transition_nas_dream_results.csv   (all NAS transition DE)
#   - transition_programs.csv            (gene × transition with tau)
#   - transition_pathway_enrichment.csv  (fgsea per transition)
#   - transition_summary.csv             (per-transition DEG counts)
#
# DEG THRESHOLD SYSTEM (Two-Tier):
#   Tier 1 (Primary): padj < 0.05, |logFC| > 0.5 — main dream DEGs (Script 05b; migrated 0.3 -> 0.5)
#   Tier 2 (Progression): padj < 0.05, no LFC filter — binary and adjacent contrasts
#   Rationale: Binary contrasts pool multiple stages, diluting per-gene fold changes.
#              Adjacent transitions have lower N (200-400 vs 1,444), making LFC estimates
#              noisier. Cross-contrast comparisons use rank-based enrichment (fgsea)
#              to avoid confounding power with biology.
#   See: figures/supplementary/sensitivity/figS_deg_threshold_landscape.pdf
#
# SLURM: cpu, 16 CPUs, 120G RAM, 48h
# Env:   micromamba activate rnaseq
#
# Usage:
#   sbatch --job-name=stg115_transitions \
#          --partition=cpu --cpus-per-task=16 --mem=120G --time=48:00:00 \
#          --output=logs/115_transitions_%j.out \
#          --error=logs/115_transitions_%j.err \
#          --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && \
#                  micromamba activate rnaseq && \
#                  cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \
#                  Rscript 115_transition_programs.R'"
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
})

# Force injection into lme4 namespace BEFORE loading variancePartition
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}

suppressPackageStartupMessages({
  library(variancePartition)
  library(BiocParallel)
  library(fgsea)
  library(msigdbr)
})

set.seed(42)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results")
OUTDIR <- file.path(RDIR, "progression")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(INT, "scripts/logs"), showWarnings = FALSE, recursive = TRUE)

NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
register(MulticoreParam(NCORES))

cat("=== 115: Stage-Transition Gene Programs ===\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("Cores:", NCORES, "\n\n")

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
cat("Loading merged DGE...\n")
dge <- readRDS(file.path(RDIR, "integration/merged_dge.rds"))
cat("  DGE:", nrow(dge), "genes x", ncol(dge), "samples\n")

cat("Loading staging metadata...\n")
meta <- fread(file.path(RDIR, "staging_classifier/modeling_metadata.csv"))
cat("  Metadata:", nrow(meta), "samples\n")

# Match metadata to DGE
meta <- meta[match(colnames(dge), meta$sample_id), ]
stopifnot(all(meta$sample_id == colnames(dge)))
cat("  Matched:", nrow(meta), "samples\n\n")

# ---------------------------------------------------------------------------
# 2. Define transition contrasts
# ---------------------------------------------------------------------------

# --- Fibrosis transitions (F0→F1→F2→F3→F4) ---
# Samples with fibrosis staging
fib_samples <- !is.na(meta$fibrosis_stage)
cat("Fibrosis-staged samples:", sum(fib_samples), "\n")
cat("Per-stage:\n")
print(table(meta$fibrosis_stage[fib_samples]))

fib_transitions <- list(
  "F0_to_F1" = c("0", "1"),
  "F1_to_F2" = c("1", "2"),
  "F2_to_F3" = c("2", "3"),
  "F3_to_F4" = c("3", "4")
)

# --- NAS group transitions ---
# Group NAS into biologically meaningful categories for sufficient power
# NAS 0-1 (78): Minimal/no disease
# NAS 2-4 (306): Steatosis/mild NAFL
# NAS 5 (139): Borderline NASH
# NAS 6-8 (137): Definite NASH
meta$nas_group_transition <- NA_character_
meta$nas_group_transition[!is.na(meta$nas_score) & meta$nas_score <= 1] <- "NAS_0_1"
meta$nas_group_transition[!is.na(meta$nas_score) & meta$nas_score >= 2 & meta$nas_score <= 4] <- "NAS_2_4"
meta$nas_group_transition[!is.na(meta$nas_score) & meta$nas_score == 5] <- "NAS_5"
meta$nas_group_transition[!is.na(meta$nas_score) & meta$nas_score >= 6] <- "NAS_6_8"

nas_samples <- !is.na(meta$nas_group_transition)
cat("\nNAS-staged samples:", sum(nas_samples), "\n")
cat("Per-group:\n")
print(table(meta$nas_group_transition[nas_samples]))

nas_transitions <- list(
  "NAS01_to_NAS24" = c("NAS_0_1", "NAS_2_4"),
  "NAS24_to_NAS5"  = c("NAS_2_4", "NAS_5"),
  "NAS5_to_NAS68"  = c("NAS_5",   "NAS_6_8")
)

# Also define individual NAS component transitions (steatosis, inflammation, ballooning)
# These provide finer-grained biology within NAS
component_cols <- c("steatosis_grade", "lobular_inflammation_grade", "ballooning_grade")

# ---------------------------------------------------------------------------
# 3. Run dream for each transition
# ---------------------------------------------------------------------------

run_transition_dream <- function(dge_full, meta_full, sample_mask, stage_col,
                                  stage_from, stage_to, transition_name) {
  # Select samples at either stage
  in_transition <- sample_mask &
    (as.character(meta_full[[stage_col]]) == stage_from |
     as.character(meta_full[[stage_col]]) == stage_to)
  if (sum(in_transition) < 20) {
    cat("  SKIP", transition_name, "- only", sum(in_transition), "samples\n")
    return(NULL)
  }

  dge_sub <- dge_full[, in_transition]
  meta_sub <- meta_full[in_transition, ]

  # Create binary contrast: 1 = higher stage, 0 = lower stage
  meta_sub$transition_group <- factor(
    ifelse(as.character(meta_sub[[stage_col]]) == stage_to, "Higher", "Lower"),
    levels = c("Lower", "Higher")
  )

  cat("  ", transition_name, ":", table(meta_sub$transition_group), "\n")

  # Determine which datasets are present in this subset
  ds_tab <- table(meta_sub$dataset)
  ds_present <- names(ds_tab[ds_tab >= 3])
  if (length(ds_present) < 2) {
    cat("    Only", length(ds_present), "dataset(s) — using fixed effect\n")
    # Fall back to fixed-effect model if <2 datasets
    meta_sub$dataset_f <- factor(meta_sub$dataset)
    info <- data.frame(
      transition_group = meta_sub$transition_group,
      dataset_f = meta_sub$dataset_f,
      stringsAsFactors = FALSE
    )
    rownames(info) <- meta_sub$sample_id
    form <- ~ transition_group + dataset_f
    use_random <- FALSE
  } else {
    # Filter to datasets with >= 3 samples
    keep_ds <- meta_sub$dataset %in% ds_present
    dge_sub <- dge_sub[, keep_ds]
    meta_sub <- meta_sub[keep_ds, ]
    meta_sub$transition_group <- droplevels(meta_sub$transition_group)

    info <- data.frame(
      transition_group = meta_sub$transition_group,
      dataset = factor(meta_sub$dataset),
      stringsAsFactors = FALSE
    )
    rownames(info) <- meta_sub$sample_id
    form <- ~ transition_group + (1|dataset)
    use_random <- TRUE
  }

  # Filter lowly expressed genes
  dge_sub <- dge_sub[filterByExpr(dge_sub, group = info$transition_group), ]
  dge_sub <- calcNormFactors(dge_sub, method = "TMM")
  cat("    Genes after filtering:", nrow(dge_sub), "\n")
  cat("    Samples:", ncol(dge_sub), "\n")

  # Run dream
  vobjDream <- tryCatch(
    voomWithDreamWeights(dge_sub, form, info, BPPARAM = MulticoreParam(NCORES)),
    error = function(e) {
      cat("    voomWithDreamWeights failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  if (is.null(vobjDream)) return(NULL)

  fit <- tryCatch(
    dream(vobjDream, form, info, BPPARAM = MulticoreParam(NCORES)),
    error = function(e) {
      cat("    dream failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  if (is.null(fit)) return(NULL)

  # Extract results for the transition contrast
  coef_name <- "transition_groupHigher"
  if (!coef_name %in% colnames(fit$coefficients)) {
    cat("    Coefficient not found. Available:", paste(colnames(fit$coefficients), collapse=", "), "\n")
    return(NULL)
  }

  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res$transition <- transition_name
  res$n_lower <- sum(meta_sub$transition_group == "Lower")
  res$n_higher <- sum(meta_sub$transition_group == "Higher")
  res$n_datasets <- length(ds_present)

  cat("    DEGs (padj<0.05):", sum(res$adj.P.Val < 0.05, na.rm = TRUE),
      "(up:", sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE),
      "down:", sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE), ")\n")

  return(as.data.table(res))
}

# --- Run all fibrosis transitions ---
cat("\n=== Fibrosis Transitions ===\n")
fib_results <- list()
for (tname in names(fib_transitions)) {
  stages <- fib_transitions[[tname]]
  fib_results[[tname]] <- run_transition_dream(
    dge, meta, fib_samples, "fibrosis_stage",
    stages[1], stages[2], tname
  )
}
fib_dt <- rbindlist(fib_results[!sapply(fib_results, is.null)], fill = TRUE)
cat("\nTotal fibrosis transition results:", nrow(fib_dt), "gene-transition pairs\n")

# --- Run all NAS group transitions ---
cat("\n=== NAS Group Transitions ===\n")
nas_results <- list()
for (tname in names(nas_transitions)) {
  stages <- nas_transitions[[tname]]
  nas_results[[tname]] <- run_transition_dream(
    dge, meta, nas_samples, "nas_group_transition",
    stages[1], stages[2], tname
  )
}
nas_dt <- rbindlist(nas_results[!sapply(nas_results, is.null)], fill = TRUE)
cat("\nTotal NAS transition results:", nrow(nas_dt), "gene-transition pairs\n")

# ---------------------------------------------------------------------------
# 4. Compute transition specificity index (tau)
# ---------------------------------------------------------------------------
cat("\n=== Computing Transition Specificity Index ===\n")

compute_tau <- function(dt, transitions) {
  # For each gene, compute tau across transitions
  # tau = (N - sum(x_i / x_max)) / (N - 1), where x_i = |t-stat| per transition
  # tau ranges 0 (ubiquitous) to 1 (transition-specific)

  # Pivot to gene × transition matrix of absolute t-statistics
  genes <- unique(dt$gene)
  trans <- unique(dt$transition)
  N <- length(trans)

  # Build t-stat matrix
  tmat <- matrix(0, nrow = length(genes), ncol = N,
                 dimnames = list(genes, trans))
  for (i in seq_len(nrow(dt))) {
    g <- dt$gene[i]
    tr <- dt$transition[i]
    tmat[g, tr] <- abs(dt$t[i])
  }

  # Compute tau per gene
  tau <- apply(tmat, 1, function(x) {
    xmax <- max(x, na.rm = TRUE)
    if (xmax == 0) return(0)
    (N - sum(x / xmax)) / (N - 1)
  })

  # Identify the "peak transition" per gene (where |t| is largest)
  peak_transition <- trans[apply(tmat, 1, which.max)]

  # Get the logFC and padj at peak transition
  peak_logfc <- sapply(seq_along(genes), function(i) {
    row <- dt[gene == genes[i] & transition == peak_transition[i]]
    if (nrow(row) > 0) row$logFC[1] else NA_real_
  })
  peak_padj <- sapply(seq_along(genes), function(i) {
    row <- dt[gene == genes[i] & transition == peak_transition[i]]
    if (nrow(row) > 0) row$adj.P.Val[1] else NA_real_
  })

  data.table(
    gene = genes,
    tau = tau,
    peak_transition = peak_transition,
    peak_logFC = peak_logfc,
    peak_padj = peak_padj,
    max_abs_t = apply(tmat, 1, max)
  )
}

fib_tau <- compute_tau(fib_dt, names(fib_transitions))
nas_tau <- compute_tau(nas_dt, names(nas_transitions))

cat("Fibrosis tau distribution:\n")
cat("  tau > 0.8 (transition-unique):", sum(fib_tau$tau > 0.8), "\n")
cat("  tau > 0.6 (moderately specific):", sum(fib_tau$tau > 0.6), "\n")
cat("  tau < 0.3 (ubiquitous):", sum(fib_tau$tau < 0.3), "\n")

cat("\nNAS tau distribution:\n")
cat("  tau > 0.8 (transition-unique):", sum(nas_tau$tau > 0.8), "\n")
cat("  tau > 0.6 (moderately specific):", sum(nas_tau$tau > 0.6), "\n")
cat("  tau < 0.3 (ubiquitous):", sum(nas_tau$tau < 0.3), "\n")

# Per-transition DEG counts
cat("\nPer-transition summary (padj < 0.05):\n")
for (tr in unique(fib_dt$transition)) {
  sub <- fib_dt[transition == tr]
  n_sig <- sum(sub$adj.P.Val < 0.05, na.rm = TRUE)
  cat("  ", tr, ":", n_sig, "DEGs\n")
}
for (tr in unique(nas_dt$transition)) {
  sub <- nas_dt[transition == tr]
  n_sig <- sum(sub$adj.P.Val < 0.05, na.rm = TRUE)
  cat("  ", tr, ":", n_sig, "DEGs\n")
}

# ---------------------------------------------------------------------------
# 4b. SAVE intermediate results before pathway enrichment
# ---------------------------------------------------------------------------
cat("\n=== Saving Intermediate Results (pre-enrichment) ===\n")
fwrite(fib_dt, file.path(OUTDIR, "transition_fib_dream_results.csv"))
cat("  Saved transition_fib_dream_results.csv (", nrow(fib_dt), "rows)\n")
fwrite(nas_dt, file.path(OUTDIR, "transition_nas_dream_results.csv"))
cat("  Saved transition_nas_dream_results.csv (", nrow(nas_dt), "rows)\n")

# Combine tau with dream results
all_dt <- rbind(fib_dt, nas_dt, fill = TRUE)
all_tau <- rbind(fib_tau, nas_tau, fill = TRUE)
fwrite(all_tau, file.path(OUTDIR, "transition_tau_index.csv"))
cat("  Saved transition_tau_index.csv (", nrow(all_tau), "rows)\n")

# Transition programs: merge tau with per-transition DE
# NOTE: merge within each stage type (fibrosis / NAS) to avoid Cartesian join —
# the same gene can appear in both fib_tau and nas_tau with different tau values.
tau_cols <- intersect(c("gene", "tau", "peak_transition", "max_transition", "max_abs_t", "peak_logFC", "peak_padj"),
                      names(fib_tau))
tp_fib <- merge(fib_dt, fib_tau[, ..tau_cols], by = "gene", all.x = TRUE)
tp_nas <- merge(nas_dt, nas_tau[, ..tau_cols], by = "gene", all.x = TRUE)
tp <- rbind(tp_fib, tp_nas, fill = TRUE)
fwrite(tp, file.path(OUTDIR, "transition_programs.csv"))
cat("  Saved transition_programs.csv (", nrow(tp), "rows)\n")

# ---------------------------------------------------------------------------
# 5. Pathway enrichment per transition (fgsea with Hallmark + KEGG)
# ---------------------------------------------------------------------------
cat("\n=== Pathway Enrichment per Transition ===\n")

# Get gene sets
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
kegg <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
  error = function(e) msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
)
pathways <- rbind(hallmark[, c("gs_name", "gene_symbol")],
                  kegg[, c("gs_name", "gene_symbol")])
pathway_list <- split(pathways$gene_symbol, pathways$gs_name)

# Run fgsea per transition
run_fgsea_transition <- function(dt, transition_name) {
  sub <- dt[transition == transition_name]
  # Use t-statistic as ranking metric (signed, captures direction + significance)
  ranks <- setNames(sub$t, sub$gene)
  ranks <- ranks[!is.na(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)

  res <- fgsea(pathways = pathway_list, stats = ranks,
               minSize = 15, maxSize = 500, nPermSimple = 10000)
  res$transition <- transition_name
  return(as.data.table(res))
}

# Fibrosis pathways
fib_pathways <- rbindlist(lapply(unique(fib_dt$transition), function(tr) {
  run_fgsea_transition(fib_dt, tr)
}), fill = TRUE)

# NAS pathways
nas_pathways <- rbindlist(lapply(unique(nas_dt$transition), function(tr) {
  run_fgsea_transition(nas_dt, tr)
}), fill = TRUE)

all_pathways <- rbind(fib_pathways, nas_pathways)
# Remove leadingEdge list column for CSV export
all_pathways$leadingEdge <- sapply(all_pathways$leadingEdge, paste, collapse = ";")

cat("Pathway enrichment: total", nrow(all_pathways), "pathway-transition tests\n")
cat("Significant (padj < 0.05):", sum(all_pathways$padj < 0.05, na.rm = TRUE), "\n")

# Top pathways per transition
cat("\nTop pathways per fibrosis transition:\n")
for (tr in unique(fib_pathways$transition)) {
  sub <- fib_pathways[transition == tr & padj < 0.05][order(pval)][1:min(5, .N)]
  if (nrow(sub) > 0) {
    cat("  ", tr, ":\n")
    for (i in seq_len(nrow(sub))) {
      cat("    ", sub$pathway[i], " (NES=", round(sub$NES[i], 2),
          ", padj=", formatC(sub$padj[i], format = "e", digits = 1), ")\n")
    }
  }
}

# ---------------------------------------------------------------------------
# 6. Combine and save results
# ---------------------------------------------------------------------------
cat("\n=== Saving Results ===\n")

# Combine fibrosis and NAS tau into unified transition program table
fib_tau$stage_type <- "fibrosis"
nas_tau$stage_type <- "NAS"
combined_tau <- rbind(fib_tau, nas_tau)

# Add transition-specific flags
combined_tau$is_transition_unique <- combined_tau$tau > 0.8 & combined_tau$peak_padj < 0.05
combined_tau$is_moderately_specific <- combined_tau$tau > 0.6 & combined_tau$peak_padj < 0.05

# Save all outputs
fwrite(fib_dt, file.path(OUTDIR, "transition_fib_dream_results.csv"))
cat("  Saved fibrosis transition dream results:", nrow(fib_dt), "rows\n")

fwrite(nas_dt, file.path(OUTDIR, "transition_nas_dream_results.csv"))
cat("  Saved NAS transition dream results:", nrow(nas_dt), "rows\n")

fwrite(combined_tau, file.path(OUTDIR, "transition_programs.csv"))
cat("  Saved transition programs:", nrow(combined_tau), "genes\n")

fwrite(all_pathways, file.path(OUTDIR, "transition_pathway_enrichment.csv"))
cat("  Saved pathway enrichment:", nrow(all_pathways), "rows\n")

# Summary table
summary_dt <- rbind(
  fib_dt[, .(n_deg_05 = sum(adj.P.Val < 0.05, na.rm=T),
             n_deg_05_lfc02 = sum(adj.P.Val < 0.05 & abs(logFC) > 0.2, na.rm=T),
             n_deg_05_lfc03 = sum(adj.P.Val < 0.05 & abs(logFC) > 0.3, na.rm=T),
             n_deg_05_lfc05 = sum(adj.P.Val < 0.05 & abs(logFC) > 0.5, na.rm=T),  # New canonical (LOO-CV stability)
             n_deg_01 = sum(adj.P.Val < 0.1, na.rm=T),
             n_up = sum(adj.P.Val < 0.05 & logFC > 0, na.rm=T),
             n_down = sum(adj.P.Val < 0.05 & logFC < 0, na.rm=T),
             n_samples = n_lower[1] + n_higher[1],
             n_datasets = n_datasets[1],
             stage_type = "fibrosis"),
         by = transition],
  nas_dt[, .(n_deg_05 = sum(adj.P.Val < 0.05, na.rm=T),
             n_deg_05_lfc02 = sum(adj.P.Val < 0.05 & abs(logFC) > 0.2, na.rm=T),
             n_deg_05_lfc03 = sum(adj.P.Val < 0.05 & abs(logFC) > 0.3, na.rm=T),
             n_deg_05_lfc05 = sum(adj.P.Val < 0.05 & abs(logFC) > 0.5, na.rm=T),  # New canonical (LOO-CV stability)
             n_deg_01 = sum(adj.P.Val < 0.1, na.rm=T),
             n_up = sum(adj.P.Val < 0.05 & logFC > 0, na.rm=T),
             n_down = sum(adj.P.Val < 0.05 & logFC < 0, na.rm=T),
             n_samples = n_lower[1] + n_higher[1],
             n_datasets = n_datasets[1],
             stage_type = "NAS"),
         by = transition]
)
fwrite(summary_dt, file.path(OUTDIR, "transition_summary.csv"))
cat("  Saved summary:", nrow(summary_dt), "transitions\n")

# Print transition-unique genes per transition
cat("\n=== Transition-Unique Genes (tau > 0.8, padj < 0.05) ===\n")
for (tr in unique(combined_tau$peak_transition)) {
  n <- sum(combined_tau$peak_transition == tr & combined_tau$is_transition_unique, na.rm=TRUE)
  cat("  ", tr, ":", n, "unique genes\n")
}

cat("\n=== F2→F3 Deep Dive (the critical reversibility boundary) ===\n")
f2f3_genes <- combined_tau[stage_type == "fibrosis" & peak_transition == "F2_to_F3" & is_transition_unique == TRUE][order(-max_abs_t)]
cat("  Transition-unique genes:", nrow(f2f3_genes), "\n")
if (nrow(f2f3_genes) > 0) {
  cat("  Top 20:\n")
  print(head(f2f3_genes[, .(gene, tau, peak_logFC, peak_padj, max_abs_t)], 20))
}

cat("\n=== 115: COMPLETE ===\n")
cat("Finished:", as.character(Sys.time()), "\n")
