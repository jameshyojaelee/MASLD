#!/usr/bin/env Rscript
# 119_celltype_resolved_progression.R
# ---------------------------------------------------------------------------
# Cell-Type-Resolved Progression Programs:
# For each cell type, identify which genes change at each disease transition.
#
# This answers: "Which cell type drives each stage transition?"
#   - F0→F1: Expect hepatocyte steatosis + macrophage inflammation
#   - F2→F3: Expect stellate cell fibrogenesis
#   - F3→F4: Expect hepatocyte death + endothelial remodeling
#
# Uses BayesPrism cell-type-specific expression from Script 118.
# For each cell type × transition:
#   1. Wilcoxon test between adjacent stages
#   2. Compute transition specificity (which cell type drives each transition)
#   3. "Attribution" analysis: intrinsic vs composition-driven change
#
# Input:
#   - results/progression/cibersortx_celltype_expression/bayesprism_*.csv.gz
#   - results/staging_classifier/modeling_metadata.csv
#   - results/progression/transition_programs.csv (from Script 115)
#
# Output (to results/progression/):
#   - celltype_transition_programs.csv  (gene × celltype × transition)
#   - celltype_transition_summary.csv   (per celltype × transition DEG counts)
#   - celltype_attribution.csv          (fraction of change per cell type)
#
# SLURM: cpu, 16 CPUs, 120G RAM, 48h
# Env:   micromamba activate rnaseq
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

cat("=== 119: Cell-Type-Resolved Progression Programs ===\n")
cat(sprintf("Started: %s\n", Sys.time()))

# ── Paths ──────────────────────────────────────────────────────────────────
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INTEG <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
BP_DIR <- file.path(INTEG, "results/progression/cibersortx_celltype_expression")
OUTDIR <- file.path(INTEG, "results/progression")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

n_cores <- min(parallel::detectCores(), 16)

# ── 1. Load metadata ──────────────────────────────────────────────────────
cat("\nLoading metadata...\n")
meta <- fread(file.path(INTEG, "results/staging_classifier/modeling_metadata.csv"))
cat(sprintf("  Metadata: %d samples\n", nrow(meta)))

# ── 2. Define transitions ─────────────────────────────────────────────────
fib_transitions <- list(
  F0_to_F1 = c(0, 1),
  F1_to_F2 = c(1, 2),
  F2_to_F3 = c(2, 3),
  F3_to_F4 = c(3, 4)
)

# NAS transitions (grouped)
nas_transitions <- list(
  NAS_0_1_to_NAS_2_4 = list(from = c(0, 1), to = c(2, 3, 4)),
  NAS_2_4_to_NAS_5   = list(from = c(2, 3, 4), to = c(5)),
  NAS_5_to_NAS_6_8   = list(from = c(5), to = c(6, 7, 8))
)

# ── 3. Find available cell-type expression files ──────────────────────────
cat("\nFinding BayesPrism cell-type expression files...\n")
bp_files <- list.files(BP_DIR, pattern = "^bayesprism_.*\\.csv\\.gz$", full.names = TRUE)
bp_files <- bp_files[!grepl("proportions|summary|result", basename(bp_files))]

if (length(bp_files) == 0) {
  stop("No BayesPrism cell-type expression files found. Run Script 118 first.")
}

cell_type_names <- gsub("bayesprism_|\\.csv\\.gz", "", basename(bp_files))
cat(sprintf("  Found %d cell types: %s\n", length(cell_type_names),
            paste(cell_type_names, collapse = ", ")))

# ── 4. Load proportions ───────────────────────────────────────────────────
prop_file <- file.path(BP_DIR, "bayesprism_proportions.csv")
if (file.exists(prop_file)) {
  proportions <- fread(prop_file)
  cat(sprintf("  Proportions: %d samples x %d types\n",
              nrow(proportions), ncol(proportions) - 1))

  # Compute median proportion per cell type (across ALL samples)
  prop_cols <- setdiff(names(proportions), "sample_id")
  median_proportions <- sapply(prop_cols, function(ct) {
    median(proportions[[ct]], na.rm = TRUE)
  })
  cat("\n  Median proportions per cell type:\n")
  for (ct in sort(names(median_proportions))) {
    flag <- if (median_proportions[ct] < 0.01) " ** LOW (<1%)" else ""
    cat(sprintf("    %-15s: %.4f (%.2f%%)%s\n", ct,
                median_proportions[ct], median_proportions[ct] * 100, flag))
  }

  # Flag cell types with median proportion < 1% as low confidence
  LOW_PROP_THRESHOLD <- 0.01
  low_confidence_celltypes <- names(median_proportions)[median_proportions < LOW_PROP_THRESHOLD]
  cat(sprintf("\n  Low-confidence cell types (median < %.0f%%): %s\n",
              LOW_PROP_THRESHOLD * 100,
              if (length(low_confidence_celltypes) > 0)
                paste(low_confidence_celltypes, collapse = ", ")
              else "none"))
} else {
  proportions <- NULL
  median_proportions <- NULL
  low_confidence_celltypes <- character(0)
  cat("  No proportions file found (will skip attribution)\n")
}

# ── 5. Per-cell-type transition DE ─────────────────────────────────────────
cat("\n=== Cell-Type-Resolved Transition DE ===\n")

# Minimum effect size filter to prevent BayesPrism deconvolution noise
# from producing biologically implausible DEG counts (e.g., 99.7% DEGs
# for rare cell types). With N=300+ samples, Wilcoxon detects even tiny
# noise differences as significant — requiring |LFC| > 0.25 ensures
# only biologically meaningful expression changes are called DEGs.
MIN_LFC <- 0.25
cat(sprintf("  Effect size filter: |LFC| > %.2f in addition to padj < 0.1\n", MIN_LFC))

run_celltype_transition_de <- function(ct_name, ct_file) {
  cat(sprintf("\n--- %s ---\n", ct_name))

  # Load cell-type-specific expression
  ct_expr <- fread(ct_file)
  sample_ids <- ct_expr$sample_id
  ct_expr <- as.matrix(ct_expr[, -1, with = FALSE])
  rownames(ct_expr) <- sample_ids
  cat(sprintf("  Expression: %d samples x %d genes\n", nrow(ct_expr), ncol(ct_expr)))

  # Align with metadata
  common <- intersect(rownames(ct_expr), meta$sample_id)
  ct_expr <- ct_expr[common, , drop = FALSE]
  meta_aligned <- meta[match(common, meta$sample_id), ]

  results <- list()

  # --- Fibrosis transitions ---
  for (trans_name in names(fib_transitions)) {
    stages <- fib_transitions[[trans_name]]
    s_from <- stages[1]
    s_to <- stages[2]

    idx_from <- which(meta_aligned$fibrosis_stage == s_from)
    idx_to <- which(meta_aligned$fibrosis_stage == s_to)

    if (length(idx_from) < 10 || length(idx_to) < 10) {
      cat(sprintf("    %s: skipped (n=%d, %d)\n", trans_name, length(idx_from), length(idx_to)))
      next
    }

    # Wilcoxon test per gene (parallelized)
    gene_pvals <- mclapply(seq_len(ncol(ct_expr)), function(g) {
      x_from <- ct_expr[idx_from, g]
      x_to <- ct_expr[idx_to, g]

      # Skip if no variance
      if (sd(c(x_from, x_to)) < 1e-10) return(c(pval = 1, lfc = 0))

      wt <- tryCatch(
        wilcox.test(x_to, x_from, alternative = "two.sided"),
        error = function(e) list(p.value = 1)
      )

      # Log fold change (mean difference in log space or pseudocount)
      mean_from <- mean(x_from)
      mean_to <- mean(x_to)
      if (mean_from > 0 && mean_to > 0) {
        lfc <- log2(mean_to / mean_from)
      } else {
        lfc <- mean_to - mean_from
      }

      c(pval = wt$p.value, lfc = lfc)
    }, mc.cores = n_cores)

    pvals <- sapply(gene_pvals, `[`, "pval")
    lfcs <- sapply(gene_pvals, `[`, "lfc")
    padj <- p.adjust(pvals, method = "BH")

    n_sig_raw <- sum(padj < 0.1, na.rm = TRUE)
    n_sig <- sum(padj < 0.1 & abs(lfcs) > MIN_LFC, na.rm = TRUE)
    cat(sprintf("    %s: %d DEGs (padj<0.1 & |LFC|>%.2f) [%d before LFC filter] [n=%d vs %d]\n",
                trans_name, n_sig, MIN_LFC, n_sig_raw, length(idx_from), length(idx_to)))

    results[[paste0(trans_name)]] <- data.table(
      gene = colnames(ct_expr),
      cell_type = ct_name,
      transition = trans_name,
      transition_type = "fibrosis",
      lfc = lfcs,
      pval = pvals,
      padj = padj,
      n_from = length(idx_from),
      n_to = length(idx_to)
    )
  }

  # --- NAS transitions ---
  for (trans_name in names(nas_transitions)) {
    spec <- nas_transitions[[trans_name]]

    idx_from <- which(meta_aligned$nas_score %in% spec$from)
    idx_to <- which(meta_aligned$nas_score %in% spec$to)

    if (length(idx_from) < 10 || length(idx_to) < 10) {
      cat(sprintf("    %s: skipped (n=%d, %d)\n", trans_name, length(idx_from), length(idx_to)))
      next
    }

    gene_pvals <- mclapply(seq_len(ncol(ct_expr)), function(g) {
      x_from <- ct_expr[idx_from, g]
      x_to <- ct_expr[idx_to, g]
      if (sd(c(x_from, x_to)) < 1e-10) return(c(pval = 1, lfc = 0))
      wt <- tryCatch(
        wilcox.test(x_to, x_from, alternative = "two.sided"),
        error = function(e) list(p.value = 1)
      )
      mean_from <- mean(x_from)
      mean_to <- mean(x_to)
      if (mean_from > 0 && mean_to > 0) {
        lfc <- log2(mean_to / mean_from)
      } else {
        lfc <- mean_to - mean_from
      }
      c(pval = wt$p.value, lfc = lfc)
    }, mc.cores = n_cores)

    pvals <- sapply(gene_pvals, `[`, "pval")
    lfcs <- sapply(gene_pvals, `[`, "lfc")
    padj <- p.adjust(pvals, method = "BH")

    n_sig_raw <- sum(padj < 0.1, na.rm = TRUE)
    n_sig <- sum(padj < 0.1 & abs(lfcs) > MIN_LFC, na.rm = TRUE)
    cat(sprintf("    %s: %d DEGs (padj<0.1 & |LFC|>%.2f) [%d before LFC filter] [n=%d vs %d]\n",
                trans_name, n_sig, MIN_LFC, n_sig_raw, length(idx_from), length(idx_to)))

    results[[trans_name]] <- data.table(
      gene = colnames(ct_expr),
      cell_type = ct_name,
      transition = trans_name,
      transition_type = "NAS",
      lfc = lfcs,
      pval = pvals,
      padj = padj,
      n_from = length(idx_from),
      n_to = length(idx_to)
    )
  }

  rbindlist(results)
}

# Run for each cell type
all_results <- list()
for (i in seq_along(bp_files)) {
  all_results[[i]] <- run_celltype_transition_de(cell_type_names[i], bp_files[i])
}
ct_de <- rbindlist(all_results)
cat(sprintf("\n  Total results: %d rows\n", nrow(ct_de)))

# ── 6. Cell-type transition summary ───────────────────────────────────────
cat("\n=== Transition Summary ===\n")
cat(sprintf("  DEG criteria: padj < 0.1 AND |LFC| > %.2f\n", MIN_LFC))

summary_dt <- ct_de[, .(
  n_degs_01 = sum(padj < 0.1 & abs(lfc) > MIN_LFC, na.rm = TRUE),
  n_degs_01_nolfc = sum(padj < 0.1, na.rm = TRUE),
  n_degs_005 = sum(padj < 0.05 & abs(lfc) > MIN_LFC, na.rm = TRUE),
  n_up = sum(padj < 0.1 & lfc > MIN_LFC, na.rm = TRUE),
  n_down = sum(padj < 0.1 & lfc < -MIN_LFC, na.rm = TRUE),
  mean_abs_lfc = mean(abs(lfc[padj < 0.1 & abs(lfc) > MIN_LFC]), na.rm = TRUE),
  n_tested = .N
), by = .(cell_type, transition, transition_type)]

# Add low-confidence flag based on cell-type median proportion
summary_dt[, low_confidence := cell_type %in% low_confidence_celltypes]
summary_dt[, pct_deg := round(n_degs_01 / n_tested * 100, 1)]

# Print nicely
for (trans in unique(summary_dt$transition)) {
  sub <- summary_dt[transition == trans]
  cat(sprintf("\n  %s:\n", trans))
  for (j in seq_len(nrow(sub))) {
    flag <- if (sub$low_confidence[j]) " [LOW CONFIDENCE]" else ""
    cat(sprintf("    %-15s: %4d DEGs (%4.1f%%) (up=%d, down=%d)%s\n",
                sub$cell_type[j], sub$n_degs_01[j], sub$pct_deg[j],
                sub$n_up[j], sub$n_down[j], flag))
  }
}

# ── 7. Cell-type attribution per transition ───────────────────────────────
cat("\n=== Cell-Type Attribution ===\n")
cat("(Which cell type contributes most to each transition?)\n")
cat("  NOTE: Attribution computed from FILTERED DEGs (padj<0.1 & |LFC|>0.25)\n")

# Compute attribution only from high-confidence cell types where possible
attribution <- summary_dt[, {
  hc <- .SD[low_confidence == FALSE]
  if (nrow(hc) > 0 && sum(hc$n_degs_01) > 0) {
    # Use only high-confidence cell types for attribution
    .(total_degs = sum(n_degs_01),
      total_degs_hc = sum(hc$n_degs_01),
      max_celltype = hc$cell_type[which.max(hc$n_degs_01)],
      max_degs = max(hc$n_degs_01),
      fraction_from_max = max(hc$n_degs_01) / sum(hc$n_degs_01),
      n_celltypes = .N,
      n_celltypes_hc = nrow(hc))
  } else {
    # Fallback: use all cell types (but flag as all-low-confidence)
    .(total_degs = sum(n_degs_01),
      total_degs_hc = 0L,
      max_celltype = cell_type[which.max(n_degs_01)],
      max_degs = max(n_degs_01),
      fraction_from_max = if (sum(n_degs_01) > 0) max(n_degs_01) / sum(n_degs_01) else 0,
      n_celltypes = .N,
      n_celltypes_hc = 0L)
  }
}, by = .(transition)]

for (i in seq_len(nrow(attribution))) {
  cat(sprintf("  %s: driven by %s (%d/%d DEGs, %.0f%%; %d/%d high-confidence cell types)\n",
              attribution$transition[i],
              attribution$max_celltype[i],
              attribution$max_degs[i],
              attribution$total_degs[i],
              attribution$fraction_from_max[i] * 100,
              attribution$n_celltypes_hc[i],
              attribution$n_celltypes[i]))
}

# ── 8. Composition vs intrinsic attribution ───────────────────────────────
if (!is.null(proportions)) {
  cat("\n=== Composition vs Intrinsic Change ===\n")

  comp_attr <- list()
  for (trans_name in names(fib_transitions)) {
    stages <- fib_transitions[[trans_name]]
    s_from <- stages[1]
    s_to <- stages[2]

    idx_from <- which(meta$fibrosis_stage == s_from)
    idx_to <- which(meta$fibrosis_stage == s_to)

    if (length(idx_from) < 10 || length(idx_to) < 10) next

    sid_from <- meta$sample_id[idx_from]
    sid_to <- meta$sample_id[idx_to]

    prop_from <- proportions[match(sid_from, proportions$sample_id),
                             -1, with = FALSE]
    prop_to <- proportions[match(sid_to, proportions$sample_id),
                           -1, with = FALSE]

    for (ct in names(prop_from)) {
      mean_from <- mean(prop_from[[ct]], na.rm = TRUE)
      mean_to <- mean(prop_to[[ct]], na.rm = TRUE)
      pct_change <- (mean_to - mean_from) / max(mean_from, 1e-5) * 100

      if (mean_from > 0.01 || mean_to > 0.01) {
        wt <- tryCatch(
          wilcox.test(prop_to[[ct]], prop_from[[ct]])$p.value,
          error = function(e) 1
        )
      } else {
        wt <- 1
      }

      comp_attr[[paste(trans_name, ct)]] <- data.frame(
        transition = trans_name,
        cell_type = ct,
        mean_prop_from = mean_from,
        mean_prop_to = mean_to,
        pct_change = pct_change,
        comp_pval = wt,
        stringsAsFactors = FALSE
      )
    }
  }

  comp_df <- do.call(rbind, comp_attr)
  comp_df$comp_padj <- p.adjust(comp_df$comp_pval, method = "BH")

  sig_comp <- comp_df[comp_df$comp_padj < 0.05 & abs(comp_df$pct_change) > 10, ]
  if (nrow(sig_comp) > 0) {
    cat("  Significant composition changes:\n")
    sig_comp <- sig_comp[order(sig_comp$comp_padj), ]
    for (i in seq_len(min(nrow(sig_comp), 10))) {
      cat(sprintf("    %s → %s: %.1f%% → %.1f%% (%+.1f%%, padj=%.2e)\n",
                  sig_comp$transition[i], sig_comp$cell_type[i],
                  sig_comp$mean_prop_from[i] * 100, sig_comp$mean_prop_to[i] * 100,
                  sig_comp$pct_change[i], sig_comp$comp_padj[i]))
    }
  }

  fwrite(comp_df, file.path(OUTDIR, "celltype_composition_changes.csv"))
  cat("  Saved celltype_composition_changes.csv\n")
}

# ── 9. Save results ───────────────────────────────────────────────────────
cat("\n=== Saving Results ===\n")

fwrite(ct_de, file.path(OUTDIR, "celltype_transition_programs.csv"))
cat(sprintf("  Saved celltype_transition_programs.csv (%d rows)\n", nrow(ct_de)))

fwrite(summary_dt, file.path(OUTDIR, "celltype_transition_summary.csv"))
cat("  Saved celltype_transition_summary.csv\n")

fwrite(attribution, file.path(OUTDIR, "celltype_attribution.csv"))
cat("  Saved celltype_attribution.csv\n")

cat(sprintf("\n=== 119: COMPLETE (%s) ===\n", Sys.time()))
