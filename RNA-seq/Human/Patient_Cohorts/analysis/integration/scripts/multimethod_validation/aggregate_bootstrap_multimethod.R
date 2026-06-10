#!/usr/bin/env Rscript
# aggregate_bootstrap_multimethod.R
# ===========================================================================
# Aggregate the per-iteration multi-method bootstrap outputs into:
#   1. selection_frequency.csv  — per gene per method, selection frequency =
#        mean(padj < 0.1) over the iters where that gene was ELIGIBLE (i.e.
#        appeared in that iter's output for that method). Eligibility is
#        iter-dependent (filterByExpr / DESeq2 prefilter / metafor min_K
#        coverage all vary by subsample), so we weight by the count of eligible
#        iters per (gene, method), NOT the raw total iteration count.
#   2. stability_summary.csv    — per method: n iters, n genes freq>0.5,
#        freq>0.9, and a summary of per-gene logFC CV (sd/|mean| across iters).
#   3. method_concordance.csv   — pairwise Spearman of selection frequencies
#        and Jaccard of the freq>0.5 gene sets between methods.
#
# Reads:  results/integration/multimethod_validation/bootstrap/iter_*.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RDIR    <- file.path(PROJECT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
BOOT_DIR <- file.path(RDIR, "multimethod_validation/bootstrap")

PADJ_THRESH <- 0.1  # selection threshold (matches LOO/atlas exploratory padj)

# --- Read all iteration files -----------------------------------------------
# Match default dream+deseq2 files (iter_{N}.csv) AND any single-method suffixed
# run (e.g. iter_metafor_{N}.csv) so the 3-way merge happens automatically.
files <- list.files(BOOT_DIR, pattern = "^iter.*_\\d+\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No iter_*.csv files in ", BOOT_DIR)
cat("Reading", length(files), "iteration files from", BOOT_DIR, "\n")

dt <- rbindlist(lapply(files, fread), fill = TRUE)
# Empty-guard files contribute 0 rows; drop any malformed/empty reads.
dt <- dt[!is.na(gene) & nzchar(gene)]
if (nrow(dt) == 0) stop("All iteration files were empty — nothing to aggregate")

methods <- sort(unique(dt$method))
n_iters_total <- uniqueN(dt$iter)
cat("Methods:", paste(methods, collapse = ", "),
    "| total iters with data:", n_iters_total, "\n")

# ELIGIBLE = CALLABLE (R3 P1-2): a (gene, method, iter) counts as eligible only
# if the method actually CALLED it (non-NA padj). DESeq2's independent filtering
# sets padj=NA for low-count genes every iter; counting those as
# eligible-but-not-selected would deflate DESeq2's selection frequency relative
# to dream (which has no NA-padj). Restricting the denominator to callable rows
# makes the per-method frequency apples-to-apples.
dt <- dt[!is.na(padj)]                                  # drop non-callable rows
dt[, selected := as.integer(padj < PADJ_THRESH)]

# ===========================================================================
# 1. Selection frequency (over CALLABLE iters only)
# ===========================================================================
freq_long <- dt[, .(
  freq       = mean(selected),          # over callable (non-NA-padj) iters
  n_eligible = .N                       # number of iters this gene was callable
), by = .(method, gene)]

# Wide: gene + {method}_freq + {method}_n_eligible
freq_wide <- dcast(freq_long, gene ~ method,
                   value.var = c("freq", "n_eligible"), fill = NA)
# dcast names columns freq_<method> / n_eligible_<method>; rename to the
# requested {method}_freq / {method}_n_eligible layout.
for (m in methods) {
  if (paste0("freq_", m) %in% names(freq_wide))
    setnames(freq_wide, paste0("freq_", m), paste0(m, "_freq"))
  if (paste0("n_eligible_", m) %in% names(freq_wide))
    setnames(freq_wide, paste0("n_eligible_", m), paste0(m, "_n_eligible"))
}
setorder(freq_wide, gene)
fwrite(freq_wide, file.path(BOOT_DIR, "selection_frequency.csv"))
cat("Wrote selection_frequency.csv:", nrow(freq_wide), "genes\n")

# ===========================================================================
# 2. Stability summary (per method)
# ===========================================================================
# Per-gene logFC CV across iters where the gene appeared (need >= 2 iters).
cv_long <- dt[, {
  lfc <- logFC[is.finite(logFC)]
  cv  <- if (length(lfc) >= 2 && abs(mean(lfc)) > 0)
           sd(lfc) / abs(mean(lfc)) else NA_real_
  .(logFC_cv = cv)
}, by = .(method, gene)]

stability <- merge(
  freq_long[, .(
    n_genes_total = .N,
    n_freq_gt_0.5 = sum(freq > 0.5),
    n_freq_gt_0.9 = sum(freq > 0.9)
  ), by = method],
  cv_long[, .(
    mean_logFC_cv   = mean(logFC_cv, na.rm = TRUE),
    median_logFC_cv = median(logFC_cv, na.rm = TRUE),
    n_genes_cv      = sum(is.finite(logFC_cv))
  ), by = method],
  by = "method", all = TRUE)

# Iteration count per method (eligibility makes this potentially method-specific).
n_iters_by_method <- dt[, .(n_iters = uniqueN(iter)), by = method]
stability <- merge(n_iters_by_method, stability, by = "method", all = TRUE)
setcolorder(stability, c("method", "n_iters", "n_genes_total",
                         "n_freq_gt_0.5", "n_freq_gt_0.9",
                         "mean_logFC_cv", "median_logFC_cv", "n_genes_cv"))
setorder(stability, method)
fwrite(stability, file.path(BOOT_DIR, "stability_summary.csv"))
cat("Wrote stability_summary.csv:", nrow(stability), "methods\n")

# ===========================================================================
# 3. Pairwise method concordance (Spearman of freqs + Jaccard of freq>0.5 sets)
# ===========================================================================
concordance <- data.table(method_a = character(), method_b = character(),
                          spearman_freq = numeric(), n_common_genes = integer(),
                          jaccard_freq_gt_0.5 = numeric(),
                          n_a_gt_0.5 = integer(), n_b_gt_0.5 = integer())

if (length(methods) >= 2) {
  pairs <- combn(methods, 2, simplify = FALSE)
  rows <- lapply(pairs, function(p) {
    ma <- p[1]; mb <- p[2]
    fa <- freq_long[method == ma, .(gene, freq_a = freq)]
    fb <- freq_long[method == mb, .(gene, freq_b = freq)]
    # Spearman over genes eligible for BOTH methods (common universe).
    m  <- merge(fa, fb, by = "gene")
    sp <- if (nrow(m) >= 3)
            suppressWarnings(cor(m$freq_a, m$freq_b, method = "spearman",
                                 use = "pairwise.complete.obs"))
          else NA_real_
    set_a <- freq_long[method == ma & freq > 0.5, gene]
    set_b <- freq_long[method == mb & freq > 0.5, gene]
    un    <- length(union(set_a, set_b))
    jac   <- if (un > 0) length(intersect(set_a, set_b)) / un else NA_real_
    data.table(method_a = ma, method_b = mb,
               spearman_freq = round(sp, 4), n_common_genes = nrow(m),
               jaccard_freq_gt_0.5 = round(jac, 4),
               n_a_gt_0.5 = length(set_a), n_b_gt_0.5 = length(set_b))
  })
  concordance <- rbindlist(rows)
} else {
  cat("Only", length(methods), "method present — concordance table will be empty\n")
}
fwrite(concordance, file.path(BOOT_DIR, "method_concordance.csv"))
cat("Wrote method_concordance.csv:", nrow(concordance), "method pairs\n")

cat("Aggregation complete.\n")
