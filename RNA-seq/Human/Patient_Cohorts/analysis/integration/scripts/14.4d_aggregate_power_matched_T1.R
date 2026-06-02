#!/usr/bin/env Rscript
# 14.4d_aggregate_power_matched_T1.R
# ---------------------------------------------------------------------------
# Aggregate per-iteration outputs from 14.4d_sex_power_matched_T1.R into
# canonical T1 outputs:
#   - per_gene_class_freq.csv  : gene × class frequency table (modal class +
#                                fraction of iterations in each class)
#   - summary.csv              : overall counts + modal class summary
#   - iter_metrics.csv         : per-iteration sample counts + sex_class totals
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ITER_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_power_matched/iterations"
OUT_DIR  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_power_matched"

# ----- iter_metrics ---------------------------------------------------------
metric_files <- list.files(ITER_DIR, pattern = "^pm_metrics_i\\d+\\.csv$", full.names = TRUE)
cat("Found", length(metric_files), "metrics files\n")
if (length(metric_files) == 0) stop("No pm_metrics files found in ", ITER_DIR)

iter_metrics <- rbindlist(lapply(metric_files, fread), fill = TRUE)
setorder(iter_metrics, iter)
fwrite(iter_metrics, file.path(OUT_DIR, "iter_metrics.csv"))
cat("Wrote iter_metrics.csv (", nrow(iter_metrics), " rows)\n", sep = "")

# ----- per-gene class frequency --------------------------------------------
cls_files <- list.files(ITER_DIR, pattern = "^pm_class_i\\d+\\.csv$", full.names = TRUE)
cat("Found", length(cls_files), "class files\n")
if (length(cls_files) == 0) stop("No pm_class files found in ", ITER_DIR)

B <- length(cls_files)
all_cls <- rbindlist(lapply(cls_files, function(f) {
  d <- fread(f, select = c("gene", "sex_class"))
  d[, iter := sub(".*_i(\\d+)\\.csv$", "\\1", basename(f))]
  d
}), fill = TRUE)

# wide: gene × sex_class counts
freq <- dcast(all_cls, gene ~ sex_class, fun.aggregate = length, value.var = "iter")
# ensure all 5 columns present
for (col in c("Concordant", "Female_biased", "Male_biased", "Divergent", "Unclassified")) {
  if (!col %in% names(freq)) freq[[col]] <- 0L
}

freq[, n_iter_observed := Concordant + Female_biased + Male_biased + Divergent + Unclassified]
freq[, frac_concordant    := Concordant    / n_iter_observed]
freq[, frac_female_biased := Female_biased / n_iter_observed]
freq[, frac_male_biased   := Male_biased   / n_iter_observed]
freq[, frac_divergent     := Divergent     / n_iter_observed]
freq[, frac_unclassified  := Unclassified  / n_iter_observed]

# modal class (most common across iterations)
class_cols <- c("Concordant", "Female_biased", "Male_biased", "Divergent", "Unclassified")
freq[, modal_class := class_cols[max.col(.SD, ties.method = "first")], .SDcols = class_cols]
freq[, modal_frac  := pmax(frac_concordant, frac_female_biased, frac_male_biased,
                            frac_divergent, frac_unclassified)]

# class-switching = number of distinct classes observed across iterations
freq[, n_distinct_classes := rowSums(.SD > 0L), .SDcols = class_cols]

setcolorder(freq, c("gene", "modal_class", "modal_frac", "n_distinct_classes",
                    "n_iter_observed",
                    class_cols,
                    "frac_concordant", "frac_female_biased", "frac_male_biased",
                    "frac_divergent", "frac_unclassified"))
fwrite(freq, file.path(OUT_DIR, "per_gene_class_freq.csv"))
cat("Wrote per_gene_class_freq.csv (", nrow(freq), " genes)\n", sep = "")

# ----- overall summary ------------------------------------------------------
summary_dt <- data.table(
  B_iter                 = B,
  n_genes                = nrow(freq),
  n_modal_concordant     = freq[modal_class == "Concordant", .N],
  n_modal_female_biased  = freq[modal_class == "Female_biased", .N],
  n_modal_male_biased    = freq[modal_class == "Male_biased", .N],
  n_modal_divergent      = freq[modal_class == "Divergent", .N],
  n_modal_unclassified   = freq[modal_class == "Unclassified", .N],
  # stable = same class in >= 80% of iterations
  n_stable_modal_ge80    = freq[modal_frac >= 0.8, .N],
  # high-switching = 3+ distinct classes
  n_high_switching_3plus = freq[n_distinct_classes >= 3, .N],
  mean_dimorphic_per_iter = mean(iter_metrics$n_dimorphic, na.rm = TRUE),
  median_dimorphic_per_iter = median(iter_metrics$n_dimorphic, na.rm = TRUE),
  mean_F_biased_per_iter  = mean(iter_metrics$n_female_biased, na.rm = TRUE),
  mean_M_biased_per_iter  = mean(iter_metrics$n_male_biased, na.rm = TRUE),
  FM_ratio_per_iter_mean  = mean(iter_metrics$n_female_biased / pmax(iter_metrics$n_male_biased, 1), na.rm = TRUE)
)
fwrite(summary_dt, file.path(OUT_DIR, "summary.csv"))
cat("Wrote summary.csv\n")

print(summary_dt)
cat("\nDone (T1 aggregation).\n")
