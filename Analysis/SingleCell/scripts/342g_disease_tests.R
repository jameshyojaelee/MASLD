#!/usr/bin/env Rscript
# 342g_disease_tests.R - Phase 3: disease tests on snRNA + bulk polyploid signature scores.
#
# snRNA: mixed-effects model across 5 disease cohorts (cohort random intercept)
# bulk: ordinal F0-F4 trend + NAS continuous test, mixed-effects with cohort random intercept

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(patchwork)
})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
setwd(ROOT)
OUT_DIR <- "Analysis/SingleCell/results_gpu_v2/ploidy"
FIG_DIR <- "figures/supplementary"

# ---------- snRNA path ----------
sn <- fread(file.path(OUT_DIR, "signature_scores_snRNA_sample_summary.csv"))
cat("snRNA samples:", nrow(sn), "\n")

# Annotate condition as binary
sn[, condition_binary := fcase(
  grepl("Healthy|HEALTHY|Normal|NORMAL", condition), "control",
  default = "disease"
)]
# Reference cohorts (Liver_Atlas, GSE136103) - keep but flag as reference
sn[, is_reference := dataset %in% c("Liver_Atlas", "GSE136103")]

# Merge GSE244832 donor pairing to get NORMAL/MASL/MASH per sample
pair <- fread("data/GSE244832/metadata/donor_pairing.csv")
pair_long <- pair[, .(rna_srr = unlist(strsplit(rna_srrs, ";"))),
                  by = .(donor_id, condition_3lvl = condition)]
sn[dataset == "GSE244832", condition_3lvl := pair_long$condition_3lvl[match(sample, pair_long$rna_srr)]]
# Update condition_binary for GSE244832 using donor_pairing labels
sn[dataset == "GSE244832" & condition_3lvl == "NORMAL", condition_binary := "control"]
sn[dataset == "GSE244832" & condition_3lvl %in% c("MASL", "MASH"), condition_binary := "disease"]

sig_cols <- grep("^sig_|^polyploid_index", names(sn), value = TRUE)
cat("\nsignature columns:", paste(sig_cols, collapse = ", "), "\n")

cat("\n=== snRNA condition breakdown ===\n")
print(table(sn$condition_binary, sn$dataset))
print(table(sn$condition_3lvl, useNA = "ifany"))

# Filter to non-reference disease cohorts
sn_test <- sn[!is_reference & !is.na(condition_binary)]
cat("\nsnRNA test set (disease cohorts):", nrow(sn_test), "samples\n")

# Mixed-effects: signature ~ condition + (1|dataset)
mixed_results <- list()
for (sc in sig_cols) {
  if (sum(!is.na(sn_test[[sc]])) < 10) next
  dat <- sn_test[!is.na(get(sc))]
  dat$y <- dat[[sc]]
  m <- tryCatch(
    lmer(y ~ condition_binary + (1 | dataset), data = dat),
    error = function(e) NULL
  )
  if (is.null(m)) next
  fe <- summary(m)$coefficients
  if (!"condition_binarydisease" %in% rownames(fe)) next
  est <- fe["condition_binarydisease", "Estimate"]
  se  <- fe["condition_binarydisease", "Std. Error"]
  p   <- fe["condition_binarydisease", "Pr(>|t|)"]
  # also pool counts
  n_ctrl <- sum(dat$condition_binary == "control")
  n_dis  <- sum(dat$condition_binary == "disease")
  mixed_results[[sc]] <- data.table(
    signature = sc, n_control = n_ctrl, n_disease = n_dis,
    mean_control = mean(dat$y[dat$condition_binary == "control"]),
    mean_disease = mean(dat$y[dat$condition_binary == "disease"]),
    estimate = est, std_error = se, p_value = p
  )
}
mixed_dt <- rbindlist(mixed_results)
mixed_dt[, padj := p.adjust(p_value, method = "BH")]
mixed_dt <- mixed_dt[order(p_value)]
fwrite(mixed_dt, file.path(OUT_DIR, "snRNA_disease_mixed_effects.csv"))
cat("\n=== snRNA mixed-effects: signature ~ condition + (1|dataset) ===\n")
print(mixed_dt)

# Per-cohort Wilcoxon (excluding controls-only cohort GSE185477)
cohort_test <- list()
for (sc in sig_cols) {
  for (ds in unique(sn_test$dataset)) {
    sub <- sn_test[dataset == ds & !is.na(get(sc)) & !is.na(condition_binary)]
    n_c <- sum(sub$condition_binary == "control"); n_d <- sum(sub$condition_binary == "disease")
    if (n_c < 2 || n_d < 2) next
    w <- tryCatch(wilcox.test(sub[[sc]] ~ sub$condition_binary, exact = FALSE), error = function(e) NULL)
    if (is.null(w)) next
    cohort_test[[paste(ds, sc, sep = "::")]] <- data.table(
      cohort = ds, signature = sc, n_control = n_c, n_disease = n_d,
      mean_control = mean(sub[[sc]][sub$condition_binary == "control"]),
      mean_disease = mean(sub[[sc]][sub$condition_binary == "disease"]),
      direction = sign(mean(sub[[sc]][sub$condition_binary == "disease"]) -
                       mean(sub[[sc]][sub$condition_binary == "control"])),
      wilcox_p = w$p.value
    )
  }
}
cohort_dt <- rbindlist(cohort_test)
fwrite(cohort_dt, file.path(OUT_DIR, "snRNA_per_cohort_wilcoxon.csv"))
cat("\n=== Per-cohort Wilcoxon (top 20 by p) ===\n")
print(cohort_dt[order(wilcox_p)][1:20])

# 3-level test within GSE244832 (NORMAL/MASL/MASH)
gse <- sn_test[dataset == "GSE244832" & !is.na(condition_3lvl)]
cat("\nGSE244832 3-level n:", table(gse$condition_3lvl), "\n")
gse_3lvl <- list()
for (sc in sig_cols) {
  if (sum(!is.na(gse[[sc]])) < 5) next
  k <- kruskal.test(gse[[sc]] ~ factor(gse$condition_3lvl, levels = c("NORMAL","MASL","MASH")))
  gse_3lvl[[sc]] <- data.table(
    signature = sc,
    mean_NORMAL = mean(gse[[sc]][gse$condition_3lvl == "NORMAL"], na.rm = TRUE),
    mean_MASL   = mean(gse[[sc]][gse$condition_3lvl == "MASL"], na.rm = TRUE),
    mean_MASH   = mean(gse[[sc]][gse$condition_3lvl == "MASH"], na.rm = TRUE),
    kruskal_p   = k$p.value
  )
}
gse_dt <- rbindlist(gse_3lvl)
fwrite(gse_dt, file.path(OUT_DIR, "snRNA_GSE244832_3level.csv"))
cat("\n=== GSE244832 3-level Kruskal-Wallis ===\n")
print(gse_dt[order(kruskal_p)])

# ---------- Bulk path ----------
bulk_path <- file.path(OUT_DIR, "signature_scores_bulk_persample.csv")
if (file.exists(bulk_path)) {
  bulk <- fread(bulk_path)
  cat("\n\n=== BULK ===\n")
  cat("samples:", nrow(bulk), "\n")
  bulk_sig_cols <- grep("^richter|^katsuda|^yin|^consensus|^polyploid_index", names(bulk), value = TRUE)

  # Fibrosis ordinal trend: kendall on fibrosis_stage continuous
  fib_results <- list()
  for (sc in bulk_sig_cols) {
    sub <- bulk[!is.na(fibrosis_stage) & !is.na(get(sc))]
    if (nrow(sub) < 30) next
    k <- cor.test(sub$fibrosis_stage, sub[[sc]], method = "kendall")
    ct <- cor.test(sub$nas_score, sub[[sc]], method = "kendall")
    # Mixed-effects continuous
    m <- tryCatch(lmer(sub[[sc]] ~ fibrosis_stage + (1 | dataset), data = sub),
                  error = function(e) NULL)
    p_fib_mixed <- if (!is.null(m)) summary(m)$coefficients["fibrosis_stage", "Pr(>|t|)"] else NA
    fib_results[[sc]] <- data.table(
      signature = sc, n = nrow(sub),
      kendall_fib_tau = unname(k$estimate), kendall_fib_p = k$p.value,
      kendall_nas_tau = unname(ct$estimate), kendall_nas_p = ct$p.value,
      mixed_fib_p = p_fib_mixed
    )
  }
  fib_dt <- rbindlist(fib_results)
  fib_dt[, padj := p.adjust(kendall_fib_p, method = "BH")]
  fib_dt <- fib_dt[order(kendall_fib_p)]
  fwrite(fib_dt, file.path(OUT_DIR, "bulk_fibrosis_trend.csv"))
  cat("\n=== Bulk fibrosis F0-F4 + NAS trend ===\n")
  print(fib_dt)

  # Per-stage means
  bulk_stage_means <- bulk[!is.na(fibrosis_stage), lapply(.SD, mean, na.rm=TRUE),
                          by = fibrosis_stage, .SDcols = bulk_sig_cols][order(fibrosis_stage)]
  fwrite(bulk_stage_means, file.path(OUT_DIR, "bulk_per_stage_means.csv"))
} else {
  cat("\nBULK SCORING NOT YET COMPLETE; skipping bulk tests\n")
}

cat("\n[", format(Sys.time()), "] DONE Phase 3 tests\n", sep="")

# ---------- Generate figure ----------
PLOIDY_FIG_DIR <- file.path(FIG_DIR, "ploidy_analysis")
dir.create(PLOIDY_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
sn_test[, condition_binary := factor(condition_binary, levels = c("control", "disease"))]
CTRL_GRAY <- "#9E9E9E"; MASL_ORANGE <- "#F57C00"; DISEASE_RED <- "#C62828"
sig_titles <- c(polyploid_index_richter = "Richter polyploid index",
                polyploid_index_katsuda = "Katsuda polyploid index",
                sig_yin_all = "Yin polyploid score")

make_pooled_plot <- function(sc) {
  dat <- sn_test[!is.na(get(sc))]
  pval <- mixed_dt[signature == sc, p_value]
  pval_str <- if (length(pval) == 1 && !is.na(pval)) sprintf("mixed-effects p = %.2f", pval) else ""
  ggplot(dat, aes(x = condition_binary, y = get(sc))) +
    geom_jitter(aes(color = condition_binary), width = 0.15, height = 0, alpha = 0.6, size = 1.8) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linewidth = 0.4) +
    scale_color_manual(values = c(control = CTRL_GRAY, disease = DISEASE_RED), guide = "none") +
    annotate("text", x = 1.5, y = Inf, label = pval_str, hjust = 0.5, vjust = 1.5, size = 3) +
    labs(title = paste0(sig_titles[sc], " (5-cohort pooled)"), x = NULL, y = sc) +
    theme_classic()
}

gse <- sn_test[dataset == "GSE244832" & !is.na(condition_3lvl)]
gse[, condition_3lvl := factor(condition_3lvl, levels = c("NORMAL", "MASL", "MASH"))]
make_3lvl_plot <- function(sc) {
  dat <- gse[!is.na(get(sc))]
  pval <- gse_dt[signature == sc, kruskal_p]
  pval_str <- if (length(pval) == 1 && !is.na(pval)) sprintf("Kruskal-Wallis p = %.3f", pval) else ""
  ggplot(dat, aes(x = condition_3lvl, y = get(sc))) +
    geom_jitter(aes(color = condition_3lvl), width = 0.15, height = 0, alpha = 0.7, size = 2) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linewidth = 0.4) +
    scale_color_manual(values = c(NORMAL = CTRL_GRAY, MASL = MASL_ORANGE, MASH = DISEASE_RED), guide = "none") +
    annotate("text", x = 2, y = Inf, label = pval_str, hjust = 0.5, vjust = 1.5, size = 3) +
    labs(title = paste0(sig_titles[sc], " (GSE244832 only)"), x = NULL, y = sc) +
    theme_classic()
}

sigs <- c("polyploid_index_richter", "polyploid_index_katsuda", "sig_yin_all")
row1 <- lapply(sigs, make_pooled_plot)
row2 <- lapply(sigs, make_3lvl_plot)
fig <- (row1[[1]] | row1[[2]] | row1[[3]]) /
       (row2[[1]] | row2[[2]] | row2[[3]]) +
       plot_annotation(tag_levels = "A")
ggsave(file.path(PLOIDY_FIG_DIR, "02_polyploid_signature_snRNA.pdf"),
       fig, width = 11, height = 8)
cat("Wrote figure:", file.path(PLOIDY_FIG_DIR, "02_polyploid_signature_snRNA.pdf"), "\n")
