#!/usr/bin/env Rscript
# 342h_bulk_hep_adjusted.R - Adjust bulk polyploid signatures for hepatocyte fraction
# (the major confound: as fibrosis progresses, hepatocyte fraction drops).

suppressPackageStartupMessages({
  library(data.table)
  library(lme4); library(lmerTest)
})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
setwd(ROOT)
OUT_DIR <- "Analysis/SingleCell/results_gpu_v2/ploidy"

bulk <- fread(file.path(OUT_DIR, "signature_scores_bulk_persample.csv"))
ct   <- fread("RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv")
cat("bulk n:", nrow(bulk), "  hep-fraction n:", nrow(ct), "\n")

# Merge
ct <- ct[, .(sample_id, hepatocyte_fraction = Hepatocytes)]
bulk2 <- merge(bulk, ct, by = "sample_id")
cat("merged:", nrow(bulk2), "\n")
cat("hepatocyte fraction range:", round(range(bulk2$hepatocyte_fraction, na.rm=TRUE), 3), "\n")
cat("hepatocyte fraction by fibrosis stage:\n")
print(bulk2[!is.na(fibrosis_stage), .(mean_hep = mean(hepatocyte_fraction, na.rm=TRUE),
                                       n = .N), by=fibrosis_stage][order(fibrosis_stage)])

sig_cols <- grep("^richter|^katsuda|^yin|^consensus|^polyploid_index", names(bulk2), value=TRUE)

# For each signature: test fibrosis_stage effect controlling for hepatocyte_fraction
results <- list()
for (sc in sig_cols) {
  sub <- bulk2[!is.na(fibrosis_stage) & !is.na(get(sc)) & !is.na(hepatocyte_fraction)]
  if (nrow(sub) < 30) next
  # Unadjusted
  m_u <- lm(sub[[sc]] ~ sub$fibrosis_stage)
  p_unadj <- summary(m_u)$coefficients["sub$fibrosis_stage", "Pr(>|t|)"]
  est_unadj <- summary(m_u)$coefficients["sub$fibrosis_stage", "Estimate"]
  # Adjusted for hepatocyte fraction
  m_a <- lm(sub[[sc]] ~ sub$fibrosis_stage + sub$hepatocyte_fraction)
  p_adj <- summary(m_a)$coefficients["sub$fibrosis_stage", "Pr(>|t|)"]
  est_adj <- summary(m_a)$coefficients["sub$fibrosis_stage", "Estimate"]
  # Mixed-effects with cohort RE
  m_m <- tryCatch(lmer(sub[[sc]] ~ fibrosis_stage + hepatocyte_fraction + (1|dataset), data=sub),
                  error = function(e) NULL)
  p_mixed <- if (!is.null(m_m)) summary(m_m)$coefficients["fibrosis_stage", "Pr(>|t|)"] else NA
  est_mixed <- if (!is.null(m_m)) summary(m_m)$coefficients["fibrosis_stage", "Estimate"] else NA
  results[[sc]] <- data.table(
    signature = sc, n = nrow(sub),
    est_unadj = est_unadj, p_unadj = p_unadj,
    est_hep_adj = est_adj, p_hep_adj = p_adj,
    est_mixed = est_mixed, p_mixed = p_mixed
  )
}
res_dt <- rbindlist(results)
res_dt[, padj_hep := p.adjust(p_hep_adj, method = "BH")]
res_dt[, padj_mixed := p.adjust(p_mixed, method = "BH")]
res_dt <- res_dt[order(p_hep_adj)]
fwrite(res_dt, file.path(OUT_DIR, "bulk_fibrosis_hep_adjusted.csv"))

cat("\n=== Bulk fibrosis effect: unadjusted vs hepatocyte-fraction-adjusted ===\n")
print(res_dt)

cat("\n=== INTERPRETATION ===\n")
cat("If p_hep_adj is much larger than p_unadj -> the original effect was driven by hepatocyte fraction loss\n")
cat("If p_hep_adj remains significant -> there is genuine within-hepatocyte ploidy signature change\n")
