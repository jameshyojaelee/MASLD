#!/usr/bin/env Rscript
# figS_multimethod_bootstrap.R
# Panel C (bootstrap stability) for the multi-method validation, 3-way:
#   panelC_selfreq_concordance.pdf  - pairwise selection-frequency concordance (hexbin facets + rho)
#   panelC2_stability_counts.pdf    - per-method # stably-selected genes (freq>0.5, >0.9)
# 500 cohort-stratified bootstrap resamples; "selection frequency" = fraction of
# resamples a gene is called (padj<0.1). Shows dream~DESeq2 highly concordant,
# metafor conservative (fewer stable genes) but equally stable effect sizes.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
B   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/multimethod_validation/bootstrap")
OUT <- file.path(BASE, "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

sf  <- fread(file.path(B, "selection_frequency.csv"))
mcc <- fread(file.path(B, "method_concordance.csv"))
mcol <- c(dream = "#1b9e77", deseq2 = "#66a61e", metafor = "#1f78b4")
pretty <- c(dream = "dream", deseq2 = "DESeq2", metafor = "metafor")

# ---- Panel C1: pairwise selection-frequency hexbin facets --------------------
pairs <- list(c("dream","deseq2"), c("dream","metafor"), c("deseq2","metafor"))
dl <- rbindlist(lapply(pairs, function(p) {
  a <- paste0(p[1],"_freq"); b <- paste0(p[2],"_freq")
  d <- sf[!is.na(get(a)) & !is.na(get(b)), .(x = get(a), y = get(b))]
  rho <- mcc[(method_a==p[1]&method_b==p[2])|(method_a==p[2]&method_b==p[1]), spearman_freq][1]
  d[, pair := sprintf("%s vs %s  (rho=%.2f)", pretty[p[1]], pretty[p[2]], rho)]
  d[, xlab := pretty[p[1]]]; d[, ylab := pretty[p[2]]]; d
}))
p1 <- ggplot(dl, aes(x, y)) +
  geom_abline(slope = 1, linetype = "dashed", colour = "#D6604D", linewidth = 0.3) +
  geom_hex(bins = 40) +
  scale_fill_viridis_c(trans = "log10", name = "genes", option = "mako", direction = -1) +
  facet_wrap(~ pair, nrow = 1) +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "selection frequency (method A)", y = "selection frequency (method B)",
       title = "Bootstrap selection-frequency concordance (500 resamples)",
       subtitle = "dream~DESeq2 tight (rho=0.88); metafor concordant in rank but conservative (rho 0.60-0.67)") +
  theme_masld(base_size = 7) +
  theme(plot.subtitle = element_text(size = 6, colour = "grey40"),
        legend.position = "right", legend.key.width = unit(0.22,"cm"),
        strip.text = element_text(size = 6.5))
# panelC_selfreq_concordance.pdf retired per user request 2026-06-05 — no longer generated.
# (p1 retained above for provenance; intentionally not written.)
# ggsave(file.path(OUT, "panelC_selfreq_concordance.pdf"), p1,
#        width = 8.0, height = 3.2, useDingbats = FALSE)
# cat("Wrote panelC_selfreq_concordance.pdf\n")

# ---- Panel C2: per-method stably-selected gene counts ------------------------
ss <- fread(file.path(B, "stability_summary.csv"))
sl <- melt(ss[, .(method, `freq>0.5` = n_freq_gt_0.5, `freq>0.9` = n_freq_gt_0.9)],
           id.vars = "method", variable.name = "threshold", value.name = "n_genes")
sl[, method := factor(method, levels = c("dream","deseq2","metafor"))]
p2 <- ggplot(sl, aes(threshold, n_genes, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  geom_text(aes(label = scales::comma(n_genes)), position = position_dodge(0.8),
            vjust = -0.3, size = 2.3) +
  scale_fill_manual(values = mcol, labels = pretty, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = scales::comma) +
  labs(x = "bootstrap selection-frequency threshold", y = "stably-selected genes",
       title = "Bootstrap stability: stably-selected genes per method",
       subtitle = "metafor selects far fewer (conservative) but its effect sizes are as stable (median logFC CV 0.51 vs dream 0.51, DESeq2 0.61)") +
  theme_masld(base_size = 7) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 5.6, colour = "grey40"))
ggsave(file.path(OUT, "panelC2_stability_counts.pdf"), p2,
       width = 5.2, height = 3.4, useDingbats = FALSE)
cat("Wrote panelC2_stability_counts.pdf\n")
print(ss[, .(method, n_freq_gt_0.5, n_freq_gt_0.9, median_logFC_cv = round(median_logFC_cv,3))])
