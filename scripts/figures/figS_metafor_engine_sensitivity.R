#!/usr/bin/env Rscript
# figS_metafor_engine_sensitivity.R
# Compact sensitivity panel: does the metafor two-stage meta-analysis depend on
# the per-study (Stage-1) DE engine? Plots meta-log2FC from limma-voom->metafor
# vs DESeq2->metafor across the common genes. Near-identity => we use ONE engine
# (limma-voom) as "metafor" in the dream/DESeq2/metafor comparison.
# Output: figures/supplementary/figS_methods_validation/multimethod_validation/panels/
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
SENS <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/multimethod_validation/sensitivity")
OUT  <- file.path(BASE, "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

m <- fread(file.path(SENS, "deseq2_vs_limma_metafor_comparison.csv"))
m <- m[is.finite(lv_logFC) & is.finite(ds_logFC)]
rho <- cor(m$lv_logFC, m$ds_logFC, method = "spearman")
pr  <- cor(m$lv_logFC, m$ds_logFC, method = "pearson")
# direction concordance among genes Tier-1 in either engine
lv <- m[!is.na(lv_padj) & lv_padj < 0.05 & abs(lv_logFC) > 0.5, gene]
ds <- m[!is.na(ds_padj) & ds_padj < 0.05 & abs(ds_logFC) > 0.5, gene]
shared <- intersect(lv, ds)
dir_conc <- mean(sign(m[gene %in% shared, lv_logFC]) == sign(m[gene %in% shared, ds_logFC]))
lim <- quantile(abs(c(m$lv_logFC, m$ds_logFC)), 0.999)

p <- ggplot(m, aes(lv_logFC, ds_logFC)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "#D6604D", linewidth = 0.4) +
  geom_hex(bins = 70) +
  scale_fill_viridis_c(trans = "log10", name = "genes", option = "mako", direction = -1) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "meta-log2FC  (limma-voom -> metafor)",
       y = "meta-log2FC  (DESeq2 -> metafor)") +
  annotate("text", x = -lim * 0.95, y = lim * 0.9, hjust = 0,
           label = sprintf("rho = %.3f", rho), size = 6/ggplot2::.pt, fontface = "plain", colour = "black") +
  theme_masld(base_size = 6) +
  theme(legend.position = "right", legend.key.width = unit(0.25, "cm"))

message(sprintf("[caption] metafor is robust to the per-study DE engine: %s common genes | Spearman rho = %.3f, Pearson = %.3f | direction concordance %.0f%% (shared DEGs); per-cohort Stage-1 rho = 0.94-0.97 -> use limma-voom as the single metafor arm",
                format(nrow(m), big.mark = ","), rho, pr, 100 * dir_conc))

ggsave(file.path(OUT, "panelF_metafor_engine_sensitivity.pdf"), p,
       width = 4.8, height = 4.4, useDingbats = FALSE)
cat(sprintf("Wrote panelF_metafor_engine_sensitivity.pdf | rho=%.3f pearson=%.3f dir=%.3f n=%d\n",
            rho, pr, dir_conc, nrow(m)))
