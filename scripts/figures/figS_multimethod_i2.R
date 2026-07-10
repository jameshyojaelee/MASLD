#!/usr/bin/env Rscript
# figS_multimethod_i2.R
# Panel E — cross-cohort heterogeneity (I^2) distribution for the metafor arm.
# Single compact histogram with Higgins reference bands; frames the high
# between-cohort heterogeneity at K=5 that motivates the dream/DESeq2/metafor
# comparison (and why metafor's RE meta-analysis is underpowered here).
# Per the investigation: report the I^2 DISTRIBUTION, not per-gene values.
# Output: figures/supplementary/figS_methods_validation/multimethod_validation/panels/
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
INT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT <- file.path(BASE, "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

m   <- fread(file.path(INT, "meta_analysis_results.csv"))   # canonical metafor (REML)
i2  <- m$meta_I2[is.finite(m$meta_I2)]
med <- median(i2); p50 <- 100 * mean(i2 > 50); p75 <- 100 * mean(i2 > 75)

p <- ggplot(data.table(I2 = i2), aes(I2)) +
  geom_histogram(bins = 50, fill = "#5E4FA2", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = c(25, 50, 75), linetype = "dashed",
             colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = med, colour = "#D6604D", linewidth = 0.5) +
  annotate("text", x = c(12.5, 37.5, 62.5, 87.5), y = Inf,
           label = c("low", "moderate", "substantial", "considerable"),
           vjust = 1.6, size = GEOM_TEXT_6PT, colour = "grey45") +
  annotate("text", x = med, y = Inf, label = sprintf("median %.0f%%", med),
           vjust = 3.2, hjust = -0.08, size = GEOM_TEXT_6PT, colour = "#D6604D", fontface = "plain") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression(I^2 ~ "(%)"), y = "Genes") +
  theme_masld(base_size = 6)

message(sprintf("[caption] Cross-cohort heterogeneity at K=5 | %s genes | median I2 %.0f%% | %.0f%% >50%% | %.0f%% >75%% (genuine, not estimand)",
                format(length(i2), big.mark = ","), med, p50, p75))
ggsave(file.path(OUT, "i2_distribution.pdf"), p,
       width = 4.6, height = 3.0, useDingbats = FALSE)
cat(sprintf("Wrote i2_distribution.pdf | n=%d median=%.1f%% >50=%.1f%% >75=%.1f%%\n",
            length(i2), med, p50, p75))
