#!/usr/bin/env Rscript
##############################################################################
# figS_robust_batch.R  (new 2026-05-08)
# Pillar D condensed: variance partition (left) + dream vs dream+SVA (right).
# KEY MESSAGE: dream's cohort RE collapses cohort variance from ~65% (raw) to
#              ~0% (residual); adding SVA on top barely changes any logFC.
##############################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
AUDIT   <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
INT     <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(FIG_SUPP, "robustness")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(OUT_DIR, "figS_robust_batch.pdf")

theme_robust <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
    axis.line   = element_line(colour = "black", linewidth = 0.3),
    axis.ticks  = element_line(colour = "black", linewidth = 0.3),
    legend.position = "top",
    legend.title  = element_blank(),
    legend.margin = margin(b = -3),
    plot.title    = element_text(face = "bold", size = 11,
                                 margin = margin(b = 2)),
    plot.subtitle = element_text(size = 9, colour = "grey30",
                                 margin = margin(b = 6)),
    plot.tag      = element_text(face = "bold", size = 12),
    plot.margin   = margin(8, 10, 8, 10)
  )

# --- Panel A: variance partition raw vs dream residual ---
raw  <- fread(file.path(INT,   "variance_partition.csv"))
resd <- fread(file.path(AUDIT, "pillar_D_residual_varpart.csv"))
melt_long <- function(dt, lbl) {
  cols <- intersect(c("dataset", "condition", "group_binary",
                      "sex", "inferred_sex", "Residuals"), names(dt))
  m <- melt(dt[, c("gene", cols), with = FALSE],
            id.vars = "gene", variable.name = "component",
            value.name = "var_frac")
  m[, partition := lbl]
  m[]
}
both <- rbind(melt_long(raw, "raw expression"),
              melt_long(resd, "dream residuals"), fill = TRUE)
both[, component_lbl := fcase(
  component %in% c("dataset"), "Cohort",
  component %in% c("condition", "group_binary"), "Disease",
  component %in% c("sex", "inferred_sex"), "Sex",
  component == "Residuals", "Residual",
  default = as.character(component))]
both[, component_lbl := factor(component_lbl,
                               levels = c("Cohort", "Disease", "Sex", "Residual"))]
both[, partition := factor(partition,
                           levels = c("raw expression", "dream residuals"))]
both <- both[is.finite(var_frac)]

# Median dots overlaid on the violins (variancePartition convention)
agg <- both[, .(median_frac = median(var_frac, na.rm = TRUE)),
            by = .(partition, component_lbl)]

PAL_PART <- c("raw expression" = "#999999", "dream residuals" = "#0072B2")

p_a <- ggplot(both, aes(x = component_lbl, y = var_frac,
                        fill = partition, colour = partition)) +
  geom_violin(position = position_dodge(width = 0.78),
              scale = "width", trim = TRUE,
              alpha = 0.55, linewidth = 0.25, width = 0.78) +
  geom_point(data = agg,
             aes(y = median_frac, group = partition),
             position = position_dodge(width = 0.78),
             shape = 21, size = 2.6, fill = "white",
             stroke = 0.6, show.legend = FALSE) +
  geom_text(data = agg,
            aes(y = median_frac,
                label = sprintf("%.0f%%", 100 * median_frac),
                group = partition),
            position = position_dodge(width = 0.78),
            vjust = -1.0, size = 2.8, colour = "grey15",
            show.legend = FALSE) +
  scale_fill_manual(values = PAL_PART) +
  scale_colour_manual(values = PAL_PART) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.10)),
                     limits = c(0, 1)) +
  labs(tag = "A",
       title = "Variance partition",
       subtitle = "Per-gene variance fraction (violins); white dots = medians",
       x = NULL, y = "Variance fraction") +
  theme_robust

# --- Panel B: dream vs dream+SVA logFC ---
sva <- fread(file.path(AUDIT, "pillar_D_sva_concordance.csv"))
sva_summary <- fread(file.path(AUDIT, "pillar_D_sva_summary.csv"))

set.seed(42)
sva[, panel_class := factor(fifelse(sig_primary | sig_sva, "DEG", "Non-DEG"),
                            levels = c("Non-DEG", "DEG"))]
sva_plot <- rbind(sva[panel_class == "DEG"],
                  sva[panel_class == "Non-DEG"][sample(.N, min(2000, .N))])

p_b <- ggplot(sva_plot, aes(x = logFC_primary, y = logFC_sva,
                            colour = panel_class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey50", linewidth = 0.3) +
  geom_point(alpha = 0.5, size = 0.5) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
           label = sprintf("rho = %.3f\nJaccard = %.3f\nn_SV = %d",
                           sva_summary$rho_logFC[1],
                           sva_summary$jaccard[1],
                           sva_summary$n_sv[1]),
           size = 3.2, colour = "grey20") +
  scale_colour_manual(values = c("Non-DEG" = "grey70", "DEG" = "#0072B2")) +
  labs(tag = "B",
       title = "dream vs dream + SVA",
       subtitle = "Surrogate variables on top of dream barely shift any logFC",
       x = "log2FC (dream)", y = "log2FC (dream + SVA)") +
  theme_robust +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

composite <- p_a + p_b + plot_layout(ncol = 2)
ggsave(OUT_PDF, composite, width = 9, height = 4, device = cairo_pdf)
cat("Saved:", OUT_PDF, "\n")
