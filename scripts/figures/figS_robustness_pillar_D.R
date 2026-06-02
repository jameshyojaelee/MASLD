#!/usr/bin/env Rscript
# figS_robustness_pillar_D.R
# Pillar D — Batch-RE diagnostics (residual variance partition + SVA).
# Panels:
#   D1: Variance fraction comparison: raw v$E vs dream-residual (cohort, disease, sex)
#   D2: dream vs dream+SVA logFC scatter (Jaccard + ρ overlay)
#   D3: SVA concordance: bar of canonical-DEG concordance status

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

AUDIT <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
INT   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(BASE, "figures/supplementary/robustness")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- D1: Raw vs residual variance partition ---
raw  <- fread(file.path(INT, "variance_partition.csv"))
resd <- fread(file.path(AUDIT, "pillar_D_residual_varpart.csv"))

melt_long <- function(dt, lbl) {
  cols <- intersect(c("dataset", "condition", "group_binary", "sex", "inferred_sex", "Residuals"), names(dt))
  m <- melt(dt[, c("gene", cols), with = FALSE],
            id.vars = "gene", variable.name = "component", value.name = "var_frac")
  m[, partition := lbl]; m
}
raw_l  <- melt_long(raw,  "raw expression")
resd_l <- melt_long(resd, "dream residuals")
both <- rbind(raw_l, resd_l, fill = TRUE)

# Map component names to a common axis
both[, component_lbl := fcase(
  component %in% c("dataset"), "Cohort",
  component %in% c("condition", "group_binary"), "Disease",
  component %in% c("sex", "inferred_sex"), "Sex",
  component == "Residuals", "Residual",
  default = as.character(component))]
agg <- both[, .(median_frac = median(var_frac, na.rm = TRUE),
                q90_frac    = quantile(var_frac, 0.9, na.rm = TRUE)),
            by = .(partition, component_lbl)]

p1 <- ggplot(agg, aes(x = component_lbl, y = median_frac, fill = partition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = median_frac, ymax = q90_frac, group = partition),
                position = position_dodge(width = 0.8), width = 0.2, colour = "grey40") +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(title = "D1 — Variance fractions: raw vs dream residual",
       subtitle = "Bars = median across genes; error bars to 90th percentile.",
       x = "Variance component", y = "Fraction of total variance", fill = "") +
  theme_pub() + theme(axis.text.x = element_text(angle = 0))

# --- D2: dream vs dream+SVA logFC scatter ---
sva <- fread(file.path(AUDIT, "pillar_D_sva_concordance.csv"))
sva_summary <- fread(file.path(AUDIT, "pillar_D_sva_summary.csv"))

p2 <- ggplot(sva, aes(x = logFC_primary, y = logFC_sva,
                       colour = sva_concordant)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf("Jaccard = %.3f\nρ = %.3f\nn_SV = %d",
                           sva_summary$jaccard[1], sva_summary$rho_logFC[1],
                           sva_summary$n_sv[1])) +
  labs(title = "D2 — dream vs dream+SVA logFC",
       x = "logFC (primary dream)", y = "logFC (dream + SVA)",
       colour = "Concordant DEG") +
  theme_pub()

# --- D3: Concordance bar ---
sva[, status := fcase(
  sig_primary & sig_sva & sva_concordant, "Both & same direction",
  sig_primary & sig_sva & !sva_concordant, "Both, opposite",
  sig_primary & !sig_sva, "Primary only",
  !sig_primary & sig_sva, "SVA only",
  default = "Neither")]
status_summary <- sva[, .N, by = status][order(-N)]
p3 <- ggplot(status_summary, aes(x = reorder(status, N), y = N, fill = status)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(title = "D3 — DEG concordance: dream vs dream+SVA",
       x = "", y = "Genes", fill = "") +
  theme_pub() + theme(legend.position = "none")

combined <- (p1 / (p2 | p3)) + plot_layout(heights = c(1, 1.1))
out_pdf <- file.path(OUT_DIR, "figS_robustness_pillar_D.pdf")
ggsave(out_pdf, combined, width = 12, height = 9, device = cairo_pdf)
cat("Saved:", out_pdf, "\n")
