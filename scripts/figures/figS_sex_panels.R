#!/usr/bin/env Rscript
# figS_sex_panels.R — Sex-modulation supplementary panels.
#
# Generates three panels in figures/supplementary/figS_sex_dimorphism/:
#   1. figS_sex_scatter.pdf            (β_F vs β_M, named candidates highlighted)
#   2. figS_sex_class_distribution.pdf (bar of cross-pillar class counts)
#   3. figS_sex_beta_int_violin.pdf    (β_int distribution per class)
#
# All three share a single 6-class color palette that's consistent across the
# whole sex supp directory (Female-biased = rose / Male-biased = navy /
# Divergent = yellow / Suggestive-other = pale green / Power-limited = amber
# / Uncertain = gray). Data source: load_sex_classification() reads
# sex_deg_classification_v3.csv which carries the cross-pillar consensus
# columns via the 06i shim.
#
# Extracted from the archived fig2_disease_progression.R sex SUPP block
# (2026-05-17 redesign) so the panels can be iterated without touching the
# main-figure script.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS_SEX_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Single canonical 6-class palette shared across all sex panels --------
class_levels <- c("Female-biased", "Male-biased", "Divergent",
                  "Suggestive (other)", "Power-limited", "Uncertain")
class_colors <- c(
  "Female-biased"      = masld_colors$female,   # rose #AD1457
  "Male-biased"        = masld_colors$male,     # navy #1A237E
  "Divergent"          = "#FFC107",             # yellow
  "Suggestive (other)" = "#A5D6A7",             # pale green
  "Power-limited"      = "#FFB300",             # amber
  "Uncertain"          = "#CFD8DC"              # cool gray
)

# Map raw class_v6_combined → display class.
combined_to_display <- function(x) {
  fcase(
    x == "Suggestive_F_biased",  "Female-biased",
    x == "Suggestive_M_biased",  "Male-biased",
    x == "Suggestive_Divergent", "Divergent",
    x == "Suggestive_NA",        "Suggestive (other)",
    x == "Power_limited",        "Power-limited",
    default = "Uncertain")
}

# ---- Load data ------------------------------------------------------------
sex <- as.data.table(load_sex_classification())
if (!"class_v6_consensus" %in% names(sex) || all(is.na(sex$class_v6_consensus))) {
  stop("Cross-pillar consensus columns missing — has 06i shim run?")
}
if (!"symbol" %in% names(sex) && "gene_symbol" %in% names(sex))
  sex[, symbol := gene_symbol]

# β_int comes from the dream M2 interaction coefficient
beta_int_col <- intersect(c("beta_interaction", "beta_int", "beta_int_v6"),
                          names(sex))[1]
if (is.na(beta_int_col))
  stop("No β_int column found; expected one of beta_interaction / beta_int")
setnames(sex, beta_int_col, "beta_int_value")

sex[, class_display := combined_to_display(class_v6_combined)]
sex[, class_display := factor(class_display, levels = class_levels)]

# Named candidate pool (for the scatter labels): named v6 Suggestive
# Female-biased + Male-biased + top 6 Divergent by |β_F − β_M|.
named <- rbindlist(list(
  sex[class_display == "Female-biased" & !grepl("^ENS", symbol)],
  sex[class_display == "Male-biased"   & !grepl("^ENS", symbol)],
  head(sex[class_display == "Divergent" & !grepl("^ENS", symbol)][
        order(-abs(beta_F - beta_M))], 6)
), fill = TRUE)
named <- named[!is.na(beta_F) & !is.na(beta_M)]

# ============================================================================
# Panel 1 — β_F vs β_M scatter
# ============================================================================
sex[, plot_z := fcase(
  class_display == "Uncertain",          0L,
  class_display == "Suggestive (other)", 1L,
  class_display == "Power-limited",      2L,
  class_display == "Divergent",          3L,
  class_display %in% c("Female-biased", "Male-biased"), 4L)]
setorder(sex, plot_z)

sc <- sex[!is.na(beta_F) & !is.na(beta_M)]
lim_v6 <- max(quantile(abs(c(sc$beta_F, sc$beta_M)), 0.99, na.rm = TRUE), 0.5)

display_counts <- sc[, .N, by = class_display]
count_lookup <- setNames(display_counts$N, as.character(display_counts$class_display))
get_n <- function(k) format(if (!is.null(count_lookup[[k]])) count_lookup[[k]] else 0L,
                             big.mark = ",")

# Legend: only the four substantive sex-modulation categories
legend_classes <- c("Female-biased", "Male-biased", "Divergent", "Uncertain")
legend_classes <- intersect(legend_classes,
                             as.character(unique(sc$class_display)))
legend_labels <- vapply(legend_classes,
                         function(k) sprintf("%s (n=%s)", k, get_n(k)),
                         character(1))

n_strong <- sum(sex$class_v6_consensus == "Strong",   na.rm = TRUE)
n_mod    <- sum(sex$class_v6_consensus == "Moderate", na.rm = TRUE)
null_note <- if (n_strong == 0 && n_mod == 0)
  "Strong = 0 and Moderate = 0 at the strict cross-pillar gate." else ""

p_scatter <- ggplot(sc, aes(x = beta_M, y = beta_F, color = class_display)) +
  geom_abline(slope =  1, intercept = 0, linetype = "dashed",
              linewidth = 0.3, color = "gray50") +
  geom_abline(slope = -1, intercept = 0, linetype = "dotted",
              linewidth = 0.25, color = "gray70") +
  rasterize_layer(geom_point(size = 0.32, alpha = 0.55, shape = 16)) +
  # Named candidates: larger filled circle + white border so the eye picks
  # them out without inventing a second color scale.
  geom_point(data = named,
             aes(x = beta_M, y = beta_F, color = class_display),
             size = 1.8, shape = 21, stroke = 0.5, fill = "white",
             show.legend = FALSE) +
  geom_point(data = named,
             aes(x = beta_M, y = beta_F, color = class_display),
             size = 1.6, shape = 16, show.legend = FALSE) +
  geom_label_repel(data = named,
                   aes(x = beta_M, y = beta_F, label = symbol),
                   inherit.aes = FALSE, color = "gray20",
                   size = 1.5, max.overlaps = 50,
                   label.padding = 0.10, box.padding = 0.4,
                   segment.size = 0.15, fill = "white",
                   alpha = 0.92, show.legend = FALSE) +
  scale_color_manual(values = class_colors,
                      breaks = legend_classes,
                      labels = legend_labels,
                      name = "Cross-pillar class",
                      drop = TRUE) +
  coord_equal(xlim = c(-lim_v6, lim_v6), ylim = c(-lim_v6, lim_v6)) +
  labs(x = expression("Male disease effect "*beta[M]),
       y = expression("Female disease effect "*beta[F]),
       title = expression("Cross-pillar sex-modulation: " * beta[F] * " vs " * beta[M])) +
  theme_masld() + theme_pub() +
  guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 1), ncol = 1)) +
  theme(legend.text  = element_text(size = 5.5),
        legend.title = element_text(size = 6.5),
        legend.key.size = unit(6, "pt"),
        legend.spacing.y = unit(2, "pt"))

save_fig(p_scatter, file.path(OUT_DIR, "figS_sex_scatter.pdf"),
         width = fig_half_width, height = 3.5, dpi = 300)
cat("Wrote: figS_sex_scatter.pdf\n")

# ============================================================================
# Panel 2 — class count bar (drops the Uncertain backdrop)
# ============================================================================
focal <- c("Female-biased", "Male-biased", "Divergent",
           "Suggestive (other)", "Power-limited")
cnt <- sex[, .N, by = class_display]
uncertain_n <- cnt[class_display == "Uncertain", N]
if (length(uncertain_n) == 0) uncertain_n <- 0L
cnt <- cnt[class_display %in% focal]
cnt[, class_display := factor(class_display, levels = rev(focal))]

p_class <- ggplot(cnt, aes(y = class_display, x = N, fill = class_display)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(aes(label = format(N, big.mark = ",")),
            hjust = -0.2, size = 2.4, fontface = "plain", color = "gray20") +
  scale_fill_manual(values = class_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(y = NULL, x = "Genes (cross-pillar consensus)",
       title = "Sex-modulation class counts",
       subtitle = sprintf("0 Strong / 0 Moderate at the strict gate.\n%s Uncertain genes not shown.",
                          format(uncertain_n, big.mark = ","))) +
  theme_masld() + theme_pub() +
  theme(plot.subtitle = element_text(size = 6.5, color = "gray30",
                                       lineheight = 1.05))

save_fig(p_class, file.path(OUT_DIR, "figS_sex_class_distribution.pdf"),
         width = fig_half_width, height = 3.5, dpi = 300)
cat("Wrote: figS_sex_class_distribution.pdf\n")

# ============================================================================
# Panel 3 — β_int violin per class (with handling for n=1 classes)
# ============================================================================
vio <- sex[!is.na(beta_int_value)]
focal_violin <- c("Female-biased", "Male-biased", "Divergent",
                  "Suggestive (other)", "Power-limited", "Uncertain")
vio[, class_display := factor(class_display, levels = rev(focal_violin))]

# Split: classes with ≥ 3 points get a violin; classes with < 3 get a point
class_n <- vio[, .N, by = class_display]
small_classes <- class_n[N < 3, as.character(class_display)]
big_classes   <- class_n[N >= 3, as.character(class_display)]

vio_big   <- vio[as.character(class_display) %in% big_classes]
vio_small <- vio[as.character(class_display) %in% small_classes]

# n=1 (Power-limited) label tag — show count
vio[, lab_text := sprintf("n=%d", .N), by = class_display]
class_n_labels <- unique(vio[, .(class_display, lab_text)])

p_violin <- ggplot(vio, aes(y = class_display, x = beta_int_value,
                             fill = class_display)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted",
             color = "gray50", linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = "solid",
             color = "gray70", linewidth = 0.2) +
  # Violins only for classes with ≥3 points
  geom_violin(data = vio_big, scale = "width", linewidth = 0.25,
              color = "gray30", alpha = 0.85, trim = TRUE) +
  # Big circles for small-N classes (so n=1 Power-limited is visible)
  geom_point(data = vio_small, size = 1.8, shape = 21, stroke = 0.5,
             color = "gray20") +
  # Jittered dots for the named-candidate sex-biased classes
  geom_point(data = vio[as.character(class_display) %in%
                          c("Female-biased", "Male-biased", "Divergent")],
             position = position_jitter(height = 0.20, seed = 42),
             size = 0.45, alpha = 0.7, color = "gray20", shape = 16) +
  # n-count tag at the right edge of each row
  geom_text(data = class_n_labels,
            aes(y = class_display, x = 2.0, label = lab_text),
            color = "gray35", size = 1.8, hjust = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = class_colors, guide = "none") +
  coord_cartesian(xlim = c(-2.0, 2.1)) +
  labs(y = NULL,
       x = expression("Sex × disease interaction "*beta[int]),
       title = expression(beta[int] * " distribution by cross-pillar class"),
       subtitle = expression("Dotted lines mark the Tier-C magnitude gate |"*beta[int]*"| = 0.5.")) +
  theme_masld() + theme_pub() +
  theme(plot.subtitle = element_text(size = 6, color = "gray30",
                                       lineheight = 1.05))

save_fig(p_violin, file.path(OUT_DIR, "figS_sex_beta_int_violin.pdf"),
         width = fig_half_width, height = 3.0, dpi = 300)
cat("Wrote: figS_sex_beta_int_violin.pdf\n")

cat("\nAll 3 sex supp panels written to:", OUT_DIR, "\n")
