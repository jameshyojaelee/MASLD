#!/usr/bin/env Rscript

# Release-driven panels for the complementary evidence-class manuscript.
# Outputs PDF only under the canonical Figure 1 and Figure 4 directories.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- Sys.getenv("MANUSCRIPT_RELEASE_ID", "2026-07-15-r2")
RELEASE_DIR <- file.path(BASE, "RNA-seq/results/manuscript_release", RELEASE_ID)
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

required <- file.path(RELEASE_DIR, c(
  "evidence_class_table.tsv",
  "orthogonality_audit.tsv",
  "evidence_class_validation_summary.tsv",
  "evidence_class_validation_adjusted.tsv"
))
if (any(!file.exists(required))) {
  stop("Missing manuscript-release output(s): ", paste(basename(required[!file.exists(required)]), collapse = ", "))
}

classes <- fread(required[1])
orth <- fread(required[2])
val <- fread(required[3])
adjusted <- fread(required[4])

class_colors <- c(
  neither = "#D7D9DA",
  genetic_only = "#1565C0",
  disease_state_only = "#C9265E",
  convergent = "#00695C"
)
class_labels <- c(
  neither = "Neither",
  genetic_only = "Genetic only",
  disease_state_only = "Disease-state only",
  convergent = "Convergent"
)

# Figure 1B: two conditional views of the same overlap. This avoids a Venn area
# illusion while showing the biologically relevant denominator in each map.
primary <- orth[
  trait_scope == "all" & genetic_definition == "susie" & pp4_threshold == 0.5
][1]
overlap <- primary$n_overlap
view <- rbind(
  data.table(
    map = "Primary SuSiE genes",
    segment = c("convergent", "genetic_only"),
    n = c(overlap, primary$n_genetic_joint - overlap)
  ),
  data.table(
    map = "Canonical DEGs\n(interval-null FDR gate)",
    segment = c("convergent", "disease_state_only"),
    n = c(overlap, primary$n_treat_deg_joint - overlap)
  )
)
view[, fraction := n / sum(n), by = map]
view[, label := sprintf("%s\n%s", comma(n), percent(fraction, accuracy = 0.1))]
view[, map := factor(map, levels = c("Canonical DEGs\n(interval-null FDR gate)", "Primary SuSiE genes"))]

p1b <- ggplot(view, aes(x = fraction, y = map, fill = segment)) +
  geom_col(width = 0.56, color = "white", linewidth = 0.35) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            size = 6 / .pt, color = "white", lineheight = 0.9) +
  scale_fill_manual(values = class_colors, labels = class_labels, name = NULL) +
  scale_x_continuous(labels = percent, expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Fraction of jointly tested map", y = NULL) +
  theme_masld(base_size = 6) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(0.38, "cm")
  )

ggsave(file.path(FIG1_DIR, "fig1b_complementary_maps.pdf"), p1b,
       width = 4.8, height = 2.25, device = cairo_pdf, bg = "white")

# Figure 1C: composition of the 473 primary SuSiE genes by phenotype source.
trait <- classes[genetic_confidence == "susie", .N, by = genetic_trait_scope]
trait <- trait[genetic_trait_scope %in% c("direct_disease", "enzyme", "both")]
trait[, label := fcase(
  genetic_trait_scope == "direct_disease", "Direct disease/PDFF only",
  genetic_trait_scope == "enzyme", "Liver enzyme only",
  default = "Both"
)]
trait[, label := factor(label, levels = c("Liver enzyme only", "Both", "Direct disease/PDFF only"))]

p1c <- ggplot(trait, aes(x = N, y = label, fill = label)) +
  geom_col(width = 0.58, show.legend = FALSE) +
  geom_text(aes(label = comma(N)), hjust = -0.15, size = 6 / .pt) +
  scale_fill_manual(values = c(
    "Liver enzyme only" = "#F4A674",
    "Both" = "#00695C",
    "Direct disease/PDFF only" = "#1565C0"
  )) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Primary SuSiE genes", y = NULL) +
  theme_masld(base_size = 6) +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(FIG1_DIR, "fig1c_genetic_trait_scope.pdf"), p1c,
       width = 4.0, height = 2.25, device = cairo_pdf, bg = "white")

# Figure 4A: assay-specific positive rates. Denominators are printed explicitly.
val[, class_label := factor(class_labels[primary_evidence_class],
                            levels = unname(class_labels[c("neither", "genetic_only",
                              "disease_state_only", "convergent")]))]
val[, endpoint_label := factor(endpoint,
  levels = c("single_cell", "spatial_svg", "proteomics"),
  labels = c("Single-cell", "Spatial SVG", "Proteomics"))]
val[, rate_label := sprintf("%s/%s", comma(n_positive), comma(n_tested))]

p4a <- ggplot(val, aes(x = positive_rate, y = endpoint_label,
                       color = class_label, size = n_tested, group = class_label)) +
  geom_point(alpha = 0.9, position = position_dodge(width = 0.56)) +
  geom_text(aes(label = rate_label), position = position_dodge(width = 0.56),
            vjust = -1.0, size = 5 / .pt, color = "grey25", show.legend = FALSE) +
  scale_color_manual(values = setNames(class_colors, class_labels), name = NULL) +
  scale_size_continuous(range = c(2.0, 5.2), guide = "none") +
  scale_x_continuous(labels = percent, limits = c(0, 0.42),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Positive rate within assay-specific background", y = NULL) +
  theme_masld(base_size = 6) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

ggsave(file.path(FIG4_DIR, "fig4a_evidence_class_positive_rates.pdf"), p4a,
       width = 5.2, height = 3.0, device = cairo_pdf, bg = "white")

# Figure 4B: adjusted effects versus neither. Confidence intervals expose the
# uncertainty of the small convergent class rather than encoding rank order.
adjusted[, class_key := sub("^primary_evidence_class", "", term)]
adjusted[, class_label := factor(class_labels[class_key],
  levels = c("Genetic only", "Disease-state only", "Convergent"))]
adjusted[, endpoint_label := factor(endpoint,
  levels = c("single_cell", "spatial_svg", "proteomics"),
  labels = c("Single-cell", "Spatial SVG", "Proteomics"))]

p4b <- ggplot(adjusted, aes(x = odds_ratio, y = endpoint_label,
                            color = class_label)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55", linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                position = position_dodge(width = 0.42), width = 0, linewidth = 0.45) +
  geom_point(position = position_dodge(width = 0.42), size = 1.8) +
  scale_color_manual(values = setNames(class_colors[c("genetic_only", "disease_state_only", "convergent")],
                                       class_labels[c("genetic_only", "disease_state_only", "convergent")]),
                     name = NULL) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50)) +
  labs(x = "Adjusted odds ratio versus neither (log scale)", y = NULL) +
  theme_masld(base_size = 6) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

ggsave(file.path(FIG4_DIR, "fig4b_evidence_class_adjusted_or.pdf"), p4b,
       width = 5.2, height = 3.0, device = cairo_pdf, bg = "white")

# Figure 4C: convergent-versus-single-map comparisons -- the honest "no universal
# convergence advantage" result. The convergent class is not detectably superior to
# disease-state-only genes in any modality (all BH q>0.05); it exceeds genetic-only
# only weakly. Denominators (small convergent n) are printed so the null reads as an
# underpowered honest-negative, not a proven equivalence.
pair <- fread(file.path(RELEASE_DIR, "evidence_class_validation_pairwise.tsv"))
pair[, endpoint_label := factor(endpoint,
  levels = c("single_cell", "spatial_svg", "proteomics"),
  labels = c("Single-cell", "Spatial SVG", "Proteomics"))]
pair[, comp_class := fifelse(comparison == "convergent_vs_disease_state_only",
                             "disease_state_only", "genetic_only")]
pair[, point_label := sprintf("q=%.2f (n=%d vs %d)", q_value, n_convergent, n_comparator)]

p4c <- ggplot(pair, aes(x = odds_ratio, y = endpoint_label, color = comp_class)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55", linewidth = 0.35) +
  geom_point(position = position_dodge(width = 0.55), size = 1.9) +
  geom_text(aes(label = point_label), position = position_dodge(width = 0.55),
            hjust = -0.12, size = 5 / .pt, color = "grey25", show.legend = FALSE) +
  scale_color_manual(
    values = c(genetic_only = class_colors[["genetic_only"]],
               disease_state_only = class_colors[["disease_state_only"]]),
    labels = c(genetic_only = "vs genetic-only",
               disease_state_only = "vs disease-state-only"),
    name = NULL) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10), limits = c(0.5, 30),
                expand = expansion(mult = c(0.02, 0))) +
  labs(x = "Convergent-class odds ratio vs single-map class (log scale)", y = NULL) +
  theme_masld(base_size = 6) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

ggsave(file.path(FIG4_DIR, "fig4c_convergence_comparison.pdf"), p4c,
       width = 5.2, height = 3.0, device = cairo_pdf, bg = "white")

cat(sprintf("[release %s] Figure 1/4 evidence-class panels written\n", RELEASE_ID))
cat(sprintf("  SuSiE jointly tested: %d; overlap: %d; non-DE: %.1f%%\n",
            primary$n_genetic_joint, primary$n_overlap, primary$pct_genetic_not_deg))
