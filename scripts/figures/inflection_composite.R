#!/usr/bin/env Rscript
# ============================================================================
# inflection_composite.R
# Fig 2 panel b — F1-F3 inflection composite (NEW; closes outline gap)
#
# Two vertically stacked tracks sharing x = 4 CRN transitions:
#   Track 1 (top):    Bayesian breakpoint posterior (P1-P6 NMF programs)
#                      across F0->F1 / F1->F2 / F2->F3 / F3->F4. Dashed line
#                      at 0.90 evidence threshold. All 6 programs stay below
#                      threshold -- no single discrete switch.
#   Track 2 (bottom): NF-kB / SMAD3 TF activity score escalation
#                      (NFKB1, RELA, SMAD3, STAT1). decoupleR norm_wmean
#                      score with p-value annotation.
#
# The "F1-F3 inflection" is visible as the steepest TF-activity rise between
# F0->F1 and F1->F2/F2->F3 transitions, while breakpoint posteriors show NO
# single-stage anchor (max 0.576 for P4@F2-F3 < 0.90 threshold).
#
# Sources:
#   - RNA-seq/results/granular_staging/nmf_changepoint_posterior.csv
#   - RNA-seq/results/stratified_causal/transition_tf_activity.csv
#
# Output: figures/main/fig3_RNAseq/panels/fig3g_inflection_composite.pdf
#   sized 90 x 65 mm.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
OUT_PDF   <- file.path(PANEL_DIR, "fig3g_inflection_composite.pdf")

# Shared x-axis levels
trans_levels <- c("F0→F1", "F1→F2", "F2→F3", "F3→F4")

# ---------------------------------------------------------------------------
# Track 1: NMF changepoint posterior (P1..P6 across 4 transitions)
# ---------------------------------------------------------------------------
cp <- fread(file.path(BASE,
  "RNA-seq/results/granular_staging/nmf_changepoint_posterior.csv"))
cp[, transition_id := factor(sub("_", "→", transition), levels = trans_levels)]
cp[, program := factor(program, levels = paste0("P", 1:6))]
prog_pal <- c(
  P1 = masld_colors$nash,         # Pro-inflammatory
  P2 = "#7B1FA2",                 # Innate-immune
  P3 = masld_colors$control,      # Parenchymal (control gray)
  P4 = "#42A5F5",                 # lncRNA
  P5 = masld_colors$nafl,         # Hepatic-metabolic
  P6 = masld_colors$fibrosis      # Stellate-myofibroblast
)
prog_label <- c(P1 = "Pro-inflam.", P2 = "Innate-imm.", P3 = "Parenchymal",
                P4 = "lncRNA", P5 = "Hep-metab.", P6 = "Stellate")

p_top <- ggplot(cp, aes(x = transition_id, y = posterior_changepoint,
                         fill = program)) +
  geom_hline(yintercept = 0.90, color = "grey30",
             linetype = "dashed", linewidth = 0.35) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7,
           color = "grey25", linewidth = 0.12) +
  annotate("text", x = 1, y = 0.93, hjust = 0,
           label = "0.90 evidence threshold (no single switch passes)",
           size = 1.7, color = "grey25", fontface = "italic") +
  scale_fill_manual(values = prog_pal, labels = prog_label,
                    name = NULL) +
  scale_y_continuous(name = "Bayesian breakpoint posterior",
                     limits = c(0, 1), expand = c(0, 0),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = NULL,
       title = "F1–F3 inflection composite",
       subtitle = "No single-transition breakpoint (max P=0.58 < 0.90); NF-κB/RELA peak in F1→F2 → F2→F3 inflection window") +
  theme_masld(base_size = 7) +
  theme(
    plot.title    = element_text(size = 7.5, face = "bold"),
    plot.subtitle = element_text(size = 5.7, color = "grey35"),
    legend.position = "right",
    legend.text   = element_text(size = 5),
    legend.key.size = unit(0.18, "cm"),
    axis.text.x   = element_blank(),
    axis.title.x  = element_blank(),
    plot.margin   = margin(2, 2, 0, 2)
  )

# ---------------------------------------------------------------------------
# Track 2: NF-kB / SMAD3 TF activity escalation
# ---------------------------------------------------------------------------
tf <- fread(file.path(BASE,
  "RNA-seq/results/stratified_causal/transition_tf_activity.csv"))
tf_of_interest <- c("NFKB1", "RELA", "SMAD3", "STAT1")
tf_use <- tf[tf %in% tf_of_interest]
setnames(tf_use, "tf", "tf_name")
tf_use[, transition_id := factor(sub("_to_", "→", transition),
                                  levels = trans_levels)]
tf_use[, tf_name := factor(tf_name, levels = tf_of_interest)]
tf_use[, sig := p_value < 0.05]

tf_pal <- c(
  NFKB1 = masld_colors$nash,        # canonical NF-kB
  RELA  = masld_colors$fibrosis,    # NF-kB family
  SMAD3 = "#1565C0",                # TGF-beta canonical
  STAT1 = "#7B1FA2"                 # interferon signaling
)

p_bot <- ggplot(tf_use, aes(x = transition_id, y = score,
                              color = tf_name, group = tf_name)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.2) +
  geom_line(linewidth = 0.6) +
  geom_point(aes(shape = sig), size = 1.5, stroke = 0.35,
             fill = "white") +
  scale_color_manual(values = tf_pal, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "p<0.05", `FALSE` = "n.s."),
                     name = NULL) +
  scale_y_continuous(name = "TF activity (decoupleR z)",
                     breaks = c(0, 3, 6, 9, 12)) +
  labs(x = NULL) +
  theme_masld(base_size = 7) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.18, "cm"),
    axis.text.x = element_text(size = 6.5, face = "bold"),
    plot.margin = margin(0, 2, 2, 2)
  ) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

panel <- p_top / p_bot + plot_layout(heights = c(1, 1.05))

ggsave(OUT_PDF, panel,
       width  = 90 / 25.4,
       height = 65 / 25.4,
       units  = "in",
       device = cairo_pdf)

fwrite(cp, file.path(DATA_DIR, "inflection_changepoint.csv"))
fwrite(tf_use, file.path(DATA_DIR, "inflection_tf_activity.csv"))
cat(sprintf("[saved] %s\n", OUT_PDF))
