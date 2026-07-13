#!/usr/bin/env Rscript
# KEY MESSAGE: The regulatory drivers of fibrosis progression activate in a
# SEQUENCE — inflammatory NF-kB (NFKB1/RELA) peaks earliest at F1->F2, the
# fibrogenic driver SMAD3 peaks one transition later at F2->F3 — and NO single
# transition behaves as a discrete switch (every changepoint posterior stays
# below 0.90). This is the regulatory-driver layer beneath the monotonic
# phenotype cascade in fig3e; it is NOT itself monotonic (the "F1-F3 inflection"
# framing is retired — the TF activity peaks-and-declines, a relay, not a switch).
# ============================================================================
# inflection_composite.R  — Fig 3G (inflammatory->fibrogenic driver relay)
#
# Two vertically stacked tracks sharing x = 4 CRN transitions:
#   Track 1 (top):    Changepoint posterior (P1-P6 NMF programs) across
#                      F0->F1 / F1->F2 / F2->F3 / F3->F4. Dashed line at the
#                      0.90 evidence threshold. NO program clears it (max 0.58)
#                      -> no single discrete switch stage. (BIC-discrete
#                      fallback; mcp unavailable — stated in caption.)
#   Track 2 (bottom): TF activity relay (decoupleR norm_wmean z). Inflammatory
#                      NF-kB drivers NFKB1/RELA peak at F1->F2; fibrogenic SMAD3
#                      peaks at F2->F3; STAT1 (interferon) is non-monotonic.
#
# Sources:
#   - RNA-seq/results/granular_staging/nmf_changepoint_posterior.csv
#   - RNA-seq/results/stratified_causal/transition_tf_activity.csv
#
# Output: figures/main/fig3_RNAseq/panels/fig3h_inflection_composite.pdf
#   (filename kept for back-compat; the panel no longer says "inflection").
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
OUT_PDF   <- file.path(PANEL_DIR, "fig3h_inflection_composite.pdf")

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
  annotate("text", x = 1, y = 0.90, hjust = 0, vjust = -0.5,
           label = "0.90", size = 1.9, color = "grey30") +
  scale_fill_manual(values = prog_pal, labels = prog_label,
                    name = NULL) +
  scale_y_continuous(name = "Changepoint\nposterior",
                     limits = c(0, 1), expand = c(0, 0),
                     breaks = c(0, 0.5, 1.0)) +
  labs(x = NULL) +
  theme_masld(base_size = 7) +
  theme(
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

# Two "waves": inflammatory NF-kB (NFKB1/RELA, magenta) crests at F1->F2,
# fibrogenic SMAD3 (teal) crests at F2->F3; STAT1 (interferon) is de-emphasised
# grey/thin because it is non-monotonic and mostly n.s. mid-trajectory.
tf_pal <- c(
  NFKB1 = masld_colors$nash,        # #C9265E inflammatory wave
  RELA  = masld_colors$fibrosis,    # #A01753 inflammatory wave (NF-kB family)
  SMAD3 = "#00695C",                # teal — fibrogenic (TGF-beta) wave
  STAT1 = "#9E9E9E"                 # grey — interferon (secondary, non-monotonic)
)
tf_lw <- c(NFKB1 = 0.75, RELA = 0.75, SMAD3 = 0.75, STAT1 = 0.45)

# Label each wave at its crest to make the inflammatory->fibrogenic hand-off
# legible (gene labels = data labels, italic + BLACK per house style).
peak_lab <- tf_use[tf_name %in% c("NFKB1", "SMAD3"),
                   .SD[which.max(score)], by = tf_name]
peak_lab[, vj := ifelse(tf_name == "NFKB1", -0.9, 1.9)]

p_bot <- ggplot(tf_use, aes(x = transition_id, y = score,
                              color = tf_name, group = tf_name)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.2) +
  geom_line(aes(linewidth = tf_name)) +
  geom_point(aes(shape = sig), size = 1.5, stroke = 0.35,
             fill = "white") +
  geom_text(data = peak_lab, aes(label = tf_name, vjust = vj),
            size = 2.1, fontface = "italic", color = "black",
            show.legend = FALSE) +
  scale_color_manual(values = tf_pal, name = NULL,
                     labels = c(expression(italic("NFKB1")),
                                expression(italic("RELA")),
                                expression(italic("SMAD3")),
                                expression(italic("STAT1")))) +
  scale_linewidth_manual(values = tf_lw, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "p<0.05", `FALSE` = "n.s."),
                     name = NULL) +
  scale_y_continuous(name = "TF activity\n(decoupleR z)",
                     breaks = c(0, 3, 6, 9),
                     expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = "Fibrosis transition") +
  theme_masld(base_size = 7) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.18, "cm"),
    axis.text.x = element_text(size = 6.5, face = "bold"),
    axis.title.x = element_text(size = 7),
    plot.margin = margin(0, 2, 2, 2)
  ) +
  guides(color = guide_legend(order = 1,
                              override.aes = list(linewidth = 0.75)),
         shape = guide_legend(order = 2))

panel <- p_top / p_bot +
  plot_layout(heights = c(0.9, 1.15)) +
  plot_annotation(
    title = "Regulatory drivers activate sequentially, not at a discrete switch",
    theme = theme(plot.title = element_text(size = 8, face = "bold")))

ggsave(OUT_PDF, panel,
       width  = 92 / 25.4,
       height = 70 / 25.4,
       units  = "in",
       device = cairo_pdf)

fwrite(cp, file.path(DATA_DIR, "inflection_changepoint.csv"))
fwrite(tf_use, file.path(DATA_DIR, "inflection_tf_activity.csv"))

# Descriptive text -> figure CAPTION (house style: no explanatory text in-plot).
message(sprintf(paste0(
  "[Fig 3G caption] Inflammatory->fibrogenic regulatory-driver relay across the four ",
  "CRN fibrosis transitions (F0->F1 ... F3->F4).\n",
  "(top) Changepoint posterior for each NMF program P1-P6; no program's posterior ",
  "exceeds the 0.90 evidence threshold (max %.2f), so no single transition behaves ",
  "as a discrete switch (BIC-discrete fallback, mcp unavailable).\n",
  "(bottom) TF activity (decoupleR norm_wmean z) per transition: inflammatory NF-kB ",
  "drivers NFKB1/RELA peak earliest at F1->F2 (%.1f / %.1f), the fibrogenic driver ",
  "SMAD3 peaks one transition later at F2->F3 (%.1f), and STAT1 (interferon) is ",
  "non-monotonic. Filled points p<0.05, open n.s. The drivers hand off in sequence ",
  "rather than switching at one stage."),
  max(cp$posterior_changepoint, na.rm = TRUE),
  tf_use[tf_name == "NFKB1" & transition_id == "F1→F2", score],
  tf_use[tf_name == "RELA"  & transition_id == "F1→F2", score],
  tf_use[tf_name == "SMAD3" & transition_id == "F2→F3", score]))
cat(sprintf("[saved] %s\n", OUT_PDF))
