#!/usr/bin/env Rscript
##############################################################################
# figS_convergence_evidence.R
# Supplementary figure for Script 46d (unsupervised convergence evidence score).
# 4 panels:
#   (a) Tier breakdown bar plot
#   (b) Concordance state distribution stacked by tier
#   (c) Protective-LOF top-30 gene table
#   (d) Held-out panel recovery curves (46d vs 46b)
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

ME <- file.path(BASE, "RNA-seq/results/multi_evidence")
FIGDIR <- FIGS_CONV_EVID_DIR

ev    <- fread(file.path(ME, "convergence_evidence.csv"))
plof  <- fread(file.path(ME, "convergence_evidence_genetic_down_coherent.csv"))
bench <- if (file.exists(file.path(ME, "convergence_evidence_benchmark.csv")))
           fread(file.path(ME, "convergence_evidence_benchmark.csv")) else NULL
vs46b <- if (file.exists(file.path(ME, "convergence_evidence_vs_46b.csv")))
           fread(file.path(ME, "convergence_evidence_vs_46b.csv")) else NULL

# ============================================================================
# Panel (a) — Tier breakdown
# ============================================================================
tier_tbl <- ev[excluded_from_ranking == FALSE, .N, by = tier]
tier_tbl[, tier := factor(tier,
  levels = c("1_Genetic_validated","2_Convergent","3_Suggestive","4_Weak"),
  labels = c("1: Genetic\nvalidated","2: Convergent\nmulti-modal",
             "3: Suggestive","4: Weak"))]
tier_colors <- c("1: Genetic\nvalidated" = "#880E4F",
                 "2: Convergent\nmulti-modal" = "#E91E63",
                 "3: Suggestive" = "#F48FB1",
                 "4: Weak" = "#BDBDBD")

message(sprintf("[caption] Convergence evidence tier distribution. Total ranked: %s  (excluded confounders: %d)",
                comma(sum(tier_tbl$N)),
                sum(ev$excluded_from_ranking)))

pa <- ggplot(tier_tbl, aes(x = tier, y = N, fill = tier)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = comma(N)), vjust = -0.3, size = GEOM_TEXT_6PT, color = "gray15") +
  scale_fill_manual(values = tier_colors, guide = "none") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Gene count") +
  theme_masld() +
  theme(axis.text.x = element_text(size = 6))

# ============================================================================
# Panel (b) — Concordance state by tier (stacked %)
# ============================================================================
state_tbl <- ev[excluded_from_ranking == FALSE,
  .N, by = .(tier, concordance_state)]
state_tbl[, concordance_state := factor(concordance_state,
  levels = c("Concordant-up","Concordant-down","Protective-LOF",
             "Conflicted","Insufficient-evidence"))]
state_tbl[, tier := factor(tier,
  levels = c("1_Genetic_validated","2_Convergent","3_Suggestive","4_Weak"),
  labels = c("Tier 1","Tier 2","Tier 3","Tier 4"))]

state_colors <- c("Concordant-up"        = "#C2185B",
                  "Concordant-down"      = "#1565C0",
                  "Protective-LOF"       = "#880E4F",
                  "Conflicted"           = "#F48FB1",
                  "Insufficient-evidence"= "#BDBDBD")

pb <- ggplot(state_tbl, aes(x = tier, y = N, fill = concordance_state)) +
  geom_col(position = "fill", width = 0.7) +
  scale_fill_manual(values = state_colors, name = NULL) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(x = NULL, y = "Share of tier") +
  theme_masld() +
  theme(legend.position = "bottom", legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(nrow = 2))

# ============================================================================
# Panel (c) — Protective-LOF top 30 (table)
# ============================================================================
plof_top <- plof[order(-convergence_score)][1:min(30, .N)]
plof_top[, label := sprintf(
  "%-9s  PP4=%.2f  n_mod=%d  stage=%s",
  human_symbol,
  ifelse(is.na(coloc_best_pp4_S2a),0,coloc_best_pp4_S2a),
  n_modalities_active,
  ifelse(is.na(dominant_stage_S1),"NA",substr(dominant_stage_S1,1,14)))]
plof_top[, y := seq(.N, 1)]

message(sprintf("[caption] Protective-LOF candidates (top %d of %d). Genetic + expression-down signature (inhibitor-target pattern).",
                nrow(plof_top), nrow(plof)))

pc <- ggplot(plof_top, aes(x = 0, y = y)) +
  geom_text(aes(label = label), hjust = 0, size = GEOM_TEXT_6PT,
            family = "mono", color = "gray15") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, nrow(plof_top) + 1), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_masld() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

# ============================================================================
# Panel (d) — Held-out panel recovery curves
# ============================================================================
if (!is.null(bench) && nrow(bench) > 0) {
  message("[caption] Held-out panel recovery. AUROC for distinguishing panel members vs rest (higher = better).")
  bench_dt <- melt(bench,
    id.vars = "panel",
    measure.vars = c("auroc_46d","auroc_46b"),
    variable.name = "model", value.name = "auroc")
  bench_dt[, model := factor(gsub("auroc_","",model),
    levels = c("46d","46b"),
    labels = c("46d (unsupervised)","46b (archetype)"))]

  pd <- ggplot(bench_dt, aes(x = panel, y = auroc, fill = model)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               linewidth = 0.3, color = "gray55") +
    geom_hline(yintercept = 0.7, linetype = "dotted",
               linewidth = 0.3, color = "gray45") +
    scale_fill_manual(values = c("46d (unsupervised)" = "#880E4F",
                                  "46b (archetype)"     = "#42A5F5"),
                      name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "AUROC") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          legend.position = "bottom",
          legend.key.size = unit(0.25, "cm"))
} else {
  pd <- placeholder("Panel benchmark CSV missing — curate data/published_gene_panels/")
}

# ============================================================================
# Assemble (4 panels — the 46b-vs-46d scatter was retired with 46b, 2026-04-23)
# ============================================================================
composite <- (pa | pb) / (pc | pd) +
  plot_layout(heights = c(1, 1.3)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 6, face = "plain"))

save_fig(composite, file.path(FIGDIR, "figS_convergence_evidence.pdf"),
         width = fig_full_width, height = 9)
message("Saved: ", file.path(FIGDIR, "figS_convergence_evidence.pdf"))
