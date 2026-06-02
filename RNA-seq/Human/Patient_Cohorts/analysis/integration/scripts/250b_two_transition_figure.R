# 250b_two_transition_figure.R — Multi-step cascade figure for supplementary.
#
# KEY MESSAGE: All four CRN transitions (F0→F1, F1→F2, F2→F3, F3→F4) show substantial
# transcriptomic shifts; F1→F2 is the smallest by both effect-size and DEG count.
# Multi-step cellular cascade rather than single switch — four-panel evidence:
# effect size per transition, unique pathway programs, NMF program activity per
# stage (monotonic, no anchors per M6), and cell-type proportion shifts (neutrophil
# bar shown as regularized log-ratio per M5).
#
# Composes 4 panels following the project FIGURE_GUIDELINES.md (theme_masld + theme_pub,
# magenta/pink/blue palette, colorblind-safe). Run AFTER 250 + 250c complete.
#
# Outputs (vector PDF, no PNG):
#   figures/supplementary/figS_granular_staging/figS_granular_two_transitions.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                           "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(PROJECT_ROOT, "scripts/figures/publication_theme.R"))
source(file.path(PROJECT_ROOT, "scripts/figures/load_figure_data.R"))

DAT <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")

# ---------------------------------------------------------------------------
# Project palette anchors
# ---------------------------------------------------------------------------
# Transitions: F0→F1 (early-blue), F1→F2 (magenta = adaptive immune onset),
# F2→F3 (dark magenta = fibrotic transition), F3→F4 (violet = late).
trans_colors <- c(
  "F0→F1" = "#42A5F5",
  "F1→F2" = masld_colors$mash,
  "F2→F3" = masld_colors$fibrosis,
  "F3→F4" = "#7B1FA2"
)

# NMF program palette — emphasise P1 / P6, mute the rest
prog_colors <- c(
  P1 = masld_colors$mash,         # Pro-inflammatory (peaks F1→F2)
  P2 = "#9E9E9E",                 # Innate-immune
  P3 = "#BDBDBD",                 # Parenchymal
  P4 = "#90CAF9",                 # lncRNA
  P5 = "#42A5F5",                 # Hepatic-metabolic
  P6 = masld_colors$fibrosis      # Stellate-myofibroblast (peaks F2→F3)
)

# Pretty contrast labels
pretty_contrast <- function(x) {
  x <- gsub("F0F1", "F0→F1", x)
  x <- gsub("F1F2", "F1→F2", x)
  x <- gsub("F2F3", "F2→F3", x)
  x <- gsub("F3F4", "F3→F4", x)
  x <- gsub("F0vF1", "F0→F1", x)
  x <- gsub("F1vF2", "F1→F2", x)
  x <- gsub("F2vF3", "F2→F3", x)
  x <- gsub("F3vF4", "F3→F4", x)
  x
}

# ---------------------------------------------------------------------------
# Panel A — Effect size per transition (lollipop chart)
# ---------------------------------------------------------------------------
eff <- read.csv(file.path(DAT, "transition_effect_sizes.csv"))
eff$transition <- pretty_contrast(eff$contrast)
eff$transition <- factor(eff$transition, levels = names(trans_colors))
eff <- eff[!is.na(eff$transition), ]

panel_A <- ggplot(eff,
                  aes(x = transition, y = median_abs_logFC_top1000, color = transition)) +
  geom_segment(aes(xend = transition, y = 0, yend = median_abs_logFC_top1000),
               linewidth = 0.6) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%d DEG", n_padj_05)),
            vjust = -1.0, size = PUB_GEOM_TEXT, color = "gray25") +
  scale_color_manual(values = trans_colors, drop = FALSE, guide = "none") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.20))) +
  labs(title = "Mean shift per CRN transition",
       subtitle = "Median |log2FC| of top 1000 DEGs",
       x = NULL, y = "Median |log2FC|") +
  theme_masld() + theme_pub()

# ---------------------------------------------------------------------------
# Panel B — Unique pathways per transition (horizontal lollipop)
# ---------------------------------------------------------------------------
pwy_path <- file.path(DAT, "two_transition_pathway_decomposition.csv")
panel_B <- ggplot() + theme_void()
if (file.exists(pwy_path)) {
  pwy <- read.csv(pwy_path, stringsAsFactors = FALSE)
  pwy_sig <- pwy[!is.na(pwy$padj) & pwy$padj < 0.05 &
                 pwy$transition %in% c("F1F2", "F2F3"), ]
  f1f2 <- pwy_sig[pwy_sig$transition == "F1F2" & pwy_sig$NES > 0, ]
  f2f3 <- pwy_sig[pwy_sig$transition == "F2F3" & pwy_sig$NES > 0, ]
  f1f2_only <- f1f2[!f1f2$pathway %in% f2f3$pathway, ]
  f2f3_only <- f2f3[!f2f3$pathway %in% f1f2$pathway, ]
  shared <- intersect(f1f2$pathway, f2f3$pathway)

  top_f1f2 <- head(f1f2_only[order(-f1f2_only$NES), ], 6)
  if (nrow(top_f1f2) > 0) top_f1f2$transition <- "F1→F2"
  top_f2f3 <- head(f2f3_only[order(-f2f3_only$NES), ], 4)
  if (nrow(top_f2f3) > 0) top_f2f3$transition <- "F2→F3"
  bar_df <- rbind(top_f1f2, top_f2f3)

  if (nrow(bar_df) > 0) {
    bar_df$label <- gsub("HALLMARK_|REACTOME_", "", bar_df$pathway)
    bar_df$label <- gsub("_", " ", bar_df$label)
    bar_df$label <- tools::toTitleCase(tolower(bar_df$label))
    bar_df$transition <- factor(bar_df$transition, levels = names(trans_colors))
    bar_df <- bar_df[order(bar_df$transition, -bar_df$NES), ]
    bar_df$label <- factor(bar_df$label, levels = rev(bar_df$label))

    panel_B <- ggplot(bar_df,
                      aes(x = label, y = NES, color = transition)) +
      geom_segment(aes(xend = label, y = 0, yend = NES), linewidth = 0.5) +
      geom_point(size = 2.0) +
      coord_flip() +
      scale_color_manual(values = trans_colors, drop = FALSE,
                         name = NULL, breaks = c("F1→F2", "F2→F3")) +
      labs(title = "Transition-unique pathways",
           subtitle = sprintf("padj<0.05; %d shared between F1→F2 and F2→F3", length(shared)),
           x = NULL, y = "Normalized enrichment score") +
      theme_masld() + theme_pub() +
      theme(legend.position = "top",
            legend.justification = "left",
            axis.text.y = element_text(size = PUB_AXIS_TEXT))
  }
}

# ---------------------------------------------------------------------------
# Panel C — NMF k=6 program activity per stage
# ---------------------------------------------------------------------------
nmf_path <- file.path(DAT, "two_transition_nmf_program_by_stage.csv")
panel_C <- ggplot() + theme_void()
if (file.exists(nmf_path)) {
  nmf <- read.csv(nmf_path, stringsAsFactors = FALSE)
  nmf$stage_lab <- factor(paste0("F", nmf$stage), levels = paste0("F", 0:4))
  nmf$program <- factor(nmf$program, levels = paste0("P", 1:6))
  nmf$emphasis <- ifelse(nmf$program %in% c("P1", "P6"), "anchor", "context")

  panel_C <- ggplot(nmf, aes(x = stage_lab, y = mean_score, group = program,
                              color = program, alpha = emphasis,
                              linewidth = emphasis)) +
    geom_line() +
    geom_point(aes(size = emphasis)) +
    annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
             fill = masld_colors$mash, alpha = 0.05) +
    annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
             fill = masld_colors$fibrosis, alpha = 0.05) +
    scale_color_manual(values = prog_colors, name = NULL) +
    scale_alpha_manual(values = c(anchor = 1, context = 0.45), guide = "none") +
    scale_linewidth_manual(values = c(anchor = 0.9, context = 0.4), guide = "none") +
    scale_size_manual(values = c(anchor = 1.8, context = 1.1), guide = "none") +
    labs(title = "NMF k=6 program activity per stage",
         subtitle = "P1 (Pro-inflammatory) and P6 (Stellate-myofibroblast) follow monotonic stage trajectories",
         x = "CRN stage", y = "Mean program score") +
    theme_masld() + theme_pub() +
    theme(legend.position = "right",
          legend.key.size = PUB_LEGEND_KEY)
}

# ---------------------------------------------------------------------------
# Panel D — Cell-type proportion shifts per transition (forest plot)
# ---------------------------------------------------------------------------
ct_path <- file.path(DAT, "two_transition_celltype_by_stage.csv")
panel_D <- ggplot() + theme_void()
if (file.exists(ct_path)) {
  ct <- read.csv(ct_path, stringsAsFactors = FALSE)
  ct_sig <- ct[!is.na(ct$padj) & ct$padj < 0.05 &
               ct$transition %in% c("F1vF2", "F2vF3"), ]
  if (nrow(ct_sig) > 0) {
    ct_sig$transition_l <- pretty_contrast(ct_sig$transition)
    ct_sig$transition_l <- factor(ct_sig$transition_l, levels = names(trans_colors))
    ct_sig$celltype <- gsub("\\.", " ", ct_sig$celltype)
    # M5 verdict (2026-05-09): regularize neutrophil log2FC to avoid divide-by-near-zero
    # storm framing. Original Neutrophils@F2vF3 = +10.39 reflects mean_F2 = 3.83e-8
    # (essentially 0%). Replace with regularized log-ratio at epsilon = 1e-4 = +3.90.
    eps <- 1e-4
    is_neut_f2f3 <- ct_sig$celltype == "Neutrophils" & ct_sig$transition == "F2vF3"
    if (any(is_neut_f2f3)) {
      ct_sig$log2FC_proportion[is_neut_f2f3] <- log2((ct_sig$mean_b[is_neut_f2f3] + eps) /
                                                     (ct_sig$mean_a[is_neut_f2f3] + eps))
    }
    # Top 6 absolute change per transition
    top_per_t <- ct_sig %>%
      group_by(transition_l) %>%
      arrange(desc(abs(log2FC_proportion))) %>%
      slice_head(n = 6) %>%
      ungroup() %>%
      arrange(transition_l, log2FC_proportion) %>%
      mutate(celltype = factor(celltype, levels = unique(celltype)))

    panel_D <- ggplot(top_per_t,
                      aes(x = log2FC_proportion, y = celltype, color = transition_l)) +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "gray60") +
      geom_segment(aes(x = 0, xend = log2FC_proportion,
                       y = celltype, yend = celltype),
                   linewidth = 0.5) +
      geom_point(size = 2.0) +
      facet_grid(transition_l ~ ., scales = "free_y", space = "free_y") +
      scale_color_manual(values = trans_colors, drop = FALSE, guide = "none") +
      labs(title = "Cell-type proportion shifts",
           subtitle = sprintf("MuSiC deconvolution; padj<0.05; Neutrophil@F2vF3 regularized at ε=%.0e (M5)", eps),
           x = "log2FC (cell-type proportion)", y = NULL) +
      theme_masld() + theme_pub() +
      theme(strip.text.y = element_text(face = "bold", angle = 0,
                                         color = "gray20"),
            axis.text.y = element_text(size = PUB_AXIS_TEXT))
  }
}

# ---------------------------------------------------------------------------
# Compose
# ---------------------------------------------------------------------------
combined <- (panel_A | panel_B) / (panel_C | panel_D) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 9)))

pdf_path <- file.path(FIGS_GRANULAR_DIR, "figS_granular_two_transitions.pdf")
ggsave(pdf_path, combined, width = 7.5, height = 7.0, units = "in",
       device = cairo_pdf)
cat("Wrote", pdf_path, "\n")
