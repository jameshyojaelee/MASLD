#!/usr/bin/env Rscript
##############################################################################
# figS_sex_stage_attribution.R
#
# Stage-attribution cross-tab figure for the 4-contrast sex × disease
# extension (Disease-vs-Ctrl + MASH-vs-MASL + MASL-vs-Ctrl + MASH-vs-Ctrl).
#
# Panels:
#   a  Tag count bar (where in the trajectory does sex bias act?)
#   b  Cross-tab: stage_attribution_tag × Disease-vs-Ctrl assigned_class
#      (does the existing D-vs-Ctrl 4-class scheme see what stage-resolved
#      tags see? i.e., do "Suggestive F_only" genes map to pan-stage or to
#      a specific transition?)
#   c  β_int(D) vs β_int(P) scatter, points colored by tag — visual of the
#      D vs progression decomposition
#   d  Top 6 genes per non-trivial tag (table panel)
#
# Inputs: results/integration/sex_v3/stage_attribution_table.csv (Script 20).
# Output: figures/supplementary/figS_sex_dimorphism/figS_sex_stage_attribution.pdf
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

FIGDIR <- FIGS_SEX_DIR
OUT    <- file.path(FIGDIR, "figS_sex_stage_attribution.pdf")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Input — Script 20 master table (4,724-row Tier-1 universe joined across
# 4 contrasts, with stage_attribution_tag column)
# ---------------------------------------------------------------------------
TAB <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration",
  "sex_v3/stage_attribution_table.csv")
stopifnot(file.exists(TAB))
dt <- fread(TAB)

# Defensive coerce — gene + tag must be character, betas numeric
dt[, gene := as.character(gene)]
dt[, stage_attribution_tag := as.character(stage_attribution_tag)]
for (col in grep("^beta_int_", names(dt), value = TRUE)) dt[[col]] <- as.numeric(dt[[col]])
for (col in grep("^padj_int_",  names(dt), value = TRUE)) dt[[col]] <- as.numeric(dt[[col]])

# ---------------------------------------------------------------------------
# Tag palette — ordered by biological narrative
#   pan-stage          (sex bias persists across trajectory)
#   D-Suggestive-only  (pooled-D bias only — onset, not progression)
#   progression-only   (sex bias appears at MASH-vs-MASL severity)
#   early-onset-direction (descriptive only; MASL-vs-Ctrl direction)
#   late-onset-direction  (descriptive only; MASH-vs-Ctrl direction)
#   discordant         (sign flips across contrasts — likely artifact)
#   no-signal          (Tier-1 disease DEG with no detectable sex modulation)
# ---------------------------------------------------------------------------
tag_levels <- c("pan-stage",
                "both-stage",
                "D-Suggestive-only",
                "progression-only",
                "early-onset-direction",
                "late-onset-direction",
                "discordant",
                "no-signal")
tag_palette <- c(
  "pan-stage"             = "#AD1457",   # rose magenta — stable bias
  "both-stage"            = "#E91E63",   # softer magenta — D+P Suggestive, E/L not concordant
  "D-Suggestive-only"     = "#F4A674",   # peach — pooled-D only
  "progression-only"      = "#C9265E",   # deep magenta — severity
  "early-onset-direction" = "#cf9ced",   # lavender — descriptive
  "late-onset-direction"  = "#6A1B9A",   # purple — descriptive
  "discordant"            = "#1A237E",   # navy — artifact
  "no-signal"             = "#9E9E9E"    # gray — null
)
dt[, tag := factor(stage_attribution_tag, levels = tag_levels)]

# ---------------------------------------------------------------------------
# Panel A: tag distribution
# ---------------------------------------------------------------------------
message("Panel A: tag bar ...")
tag_counts <- dt[, .N, by = tag][order(match(tag, tag_levels))]
tag_counts[, pct := 100 * N / sum(N)]

pA <- ggplot(tag_counts, aes(x = tag, y = N, fill = tag)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", N, pct)),
            vjust = -0.3, size = PUB_GEOM_TEXT) +
  scale_fill_manual(values = tag_palette, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Tier-1 DEGs", title = "Stage-attribution tag distribution",
       subtitle = "All 4,724 Tier-1 disease DEGs classified by where sex bias acts") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))

# ---------------------------------------------------------------------------
# Panel B: cross-tab tag × Disease-vs-Ctrl assigned_class
# Pulls the v3 4-class table from sex_v3/sex_deg_classification_v3.csv so
# we can show how the stage-resolved tags map onto the existing 4-class
# scheme (F_only / M_only / divergent / concordant / uncertain).
# ---------------------------------------------------------------------------
message("Panel B: tag x D-vs-Ctrl class cross-tab ...")
V3_CSV <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration",
  "sex_v3/sex_deg_classification_v3.csv")
v3 <- if (file.exists(V3_CSV)) fread(V3_CSV) else NULL

if (!is.null(v3) && "assigned_class_gated" %in% names(v3)) {
  v3_small <- v3[, .(gene, d_class = assigned_class_gated)]
  ct <- merge(dt[, .(gene, tag)], v3_small, by = "gene", all.x = TRUE)
  ct[is.na(d_class), d_class := "not_tested"]
  d_class_levels <- c("F_only", "M_only", "divergent",
                       "concordant", "uncertain", "not_tested")
  ct[, d_class := factor(d_class, levels = d_class_levels)]
  agg <- ct[, .N, by = .(tag, d_class)]
  agg[, pct_within_tag := 100 * N / sum(N), by = tag]

  pB <- ggplot(agg, aes(x = tag, y = d_class, fill = pct_within_tag)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(N > 0, as.character(N), "")),
              size = PUB_GEOM_TEXT) +
    scale_fill_gradient(low = "#F4F4F4", high = "#AD1457",
                         name = "% within tag",
                         labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = "Disease-vs-Ctrl class (v3)",
         title = "Stage-attribution tag x Disease-vs-Ctrl 4-class scheme",
         subtitle = "Cell counts = genes; color = column %") +
    theme_masld() + theme_pub() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
} else {
  message("  v3 D-vs-Ctrl class CSV missing; Panel B falls back to placeholder")
  pB <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5,
             label = "v3 D-vs-Ctrl class CSV missing",
             size = 2.4) +
    labs(title = "Panel B unavailable")
}

# ---------------------------------------------------------------------------
# Panel C: beta_int(D) vs beta_int(P) scatter, colored by tag
# ---------------------------------------------------------------------------
message("Panel C: beta_int D vs P scatter ...")
need_cols <- c("beta_int_D", "beta_int_P")
if (all(need_cols %in% names(dt))) {
  scatter_dt <- dt[!is.na(beta_int_D) & !is.na(beta_int_P)]
  pC <- ggplot(scatter_dt, aes(x = beta_int_D, y = beta_int_P, color = tag)) +
    geom_hline(yintercept = 0, color = "gray70", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "gray70", linewidth = 0.2) +
    geom_abline(slope = 1, intercept = 0, color = "gray80",
                linetype = "dashed", linewidth = 0.2) +
    rasterize_layer(geom_point(size = 0.25, alpha = 0.6)) +
    scale_color_manual(values = tag_palette,
                       guide = guide_legend(override.aes = list(size = 1.4))) +
    labs(x = expression(beta[int] * " (Disease vs Ctrl)"),
         y = expression(beta[int] * " (MASH vs MASL)"),
         color = "Tag",
         title = "Sex x disease beta_int: pooled-D vs severity",
         subtitle = "Diagonal = perfectly stable bias; off-axis = stage-specific") +
    theme_masld() + theme_pub()
} else {
  pC <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5,
             label = "beta_int_D / beta_int_P column missing", size = 2.4) +
    labs(title = "Panel C unavailable")
}

# ---------------------------------------------------------------------------
# Panel D: top 6 genes per non-trivial tag (table panel via geom_text)
# ---------------------------------------------------------------------------
message("Panel D: top-genes table ...")
non_trivial <- setdiff(tag_levels, c("no-signal", "discordant"))
# Rank within each tag by max |beta_int| across the 4 contrasts
beta_cols <- grep("^beta_int_", names(dt), value = TRUE)
dt[, max_abs_beta := pmax(abs(beta_int_D), abs(beta_int_E),
                           abs(beta_int_L), abs(beta_int_P), na.rm = TRUE)]
top_genes <- dt[tag %in% non_trivial][
  order(-max_abs_beta), .SD[1:6], by = tag]
top_genes <- top_genes[!is.na(gene)]
top_genes[, tag := factor(tag, levels = non_trivial)]
top_genes[, rank_in_tag := seq_len(.N), by = tag]

pD <- ggplot(top_genes,
              aes(x = tag, y = -rank_in_tag, label = gene, color = tag)) +
  geom_text(size = PUB_GEOM_TEXT + 0.2, fontface = "italic", hjust = 0.5) +
  scale_color_manual(values = tag_palette, guide = "none") +
  scale_y_continuous(breaks = NULL) +
  labs(x = NULL, y = NULL,
       title = "Top 6 genes per non-trivial tag",
       subtitle = "Ranked by max |beta_int| across the 4 contrasts") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        axis.text.y = element_blank())

# ---------------------------------------------------------------------------
# Compose + save
# ---------------------------------------------------------------------------
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "a",
                  title = "Stage-attribution of sex bias across the MASLD trajectory",
                  subtitle = "4 sex x disease contrasts: D-vs-Ctrl (formal, n=847) + MASH-vs-MASL (formal, n=687) + MASL-vs-Ctrl + MASH-vs-Ctrl (descriptive)",
                  theme = theme(plot.title    = element_text(size = PUB_TITLE,    face = "bold"),
                                plot.subtitle = element_text(size = PUB_SUBTITLE, color = "gray30")))
ggsave(OUT, fig,
       width = fig_full_width, height = 1.1 * fig_full_width, units = "in")
message("Wrote: ", OUT)
