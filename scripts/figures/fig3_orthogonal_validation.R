##############################################################################
# Fig 3 — Orthogonal validation of regulatory variants (Liang refinement)
#
# Two compact subpanels combined via patchwork (~6in x 3in total):
#   A. Roadmap E066 chromHMM overlap of disease-regulon-disrupting variants
#      Horizontal stacked bar over 7 chromHMM states (TssA, Enh, EnhG,
#      TxFlnk, TxWk, Het, Quies). Magenta for active, gray for inactive.
#      Annotation: 30/33 (90.9%) in active states; 28.8x enriched (p=4.6e-219).
#
#   B. Currin 138-donor caQTL agreement per disease-TF (top 15 by % agreement)
#      Horizontal lollipop, 4-way validated TFs (THRB/HNF4A/RORA/MLXIPL/MAX)
#      highlighted in magenta + gold outline; others in blue.
#
# Inputs:
#   GWAS/finemapping/results/gwas_atac/roadmap_disease_regulon_replication.csv
#   GWAS/finemapping/results/gwas_atac/roadmap_A8_summary.csv
#   GWAS/finemapping/results/gwas_atac/currin_caqtl_per_tf.csv
#
# Output:
#   figures/main/fig2_genetics/panels/fig3_orthogonal_validation.pdf
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(grid)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_PDF <- file.path(FIGS05_DIR, "figS05_scatac_orthogonal_validation.pdf")

GWAS_ATAC_DIR <- file.path(BASE, "GWAS/finemapping/results/gwas_atac")

# ---------------------------------------------------------------------------
# Color helpers
# ---------------------------------------------------------------------------
MAGENTA_DEEP <- "#C9265E"
MAGENTA_MID  <- "#E04C82"
MAGENTA_LITE <- "#F08AAC"
GRAY_DEEP    <- "#616161"
GRAY_MID     <- "#9E9E9E"
GRAY_LITE    <- "#CFCFCF"
BLUE         <- "#1565C0"
GOLD         <- "#FFB300"

# ---------------------------------------------------------------------------
# Subpanel A — Roadmap chromHMM state overlap
# ---------------------------------------------------------------------------
# Roadmap 15-state code -> short label mapping (E066 model used in 56h)
state_label <- c(
  "E1"  = "TssA",  "E2"  = "TssAFlnk", "E3"  = "TxFlnk",
  "E4"  = "Tx",    "E5"  = "TxWk",     "E6"  = "EnhG",
  "E7"  = "Enh",   "E8"  = "ZNF/Rpts", "E9"  = "Het",
  "E10" = "TssBiv","E11" = "BivFlnk",  "E12" = "EnhBiv",
  "E13" = "ReprPC","E14" = "ReprPCWk", "E15" = "Quies"
)

dis_reg <- fread(file.path(GWAS_ATAC_DIR, "roadmap_disease_regulon_replication.csv"))
a8      <- fread(file.path(GWAS_ATAC_DIR, "roadmap_A8_summary.csv"))

# Tabulate variants by chromHMM state
state_counts <- dis_reg[, .N, by = chromhmm_state]
state_counts[, label := state_label[chromhmm_state]]

# Requested display ordering (top->bottom = active first, then inactive)
display_states <- c("TssA", "TssAFlnk", "Enh", "EnhG", "TxFlnk", "Tx",
                    "TxWk", "Het", "Quies")
# Map any missing states to 0 so the stacked bar has consistent ordering
plot_df <- data.table(label = display_states)
plot_df <- merge(plot_df, state_counts[, .(label, N)], by = "label", all.x = TRUE)
plot_df[is.na(N), N := 0]
plot_df[, label := factor(label, levels = rev(display_states))]
plot_df[, pct := 100 * N / sum(N)]

# Active states (E066 model): TssA, TssAFlnk, EnhG, Enh (per 56h_roadmap_replication.R)
active_states <- c("TssA", "TssAFlnk", "EnhG", "Enh")
plot_df[, is_active := label %in% active_states]

# Per-state magenta/gray shades (gradient by activity strength)
state_colors <- c(
  "TssA"     = "#C9265E",  # Liang deep magenta — active promoter
  "TssAFlnk" = "#D14A7E",
  "Enh"      = "#E04C82",  # active enhancer
  "EnhG"     = "#EE7BAA",
  "TxFlnk"   = MAGENTA_LITE,
  "Tx"       = "#BFBFBF",
  "TxWk"     = GRAY_LITE,
  "Het"      = GRAY_MID,
  "Quies"    = GRAY_DEEP
)

n_active <- sum(plot_df[is_active == TRUE]$N)
n_total  <- sum(plot_df$N)
pct_active <- round(100 * n_active / n_total, 1)

# Pull headline OR/p from A8 summary
or_hep <- a8[metric == "fisher_or_hep_vs_nopeak_chromhmm_active"]$value
pv_hep <- a8[metric == "fisher_or_hep_vs_nopeak_pvalue"]$value

annot_A <- sprintf("30/33 (90.9%%) in active states\n28.8x enriched (p = 4.6e-219)")

# Force the stacked-bar to a single horizontal row (y dummy)
plot_df[, y := "Disease-regulon variants"]
plot_df[, label := factor(label, levels = display_states)]

p_A <- ggplot(plot_df, aes(x = N, y = y, fill = label)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = state_colors, name = "chromHMM state",
                    breaks = display_states,
                    guide = guide_legend(ncol = 1, reverse = FALSE,
                                         keyheight = unit(0.30, "cm"),
                                         keywidth  = unit(0.30, "cm"))) +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.02)),
                     breaks = c(0, 10, 20, 30),
                     labels = function(x) paste0(x)) +
  labs(tag = "a", x = "Variants in chromHMM state (n)", y = NULL) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(
    legend.position = "right",
    legend.justification = "top",
    legend.title = element_text(size = PUB_LEGEND_TIT, face = "bold"),
    legend.text  = element_text(size = PUB_LEGEND),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    plot.tag = element_text(size = PUB_TITLE + 1, face = "bold"),
    plot.tag.position = c(0.02, 0.97),
    panel.grid = element_blank()
  )

# ---------------------------------------------------------------------------
# Subpanel B — Currin caQTL per-TF agreement (lollipop)
# ---------------------------------------------------------------------------
caqtl <- fread(file.path(GWAS_ATAC_DIR, "currin_caqtl_per_tf.csv"))

# 4-way validated TFs (from headline: THRB, HNF4A, RORA, MLXIPL, MAX)
fourway_tfs <- c("THRB", "HNF4A", "RORA", "MLXIPL", "MAX")

# Filter: TFs with n_caqtl_tested >= 4 OR 4-way validated (always retain
# the 4-way set even if MAX has n=1) AND a non-NA agreement %. The 4-way set
# is the scientific anchor of this subpanel; the remaining 10 slots show the
# strongest disease-regulon + high-n agreement TFs as comparators.
caqtl[, is_4way := tf_name %in% fourway_tfs]

caqtl_pool <- caqtl[!is.na(pct_motif_caqtl_agree) &
                    (n_caqtl_tested >= 4 | is_4way == TRUE)]

# Rank: 4-way first (in headline order, since those are the named TFs);
# remaining slots filled by disease-regulon TFs at n>=4 ranked by % agree,
# then by other top-n high-agreement TFs.
fourway_rows <- caqtl_pool[is_4way == TRUE]
fourway_rows[, rank_4way := match(tf_name, fourway_tfs)]
setorder(fourway_rows, rank_4way)

remaining <- caqtl_pool[is_4way == FALSE & n_caqtl_tested >= 4]
# Priority: disease-regulon members first, then by % agreement, then by n
setorder(remaining, -any_disease_regulon, -pct_motif_caqtl_agree,
         -n_caqtl_tested)

n_slots_remaining <- max(0, 15 - nrow(fourway_rows))
caqtl_top <- rbindlist(list(fourway_rows[, !"rank_4way"],
                            head(remaining, n_slots_remaining)),
                       use.names = TRUE, fill = TRUE)

# Order display: descending % agreement so it reads as a true lollipop ranking
setorder(caqtl_top, -pct_motif_caqtl_agree, -n_caqtl_tested)

# Order by % agreement (top -> bottom)
caqtl_top[, tf_name := factor(tf_name, levels = rev(caqtl_top$tf_name))]
caqtl_top[, point_color := ifelse(is_4way, MAGENTA_DEEP, BLUE)]
caqtl_top[, lollipop_segment := BLUE]
caqtl_top[is_4way == TRUE, lollipop_segment := MAGENTA_DEEP]
caqtl_top[, annot := sprintf("n=%d/%d", n_motif_caqtl_agree, n_caqtl_tested)]

# Layer points and segments to expose the gold outline on the 4-way TFs
p_B <- ggplot(caqtl_top, aes(y = tf_name)) +
  geom_vline(xintercept = 50, linetype = "dashed",
             color = "gray60", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = pct_motif_caqtl_agree, yend = tf_name,
                   color = is_4way), linewidth = 0.5,
               show.legend = FALSE) +
  # Gold outline ring for 4-way TFs (slightly larger)
  geom_point(data = caqtl_top[is_4way == TRUE],
             aes(x = pct_motif_caqtl_agree),
             color = GOLD, size = 3.2, shape = 16) +
  geom_point(aes(x = pct_motif_caqtl_agree, fill = is_4way,
                 color = is_4way),
             size = 2.0, shape = 21, stroke = 0.6) +
  scale_color_manual(values = c("FALSE" = BLUE, "TRUE" = MAGENTA_DEEP)) +
  scale_fill_manual(values  = c("FALSE" = BLUE, "TRUE" = MAGENTA_DEEP)) +
  scale_x_continuous(limits = c(0, 105),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(tag = "b",
       x = "caQTL concordance (%)",
       y = NULL) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "bold", size = PUB_AXIS_TEXT,
                               color = "black"),
    axis.title.x = element_text(face = "bold"),
    plot.tag = element_text(size = PUB_TITLE + 1, face = "bold"),
    plot.tag.position = c(0.02, 0.97),
    panel.grid = element_blank()
  )

# ---------------------------------------------------------------------------
# Compose via patchwork (A: ~2.5in, B: ~3.5in -> total ~6in x 3in)
# ---------------------------------------------------------------------------
composed <- p_A + p_B + plot_layout(widths = c(2.5, 3.5))

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
ggsave(OUT_PDF, composed, width = 6.0, height = 3.0, device = cairo_pdf)
cat(sprintf("[Fig 3 orthogonal] Wrote %s\n", OUT_PDF))

# Sanity readout
cat(sprintf("[A] active variants: %d / %d (%.1f%%)\n",
            n_active, n_total, pct_active))
cat("[B] Top-15 TFs by caQTL agreement (n_caqtl_tested >= 4):\n")
print(caqtl_top[, .(tf_name, pct_motif_caqtl_agree,
                    n_caqtl_tested, n_motif_caqtl_agree,
                    is_4way, any_disease_regulon)])
