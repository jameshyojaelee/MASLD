#!/usr/bin/env Rscript
# fig1_volcano.R — Integrated dream mega-analysis volcano plot for Fig 1.
#
# Tier 1 thresholds: padj < 0.05, |log2FC| > 0.5.
# 1,918 primary DEGs from the dream mega-analysis (5 control-bearing
# cohorts, 847 samples). Cohort presentation: 9 cohorts collected,
# 1,259 QC-passing (PRJNA512027 excluded for L0/S0 batch confound).
#
# Output: figures/main/fig1_atlas_overview/panels/fig1_volcano.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG1_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

LFSR_CUT  <- 0.05   # applied to dream_lfsr (ashr local false sign rate)
LFC_CUT   <- 0.5    # applied to dream_shrunk_logFC (ashr-shrunk effect size)
N_TOP_DIR <- 8     # top-padj genes labelled per direction

# FIG1_VOLCANO_VECTOR=TRUE renders point clouds as native vector geometry and
# writes to fig1_volcano_vector.pdf (heavier file, but fully editable).
VECTOR_MODE <- identical(toupper(Sys.getenv("FIG1_VOLCANO_VECTOR", "FALSE")), "TRUE")

# Curated MASLD-relevant anchor genes — labelled if Tier 1.
CURATED <- c(
  "THRB",     # resmetirom target
  "HNF4A",    # master hepatocyte TF
  "NR1H4",    # FXR / obeticholic acid
  "PPARA",    # elafibranor / fenofibrate
  "RORA",     # cross-ancestry COLOC
  "AKR1B10",  # largest bulk effect, MASLD biomarker
  "TREM2",    # LAM macrophage marker
  "LPL",      # lipoprotein lipase, lipid metabolism
  "GSTM1",    # antioxidant
  "POSTN",    # fibrosis ECM
  "GSN",      # advanced fibrosis
  "FABP4",    # macrophage / lipid
  "SPP1",     # scar-associated macrophage
  "COL1A1",   # fibrosis
  "HKDC1"     # F1->F2 progression driver
)

message("Loading dream results...")
dream <- load_dream_results()

# ----------------------------------------------------------------------------
# Build plotting frame
# ----------------------------------------------------------------------------
volc <- dream[!is.na(dream_lfsr) & !is.na(dream_shrunk_logFC),
              .(symbol, dream_shrunk_logFC, dream_lfsr, dream_logFC, dream_padj)]
volc[, neglog10lfsr := -log10(dream_lfsr)]

# Three-class direction × significance (lfsr + shrunk LFC, canonical 2026-06-02)
volc[, status := fcase(
  dream_lfsr < LFSR_CUT &  dream_shrunk_logFC >  LFC_CUT, "Up",
  dream_lfsr < LFSR_CUT &  dream_shrunk_logFC < -LFC_CUT, "Down",
  default = "n.s."
)]
volc[, status := factor(status, levels = c("n.s.", "Down", "Up"))]
setorder(volc, status)   # n.s. plotted first, sig on top

n_up   <- sum(volc$status == "Up")
n_down <- sum(volc$status == "Down")
n_ns   <- sum(volc$status == "n.s.")
message(sprintf("DEG counts at lfsr<%.2g, |shrunk_logFC|>%.1f: %s up, %s down, %s n.s.",
                LFSR_CUT, LFC_CUT, comma(n_up), comma(n_down), comma(n_ns)))

# ----------------------------------------------------------------------------
# Label set: top-padj per direction + curated anchors
# ----------------------------------------------------------------------------
sig_up   <- volc[status == "Up"   & symbol != "" & !grepl("^ENSG", symbol)]
sig_down <- volc[status == "Down" & symbol != "" & !grepl("^ENSG", symbol)]
setorder(sig_up,   dream_lfsr)
setorder(sig_down, dream_lfsr)

top_up   <- head(sig_up,   N_TOP_DIR)
top_down <- head(sig_down, N_TOP_DIR)

curated_lbl <- volc[symbol %in% CURATED & status != "n.s." & !grepl("^ENSG", symbol)]

label_df <- unique(rbind(top_up, top_down, curated_lbl), by = "symbol")
message(sprintf("Labelling %d genes (%d top-padj per direction + %d curated anchors)",
                nrow(label_df), N_TOP_DIR, nrow(curated_lbl)))

# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
volc_colors <- c(
  "Up"   = masld_colors$up,    # #C9265E Liang magenta
  "Down" = masld_colors$down,  # #1565C0 blue
  "n.s." = "#D8D8D8"           # very light gray
)

x_lim <- max(abs(volc$dream_shrunk_logFC), na.rm = TRUE) * 1.04
y_max <- max(volc$neglog10lfsr, na.rm = TRUE)
y_lim <- y_max * 1.05

p <- ggplot(volc, aes(x = dream_shrunk_logFC, y = neglog10lfsr, color = status)) +
  # Dashed thresholds — drawn behind everything
  geom_hline(yintercept = -log10(LFSR_CUT),
             linetype = "dashed", color = "gray70", linewidth = 0.25) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT),
             linetype = "dashed", color = "gray70", linewidth = 0.25) +
  # Volcano cloud — rasterised by default (lean PDF); native vector when
  # FIG1_VOLCANO_VECTOR=TRUE
  (if (VECTOR_MODE) {
     geom_point(data = volc[status == "n.s."],
                size = 0.35, alpha = 0.45, shape = 16)
   } else {
     rasterize_layer(geom_point(data = volc[status == "n.s."],
                                size = 0.35, alpha = 0.45, shape = 16))
   }) +
  (if (VECTOR_MODE) {
     geom_point(data = volc[status != "n.s."],
                size = 0.55, alpha = 0.85, shape = 16)
   } else {
     rasterize_layer(geom_point(data = volc[status != "n.s."],
                                size = 0.55, alpha = 0.85, shape = 16))
   }) +
  # Halo around labelled points so they pop
  geom_point(data = label_df,
             aes(x = dream_shrunk_logFC, y = neglog10lfsr, fill = status),
             color = "black", shape = 21, size = 1.25,
             stroke = 0.25, inherit.aes = FALSE) +
  # Labels
  geom_text_repel(data = label_df,
                  aes(x = dream_shrunk_logFC, y = neglog10lfsr, label = symbol),
                  inherit.aes = FALSE,
                  size = 2.2, color = "black", fontface = "italic",
                  segment.size = 0.2, segment.color = "gray45",
                  box.padding = 0.4, point.padding = 0.2,
                  min.segment.length = 0,
                  max.overlaps = Inf, force = 4, seed = 42,
                  show.legend = FALSE) +
  # Per-direction count annotations
  annotate("text", x = -x_lim * 0.98, y = y_lim * 0.97,
           label = sprintf("%s ↓", comma(n_down)),
           hjust = 0, vjust = 1, size = 2.4, fontface = "bold",
           color = masld_colors$down) +
  annotate("text", x =  x_lim * 0.98, y = y_lim * 0.97,
           label = sprintf("%s ↑", comma(n_up)),
           hjust = 1, vjust = 1, size = 2.4, fontface = "bold",
           color = masld_colors$up) +
  scale_color_manual(values = volc_colors, guide = "none") +
  scale_fill_manual(values = volc_colors, guide = "none") +
  scale_x_continuous(limits = c(-x_lim, x_lim),
                     expand = expansion(mult = 0),
                     breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(0, y_lim),
                     expand = expansion(mult = c(0, 0)),
                     breaks = pretty_breaks(n = 5)) +
  labs(
    x = expression("ashr shrunk log"[2]*" fold change"),
    y = expression(-log[10]~"lfsr")
  ) +
  theme_masld(base_size = 7) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(linewidth = 0.18, color = "gray92"),
    legend.position    = "none",
    axis.title         = element_text(size = 7.5),
    plot.margin        = margin(4, 6, 2, 4)
  )

out_name <- if (VECTOR_MODE) "fig1_volcano_vector.pdf" else "fig1_volcano.pdf"
out_pdf  <- file.path(PANEL_DIR, out_name)
save_fig(p, out_pdf, width = fig_half_width * 1.15, height = 3.2)

# Companion CSV: every labelled gene (only on the rasterised default run, to
# avoid clobbering when the vector pass is rendered after edits)
if (!VECTOR_MODE) {
  fwrite(label_df[, .(symbol, dream_shrunk_logFC, dream_lfsr, dream_logFC, dream_padj, status)],
         file.path(PANEL_DIR, "fig1_volcano_labels.csv"))
}

if (file.exists(out_pdf)) {
  message(sprintf("\nOutput: %s (%s)", out_pdf,
                  utils:::format.object_size(file.size(out_pdf), "auto")))
}
