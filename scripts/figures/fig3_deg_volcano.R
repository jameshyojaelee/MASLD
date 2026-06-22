#!/usr/bin/env Rscript
# fig3_deg_volcano.R — Canonical bulk DEG volcano for Fig 3 (RNA-seq).
#
# Source: canonical_deg_results.csv (limma-voom quality-weighted, C2 design;
# canonical since 2026-06-08, superseding dream). Loaded via load_dream_results()
# which emits bulk_* column names.
#
# Tier 1 DEG definition (ashr): lfsr < 0.05 AND |shrunk_logFC| > 0.5.
# 5 control-bearing cohorts, 846 samples. Cohort presentation: 9 cohorts
# collected, 1,259 QC-passing (PRJNA512027 excluded for L0/S0 batch confound).
#
# Output: figures/main/fig3_RNAseq/panels/figs3a_deg_volcano.pdf
# (moved from fig1_atlas_overview 2026-06-12 — this is an RNA-seq DEG panel.)

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

PANEL_DIR <- file.path(FIG2_DIR, "panels")   # FIG2_DIR = figures/main/fig3_RNAseq
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

PADJ_CUT  <- 0.05   # applied to bulk_padj (BH-adjusted P); matches the manuscript Fig 3C "padj<0.05" definition
LFC_CUT   <- 0.5    # applied to bulk_logFC (raw log2 fold change)
N_TOP_DIR <- 8     # top-padj genes labelled per direction

# DEG_VOLCANO_VECTOR=TRUE renders point clouds as native vector geometry and
# writes to deg_volcano_vector.pdf (heavier file, but fully editable).
VECTOR_MODE <- identical(toupper(Sys.getenv("DEG_VOLCANO_VECTOR", "FALSE")), "TRUE")

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

message("Loading canonical DEG results (limma-voom C2)...")
dream <- load_dream_results()

# ----------------------------------------------------------------------------
# Build plotting frame
# ----------------------------------------------------------------------------
# x = raw log2 fold change, y = -log10(BH-adjusted P). padj (unlike ashr lfsr)
# never underflows to exactly 0 here (min padj = 1.5e-52), so the y-axis is
# continuous: no pile-up stripe at the top and no hollow band on the down side
# — the two artefacts the old -log10(lfsr) axis produced from 241 lfsr==0 genes.
volc <- dream[!is.na(bulk_padj) & !is.na(bulk_logFC),
              .(symbol, bulk_logFC, bulk_padj, bulk_shrunk_logFC, bulk_lfsr)]
# Defensive guard only (no padj==0 in the current data): floor any exact zero
# at the smallest positive padj so -log10 stays finite.
min_pos_padj <- min(volc[bulk_padj > 0, bulk_padj], na.rm = TRUE)
volc[bulk_padj <= 0, bulk_padj := min_pos_padj]
volc[, neglog10padj := -log10(bulk_padj)]

# Three-class direction × significance — padj<0.05 & |raw log2FC|>0.5, matching
# the manuscript Fig 3C definition (1,853 DEGs: 1,438 up / 415 down).
volc[, status := fcase(
  bulk_padj < PADJ_CUT &  bulk_logFC >  LFC_CUT, "Up",
  bulk_padj < PADJ_CUT &  bulk_logFC < -LFC_CUT, "Down",
  default = "n.s."
)]
volc[, status := factor(status, levels = c("n.s.", "Down", "Up"))]
setorder(volc, status)   # n.s. plotted first, sig on top

n_up   <- sum(volc$status == "Up")
n_down <- sum(volc$status == "Down")
n_ns   <- sum(volc$status == "n.s.")
message(sprintf("DEG counts at padj<%.2g, |log2FC|>%.1f: %s up, %s down, %s n.s.",
                PADJ_CUT, LFC_CUT, comma(n_up), comma(n_down), comma(n_ns)))

# ----------------------------------------------------------------------------
# Label set: top-padj per direction + curated anchors
# ----------------------------------------------------------------------------
sig_up   <- volc[status == "Up"   & symbol != "" & !grepl("^ENSG", symbol)]
sig_down <- volc[status == "Down" & symbol != "" & !grepl("^ENSG", symbol)]
setorder(sig_up,   bulk_padj)
setorder(sig_down, bulk_padj)

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

x_lim <- max(abs(volc$bulk_logFC), na.rm = TRUE) * 1.04
y_max <- max(volc$neglog10padj, na.rm = TRUE)
y_lim <- y_max * 1.05

p <- ggplot(volc, aes(x = bulk_logFC, y = neglog10padj, color = status)) +
  # Dashed thresholds — drawn behind everything
  geom_hline(yintercept = -log10(PADJ_CUT),
             linetype = "dashed", color = "gray70", linewidth = 0.25) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT),
             linetype = "dashed", color = "gray70", linewidth = 0.25) +
  # Volcano cloud — rasterised by default (lean PDF); native vector when
  # DEG_VOLCANO_VECTOR=TRUE
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
             aes(x = bulk_logFC, y = neglog10padj, fill = status),
             color = "black", shape = 21, size = 1.25,
             stroke = 0.25, inherit.aes = FALSE) +
  # Labels
  geom_text_repel(data = label_df,
                  aes(x = bulk_logFC, y = neglog10padj, label = symbol),
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
    x = expression("log"[2]*" fold change"),
    y = expression(-log[10]~"adjusted "*italic(P))
  ) +
  theme_masld(base_size = 7) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(linewidth = 0.18, color = "gray92"),
    legend.position    = "none",
    axis.title         = element_text(size = 7.5),
    plot.margin        = margin(4, 6, 2, 4)
  )

out_name <- if (VECTOR_MODE) "deg_volcano_vector.pdf" else "figs3a_deg_volcano.pdf"
out_pdf  <- file.path(PANEL_DIR, out_name)
save_fig(p, out_pdf, width = fig_half_width * 1.15, height = 3.2)

# Companion CSV: every labelled gene (only on the rasterised default run, to
# avoid clobbering when the vector pass is rendered after edits)
if (!VECTOR_MODE) {
  fwrite(label_df[, .(symbol, bulk_shrunk_logFC, bulk_lfsr, bulk_logFC, bulk_padj, status)],
         file.path(PANEL_DIR, "deg_volcano_labels.csv"))
}

if (file.exists(out_pdf)) {
  message(sprintf("\nOutput: %s (%s)", out_pdf,
                  utils:::format.object_size(file.size(out_pdf), "auto")))
}
