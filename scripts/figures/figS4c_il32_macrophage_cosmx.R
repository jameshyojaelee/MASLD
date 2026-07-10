#!/usr/bin/env Rscript
# ==============================================================================
# Fig S4c — CosMx single-cell proximity (Govaere 2026): IL32-hepatocyte -> macrophage.
#
# Split out of fig4h_spatial_ccc_consensus.R 2026-07-08 (user call: that panel crammed
# two genuinely distinct exhibits into one main-text callout; this one is a narrower
# single-mechanism deep-dive better suited to supplement).
#
# An INDEPENDENT cohort AND an independent single-cell-resolution technology (522k
# segmented cells, 968-gene panel). A finding Visium's 55-µm spots cannot resolve:
# IL32-high hepatocytes sit CLOSER to the macrophage compartment, and their nearest
# macrophages carry more CD74 (MIF / MHC-II receptor). Significance = slide-level
# DIRECTION concordance across 3 MASH slides (cells within a slide are
# pseudoreplicated), NOT a cell-level p.
#   pA = per-slide nearest-macrophage distance shift (IL32-high − IL32-low), µm;
#        negative = IL32-high hepatocytes closer.
#   pB = per-slide Spearman ρ between hepatocyte IL32 and mean CD74 of its 5
#        nearest macrophages.
#
# PROXIMITY-level recovery of inferred signaling, NOT cell-type-resolved proof.
# Main-text companion: fig4h_spatial_ccc_consensus.R (Visium spot-adjacency, GSE192741).
#
# Data: Analysis/Spatial/results/govaere2026/il32_colocalization/il32_coloc_per_slide.csv
# Output: figures/main/fig4_validation/panels/figS4c.pdf
# Env: rnaseq
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

ps <- fread(file.path(BASE,
  "Analysis/Spatial/results/govaere2026/il32_colocalization/il32_coloc_per_slide.csv"))
ps[, is_mash := tolower(as.character(is_mash)) == "true"]
ps <- ps[order(is_mash, slide)]
ps[, label := ifelse(is_mash, slide, paste0(slide, " (Normal)"))]
ps[, grp   := ifelse(is_mash, "MASH", "Normal")]
ps[, label := factor(label, levels = label)]
cos_col <- c("MASH" = masld_colors$up, "Normal" = masld_colors$ns)
n_mash <- sum(ps$is_mash)
n_dist <- sum(ps$is_mash & ps$delta_high_minus_low_um < 0)
n_cd74 <- sum(ps$is_mash & ps$rho_il32_vs_knn_cd74 > 0)
mean_delta <- mean(ps$delta_high_minus_low_um[ps$is_mash])
mean_rho   <- mean(ps$rho_il32_vs_knn_cd74[ps$is_mash])
cat(sprintf("[figS4c] CosMx: distance %d/%d MASH closer (mean Δ %.2f µm); nbr-CD74 %d/%d MASH positive (mean ρ %+.3f)\n",
            n_dist, n_mash, mean_delta, n_cd74, n_mash, mean_rho))

# pA — nearest-macrophage distance shift (IL32-high − IL32-low), µm; <0 = closer (BARS, not lollipop)
pA <- ggplot(ps, aes(delta_high_minus_low_um, label, fill = grp)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_col(width = 0.62, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = cos_col, guide = "none") +
  scale_x_continuous(limits = c(-1.6, 0.1), breaks = c(-1.5, -1, -0.5, 0)) +
  labs(x = "nearest-macrophage Δ (µm)", y = NULL) +
  ggtitle("Govaere CosMx · IL32-hep → macrophage") +
  theme_masld() +
  theme(plot.title = element_text(size = 6, hjust = 0, face = "plain"),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(2, 2, 2, 2))

# pB — hepatocyte IL32 vs mean CD74 of 5 nearest macrophages (BARS)
pB <- ggplot(ps, aes(rho_il32_vs_knn_cd74, label, fill = grp)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_col(width = 0.62, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = cos_col, name = NULL,
                    guide = guide_legend(keyheight = unit(0.3, "lines"))) +
  scale_x_continuous(limits = c(-0.05, 0.15), breaks = c(0, 0.05, 0.10, 0.15)) +
  labs(x = "ρ: IL32 vs neighbour CD74", y = NULL) +
  ggtitle("IL32 vs neighbour CD74") +
  theme_masld() +
  theme(plot.title = element_text(size = 6, hjust = 0, face = "plain"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = c(0.72, 0.20), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(0.28, "lines"),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(2, 2, 2, 2))

p <- pA + pB + plot_layout(widths = c(1.15, 1))

out <- file.path(FIG4_DIR, "panels", "figS4c.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
ggsave(out, p, width = 4.6, height = 1.9, device = grDevices::cairo_pdf)
cat("[figS4c] saved:", out, "\n")

message(sprintf(paste0(
  "CAPTION (Fig S4c): CosMx single-cell proximity (Govaere 2026 — an INDEPENDENT cohort and an ",
  "independent single-cell-resolution technology, 522k cells / 968-gene panel): a finding Visium's ",
  "55-µm spots cannot resolve. Left: per slide, the shift in median nearest-macrophage distance between ",
  "IL32-high and IL32-low hepatocytes (µm); negative = IL32-high hepatocytes closer, %d/%d MASH slides ",
  "(mean %.2f µm). Right: per slide, Spearman ρ between hepatocyte IL32 and the mean CD74 of its 5 ",
  "nearest macrophages, %d/%d MASH slides positive (mean %+.3f). Significance = slide-level DIRECTION ",
  "concordance across the 3 MASH slides, NOT a cell-level p (cells within a slide are pseudoreplicated); ",
  "Leuven_2 is a mixed-Normal slide (grey, not in the concordance count). The recovered MetMac is ",
  "MHC-II-high, NOT the canonical GPNMB+ lipid-associated macrophage (968-gene panel limit) — no LAM ",
  "identity claim. PROXIMITY-level recovery of inferred signaling, NOT cell-type-resolved proof. Main-",
  "text companion: Fig 4h (Visium spot-adjacency, GSE192741)."),
  n_dist, n_mash, mean_delta, n_cd74, n_mash, mean_rho))
