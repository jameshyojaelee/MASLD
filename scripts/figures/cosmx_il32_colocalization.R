#!/usr/bin/env Rscript
# KEY MESSAGE (single-cell spatial novelty -- Visium's 55-um spots cannot resolve
# this): in the Govaere CosMx single-cell spatial data, hepatocytes that express
# IL32 sit closer to the macrophage compartment, and their nearest macrophages
# carry more CD74 (the MIF / MHC-II receptor chain). The direction is concordant
# in ALL THREE MASH slides. Significance = slide-level direction concordance
# (n=3 MASH), NOT cell-level p (cells within a slide are pseudoreplicated).
#
# Panel A: per slide, median nearest-macrophage distance shift
#          (IL32-high minus IL32-low hepatocytes), in um. Negative = IL32-high
#          hepatocytes are closer. 3/3 MASH slides negative.
# Panel B: per slide, Spearman rho between hepatocyte IL32 and the mean CD74 of
#          its 5 nearest macrophages. 3/3 MASH slides positive.
#
# CAVEAT (legend): effect on raw distance is small because macrophages are
# ubiquitous; the reproducible signal is the consistent DIRECTION across slides
# and the neighbour-CD74 relationship. The recovered "MetMac" is MHC-II-high,
# NOT the canonical GPNMB+ lipid-associated macrophage (968-gene panel) -- no LAM
# claim. Leuven_2 is a mixed-Normal slide (shown faded, not in the concordance).
#
# Source: Analysis/Spatial/results/govaere2026/il32_colocalization/il32_coloc_per_slide.csv
# Output: figures/main/fig5_molecular_context/cosmx_il32_colocalization.pdf
# Env: rnaseq
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

ps <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/govaere2026/il32_colocalization/il32_coloc_per_slide.csv"),
  stringsAsFactors = FALSE)
ps$is_mash <- tolower(as.character(ps$is_mash)) == "true"
ps$il32high_closer <- tolower(as.character(ps$il32high_closer)) == "true"
ps <- ps %>% mutate(
  label = ifelse(is_mash, slide, paste0(slide, " (Normal)")),
  grp   = ifelse(is_mash, "MASH", "Normal")) %>%
  arrange(is_mash, slide) %>%
  mutate(label = factor(label, levels = label))
cols <- c(MASH = masld_colors$up, Normal = masld_colors$ns)
n_mash <- sum(ps$is_mash)
n_dist <- sum(ps$is_mash & ps$delta_high_minus_low_um < 0)
n_cd74 <- sum(ps$is_mash & ps$rho_il32_vs_knn_cd74 > 0)

# ── Panel A: nearest-macrophage distance shift (IL32-high - IL32-low) ─────────
pA <- ggplot(ps, aes(delta_high_minus_low_um, label, color = grp)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_segment(aes(x = 0, xend = delta_high_minus_low_um, yend = label), linewidth = 1.0) +
  geom_point(size = 2.8) +
  scale_color_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(-1.6, 0.25)) +
  labs(x = "Nearest-macrophage distance shift (um)", y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.title = element_text(size = 6), axis.text = element_text(size = 6, color = "black"),
        panel.grid.major.y = element_blank())

# ── Panel B: hepatocyte IL32 vs nearest-macrophage CD74 ──────────────────────
pB <- ggplot(ps, aes(rho_il32_vs_knn_cd74, label, color = grp)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_segment(aes(x = 0, xend = rho_il32_vs_knn_cd74, yend = label), linewidth = 1.0) +
  geom_point(size = 2.8) +
  scale_color_manual(values = cols, name = NULL) +
  scale_x_continuous(limits = c(-0.06, 0.17), breaks = c(0, 0.05, 0.10, 0.15)) +
  labs(x = "rho: IL32 vs neighbour CD74", y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.title = element_text(size = 6), axis.text = element_text(size = 6, color = "black"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = c(0.78, 0.22), legend.text = element_text(size = 6),
        legend.key.size = unit(0.32, "cm"))

p_out <- (pA | pB) + plot_layout(widths = c(1.18, 1))
out <- file.path(FIG4_DIR, "_supp", "cosmx_il32_colocalization.pdf")
cairo_pdf(out, width = fig_full_width * 0.98, height = 2.5)
print(p_out); dev.off()
message(sprintf("CosMx IL32: distance %d/%d MASH closer; neighbour-CD74 %d/%d MASH positive",
                n_dist, n_mash, n_cd74, n_mash))
message("Saved: ", out)
