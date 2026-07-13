#!/usr/bin/env Rscript
# =============================================================================
# Figure 05 -- Human vs Mouse DEG comparison for the Cas13 library (ashr basis)
# KEY MESSAGE: the library's human spine and mouse-confirmed tier reflect a shared
#   cross-species program; human and mouse ashr effect sizes are directionally
#   concordant. (Transcriptomic only -- no genetic/COLOC evidence.)
#
# DEG definitions = the library's own ashr instruments:
#   Human DEG : lfsr<0.05 & |shrunk_logFC| > 0.2                 (human spine)
#   Mouse DEG : lfsr<0.05 & |shrunk_logFC| > 0.5 in >=3 of 4 diets (mouse-confirmed)
# Effect axes = ashr shrunk_logFC (mouse = mean shrunk across the 4 diets).
# Ortholog-aligned (1 mouse ortholog per human, one2one-preferred) like the rebuild.
#
# Panels:
#   A  Scatter human vs mouse shrunk_logFC (DEG-in-either), colored by category
#   B  Overlap counts (human-only / shared-concordant / shared-discordant / mouse-only)
# Output: Cas13_Library_Design/figures/05_human_vs_mouse_degs.pdf
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT_DIR <- FIGS_CAS13LIB_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
strip_v <- function(x) sub("[.][0-9]+$", "", x)

INT     <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
PERDIET <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet_cas13")  # Cas13 library Western pool; decoupled from paper 4-model per_diet (2026-06-16)
ORTHO   <- file.path(BASE, "data/external/orthologs/master_ortholog_table.tsv.gz")
DIETS   <- c("MCD", "CDAHFD", "Western", "HFD")
LFSR <- 0.05; H_SHRUNK <- 0.2; M_SHRUNK <- 0.5; MIN_DIETS <- 3L

# -- human ashr ---------------------------------------------------------------
ash <- fread(file.path(INT, "canonical_deg_results.csv"), select = c("gene", "shrunk_logFC", "lfsr", "symbol"))  # canonical limma-voom QW C2 (2026-06-24: was metafor meta_results_ashr.csv)
ash[, hb := strip_v(gene)]
ash[, hDEG := !is.na(lfsr) & lfsr < LFSR & abs(shrunk_logFC) > H_SHRUNK]

# -- mouse ashr per-diet -> per-gene mean shrunk + cross-diet sig count -------
mall <- rbindlist(lapply(DIETS, function(d) {
  x <- fread(file.path(PERDIET, paste0(d, "_de_results.csv")))
  x[, .(gb = strip_v(gene), shrunk = shrunk_logFC,
        sig = !is.na(lfsr) & lfsr < LFSR & abs(shrunk_logFC) > M_SHRUNK)]
}))
mouse <- mall[, .(mouse_mean_shrunk = mean(shrunk, na.rm = TRUE),
                  n_diets_sig = sum(sig)), by = gb]
mouse[, mDEG := n_diets_sig >= MIN_DIETS]

# -- ortholog: one mouse per human (deterministic, one2one-preferred) ---------
ortho <- fread(cmd = paste0("zcat ", ORTHO),
               select = c("mouse_ensembl", "human_ensembl", "confidence_tier", "is_one2one"))
ortho[, mouse_ensembl := strip_v(mouse_ensembl)]; ortho[, human_ensembl := strip_v(human_ensembl)]
ortho <- ortho[confidence_tier %in% c("H", "M")]
ortho[, trank := match(confidence_tier, c("H", "M"))]
ortho[, o2o := ifelse(is_one2one %in% c(TRUE, "True", "TRUE", "true"), 0L, 1L)]
setorder(ortho, human_ensembl, trank, o2o, mouse_ensembl)
ortho_byhuman <- unique(ortho, by = "human_ensembl")

# -- align human + mouse ortholog ---------------------------------------------
hm <- merge(ash[, .(hb, symbol, h_shrunk = shrunk_logFC, h_lfsr = lfsr, hDEG)],
            ortho_byhuman[, .(hb = human_ensembl, mgb = mouse_ensembl)], by = "hb")
hm <- merge(hm, mouse[, .(mgb = gb, mouse_mean_shrunk, n_diets_sig, mDEG)], by = "mgb")
hm[, category := fcase(
  hDEG & mDEG & sign(h_shrunk) == sign(mouse_mean_shrunk), "Shared concordant",
  hDEG & mDEG,                                             "Shared discordant",
  hDEG & !mDEG,                                            "Human-only",
  !hDEG & mDEG,                                            "Mouse-only",
  default = "Not DEG")]

shared <- hm[hDEG & mDEG]
conc_pct <- 100 * mean(sign(shared$h_shrunk) == sign(shared$mouse_mean_shrunk))
rho <- cor(shared$h_shrunk, shared$mouse_mean_shrunk, method = "spearman")
cat(sprintf("aligned genes: %d | human DEG: %d | mouse DEG: %d | shared: %d (%.0f%% concordant, rho=%.2f)\n",
            nrow(hm), sum(hm$hDEG), sum(hm$mDEG), nrow(shared), conc_pct, rho))

cat_levels <- c("Shared concordant", "Shared discordant", "Human-only", "Mouse-only")
pal <- c("Shared concordant" = "#2E7D32", "Shared discordant" = "#C62828",
         "Human-only" = "#1565C0", "Mouse-only" = "#F57C00")
deg <- hm[category != "Not DEG"]
deg[, category := factor(category, levels = cat_levels)]

lim <- max(abs(c(deg$h_shrunk, deg$mouse_mean_shrunk)), na.rm = TRUE)
pA <- ggplot(deg, aes(h_shrunk, mouse_mean_shrunk, color = category)) +
  geom_hline(yintercept = 0, color = "gray70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "gray70", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = pal, name = NULL) +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "Human shrunk log2FC (ashr)", y = "Mouse mean shrunk log2FC (ashr, 4 diets)") +
  theme_masld() + theme_pub() +
  theme(legend.position = "right", legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6))

cnts <- deg[, .N, by = category][order(match(category, cat_levels))]
cnts[, category := factor(category, levels = rev(cat_levels))]
pB <- ggplot(cnts, aes(N, category, fill = category)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = scales::comma(N)), hjust = -0.1, size = GEOM_TEXT_6PT) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)), labels = scales::comma) +
  labs(x = "Genes", y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_text(size = 6))

fig <- pA + pB + plot_layout(widths = c(1.6, 1)) +
  plot_annotation(caption = "Human DEG: lfsr<0.05 & |shrunk log2FC|>0.2 (spine).  Mouse DEG: lfsr<0.05 & |shrunk|>0.5 in >=3 of 4 diets (mouse-confirmed).")
ggsave(file.path(OUT_DIR, "human_vs_mouse_degs.pdf"), fig,
       width = 7.2, height = 3.6, device = pdf_device)
cat("Saved: 05_human_vs_mouse_degs.pdf (ashr basis)\n")
