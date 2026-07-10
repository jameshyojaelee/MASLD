#!/usr/bin/env Rscript
# crossmodal_module_proteomics.R — SUPPLEMENTARY figure: Fig3 Hotspot single-cell modules x Fig4 liver
# proteomics (PXD051911 DIA-MS, n=58, documented histology). Reads the linkage tables written by
# Analysis/Proteomics/scripts/module_proteomics_linkage.R and renders individual panel PDFs.
#
#   A  figS_crossmod_enrichment.pdf   — fgsea NES of each disease-significant module gene set over the
#                                       liver protein moderated-t rank (which transcriptomic modules
#                                       reach the proteome), coloured by scRNA direction, * = padj<0.05.
#   B  figS_crossmod_concordance.pdf  — beta_scRNA (fig3g 3-stage slope) vs beta_protein (No_MASLD ->
#                                       MASL -> MASH ordinal slope) per module; Spearman rho + quadrants.
#   C  figS_crossmod_cascade_modules  — per-module protein score across No_MASLD/MASL/MASH, split by
#                                       scRNA up/down (the fig3e handoff, in liver protein).
#   C2 figS_crossmod_cascade_programs — the same for the bulk NMF k=6 programs (companion).
#
# Conventions: 6 pt Helvetica plain, black text, no colored fonts/titles, control/ref gray #9E9E9E,
# individual PDFs, captions via message(). Bars/heatmaps/scatters/lines only (no lollipop).
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
HAS_REPEL <- requireNamespace("ggrepel", quietly = TRUE)
set.seed(42)
FAM  <- "Helvetica"
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
RES  <- file.path(BASE, "Analysis/Proteomics/results")
SUPP <- file.path(FIG_SUPP, "figS_crossmodal_module_proteomics")
dir.create(SUPP, showWarnings = FALSE, recursive = TRUE)

CT_MAP  <- c(hepatocytes = "Hepatocytes", macrophages = "Macrophages", fibroblasts = "Fibroblasts",
             cholangiocytes = "Cholangiocytes", tcells = "T cells")
CT_COLS <- ct_palette[CT_MAP]                         # named by pretty label
DIR_COLS <- c(up = masld_colors$up, down = masld_colors$down)
GREY <- "#9E9E9E"; TXT <- 6 / ggplot2::.pt
prettyct <- function(x) factor(CT_MAP[x], levels = CT_MAP)

# ── Panel A: module enrichment in the liver proteome ──────────────────────────
enr <- fread(file.path(RES, "module_proteomics_enrichment.csv"))[cohort == "PXD051911_liver"]
enr[, `:=`(ct = prettyct(cell_type), sig = padj < 0.05,
           lab = paste0(module_name),
           lab_it = paste0(substr(module_name, 1, 34)))]
setorder(enr, cell_type, NES)
enr[, row := factor(mid, levels = mid)]
pA <- ggplot(enr, aes(NES, row, fill = scrna_dir)) +
  geom_col(width = 0.72, aes(color = sig), linewidth = 0.3) +
  geom_vline(xintercept = 0, color = GREY, linewidth = 0.3) +
  geom_text(data = enr[sig == TRUE], aes(x = NES + sign(NES) * 0.12, label = "*"),
            size = TXT * 1.4, color = "black", vjust = 0.75) +
  scale_fill_manual(values = DIR_COLS, name = "scRNA dir") +
  scale_color_manual(values = c(`TRUE` = "black", `FALSE` = NA), guide = "none") +
  scale_y_discrete(labels = setNames(enr$lab_it, enr$row)) +
  facet_grid(ct ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Protein enrichment NES (PXD051911 liver, disease vs control)", y = NULL) +
  theme_masld_compact() +
  theme(strip.text.y = element_text(angle = 0), legend.position = "top",
        axis.text.y = element_text(size = 5))
save_fig(pA, file.path(SUPP, "figS_crossmod_enrichment.pdf"),
         width = 5.4, height = 7.2)

# ── Panel B: cross-modal stage-slope concordance (centerpiece) ────────────────
cc <- fread(file.path(RES, "crossmodal_module_concordance.csv"))
cc[, ct := prettyct(cell_type)]
ctt <- suppressWarnings(cor.test(cc$scrna_beta_fig3g, cc$beta_protein_saf, method = "spearman"))
sc  <- 100 * mean(cc$sign_concordant)
lab_pts <- cc[order(-abs(scrna_beta_fig3g))][
  (!is.na(nes_padj) & nes_padj < 0.05) | seq_len(.N) <= 8]
xr <- max(abs(cc$scrna_beta_fig3g)); yr <- max(abs(cc$beta_protein_saf), na.rm = TRUE)
pB <- ggplot(cc, aes(scrna_beta_fig3g, beta_protein_saf)) +
  annotate("rect", xmin = 0, xmax = xr * 1.1, ymin = 0, ymax = yr * 1.15, fill = masld_colors$up, alpha = 0.05) +
  annotate("rect", xmin = -xr * 1.1, xmax = 0, ymin = -yr * 1.15, ymax = 0, fill = masld_colors$down, alpha = 0.05) +
  geom_hline(yintercept = 0, color = GREY, linewidth = 0.3) +
  geom_vline(xintercept = 0, color = GREY, linewidth = 0.3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.4, linetype = "22") +
  geom_point(aes(color = ct), size = 1.3, alpha = 0.9) +
  scale_color_manual(values = CT_COLS, name = "Cell type", na.value = GREY) +
  labs(x = expression("scRNA module stage slope  " * beta[scRNA] * "  (fig3g: Healthy" %->% "Steatohepatitis)"),
       y = expression("Liver-protein module stage slope  " * beta[protein] * "  (No_MASLD" %->% "MASH)")) +
  annotate("text", x = -xr, y = yr * 1.1, hjust = 0, size = TXT,
           label = sprintf("Spearman rho = %.2f, p = %.3f\nsign-concordant %.0f%% (%d/%d)",
                           ctt$estimate, ctt$p.value, sc, sum(cc$sign_concordant), nrow(cc))) +
  theme_masld_compact() + theme(legend.position = "right")
if (HAS_REPEL) pB <- pB + ggrepel::geom_text_repel(
  data = lab_pts, aes(label = name), size = TXT, color = "black",
  min.segment.length = 0, segment.size = 0.2, max.overlaps = 20, box.padding = 0.3)
# Panel B (centerpiece) is PROMOTED into the main Fig4 validation panel set as figS4m
# (per 2026-07-10 request); the enrichment (A) + cascade (C) stay in the standalone supp block.
FIG4_PANELS <- file.path(FIG4_DIR, "panels"); dir.create(FIG4_PANELS, showWarnings = FALSE, recursive = TRUE)
save_fig(pB, file.path(FIG4_PANELS, "figS4m.pdf"), width = 5.6, height = 4.2)

# ── Panel C: per-module protein cascade across saf stages ─────────────────────
STG <- c("No_MASLD", "MASL", "MASH")
casc <- fread(file.path(RES, "module_cascade_protein.csv"))
casc[, `:=`(saf = factor(saf, levels = STG), ct = prettyct(cell_type),
            dir = factor(ifelse(scrna_dir == "up", "scRNA up-modules", "scRNA down-modules"),
                         levels = c("scRNA up-modules", "scRNA down-modules")))]
pC <- ggplot(casc, aes(saf, mean_score, group = mid)) +
  geom_hline(yintercept = 0, color = GREY, linewidth = 0.3) +
  geom_line(aes(color = ct), linewidth = 0.3, alpha = 0.45) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linewidth = 0.8) +
  scale_color_manual(values = CT_COLS, name = "Cell type") +
  facet_wrap(~ dir, nrow = 1) +
  labs(x = NULL, y = "Module protein score (mean-z)") +
  theme_masld_compact() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(pC, file.path(SUPP, "figS_crossmod_cascade_modules.pdf"), width = 5.4, height = 3.0)

# NOTE: a bulk NMF k=6 program cascade companion was intentionally omitted — the canonical
# program_topgenes integer index does NOT map to the P1..P6 curated codes (program 1 = muscle,
# not P1-Inflammatory-EMT) and the raw loadings carry confounder genes. The Hotspot-module cascade
# (Panel C) is the compartment-matched, well-defined system and is sufficient for this figure.

# ── captions ──────────────────────────────────────────────────────────────────
message(sprintf(paste0(
  "\n[figS_crossmodal_module_proteomics] enrichment(A)+cascade(C) -> %s; ",
  "concordance(B) PROMOTED -> fig4_validation/panels/figS4m.pdf\n",
  "A enrichment: %d disease-significant Hotspot modules; %d with protein NES padj<0.05.\n",
  "B concordance: Spearman rho=%.2f, p=%.3f over %d modules; %.0f%% sign-concordant. ",
  "Quadrants show transcript->protein propagation vs buffering. beta_scRNA is within-cell-type ",
  "(composition-controlled); bulk beta_protein is not, so discordance = buffering OR composition shift.\n",
  "C/C2 cascade: module & NMF-program protein scores across documented No_MASLD->MASL->MASH liver stages."),
  SUPP, nrow(enr), sum(enr$padj < 0.05, na.rm = TRUE),
  ctt$estimate, ctt$p.value, nrow(cc), sc))
