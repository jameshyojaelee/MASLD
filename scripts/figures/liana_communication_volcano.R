#!/usr/bin/env Rscript
# liana_communication_volcano.R  — CANDIDATE panel (evidence-plane grammar)
# The FULL LIANA cell-cell communication landscape as a volcano, instead of the
# 3J heatmap of only the ~10 pre-selected headline pairs. Every tested
# (cell-type-pair × ligand-receptor) is a point; the headline pairs anchor it.
#
# x = coarse_Estimate_SH (Steatohepatitis-vs-healthy communication effect from the
#     coarse LMM; SIGN VERIFIED by biology — fibroblast->hepatocyte collagen-integrin
#     pairs [COL4A2->ITGA1/ITGB1] are POSITIVE = strengthened in disease).
# y = -log10(coarse_pval_SH). Marks coloured strengthened/weakened; labels black.
# Original ccc_v3_panels.R / 3J heatmap UNTOUCHED.
#
# Output: FIG4_SC_DIR/panels/figS_liana_communication_volcano.pdf (supplementary panel)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(ggrastr)
})
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG4_SC_DIR, "panels", "supplementary"); DATA_DIR <- file.path(PANEL_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
lab_size <- 6 / ggplot2::.pt

d <- fread(file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3/stage_lr_headline_v3.tsv"),
  select = c("ct_pair", "lr_pair", "coarse_Estimate_SH", "coarse_pval_SH", "headline"))
d <- d[!is.na(coarse_Estimate_SH) & !is.na(coarse_pval_SH) & coarse_pval_SH > 0]
d[, neglogp := -log10(coarse_pval_SH)]
d[, dir := fifelse(coarse_Estimate_SH > 0, "strengthened", "weakened")]
# headline = multi-gate (bootstrap/permutation), NOT the LMM p on the y-axis. Mark whether a
# headline pair also clears the LMM p<0.05 line so open circles below the line are self-explanatory.
d[, lmm_sig := coarse_pval_SH < 0.05]
cat(sprintf("[LIANA volcano] %d points, %d headline (%d LMM-sig, %d bootstrap-only)\n",
    nrow(d), sum(d$headline, na.rm=TRUE), sum(d$headline & d$lmm_sig, na.rm=TRUE),
    sum(d$headline & !d$lmm_sig, na.rm=TRUE)))

# labels for headline pairs: "LIGAND->RECEPTOR (Sender->Receiver)" — the cell-type
# circuit disambiguates L-R pairs that recur across circuits (e.g. C4BPA->BMPR2 is both
# Hep->Endo and Hep->Fib; CDH1->PTPRM is both Hep->Endo and Hep->Hep). Format matches 3J.
ct_abbr <- c("Hepatocytes"="Hep", "Endothelial cells"="Endo", "Fibroblasts"="Fib",
             "Cholangiocytes"="Chol", "Macrophages"="Mac", "T cells"="Tcell",
             "B cells"="Bcell", "Mast cells"="Mast", "Plasma cells"="Plasma")
abbr1 <- function(x) { y <- ct_abbr[x]; ifelse(is.na(y), x, unname(y)) }
lab <- d[headline == TRUE][order(-neglogp)]
lr_parts <- tstrsplit(lab$lr_pair, "__", fixed = TRUE)
lig <- lr_parts[[1]]; rec <- gsub("_", "+", lr_parts[[2]])
ct_parts <- tstrsplit(lab$ct_pair, "->", fixed = TRUE)
send <- abbr1(ct_parts[[1]]); recv <- abbr1(ct_parts[[2]])
lab[, lr := paste0(lig, "→", rec, " (", send, "→", recv, ")")]  # plain text, console print
# plotmath: gene symbols italic, cell-type circuit upright (font rule: italics ONLY on genes)
lab[, plab := paste0('italic("', lig, '")%->%italic("', rec, '")*"  (', send, '"%->%"', recv, ')"')]

COL_UP <- "#C9265E"   # strengthened in disease
COL_DN <- "#1565C0"   # weakened
COL_BG <- masld_colors$control
xr <- max(abs(d$coarse_Estimate_SH), na.rm = TRUE)

p <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "gray85") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, colour = "gray55") +
  rasterize(geom_point(data = d, aes(coarse_Estimate_SH, neglogp),
                       colour = COL_BG, size = 0.4, alpha = 0.3, shape = 16), dpi = 600) +
  # LMM-significant headliners = filled + white halo; bootstrap-only headliners = OPEN circle
  geom_point(data = d[headline == TRUE & lmm_sig == TRUE], aes(coarse_Estimate_SH, neglogp, colour = dir),
             size = 1.7, shape = 16) +
  geom_point(data = d[headline == TRUE & lmm_sig == TRUE], aes(coarse_Estimate_SH, neglogp),
             colour = "white", size = 1.7, shape = 1, stroke = 0.3) +
  geom_point(data = d[headline == TRUE & lmm_sig == FALSE], aes(coarse_Estimate_SH, neglogp, colour = dir),
             size = 1.7, shape = 1, stroke = 0.7) +
  geom_text_repel(data = lab, aes(coarse_Estimate_SH, neglogp, label = plab),
    parse = TRUE, colour = "black", size = lab_size,
    box.padding = 0.5, point.padding = 0.3, segment.size = 0.2, segment.color = "gray60",
    min.segment.length = 0, max.overlaps = Inf, seed = 42, force = 8,
    bg.color = "white", bg.r = 0.12) +
  scale_colour_manual(values = c(strengthened = COL_UP, weakened = COL_DN), guide = "none") +
  scale_x_continuous(limits = c(-xr, xr)) +
  labs(x = "Communication change in disease (Steatohepatitis vs healthy)",
       y = expression(-log[10] * " p")) +
  theme_masld_compact()

message(sprintf(paste0("CAPTION (LIANA communication volcano): all %d tested cell-type-pair x ",
  "ligand-receptor interactions (gray); x = coarse LMM Steatohepatitis-vs-healthy effect (>0 = ",
  "strengthened in disease, e.g. fibroblast->hepatocyte collagen-integrin), y = -log10(LMM p), dashed = ",
  "p<0.05. Headline pairs pass a separate bootstrap/permutation gate (Methods): FILLED = also LMM-",
  "significant (%d); OPEN circles = bootstrap-robust but below the LMM line (%d; the fibroblast->",
  "hepatocyte collagen-integrin pairs). The full landscape lets the reader judge the headline pairs ",
  "in context (cf. the 3J heatmap of only these pairs)."),
  nrow(d), sum(d$headline & d$lmm_sig, na.rm = TRUE), sum(d$headline & !d$lmm_sig, na.rm = TRUE)))

out <- file.path(PANEL_DIR, "figS_liana_communication_volcano.pdf")
save_fig(p, out, width = fig_half_width + 0.4, height = 3.0)
cat("[LIANA volcano] Saved:", out, "\n")
fwrite(d[order(coarse_pval_SH)][, .(ct_pair, lr_pair, coarse_Estimate_SH = round(coarse_Estimate_SH,3),
        coarse_pval_SH = signif(coarse_pval_SH,3), headline)],
       file.path(DATA_DIR, "figS_liana_communication_volcano.csv"))
print(lab[, .(ct_pair, lr, est = round(coarse_Estimate_SH,2), p = signif(coarse_pval_SH,2))])
