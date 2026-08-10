#!/usr/bin/env Rscript
# figS05_scatac_cross_cohort_concordance.R
# ==============================================================================
# Cross-cohort concordance companion to figS05_scatac_disease_vs_control (which is
# the underpowered n=18, 5-control GSE244832-only donor-level DA). Adds the new
# GSE281367 snATAC cohort (Zhu/Hu 2026; 6 MASH/6 NORMAL; ~6x deeper pseudobulk)
# as a SUPPLEMENTARY CORROBORATION cohort (NOT headline replication; the signal is
# one-cohort-dominated). Frame = supplementary corroboration, NOT the inflated
# naive-pooled counts (see Analysis/ATAC/Human_External/LAYER_A_GSE281367_DA_RESULTS.md).
#
# Per-peak cross-cohort logFC scatter, one facet per major cell type:
#   x = GSE244832 log2FC (existing, n=18, 5 control -> underpowered)
#   y = GSE281367 log2FC (new snATAC, n=12, balanced 6/6)
# Honest reading the panel encodes:
#   - Stellate  : 87 peaks significant in BOTH cohorts (magenta), 98% concordant
#                 direction -> genuine independent co-replication.
#   - Hepatocyte: GSE281367 detects 11,547 (blue); GSE244832 detects 0 co-sig
#                 (underpowered) but 83.4% agree in DIRECTION, r=0.49 -> directional
#                 corroboration, NOT co-significance.
#   - Mac/Chol  : underpowered (r=0.03 / 0.21) -> shown for transparency, flagged.
# Each cohort fit SEPARATELY (limma voom-QW; per-cohort dispersion). Peaks deduped
# (24% coordinate-duplicate intervals in the canonical peak set).
# Env: rnaseq.  Output: figures/main/fig4_validation/panels/figS4o.pdf (promoted 2026-07-14;
#   supplementary corroboration panel for Fig4, sits with figS4g/figS4h etc.)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrastr) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
RES  <- file.path(BASE, "Analysis/ATAC/Human_External/results")
MAG  <- "#C9265E"; BLUE <- "#1565C0"; GRAY <- "#9E9E9E"; FDRcut <- 0.05

cts  <- data.table(ct = c("hep","stellate","macrophage","cholangiocyte"),
                   lab = c("Hepatocyte","Stellate","Macrophage","Cholangiocyte"))
summ <- fread(file.path(RES, "replication_summary.csv"))
summ[, lab := cts$lab[match(cell_type, cts$ct)]]

d <- rbindlist(lapply(seq_len(nrow(cts)), function(i){
  f <- file.path(RES, sprintf("replication_peaks_%s.csv", cts$ct[i]))
  if (!file.exists(f)) return(NULL)
  x <- fread(f); x[, ct := cts$lab[i]]; x
}))
stopifnot(nrow(d) > 0)

# 3-category palette (Liang magenta/blue/gray): the replication story is
# "both" vs "new-cohort-only"; 244-only + n.s. collapse to gray.
d[, cat := "Other / n.s."]
d[FDR_281 < FDRcut, cat := "GSE281367 only (new)"]
d[FDR_281 < FDRcut & FDR_244 < FDRcut, cat := "Significant in both"]
d[, cat := factor(cat, levels = c("Significant in both","GSE281367 only (new)","Other / n.s."))]
d[, ct := factor(ct, levels = cts$lab)]

# downsample the gray cloud per facet for render size; keep ALL significant points
sig <- d[cat != "Other / n.s."]
ns  <- d[cat == "Other / n.s."][, .SD[if (.N > 6000) sample(.N, 6000) else seq_len(.N)], by = ct]
dd  <- rbind(sig, ns)

ann <- summ[, .(ct = factor(lab, levels = cts$lab),
  txt = sprintf("281 sig: %s\nco-sig: %d\ndir concord: %.0f%%\nr = %.2f%s",
    format(n_sig_GSE281367_alone, big.mark = ","), n_sig_both_replicated,
    pct_281sig_dir_concordant_in_244, logFC_r_genomewide,
    ifelse(reportable == "YES", "", "  (underpowered)")))]

cols <- c("Significant in both" = MAG, "GSE281367 only (new)" = BLUE, "Other / n.s." = GRAY)
p <- ggplot(dd, aes(logFC_244, logFC_281)) +
  geom_hline(yintercept = 0, color = GRAY, linewidth = 0.2) +
  geom_vline(xintercept = 0, color = GRAY, linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = GRAY, linewidth = 0.2) +
  rasterise(geom_point(aes(color = cat), size = 0.22, alpha = 0.5, stroke = 0), dpi = 400) +
  geom_text(data = ann, aes(label = txt), x = -Inf, y = Inf, hjust = -0.06, vjust = 1.08,
            size = 1.75, lineheight = 0.92, color = "black") +
  facet_wrap(~ ct, nrow = 1, scales = "free") +
  scale_color_manual(values = cols, name = NULL, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 1.2, alpha = 1))) +
  labs(x = "GSE244832 log2FC  (existing multiome, n=18, 5 control)",
       y = "GSE281367 log2FC\n(new snATAC, n=12, 6/6)") +
  theme_masld(6) +
  theme(legend.position = "bottom", legend.key.size = unit(0.3, "lines"),
        panel.spacing = unit(0.5, "lines"))

OUT <- file.path(FIG4_DIR, "panels", "figS4o.pdf")   # promoted from figS05 (2026-07-14)
ggsave(OUT, p, width = 7.2, height = 2.5, useDingbats = FALSE)
message(sprintf("[wrote] %s", OUT))
message("[caption] Cross-cohort concordance of donor-level MASLD-vs-control differential ",
        "accessibility across cohorts. Each peak's log2 fold-change in the existing ",
        "GSE244832 multiome cohort (x; n=18, 5 controls) vs the new GSE281367 snATAC ",
        "cohort (y; n=12, balanced). Both cohorts fit separately (limma voom quality-",
        "weights; per-cohort dispersion). Stellate cells co-replicate (87 peaks ",
        "significant in both, magenta; 98% concordant direction). Hepatocyte differential ",
        "accessibility is detected in the deeper GSE281367 cohort (11,547 peaks) and ",
        "directionally corroborated by the underpowered GSE244832 cohort (83% same ",
        "direction, r=0.49) but does not reach co-significance (0 peaks). Macrophage and ",
        "cholangiocyte are underpowered (few testable peaks; r=0.03/0.21) and shown for ",
        "transparency only. Naive count-pooling counts are not used (dispersion-inflated). ",
        "This is a supplementary corroboration layer, not headline independent replication: ",
        "the signal is carried mainly by the deeper GSE281367 cohort, with GSE244832 ",
        "providing directional (largely sub-significant) support.")
