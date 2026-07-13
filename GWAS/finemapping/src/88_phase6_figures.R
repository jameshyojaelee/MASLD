#!/usr/bin/env Rscript
# =============================================================================
# 88_phase6_figures.R — Phase-6 ToS-CLEAN seqfunc panels (house style)
# =============================================================================
# ADDITIVE. Zero-shot, open-license direction panels from src/84 (caQTL) + src/85
# (eQTL); Decima cell-type (src/76b); credible-set PIP (src/89). Supersedes the
# retired AlphaGenome-TRAINED supervised panels (ToS).
#
# House style: 6pt, theme_masld, control gray #9E9E9E, magenta = signal, captions
# via message() (never baked into art), individual PDFs, plain-language labels.
# AlphaGenome bars = ZERO-SHOT ANNOTATION only (never trained on; non-commercial).
# APPLY-ONLY firewall: annotation/benchmark only, never a scored convergence channel.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

SF      <- file.path(BASE, "GWAS/finemapping/results/seqfunc")
OUT_DIR <- file.path(BASE, "figures/supplementary/seqfunc")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
rendered <- character(0)
GRAY <- "#9E9E9E"

grpf <- function(lic) fifelse(lic %in% c("open", "open_zeroshot"), "Open model",
                       fifelse(lic == "baseline", "Baseline", "AlphaGenome (annotation)"))

# =============================================================================
# PANEL — direction_zeroshot.pdf   (accessibility vs expression, one figure)
# =============================================================================
# The KEY split in one view: sequence-to-function models call the DIRECTION of a variant's
# effect on chromatin ACCESSIBILITY (caQTL, top) but not on gene EXPRESSION
# (eQTL, bottom). Two facets, shared accuracy axis, 0.5 = chance.
# -----------------------------------------------------------------------------
ca <- fread(file.path(SF, "caqtl_zeroshot_permodel.tsv"))[level == "signal_level_clumped"]
disp_ca <- c(borzoi_atac = "Borzoi (ATAC)", borzoi_dnase = "Borzoi (DNase)",
             chrombpnet_hepg2 = "ChromBPNet (ATAC)", alphagenome_atac = "AlphaGenome (ATAC)",
             alphagenome_dnase = "AlphaGenome (DNase)",
             positional_MAF_peakdist_logistic = "MAF + peak distance")
ca <- ca[model %in% names(disp_ca)]
ca[, `:=`(label = disp_ca[model], modality = "Accessibility (caQTL)", grp = grpf(license))]

eq <- fread(file.path(SF, "eqtl_signal_diagnose_permodel.tsv"))[level == "signal_level_one_per_signal_id"]
disp_eq <- c(borzoi_logsed_liver = "Borzoi (liver)", borzoi_logsed_allliver = "Borzoi (all-liver)",
             borzoi_acc_delta = "Borzoi (accessibility)", alphagenome_gene_lfc = "AlphaGenome",
             alphagenome_gene_lfc_alltracks = "AlphaGenome (all-tracks)",
             positional_MAF_tssdist_logistic = "MAF + TSS distance")
eq <- eq[model %in% names(disp_eq)]
eq[, `:=`(label = disp_eq[model], modality = "Expression (eQTL)", grp = grpf(license))]

D <- rbind(ca[, .(label, auroc, ci_lo, ci_hi, grp, modality)],
           eq[, .(label, auroc, ci_lo, ci_hi, grp, modality)])
D[, modality := factor(modality, levels = c("Accessibility (caQTL)", "Expression (eQTL)"))]
setorder(D, modality, auroc)
D[, label := factor(label, levels = unique(label))]
D[, grp := factor(grp, levels = c("Open model", "AlphaGenome (annotation)", "Baseline"))]

pDir <- ggplot(D, aes(x = auroc, y = label, color = grp, shape = grp)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3, color = GRAY) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2, linewidth = 0.35) +
  geom_point(size = 1.7) +
  facet_grid(modality ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_color_manual(values = c("Open model" = masld_colors$up,
                                "AlphaGenome (annotation)" = GRAY, "Baseline" = GRAY), name = NULL) +
  scale_shape_manual(values = c("Open model" = 16, "AlphaGenome (annotation)" = 17,
                                "Baseline" = 15), name = NULL) +
  scale_x_continuous(limits = c(0.3, 0.82), breaks = seq(0.3, 0.8, 0.1)) +
  labs(x = "Direction-prediction accuracy (auROC)", y = NULL,
       title = "Effect direction: accessibility is learnable, expression is not") +
  theme_masld() +
  theme(legend.position = "bottom", legend.margin = margin(t = -2),
        strip.placement = "outside", panel.spacing = unit(0.15, "cm"))

message("[direction_zeroshot] CAPTION: Zero-shot sign-concordance accuracy (auROC) for the ",
        "DIRECTION of a variant's effect on chromatin accessibility (caQTL, top; Currin 2025 ",
        "bulk liver n=138, 253-signal panel) and gene expression (eQTL, bottom; Broadaway ",
        "liver-microarray, 169 signals). Dot = auROC, whisker = block-bootstrap 95% CI, dashed ",
        "line = 0.5 chance. Open-license models (magenta) call ACCESSIBILITY direction (best ",
        "Borzoi ATAC 0.749 [0.720,0.777], +0.23 over a MAF+distance baseline) but are ",
        "SUB-THRESHOLD for EXPRESSION (best 0.592; the 0.65 pre-set bar not met). AlphaGenome's ",
        "own eQTL scorer is also at chance (0.545) — so the expression result is PANEL ",
        "ASCERTAINMENT, not model failure. AlphaGenome (triangles) = zero-shot annotation only. ",
        "The two facets are NOT a paired comparison (different assay/ascertainment/n).")

save_fig(pDir, file.path(OUT_DIR, "direction_zeroshot.pdf"), width = 4.4, height = 3.6)
rendered <- c(rendered, "direction_zeroshot.pdf")

# =============================================================================
# PANEL — decima_celltype.pdf   (top cell type per risk gene — an annotation)
# =============================================================================
dec <- fread(file.path(SF, "decima_celltype.tsv"))
d5 <- dec[, .N, by = argmax_celltype][order(-N)]
pretty_ct <- c(hepatocyte = "Hepatocyte", kupffer_macrophage = "Kupffer / macrophage",
               lymphoid_T_NK_B = "Lymphoid (T/NK/B)", LSEC_endothelial = "LSEC / endothelial",
               stellate_fibroblast = "Stellate / fibroblast", cholangiocyte = "Cholangiocyte")
d5[, pretty := pretty_ct[argmax_celltype]]
d5[, is_hep := argmax_celltype == "hepatocyte"]
d5 <- d5[order(N)]
d5[, pretty := factor(pretty, levels = pretty)]

pC <- ggplot(d5, aes(x = N, y = pretty, fill = is_hep)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), hjust = -0.3, size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_manual(values = c("TRUE" = masld_colors$up, "FALSE" = GRAY), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Number of risk genes", y = NULL,
       title = "Predicted top cell type per risk gene") +
  theme_masld()

message("[decima_celltype] CAPTION: For each nominated risk gene, the liver cell type where ",
        "Decima predicts its effect is strongest (n=120). Read as a PER-GENE ANNOTATION, NOT a ",
        "compartment proportion: the aggregate split is not robust (a balanced 44-56% hepatocyte ",
        "among confident calls). The defensible evidence is the next panel + the SORT1 positive ",
        "control (-> hepatocyte). APPLY-ONLY annotation; never a scored convergence channel.")

save_fig(pC, file.path(OUT_DIR, "decima_celltype.pdf"), width = 3.4, height = 2.3)
rendered <- c(rendered, "decima_celltype.pdf")

# =============================================================================
# PANEL — decima_snatac_validation.pdf   (does the cell-type call match ATAC?)
# =============================================================================
snat <- data.table(
  grp = factor(c("Observed", "Expected by chance"), levels = c("Expected by chance", "Observed")),
  val = c(19, 14.2), lab = c("19 of 29", "~14 of 29"))
pC2 <- ggplot(snat, aes(x = val, y = grp, fill = grp)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = lab), hjust = -0.12, size = GEOM_TEXT_6PT, color = "black") +
  annotate("text", x = 19, y = 2.45, label = "permutation p = 0.017",
           size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_manual(values = c("Observed" = masld_colors$up, "Expected by chance" = GRAY),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 25), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Risk genes whose predicted cell type is accessible", y = NULL,
       title = "Predicted cell type matches measured accessibility") +
  theme_masld() +
  coord_cartesian(ylim = c(0.5, 2.6), clip = "off")

message("[decima_snatac_validation] CAPTION: The load-bearing Decima check. Of 29 nominated ",
        "variants overlapping >=1 snATAC peak, the predicted (argmax) cell type is one of the ",
        "accessible cell types for 19 (66%) — significantly above a peak-multiplicity-preserving ",
        "permutation null (~14 expected; perm-p 0.017 uniform / 0.025 conservative). SORT1 -> ",
        "hepatocyte is the positive control. A semi-independent sanity check (same cohort), NOT ",
        "disjoint validation and NOT a compartment proportion.")

save_fig(pC2, file.path(OUT_DIR, "decima_snatac_validation.pdf"), width = 3.4, height = 1.7)
rendered <- c(rendered, "decima_snatac_validation.pdf")

# =============================================================================
# PANEL — pip_weighted_mechanism.pdf   (how well each signal resolves to 1 variant)
# =============================================================================
pw <- fread(file.path(SF, "pip_weighted_mechanism_per_cs.tsv"))
pw[, conc := factor(concentration, levels = c("concentrated", "moderate", "diffuse"),
                    labels = c("Well-resolved\n(≈1 variant)", "Moderate",
                               "Diffuse\n(many variants)"))]
d6 <- pw[, .N, by = conc][order(conc)]

pD <- ggplot(d6, aes(x = conc, y = N, fill = conc)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), vjust = -0.4, size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_manual(values = c("Well-resolved\n(≈1 variant)" = "#4D4D4D",
                               "Moderate" = "#9E9E9E",
                               "Diffuse\n(many variants)" = "#D9D9D9"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Fine-mapped signals",
       title = "Most nominations resolve to a single variant") +
  theme_masld()

message("[pip_weighted_mechanism] CAPTION: The 81 nominations map to 253 fine-mapped credible ",
        "sets (a lead recurs across GWAS/trait signals). Within each (posterior PIP sums to 1), ",
        "the mechanism class is PIP-weighted over ALL member variants. 159/253 signals resolve to ",
        "essentially one variant; every signal carries >=50% posterior mass on one mechanism, and ",
        "proper PIP-weighting changes the coarse mechanism call in only 5/253 (2%) — one coding-> ",
        "regulatory correction. So the coarse mechanism is ROBUST to fine-mapping uncertainty; the ",
        "residual uncertainty is which variant in the diffuse tail. No global architecture claim ",
        "(manuscript 15.9% coverage gate). APPLY-ONLY annotation.")

save_fig(pD, file.path(OUT_DIR, "pip_weighted_mechanism.pdf"), width = 3.0, height = 2.4)
rendered <- c(rendered, "pip_weighted_mechanism.pdf")

message("\n[88_phase6_figures] rendered: ", paste(rendered, collapse = ", "))
