#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# figS5 (SUPPLEMENTARY, MOUSE) — MOCK-UP for review.
# "An approved MASH drug (semaglutide) reverses the disease transcriptome in
#  mouse models" — presented as THERAPEUTIC-TRACTABILITY CONTEXT, NOT atlas
# validation (the reversal is generic, reach-test null; and the atlas paper is
# human-only). Included only if the human-only rule is waived for a labelled
# mouse supplement.
#
#   a  reversal scatter: mouse disease meta log2FC vs semaglutide meta log2FC
#      (disease-significant orthologs); anti-diagonal = reversal
#   b  per-dataset reversal correlation (4 datasets, DESeq2 + limma) — reproducible;
#      weakest in the weight-stable CDA-HFD model
#   c  weight-INDEPENDENT reversal (CDA-HFD, GSE294630): fibrosis collagens/TIMPs
#      down, lipid genes null (Jara Nat Med 2025-concordant)
#
# Style: 6pt Helvetica plain, text black, gene symbols italic, NO lollipop, PDF.
# Data: RNA-seq/results/glp1ra/mouse_reversal/{mouse_semaglutide_reversal.csv,
#       mouse_semaglutide_reversal_meta.csv}. Env: rnaseq.
# Output: figures/supplementary/figS_glp1ra/figS5_mouse_semaglutide_reversal.pdf
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
OUT <- FIGS_GLP1RA_DIR; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
MR  <- file.path(BASE, "RNA-seq/results/glp1ra/mouse_reversal")
GEOM_TXT <- 6 / ggplot2::.pt
UP <- "#C9265E"; DOWN <- "#1565C0"; CTRL <- "#9E9E9E"
tag_theme <- theme(plot.tag = element_text(size = 6, face = "plain"))

# ── a — reversal scatter (disease vs semaglutide meta log2FC) ─────────────────
meta <- fread(file.path(MR, "mouse_semaglutide_reversal_meta.csv"))
sig  <- meta[!is.na(disease_meta_lfc) & !is.na(treatment_meta_lfc) & disease_meta_padj < 0.05]
rr   <- cor(sig$disease_meta_lfc, sig$treatment_meta_lfc, use = "complete.obs")
pctrev <- 100 * mean(sign(sig$disease_meta_lfc) != sign(sig$treatment_meta_lfc))
lab <- sig[order(-abs(disease_meta_lfc))][sign(disease_meta_lfc) != sign(treatment_meta_lfc)][1:8]
xlim <- quantile(sig$disease_meta_lfc, c(0.003, 0.997)); ylim <- quantile(sig$treatment_meta_lfc, c(0.01, 0.99))
pA <- ggplot(sig, aes(disease_meta_lfc, treatment_meta_lfc)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = DOWN, alpha = 0.05) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = DOWN, alpha = 0.05) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.25) +
  rasterize_layer(geom_point(colour = "grey55", size = 0.25, alpha = 0.5)) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.4) +
  geom_point(data = lab, colour = UP, size = 0.8) +
  geom_text_repel(data = lab, aes(label = human_symbol), size = GEOM_TXT, fontface = "italic",
                  segment.size = 0.2, min.segment.length = 0, max.overlaps = 20, colour = "black") +
  annotate("text", x = xlim[1], y = ylim[1], hjust = 0, vjust = 0, size = GEOM_TXT, colour = "black",
           label = sprintf("r = %.2f\n%.0f%% reversed", rr, pctrev)) +
  coord_cartesian(xlim = xlim, ylim = ylim) +
  labs(x = "disease log2FC (MASH vs control)", y = "semaglutide log2FC (treated vs MASH)") +
  theme_masld()

# ── b — per-dataset reversal correlation ─────────────────────────────────────
long <- fread(file.path(MR, "mouse_semaglutide_reversal.csv"))
perds <- long[!is.na(treatment_lfc) & !is.na(disease_lfc),
              .(r = cor(disease_lfc, treatment_lfc, use = "complete.obs"),
                n = .N, model_type = model_type[1]), by = dataset]
perds[, dataset := factor(dataset, levels = dataset[order(r)])]
perds[, mlab := fifelse(model_type == "weight_independent", "weight-independent\n(CDA-HFD)", "metabolic")]
pB <- ggplot(perds, aes(r, dataset, fill = mlab)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%.2f", r)), hjust = 1.2, size = GEOM_TXT, colour = "white") +
  scale_fill_manual(values = c("metabolic" = UP, "weight-independent\n(CDA-HFD)" = DOWN), name = NULL) +
  scale_x_continuous(limits = c(min(perds$r) * 1.15, 0)) +
  labs(x = "reversal correlation (disease vs semaglutide log2FC)", y = NULL) +
  theme_masld() +
  theme(legend.position = "bottom", legend.margin = margin(t = -4),
        legend.key.height = unit(0.22, "cm"))

# ── c — weight-independent CDA-HFD: fibrosis down, lipid null ─────────────────
cda <- long[dataset == "GSE294630"]
fib <- c("COL1A1", "COL1A2", "COL3A1", "TIMP1", "TIMP2", "MMP13")
lip <- c("FASN", "ACACA", "SREBF1", "PPARG")
cc <- cda[human_symbol %in% c(fib, lip),
          .(human_symbol, treatment_lfc, treatment_padj)]
cc[, cls := fifelse(human_symbol %in% fib, "fibrosis", "lipid")]
cc[, human_symbol := factor(human_symbol, levels = rev(c(fib, lip)))]
cc[, star := fifelse(treatment_padj < 0.05, "*", "")]
pC <- ggplot(cc, aes(treatment_lfc, human_symbol, fill = cls)) +
  geom_col(width = 0.66) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.25) +
  geom_text(aes(label = star, x = treatment_lfc - 0.03), hjust = 1, size = 6 / ggplot2::.pt * 1.4,
            colour = "black", vjust = 0.78) +
  scale_fill_manual(values = c("fibrosis" = DOWN, "lipid" = CTRL), name = NULL) +
  labs(x = "semaglutide log2FC (CDA-HFD)", y = NULL) +
  theme_masld() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.position = "bottom", legend.margin = margin(t = -4),
        legend.key.height = unit(0.22, "cm"))

f <- (pA | pB | pC) + plot_layout(widths = c(1.15, 1.1, 1)) +
  plot_annotation(tag_levels = "a") & tag_theme
ggsave(file.path(OUT, "figS5_mouse_semaglutide_reversal.pdf"), f,
       width = 7.4, height = 2.6, device = pdf_device)

message("figS5 mouse reversal (MOCK-UP) written: ",
        file.path(OUT, "figS5_mouse_semaglutide_reversal.pdf"))
message(sprintf("  panel a: r=%.2f among %d disease-significant orthologs; %.0f%% reversed", rr, nrow(sig), pctrev))
print(perds[, .(dataset, r = round(r, 2), n, model_type)])
message("CAPTION (MOUSE SUPPLEMENT, therapeutic-tractability context — NOT atlas validation): ",
        "In mouse MASH liver, semaglutide reverses the disease transcriptome. (a) Disease-significant ",
        "orthologs anti-correlate between the disease and semaglutide contrasts (Pearson r above). ",
        "(b) Reversal reproduces across 4 datasets and 2 methods (DESeq2 + limma), weakest in the ",
        "weight-stable CDA-HFD model. (c) In CDA-HFD (weight-independent), semaglutide suppresses ",
        "fibrosis collagens/TIMPs (* BH-FDR<0.05) with a null lipid response (Jara Nat Med 2025-concordant). ",
        "This reversal is not specific to atlas-nominated targets and is shown as therapeutic context only.")
