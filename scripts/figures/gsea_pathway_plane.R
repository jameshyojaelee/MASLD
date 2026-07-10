#!/usr/bin/env Rscript
# gsea_pathway_plane.R  — CANDIDATE panel (evidence-plane grammar; N=50 = labeled scatter)
# GSEA Hallmark landscape at the most advanced transition (F4 vs F0) as a pathway
# volcano: x = NES, y = -log10(padj), 50 Hallmark pathways, disease-relevant ones
# labeled. An alternative to the stage-escalation dotplot. Original GSEA dotplot UNTOUCHED.
# NOTE: 50 pathways is small for a "cloud" — this reads as a labeled scatter. A denser
# version would use a bigger collection (C2/GO) or the cross-species human-vs-mouse NES plane.
#
# Output: FIG2_DIR/panels/gsea_hallmark_plane.pdf (exploratory; no Fig-3 letter)

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG2_DIR, "panels"); DATA_DIR <- file.path(PANEL_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
lab_size <- 6 / ggplot2::.pt

d <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_stage_gsea.csv"))
d <- d[contrast == "F4_vs_F0" & !is.na(NES) & !is.na(padj) & padj > 0]
d[, neglogp := -log10(padj)]
d[, dir := fifelse(NES > 0, "up", "down")]
d[, short := sub("HALLMARK_", "", pathway)]

# disease-relevant Hallmark anchors (labeled)
anchors <- c("EPITHELIAL_MESENCHYMAL_TRANSITION","TNFA_SIGNALING_VIA_NFKB","INFLAMMATORY_RESPONSE",
             "IL6_JAK_STAT3_SIGNALING","HYPOXIA","ANGIOGENESIS","INTERFERON_GAMMA_RESPONSE",
             "FATTY_ACID_METABOLISM","BILE_ACID_METABOLISM","OXIDATIVE_PHOSPHORYLATION",
             "XENOBIOTIC_METABOLISM","ADIPOGENESIS")
pretty <- function(x) {
  x <- gsub("_", " ", x); x <- tolower(x)
  x <- gsub("emt|epithelial mesenchymal transition", "EMT", x)
  x <- gsub("tnfa signaling via nfkb", "TNFα/NFκB", x)
  x <- gsub("il6 jak stat3 signaling", "IL6/JAK/STAT3", x)
  x <- gsub("interferon gamma response", "IFNγ response", x)
  x <- gsub("oxidative phosphorylation", "OxPhos", x)
  tools::toTitleCase(x)
}
d[, lab := ifelse(short %in% anchors, pretty(short), NA_character_)]
lab <- d[!is.na(lab)]

COL_UP <- "#C9265E"; COL_DN <- "#1565C0"; COL_BG <- masld_colors$control
xr <- max(abs(d$NES)) * 1.05
p <- ggplot(d, aes(NES, neglogp)) +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "gray85") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, colour = "gray55") +
  geom_point(data = d[is.na(lab)], colour = COL_BG, size = 1.1, alpha = 0.7, shape = 16) +
  geom_point(data = lab, aes(colour = dir), size = 1.8, shape = 16) +
  geom_point(data = lab, colour = "white", size = 1.8, shape = 1, stroke = 0.3) +
  geom_text_repel(data = lab, aes(label = lab), colour = "black", size = lab_size,
    box.padding = 0.5, point.padding = 0.3, segment.size = 0.2, segment.color = "gray60",
    min.segment.length = 0, max.overlaps = Inf, seed = 42, force = 8,
    bg.color = "white", bg.r = 0.12) +
  scale_colour_manual(values = c(up = COL_UP, down = COL_DN), guide = "none") +
  scale_x_continuous(limits = c(-xr, xr)) +
  labs(x = "Enrichment (NES), advanced fibrosis (F4 vs F0)", y = expression(-log[10] * " padj")) +
  theme_masld_compact()

message(sprintf(paste0("CAPTION (GSEA Hallmark landscape): 50 Hallmark pathways at F4-vs-F0 ",
  "(%d significant, padj<0.05); x = NES (>0 up in advanced fibrosis), y = -log10 padj, dashed = ",
  "padj<0.05. Inflammatory/fibrogenic programs (EMT, TNFα/NFκB, inflammatory, IL6/JAK/STAT3, ",
  "hypoxia, angiogenesis) rise; metabolic programs (fatty-acid, bile-acid, OxPhos, xenobiotic) fall. ",
  "N=50 reads as a labeled scatter, not a dense cloud."), sum(d$padj < 0.05)))

out <- file.path(PANEL_DIR, "gsea_hallmark_plane.pdf")
save_fig(p, out, width = fig_half_width + 0.3, height = 3.0)
cat("[GSEA plane] Saved:", out, "\n")
fwrite(d[order(-neglogp)][, .(pathway = short, NES = round(NES,3), padj = signif(padj,3))],
       file.path(DATA_DIR, "gsea_hallmark_plane.csv"))
