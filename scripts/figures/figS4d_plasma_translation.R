#!/usr/bin/env Rscript
# ==============================================================================
# Fig S4d (DEMOTED from main 4d 2026-07-08) — LIVER -> PLASMA secreted-target translation (biomarker hook).
#
# Fig 4 spine = "physical corroboration": prioritized MASLD targets are
# corroborated across the physical compartments a therapeutic/biomarker engages.
# This panel follows one set of secreted ligands down the compartment ladder:
#     tissue mRNA (bulk RNA-seq)  ->  liver protein (DIA-MS)  ->  plasma protein (DIA-MS)
# to ask which prioritized targets leak a measurable, concordant plasma signal
# (i.e. are detectable non-invasively).
#
# SET = 18 CCC-inferred, MASLD-up, plasma-panel-detectable secreted ligands
#   (RNA-seq/results/secretome_chain/secretome_chain_triple_concordant_MASLD_up.csv).
#
# HONESTY (see caption): the CSV's `in_olink` is merely Olink-panel DETECTABILITY
# and `tissue_plasma_concordant` is tissue-scRNA x bulk-mRNA direction concordance
# (NO plasma-level DE was ever computed there). So the ONLY genuine plasma quantity
# is the independent DIA-MS plasma proteome (PXD052937), which has just 7 controls
# and measures 6/18 of these ligands, 0 of them significant. The primary,
# genuinely-significant axis is therefore the TISSUE mRNA log2FC (LEFT); the
# secreted-protein arm (RIGHT) is a directional TREND, not powered DE.
#
# NOT a CCC panel (LIANA is provenance only, kept out of Fig4's spine) and NOT a
# tile/heatmap matrix. Row-aligned lollipops.
#
# Output: figures/main/fig4_validation/panels/figS4d.pdf (supp; kept in panels/ per figS4 convention)
# Env: rnaseq
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── the 18 triple-concordant secreted ligands ────────────────────────────────
sec <- fread(file.path(BASE, "RNA-seq/results/secretome_chain/secretome_chain_triple_concordant_MASLD_up.csv"))
sec <- sec[, .(gene = ligand_complex, source_ct = primary_source_ct,
               bulk_lfc, bulk_padj)]

# ── liver protein (DIA-MS PXD051911) + plasma protein (DIA-MS PXD052937) ───────
con <- fread(file.path(BASE, "Analysis/Proteomics/results/protein_transcript_concordance_v3.csv"))
liv <- unique(con[dataset == "PXD051911", .(gene, liver_lfc = protein_logFC, liver_padj = protein_padj)], by = "gene")
pla <- fread(file.path(BASE, "Analysis/Proteomics/results/pxd052937_dep_volcano.csv"))
pla <- unique(pla[, .(gene, plasma_lfc = logFC, plasma_padj = padj)], by = "gene")

d <- Reduce(function(a, b) merge(a, b, by = "gene", all.x = TRUE), list(sec, liv, pla))
setorder(d, bulk_lfc)                                   # ascending -> highest bulk_lfc at top
gene_levels <- d$gene; ng <- nrow(d)
gi <- function(g) as.integer(factor(g, levels = gene_levels))
efy <- which(gene_levels == "EFEMP1")                   # exemplar: coherent across all 3 compartments

n_liver  <- sum(!is.na(d$liver_lfc))
n_plasma <- sum(!is.na(d$plasma_lfc))
n_plasma_sig <- sum(d$plasma_padj < 0.05, na.rm = TRUE)
cat(sprintf("[plasmaS4] %d ligands; liver measured %d; plasma measured %d (sig %d); EFEMP1 row=%d\n",
            ng, n_liver, n_plasma, n_plasma_sig, efy))

# shared log2FC scale across both compartments so attenuation is visible
XLO <- -0.35; XHI <- 1.25
tissue_col <- masld_colors$up                           # #C9265E magenta (disease-up transcript)
liver_col  <- "#2166AC"                                 # blue (liver protein)
plasma_col <- "#2CA25F"                                 # green (plasma protein — biomarker endpoint)

# ── LEFT: TISSUE mRNA log2FC (PRIMARY, genuinely significant) ─────────────────
dL <- d[, .(gene, lfc = bulk_lfc, sig = bulk_padj < 0.05, y = gi(gene))]
pL <- ggplot(dL) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_segment(aes(x = 0, xend = lfc, y = y, yend = y), linewidth = 0.4, colour = tissue_col) +
  geom_point(aes(x = lfc, y = y, shape = sig), colour = tissue_col, fill = "white", size = 1.5, stroke = 0.5) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), breaks = c(TRUE, FALSE),
                     labels = c("padj<0.05", "not significant"), name = NULL) +
  scale_y_continuous(limits = c(0.4, ng + 0.6), breaks = seq_len(ng), labels = gene_levels, expand = c(0, 0)) +
  scale_x_continuous(limits = c(XLO, XHI), breaks = c(0, 0.5, 1.0), expand = c(0, 0)) +
  labs(x = expression(log[2]~"FC"), y = NULL) +
  ggtitle("Tissue mRNA") +
  theme_masld() +
  theme(plot.title = element_text(size = 6, hjust = 0, face = "plain"),
        axis.text.y = element_text(face = "italic", size = 6, colour = "black"),
        legend.position = "bottom", legend.key.size = unit(0.3, "lines"),
        legend.text = element_text(size = 6), legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(3, "pt"), panel.grid = element_blank())

# ── RIGHT: SECRETED protein log2FC (liver + plasma DIA-MS; underpowered TREND) ─
rr <- rbindlist(list(
  d[!is.na(liver_lfc),  .(gene, compartment = "Liver protein",  lfc = liver_lfc,  padj = liver_padj,  yo =  0.18)],
  d[!is.na(plasma_lfc), .(gene, compartment = "Plasma protein", lfc = plasma_lfc, padj = plasma_padj, yo = -0.18)]))
rr[, compartment := factor(compartment, levels = c("Liver protein", "Plasma protein"))]
rr[, sig := padj < 0.05]                                # all FALSE -> every point open (underpowered)
rr[, y := gi(gene) + yo]
palR <- c("Liver protein" = liver_col, "Plasma protein" = plasma_col)
pR <- ggplot(rr) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_segment(aes(x = 0, xend = lfc, y = y, yend = y, colour = compartment), linewidth = 0.4) +
  geom_point(aes(x = lfc, y = y, colour = compartment, shape = sig), fill = "white", size = 1.5, stroke = 0.5) +
  scale_colour_manual(values = palR, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_y_continuous(limits = c(0.4, ng + 0.6), breaks = seq_len(ng), labels = NULL, expand = c(0, 0)) +
  scale_x_continuous(limits = c(XLO, XHI), breaks = c(0, 0.5, 1.0), expand = c(0, 0)) +
  labs(x = expression(log[2]~"FC"), y = NULL) +
  ggtitle("Secreted protein (DIA-MS)") +
  theme_masld() +
  theme(plot.title = element_text(size = 6, hjust = 0, face = "plain"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.key.size = unit(0.3, "lines"),
        legend.text = element_text(size = 6), legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(3, "pt"), panel.grid = element_blank())

comp <- pL + pR + plot_layout(widths = c(1, 0.92))
out <- file.path(FIG4_DIR, "panels", "figS4d.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
ggsave(out, comp, width = 4.1, height = 2.9, device = grDevices::cairo_pdf)
cat("[plasmaS4] saved:", out, "\n")

ef_bulk   <- d[gene == "EFEMP1", bulk_lfc]
ef_liver  <- d[gene == "EFEMP1", liver_lfc]
ef_plasma <- d[gene == "EFEMP1", plasma_lfc]
ef_ppadj  <- d[gene == "EFEMP1", plasma_padj]
message(sprintf(paste0(
  "CAPTION (Fig S4d, plasma translation; demoted from main 4d): Liver->plasma translation of prioritized secreted targets (biomarker hook). ",
  "%d MASLD-up secreted ligands (CCC-inferred, all detectable on the Olink plasma panel) followed down the ",
  "compartment ladder tissue mRNA -> liver protein -> plasma protein; ordered by tissue log2FC. ",
  "LEFT (PRIMARY): bulk RNA-seq log2FC (disease vs control, n=846); filled = padj<0.05, open = trend. ",
  "RIGHT: the secreted-protein readout from two independent DIA-MS proteomes - liver (PXD051911; %d/%d ligands measured) ",
  "and plasma (PXD052937; %d/%d measured). EVERY protein point is OPEN because the plasma proteome rests on only 7 ",
  "controls and reaches significance for 0/%d measured ligands - this is a directional translation TREND, not powered ",
  "differential expression. EFEMP1 is the exemplar that stays concordantly UP across all three compartments ",
  "(mRNA +%.2f, liver protein +%.2f, plasma protein +%.2f, plasma padj=%.2f, n.s.). This is a compartment-translation ",
  "claim (secreted protein reaching circulation), not a cell-cell-communication claim; LIANA/CCC is provenance only and ",
  "is deliberately NOT an axis here (the RNA-only CCC claim lives in Fig 3). Note: the source CSV's `tissue_plasma_concordant` ",
  "flag is tissue-scRNA x bulk-mRNA direction concordance gated on plasma-panel detectability - it carries NO plasma-level ",
  "quantity, so the DIA-MS plasma arm above is the only genuine circulating measurement."),
  ng, n_liver, ng, n_plasma, ng, n_plasma, ef_bulk, ef_liver, ef_plasma, ef_ppadj))
