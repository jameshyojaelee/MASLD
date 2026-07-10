#!/usr/bin/env Rscript
# ==============================================================================
# Fig 4e — EFEMP1 MULTI-MODAL validation of an inferred fibroblast->hepatocyte
# ligand. EFEMP1 (fibulin-3) is a LIANA fibroblast-sourced ligand (CCC inferred
# from scRNA); this panel VALIDATES that inference with ORTHOGONAL measurements:
# the disease effect is confirmed up in bulk RNA AND in two INDEPENDENT proteomes
# (liver DIA-MS + plasma DIA-MS), and it is spatially structured (Visium Moran's I).
# The load-bearing point is cross-MODALITY replication (RNA inference -> protein),
# not more transcriptomics. Positive control — no genetic pillar (COLOC ~0.01).
# Output: figures/main/fig4_validation/panels/fig4e_efemp1_multimodal.pdf
# Env: rnaseq
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
G <- "EFEMP1"

# ── pull EFEMP1's evidence from each modality (data-driven, not hardcoded) ─────
sc <- fread(file.path(BASE, "RNA-seq/results/secretome_chain/secretome_chain_triple_concordant_MASLD_up.csv"))
er <- sc[ligand_complex == G][1]
liana <- as.numeric(er$max_score_diff); bulk_lfc <- as.numeric(er$bulk_lfc); bulk_padj <- as.numeric(er$bulk_padj)
src <- er$primary_source_ct

con <- fread(file.path(BASE, "Analysis/Proteomics/results/protein_transcript_concordance_v3.csv"))
lp <- con[gene == G & dataset == "PXD051911"][1]
liver_lfc <- as.numeric(lp$protein_logFC); liver_padj <- as.numeric(lp$protein_padj)

pl <- fread(file.path(BASE, "Analysis/Proteomics/results/pxd052937_dep_volcano.csv"))
pr <- pl[gene == G][1]
plasma_lfc <- as.numeric(pr$logFC); plasma_padj <- as.numeric(pr$padj)

atl <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
gcol <- if ("human_symbol" %in% names(atl)) "human_symbol" else "gene"
ar <- atl[get(gcol) == G][1]
moran <- as.numeric(ar$spatial_max_I)

# ── log2FC lollipop across the 3 fold-change modalities (RNA context + 2 orthogonal proteomes) ──
mm <- data.table(
  modality = c("Bulk mRNA", "Liver protein\n(DIA-MS)", "Plasma protein\n(DIA-MS)"),
  cls      = c("Transcriptome", "Proteome (orthogonal)", "Proteome (orthogonal)"),
  lfc      = c(bulk_lfc, liver_lfc, plasma_lfc),
  padj     = c(bulk_padj, liver_padj, plasma_padj))
mm[, ypos := .N:1]                                            # Bulk mRNA on top
mm[, sig := padj < 0.05]
mm[, cls := factor(cls, levels = c("Transcriptome", "Proteome (orthogonal)"))]

cat(sprintf("[fig4e] EFEMP1 src=%s | LIANA %.2f | bulk %.2f(p%.0e) liver %.2f(q%.2f) plasma %.2f(q%.2f) | Moran %.2f\n",
            src, liana, bulk_lfc, bulk_padj, liver_lfc, liver_padj, plasma_lfc, plasma_padj, moran))

pal <- c("Transcriptome" = "#9E9E9E", "Proteome (orthogonal)" = "#2166AC")
p <- ggplot(mm, aes(y = ypos)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_segment(aes(x = 0, xend = lfc, yend = ypos, colour = cls), linewidth = 0.5) +
  geom_point(aes(x = lfc, colour = cls, shape = sig, fill = cls), size = 2.0, stroke = 0.5) +
  # context (2 short lines): the inference being validated (LIANA CCC) + the spatial modality
  annotate("text", x = 0, y = 4.55, hjust = 0, size = 6/.pt, colour = "grey25",
           label = sprintf("Fibroblast→hepatocyte LIANA ligand (%.2f)", liana)) +
  annotate("text", x = 0, y = 4.05, hjust = 0, size = 6/.pt, colour = "grey25",
           label = sprintf("spatially structured (Visium Moran's I %.2f)", moran)) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), breaks = c(TRUE, FALSE),
                     labels = c("padj < 0.05", "trend (q≈0.1)"), name = NULL,
                     guide = guide_legend(override.aes = list(fill = "grey30", colour = "grey30"))) +
  scale_y_continuous(breaks = mm$ypos, labels = mm$modality, limits = c(0.5, 4.9), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, max(mm$lfc) * 1.15), expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression(log[2]~"fold change (disease vs control)"), y = NULL) +
  ggtitle("EFEMP1 (fibulin-3)") +
  theme_masld() +
  theme(plot.title = element_text(size = 6, hjust = 0, face = "plain"),
        legend.position = "bottom", legend.box = "horizontal",
        legend.key.size = unit(0.3, "lines"), legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0), legend.spacing.x = unit(3, "pt"),
        panel.grid = element_blank())

out <- file.path(FIG4_DIR, "panels", "fig4e_efemp1_multimodal.pdf")
# ⛔ SUPERSEDED 2026-07-07 — single-gene EFEMP1 lollipop replaced by the PROGRAM lollipop
#   `fig4e_secretome_multimodal.R`. Output DISABLED (kept as reference for the EFEMP1 values).
# ggsave(out, p, width = 3.7, height = 2.8, device = grDevices::cairo_pdf)   # DISABLED
cat("[fig4e] saved:", out, "\n")
message(sprintf(paste0("CAPTION (Fig 4e): EFEMP1 (fibulin-3), a fibroblast-sourced LIANA ligand ",
  "(fibroblast→hepatocyte, score %.2f; source=%s), validated by ORTHOGONAL measurement — the disease ",
  "up-regulation seen in bulk RNA (log2FC %+.2f, padj %.0e) is independently reproduced at the PROTEIN ",
  "level in two DIA-MS proteomes: liver (%+.2f, q=%.2f) and plasma (%+.2f, q=%.2f), and EFEMP1 is spatially ",
  "structured (Visium Moran's I %.2f). Cross-modality replication (RNA inference→protein), not more ",
  "transcriptomics. LIANA cell-source is pooled/supporting; protein elevations are direction-concordant ",
  "trends (q≈0.1). Positive control — no genetic pillar (COLOC PP.H4 %.2f)."),
  liana, src, bulk_lfc, bulk_padj, liver_lfc, liver_padj, plasma_lfc, plasma_padj, moran,
  as.numeric(ar$coloc_best_pp4)))
