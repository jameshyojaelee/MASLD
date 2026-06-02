#!/usr/bin/env Rscript
# KEY MESSAGE: Three nuclear receptors — THRB (FDA-approved), RORA (cross-ancestry
# replicated), NR1H4 (Phase 3) — are all recovered through genetic colocalisation
# and transcriptomic suppression, validating the atlas pipeline on established targets.
#
# Output: figures/main/fig4_validation/panels/fig4c_nr_triad.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Data ──────────────────────────────────────────────────────────────────────
atlas <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/integration/multi_evidence_atlas_with_spatial.csv"),
  stringsAsFactors = FALSE)

# GWAS-ATAC motif counts (from manuscript; update path if CSV available)
gwas_atac <- data.frame(
  human_symbol  = c("THRB", "RORA", "NR1H4"),
  n_motif_vars  = c(10L, 11L, NA_integer_),
  in_regulon    = c(2L,  1L, NA_integer_),
  cross_anc     = c("EUR + EAS", "EUR + EAS + SAS", ""),
  drug_label    = c("FDA-approved\n(Resmetirom)", "Preclinical", "Phase 3\n(OCA)"),
  drug_color    = c("#2E7D32", "#9E9E9E", "#F9A825"),
  stringsAsFactors = FALSE
)

genes_oi <- c("THRB", "RORA", "NR1H4")
nr <- atlas %>%
  filter(human_symbol %in% genes_oi) %>%
  select(human_symbol, dream_logFC, dream_padj,
         coloc_susie_best_pp4, best_protein_logFC) %>%
  left_join(gwas_atac, by = "human_symbol") %>%
  mutate(gene = factor(human_symbol, levels = genes_oi))

gene_colors <- c(THRB = "#1565C0", RORA = "#5C6BC0", NR1H4 = "#6A1B9A")

# ── P1: mRNA LFC  (lollipop, not bar) ────────────────────────────────────────
p_expr <- ggplot(nr, aes(y = gene, x = dream_logFC, color = gene)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_segment(aes(y = gene, yend = gene, x = 0, xend = dream_logFC),
               linewidth = 0.6) +
  geom_point(size = 2.5) +
  geom_text(aes(label = sprintf("%.2f", dream_logFC),
                x = dream_logFC + ifelse(dream_logFC < 0, -0.02, 0.02)),
            hjust = ifelse(nr$dream_logFC < 0, 1, 0),
            size = PUB_GEOM_TEXT, color = "black") +
  scale_color_manual(values = gene_colors, guide = "none") +
  scale_x_continuous(limits = c(-0.55, 0.1), breaks = c(-0.4, -0.2, 0)) +
  labs(x = "mRNA log₂FC", y = NULL,
       title = "Expression (dream)") +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_text(face = "bold.italic", size = PUB_AXIS_TEXT + 1))

# ── P2: COLOC PP.H4 (dot, sorted by value) ───────────────────────────────────
p_coloc <- ggplot(nr, aes(y = gene, x = coloc_susie_best_pp4, color = gene)) +
  geom_vline(xintercept = 0.5, linewidth = 0.3, linetype = "dashed", color = "gray60") +
  geom_segment(aes(y = gene, yend = gene, x = 0, xend = coloc_susie_best_pp4),
               linewidth = 0.6) +
  geom_point(size = 2.5) +
  geom_text(aes(label = ifelse(coloc_susie_best_pp4 > 0.9,
                               sprintf("%.4f", coloc_susie_best_pp4),
                               sprintf("%.3f",  coloc_susie_best_pp4))),
            hjust = -0.2, size = PUB_GEOM_TEXT, color = "black") +
  # Cross-ancestry annotation
  geom_text(data = filter(nr, cross_anc != ""),
            aes(label = cross_anc, x = 0.02), hjust = 0,
            size = PUB_GEOM_TEXT - 0.2, color = "gray40",
            nudge_y = -0.3) +
  scale_color_manual(values = gene_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 1.15), breaks = c(0, 0.5, 1.0)) +
  labs(x = "SuSiE-COLOC PP.H4", y = NULL, title = "Genetic colocalisation") +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# ── P3: GWAS-ATAC motif disruptions (Cleveland dot) ──────────────────────────
p_atac <- ggplot(nr, aes(y = gene, x = n_motif_vars, color = gene)) +
  geom_segment(aes(y = gene, yend = gene, x = 0,
                   xend = ifelse(is.na(n_motif_vars), 0, n_motif_vars)),
               linewidth = 0.6) +
  geom_point(aes(x = ifelse(is.na(n_motif_vars), 0, n_motif_vars)), size = 2.5) +
  geom_text(aes(label = ifelse(is.na(n_motif_vars), "n/a",
                               sprintf("%d (%d in regulon)", n_motif_vars, in_regulon)),
                x = ifelse(is.na(n_motif_vars), 0.5, n_motif_vars + 0.3)),
            hjust = 0, size = PUB_GEOM_TEXT, color = "black") +
  scale_color_manual(values = gene_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 16)) +
  labs(x = "GWAS-ATAC motif disruptions", y = NULL,
       title = "Regulatory variants") +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# ── Combine (3 panels, horizontal layout) ────────────────────────────────────
p_out <- (p_expr | p_coloc | p_atac) +
  plot_layout(widths = c(1.2, 1.2, 1.4))

out <- file.path(FIG4_DIR, "panels", "fig4c_nr_triad.pdf")
pdf(out, width = fig_full_width, height = 1.6, useDingbats = FALSE)
print(p_out)
dev.off()
message("Saved: ", out)
