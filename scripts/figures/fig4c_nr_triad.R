#!/usr/bin/env Rscript
# KEY MESSAGE: Three nuclear receptors — THRB (FDA-approved), NR1H4 (clinical-
# stage), RORA (preclinical) — are recovered through genetic colocalisation and
# transcriptomic suppression without prior curation. Points are coloured by drug-
# development stage so the panel reads as "the prioritiser recovers established
# and emerging targets," not an arbitrary gene key.
#
# Output: figures/main/fig4_validation/nuclear_receptor_triad.pdf
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

# GWAS-ATAC motif-disruption counts — READ from the canonical motifbreakR table
# (same source as target_evidence_matrix.R) so the two panels cannot disagree.
# NEVER hardcode these (a stale hardcode previously gave RORA=11/NR1H4=0; the
# file gives 10/9/2 = THRB / RORA / NR1H4::RXRA).
md <- read.csv(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"),
  stringsAsFactors = FALSE)
nmot <- function(g) sum(grepl(g, md$tf_name, ignore.case = TRUE))

gwas_atac <- data.frame(
  human_symbol = c("THRB", "RORA", "NR1H4"),
  n_motif_vars = c(nmot("THRB"), nmot("RORA"), nmot("NR1H4")),
  drug_stage   = c("FDA-approved", "Preclinical", "Clinical-stage"),
  stringsAsFactors = FALSE
)

genes_oi <- c("THRB", "RORA", "NR1H4")
nr <- atlas %>%
  filter(human_symbol %in% genes_oi) %>%
  select(human_symbol, bulk_logFC, bulk_padj, coloc_susie_best_pp4) %>%
  left_join(gwas_atac, by = "human_symbol") %>%
  mutate(gene = factor(human_symbol, levels = genes_oi))

# Drug-development stage: traffic-light (approved green -> clinical amber ->
# preclinical neutral gray). This is the single colour axis across all 3 panels.
stage_levels <- c("FDA-approved", "Clinical-stage", "Preclinical")
stage_cols   <- c("FDA-approved"   = "#2E7D32",
                  "Clinical-stage" = "#F9A825",
                  "Preclinical"    = "#9E9E9E")
nr$drug_stage <- factor(nr$drug_stage, levels = stage_levels)

base_lolli <- function(df, x, xlab, xlim, xbreaks, vline = NULL, vtype = "solid",
                       show_y = FALSE) {
  p <- ggplot(df, aes(y = gene, x = .data[[x]], color = drug_stage))
  if (!is.null(vline))
    p <- p + geom_vline(xintercept = vline, linewidth = 0.3,
                        linetype = vtype, color = "gray60")
  p <- p +
    # Lollipop stems anchor at 0 (the null), so negative LFCs point LEFT from 0.
    geom_segment(aes(y = gene, yend = gene, x = 0, xend = .data[[x]]),
                 linewidth = 0.6) +
    geom_point(size = 2.6) +
    scale_color_manual(values = stage_cols, drop = FALSE, name = "Drug stage") +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    labs(x = xlab, y = NULL) +
    theme_masld() + theme_pub()
  if (show_y) {
    p <- p + theme(axis.text.y = element_text(face = "bold.italic",
                                              size = PUB_AXIS_TEXT + 1))
  } else {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  p
}

p_expr  <- base_lolli(nr, "bulk_logFC", "mRNA log2FC",
                      c(-0.55, 0.1), c(-0.4, -0.2, 0), vline = 0, show_y = TRUE)
# Axis = "COLOC PP.H4" (NOT "SuSiE-"): THRB and NR1H4 have no SuSiE posterior at
# their locus, so coloc_susie_best_pp4 carries the coloc.abf fallback for them.
p_coloc <- base_lolli(nr, "coloc_susie_best_pp4", "COLOC PP.H4",
                      c(0, 1.05), c(0, 0.5, 1.0), vline = 0.5, vtype = "dashed")
p_atac  <- base_lolli(nr, "n_motif_vars", "GWAS-ATAC motif disruptions",
                      c(0, 13), c(0, 5, 10))

# ── Combine (3 lollipops; one shared drug-stage legend at the bottom) ─────────
p_out <- (p_expr | p_coloc | p_atac) +
  plot_layout(widths = c(1.2, 1.2, 1.4), guides = "collect") &
  theme(legend.position = "bottom", legend.key.size = PUB_LEGEND_KEY)

out <- file.path(FIG4_DIR, "nuclear_receptor_triad.pdf")
pdf(out, width = fig_full_width, height = 1.85, useDingbats = FALSE)
print(p_out)
dev.off()
message("Saved: ", out)
