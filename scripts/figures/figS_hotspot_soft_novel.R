#!/usr/bin/env Rscript
# ============================================================================
# figS_hotspot_soft_novel.R
#
# KEY MESSAGE: Local-autocorrelation discovery surfaces two clinically
# meaningful modules whose gene composition is distinct from cNMF / bulk NMF /
# Hallmark / SCENIC+ / curated panels at the top-50 Jaccard level. mac__16 =
# tissue-resident Kupffer attrition; hep__25 = UPR / ERN1 stress signature.
# Both are gene-composition-novel (best Jaccard < 0.03 vs every reference
# panel) yet emerge with disease-coherent direction.
#
# A — mac__16 Kupffer-loss hero (gene weights + F-stage trajectory)
# B — hep__25 UPR / ERN1 hero (gene weights + F-stage trajectory)
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- FIGS_HOTSPOT_PANELS_DIR
DATA_DIR  <- FIGS_HOTSPOT_DATA_DIR
PREFIX    <- "soft_novel_"

HS_RES <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
META_FILE <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv")

mods <- fread(file.path(HS_RES, "all_modules.tsv"))
ds <- fread(file.path(HS_RES, "donor_scores_all.tsv"))
meta <- fread(META_FILE,
              select = c("sample", "disease_stage_coarse", "dataset",
                         "exclude_stage_analysis"))
# Backwards-compat for protocol contamination remediation
if (!"exclude_stage_analysis" %in% names(meta)) meta[, exclude_stage_analysis := FALSE]
meta <- meta[exclude_stage_analysis != TRUE]
ds <- ds[meta, on = "sample", nomatch = 0]
ds[, module_int := as.integer(sub(".*__", "", module))]
# Dataset-center per (cell_type, module, dataset) so the violins match what
# 505's lmer(score ~ stage + (1|dataset)) actually estimates. Raw pooled
# violins flip the direction via Simpson's paradox (within-dataset Spearman
# is negative for mac__16 even though pooled is positive).
ds[, score_centered := score - mean(score, na.rm = TRUE),
   by = .(cell_type, module_int, dataset)]

# Cirrhosis excluded: hepatocytes are depleted by fibrotic replacement in
# cirrhotic tissue; the 19 available cirrhosis donors (GSE202379 snRNA-seq,
# single dataset) represent a survivor-selected subpopulation not comparable
# to earlier-stage hepatocytes. 3-stage cascade is the primary analysis axis.
STAGE_LEVELS <- c("Healthy", "Steatosis", "Steatohepatitis")
STAGE_COLORS <- c(Healthy        = masld_colors$control,
                  Steatosis      = masld_colors$masl,
                  Steatohepatitis= masld_colors$mash)

make_softnovel_hero <- function(ct, mod_int, header, n_genes = 14, bar_color) {
  gf <- file.path(HS_RES, ct, "module_genes.tsv")
  g <- fread(gf)[module == mod_int][order(-weight)][1:n_genes]
  g[, gene := factor(gene, levels = rev(gene))]
  m <- mods[cell_type == ct & module == mod_int]

  p_genes <- ggplot(g, aes(weight, gene)) +
    geom_col(width = 0.7, fill = bar_color, color = "black", linewidth = 0.15) +
    labs(x = "Hotspot weight (mean local-corr Z)", y = NULL,
         subtitle = sprintf("%s\nbest-jacc=%.3f vs %s",
                            header, m$best_match_jaccard,
                            m$best_match_program)) +
    theme_masld() + theme_pub() +
    theme(axis.text.y = element_text(size = 5, face = "bold"))

  d <- ds[cell_type == ct & module_int == mod_int]
  d <- d[disease_stage_coarse %in% STAGE_LEVELS]
  d[, disease_stage_coarse := factor(disease_stage_coarse, levels = STAGE_LEVELS)]
  trend <- d[, .(med = median(score_centered, na.rm = TRUE)),
             by = disease_stage_coarse]
  p_traj <- ggplot(d, aes(disease_stage_coarse, score_centered,
                          fill = disease_stage_coarse)) +
    geom_violin(trim = TRUE, scale = "width", linewidth = 0.2, alpha = 0.85) +
    geom_line(data = trend, aes(x = disease_stage_coarse, y = med, group = 1),
              inherit.aes = FALSE, color = "gray25", linewidth = 0.35) +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.1,
                 color = "black", fill = "white", stroke = 0.3) +
    scale_fill_manual(values = STAGE_COLORS, guide = "none") +
    labs(x = NULL, y = "Module score (dataset-centered)",
         subtitle = sprintf("β=%+0.2f  q=%.1e  stab=%.2f",
                            m$disease_stage_beta, m$disease_stage_q,
                            m$stability_score)) +
    theme_masld() + theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  list(genes = p_genes, traj = p_traj, gene_data = g, donor_data = d)
}

A_parts <- make_softnovel_hero("macrophages", 16,
                               "Kupffer-loss soft-novel (mac__16)",
                               n_genes = 14,
                               bar_color = ct_palette[["Macrophages"]])
B_parts <- make_softnovel_hero("hepatocytes", 25,
                               "UPR / ERN1 soft-novel (hep__25)",
                               n_genes = 14,
                               bar_color = ct_palette[["Hepatocytes"]])

pA <- A_parts$genes | A_parts$traj
pB <- B_parts$genes | B_parts$traj

save_fig(pA, file.path(PANEL_DIR, paste0(PREFIX, "A_mac16_kupffer_loss.pdf")),
         width = fig_full_width, height = 2.5)
save_fig(pB, file.path(PANEL_DIR, paste0(PREFIX, "B_hep25_upr.pdf")),
         width = fig_full_width, height = 2.5)

composite <- pA / pB +
  patchwork::plot_annotation(tag_levels = "a")
save_fig(composite,
         file.path(FIGS_HOTSPOT_DIR, "figS_hotspot_soft_novel.pdf"),
         width = fig_full_width, height = 5.2)

fwrite(A_parts$gene_data, file.path(DATA_DIR, paste0(PREFIX, "A_mac16_genes.csv")))
fwrite(A_parts$donor_data, file.path(DATA_DIR, paste0(PREFIX, "A_mac16_donors.csv")))
fwrite(B_parts$gene_data, file.path(DATA_DIR, paste0(PREFIX, "B_hep25_genes.csv")))
fwrite(B_parts$donor_data, file.path(DATA_DIR, paste0(PREFIX, "B_hep25_donors.csv")))

cat("Wrote figS_hotspot_soft_novel: 2 panels + composite\n")
