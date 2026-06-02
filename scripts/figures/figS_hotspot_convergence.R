#!/usr/bin/env Rscript
# ============================================================================
# figS_hotspot_convergence.R
#
# KEY MESSAGE: Two Hotspot modules (HKDC1-stress hep__24 and NR1H4-loss fib__7)
# co-validate with bulk NMF, drug pipelines, and COLOC. The triple-pipeline
# table anchors the top therapy-relevant genes in three orthogonal sources.
#
# A — HKDC1 spotlight (hep__24 — gene weights + donor score by stage)
# B — NR1H4/FXR loss spotlight (fib__7 — gene weights + donor score by stage)
# C — Cross-pipeline anchor: 5 drug-relevant genes x {Hotspot, COLOC, drug}
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- FIGS_HOTSPOT_PANELS_DIR
DATA_DIR  <- FIGS_HOTSPOT_DATA_DIR
PREFIX    <- "convergence_"

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
# Dataset-center scores per (cell_type, module, dataset) — matches the
# (1|dataset) random intercept used by 505 to estimate β.
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

# Helper: spotlight a single module (gene-weight + stage-trajectory)
make_spotlight <- function(ct, mod_int, title, n_genes = 12, bar_color) {
  gf <- file.path(HS_RES, ct, "module_genes.tsv")
  g <- fread(gf)[module == mod_int][order(-weight)][1:n_genes]
  g[, gene := factor(gene, levels = rev(gene))]
  pg <- ggplot(g, aes(weight, gene)) +
    geom_col(width = 0.7, fill = bar_color, color = "black", linewidth = 0.15) +
    labs(x = "Hotspot weight", y = NULL,
         subtitle = sprintf("%s — top %d members", title, n_genes)) +
    theme_masld() + theme_pub() +
    theme(axis.text.y = element_text(size = 5, face = "bold"))

  d <- ds[cell_type == ct & module_int == mod_int]
  d <- d[disease_stage_coarse %in% STAGE_LEVELS]
  d[, disease_stage_coarse := factor(disease_stage_coarse, levels = STAGE_LEVELS)]
  m <- mods[cell_type == ct & module == mod_int]
  trend <- d[, .(med = median(score_centered, na.rm = TRUE)),
             by = disease_stage_coarse]
  pt <- ggplot(d, aes(disease_stage_coarse, score_centered,
                      fill = disease_stage_coarse)) +
    geom_violin(trim = TRUE, scale = "width",
                linewidth = 0.2, alpha = 0.85) +
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
  list(genes = pg, traj = pt)
}

A_parts <- make_spotlight("hepatocytes", 24, "HKDC1 stress (hep__24)",
                          bar_color = ct_palette[["Hepatocytes"]])
B_parts <- make_spotlight("fibroblasts", 7,  "NR1H4 / FXR loss (fib__7)",
                          bar_color = ct_palette[["Fibroblasts"]])

pA <- A_parts$genes | A_parts$traj
pB <- B_parts$genes | B_parts$traj

# ----------------------------------------------------------------------------
# Panel C — triple-pipeline cross-reference for 5 priority genes
# ----------------------------------------------------------------------------
ATLAS_ALL <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
susie_col <- intersect(c("coloc_best_pp4_polyfun", "coloc_best_susie_pp4_polyfun",
                          "coloc_best_susie_pp4"), names(ATLAS_ALL))[1]
abf_col   <- intersect(c("coloc_abf_best_pp4", "coloc_pp4"), names(ATLAS_ALL))[1]
keep <- c("human_symbol", "dream_logFC", "dream_padj", susie_col, abf_col)
ATLAS <- ATLAS_ALL[, ..keep]
setnames(ATLAS, c("human_symbol", susie_col, abf_col),
         c("gene", "coloc_susie_pp4", "coloc_abf_pp4"))

# Narrative-anchored CT-specific Hotspot module per gene (NOT atlas-wide top —
# `hotspot_top_module` is dominated by the global run which obscures lineage).
NARRATIVE <- data.table(
  gene = c("HKDC1", "NR1H4", "THRB", "PDGFRA", "RORA"),
  ct   = c("hepatocytes", "fibroblasts", "endothelial_cells", "global", "fibroblasts"),
  module_int = c(24L, 7L, 22L, 9L, 13L)
)
NARRATIVE[, module_label := sprintf("%s__%d", ct, module_int)]
# Pull β/q/stab for each from all_modules.tsv (the canonical phenotype stats)
ms <- mods[, .(ct = cell_type, module_int = module, beta = disease_stage_beta,
               q = disease_stage_q, progression_module, stability_score)]
NARRATIVE <- merge(NARRATIVE, ms, by = c("ct", "module_int"), all.x = TRUE)

# Drug — schema (verified 2026-05-17): drug, target_gene, stage, moa, in_atlas, atlas_support
DRUG <- fread(file.path(BASE, "RNA-seq/results/drug_repurposing/clinical_drug_validation_table_v2.csv"))
genes_panel <- c("HKDC1", "NR1H4", "THRB", "PDGFRA", "RORA")
drug_rows <- DRUG[target_gene %in% genes_panel,
                  .(gene = target_gene, drug, evidence = atlas_support)]
drug_summary <- drug_rows[, .(drug     = paste(unique(drug), collapse = "; "),
                              evidence = paste(unique(evidence), collapse = "/")),
                          by = gene]

# Cross-reference atlas signals
atlas_signals <- ATLAS[gene %in% genes_panel,
                       .(gene,
                         coloc_pp4 = pmax(coloc_susie_pp4, coloc_abf_pp4, na.rm = TRUE),
                         dream_logFC = dream_logFC)]
atlas_signals[, coloc_pp4 := round(ifelse(is.finite(coloc_pp4), coloc_pp4, NA), 2)]
atlas_signals[, dream_logFC := round(dream_logFC, 2)]

pickC <- merge(NARRATIVE, atlas_signals, by = "gene", all.x = TRUE)
pickC <- merge(pickC, drug_summary, by = "gene", all.x = TRUE)
pickC[is.na(drug), drug := "—"]
pickC[is.na(evidence), evidence := "—"]
pickC[, gene := factor(gene, levels = genes_panel)]
pickC <- pickC[order(gene)]
pickC[, sig_str := sprintf("β=%+0.2f\nq=%.1e", beta, q)]
pickC[, drug_str := ifelse(drug == "None" | drug == "—", "—",
                            paste0(drug, "\n(", evidence, ")"))]

# Reshape to long for tile plot
plot_long <- melt(pickC[, .(gene, module_label, sig_str,
                            coloc_pp4, dream_logFC, drug_str)],
                  id.vars = "gene", variable.name = "field", value.name = "value")
plot_long[, value_str := ifelse(is.na(value) | value == "NA", "—",
                                  as.character(value))]
plot_long[, field := factor(field, levels = c("module_label", "sig_str",
                                              "coloc_pp4", "dream_logFC",
                                              "drug_str"),
                            labels = c("Hotspot module\n(CT-specific)",
                                       "Disease β / q",
                                       "COLOC PP4\n(best)",
                                       "dream\nlogFC",
                                       "Drug pipeline\n(evidence)"))]

pC <- ggplot(plot_long, aes(field, gene)) +
  geom_tile(fill = "white", color = "gray70", linewidth = 0.3) +
  geom_text(aes(label = value_str), size = 1.7, lineheight = 0.95) +
  labs(x = NULL, y = NULL,
       subtitle = "Cross-pipeline anchor: Hotspot module + dream + COLOC + drug pipeline") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,
                                   size = PUB_AXIS_TEXT, lineheight = 1.1),
        axis.text.y = element_text(face = "bold", size = PUB_AXIS_TEXT + 1),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

# ----------------------------------------------------------------------------
# Save panels + composite
# ----------------------------------------------------------------------------
save_fig(pA, file.path(PANEL_DIR, paste0(PREFIX, "A_hkdc1_spotlight.pdf")),
         width = fig_full_width, height = 2.4)
save_fig(pB, file.path(PANEL_DIR, paste0(PREFIX, "B_nr1h4_spotlight.pdf")),
         width = fig_full_width, height = 2.4)
save_fig(pC, file.path(PANEL_DIR, paste0(PREFIX, "C_triple_pipeline.pdf")),
         width = fig_full_width, height = 1.8)

composite <- pA / pB / pC +
  patchwork::plot_layout(heights = c(1, 1, 1)) +
  patchwork::plot_annotation(tag_levels = "a")
save_fig(composite,
         file.path(FIGS_HOTSPOT_DIR, "figS_hotspot_convergence.pdf"),
         width = fig_full_width, height = 6.8)

fwrite(pickC, file.path(DATA_DIR, paste0(PREFIX, "C_triple_pipeline.csv")))
cat("Wrote figS_hotspot_convergence: 3 panels + composite\n")
