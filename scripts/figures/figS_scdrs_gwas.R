#!/usr/bin/env Rscript
# ============================================================================
# figS_scdrs_gwas.R  (v5, 2026-05-14 — heatmap + bar + disease-stage)
#
# KEY MESSAGE: Hepatocytes are the only cell type with significant scDRS
# enrichment for MASLD COLOC-causal genes. The Hep signal is highest in
# healthy livers and is partially lost in disease (identity-loss pattern).
#
# Three panels:
#   A. Heatmap — 12 CTs × 5 COLOC panels, fill = z, dot = FDR<0.05.
#   B. Bar — Hepatocyte z across the same 5 COLOC panels.
#   C. Bar — Per-cell median scDRS by disease stage (Healthy / MASL / MASH)
#      across top 4 CTs, liver-enzyme COLOC panel.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR    <- FIGS_SCDRS_DIR
DATA_DIR   <- FIGS_SCDRS_DATA_DIR
PANELS_DIR <- file.path(FIGS_SCDRS_DIR, "panels")
dir.create(PANELS_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat("[figS_scdrs_gwas] ", ..., "\n", sep = "")
save_panel <- function(p, slug, width, height) {
  out <- file.path(PANELS_DIR, paste0("figS_scdrs_gwas_", slug, ".pdf"))
  save_fig(p, out, width = width, height = height)
  log_msg("  panel saved: ", out)
}

SIG_DIR     <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures")
GWAS_ENRICH <- file.path(SIG_DIR, "scdrs_gwas_celltype_enrichment.csv")
STAGE_MED   <- file.path(SIG_DIR, "figure_prep/gwas_cell_scores_median_by_stage.csv")

BASE_SIZE <- 7
LBL_SIZE  <- BASE_SIZE / ggplot2::.pt

theme_min <- function() {
  theme_masld(base_size = BASE_SIZE) +
    theme(
      plot.title    = element_text(size = BASE_SIZE, face = "bold"),
      plot.subtitle = element_blank(),
      panel.grid    = element_blank(),
      legend.key.size = unit(0.3, "cm")
    )
}

# 5 panels: keep Pooled PP>0.5 (broad MASLD signal) + 4 disease phenotypes.
PHENO_ORDER <- c(
  "gwas_coloc_pooled_pp05"    = "Pooled PP>0.5",
  "gwas_coloc_liver_enzymes"  = "Liver enzymes",
  "gwas_coloc_nafld_specific" = "NAFLD-specific",
  "gwas_coloc_pdff"           = "PDFF (imaging)",
  "gwas_coloc_cirrhosis_hcc"  = "Cirrhosis / HCC"
)

# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------
log_msg("loading enrichment table (PolyFun, 5 panels)")
dt <- fread(GWAS_ENRICH)
dt <- dt[panel == "polyfun" & trait_class %in% names(PHENO_ORDER)]
dt[, pheno := factor(PHENO_ORDER[as.character(trait_class)],
                     levels = unname(PHENO_ORDER))]
dt[, sig := assoc_mcp < 0.05]

# Rank CTs by max z across the 5 phenotypes
ct_rank <- dt[, .(rank_z = max(assoc_mcz, na.rm = TRUE)),
              by = cell_type][order(-rank_z)]
dt[, cell_type := factor(cell_type, levels = rev(ct_rank$cell_type))]

# ----------------------------------------------------------------------------
# Panel A — CT × phenotype heatmap
# ----------------------------------------------------------------------------
log_msg("Panel A")
A_data <- copy(dt[, .(cell_type, pheno, assoc_mcz, assoc_mcp, sig)])
fwrite(A_data, file.path(DATA_DIR, "figS_scdrs_gwas_panelA_data.csv"))

zmax <- max(abs(A_data$assoc_mcz), na.rm = TRUE)

panel_A <- ggplot(A_data, aes(x = pheno, y = cell_type, fill = assoc_mcz)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_point(data = A_data[sig == TRUE],
             colour = "black", size = 0.9, shape = 1, stroke = 0.5) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0, limits = c(-zmax, zmax),
                       name = "z-score") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "A. Cell-type enrichment per COLOC panel  (○ = FDR < 0.05)") +
  theme_min() +
  theme(axis.line  = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right")

save_panel(panel_A, "A_celltype_heatmap", width = 5.2, height = 4.0)

# ----------------------------------------------------------------------------
# Panel B — Hepatocyte bar chart across panels
# ----------------------------------------------------------------------------
log_msg("Panel B")
B_data <- A_data[cell_type == "Hepatocytes"]
B_data[, pheno := factor(pheno, levels = B_data[order(assoc_mcz)]$pheno)]
fwrite(B_data, file.path(DATA_DIR, "figS_scdrs_gwas_panelB_data.csv"))

panel_B <- ggplot(B_data, aes(x = pheno, y = assoc_mcz,
                              fill = sig)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = 1.96, colour = "grey50",
             linetype = "dashed", linewidth = 0.3) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.2f%s", assoc_mcz,
                                ifelse(sig, " *", ""))),
            vjust = ifelse(B_data$assoc_mcz > 0, -0.3, 1.2),
            size = LBL_SIZE * 0.85, colour = "grey15") +
  scale_fill_manual(values = c(`FALSE` = "#CCCCCC", `TRUE` = "#C2185B"),
                    labels = c("FDR ≥ 0.05", "FDR < 0.05"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "Hepatocyte z-score",
       title = "B. Hepatocyte signal by COLOC panel") +
  theme_min() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right")

save_panel(panel_B, "B_hepatocyte_bar", width = 4.4, height = 3.4)

# ----------------------------------------------------------------------------
# Panel C — Per-cell median norm_score by disease stage (liver enzymes panel)
# Labels: Healthy / MASL / MASH (using disease_stage_coarse).
# ----------------------------------------------------------------------------
log_msg("Panel C")
TOP_CT <- c("Hepatocytes", "Macrophages", "Fibroblasts", "Cholangiocytes")
STAGE_RENAME <- c("Healthy" = "Healthy",
                  "Steatosis" = "MASL",
                  "Steatohepatitis" = "MASH")
STAGE_ORDER  <- c("Healthy", "MASL", "MASH")
STAGE_FILL   <- c("Healthy" = "#9E9E9E",
                  "MASL"    = "#F48FB1",
                  "MASH"    = "#C2185B")

cs <- fread(STAGE_MED)
cs <- cs[trait_class == "gwas_coloc_liver_enzymes" & cell_type %in% TOP_CT]
cs[, stage := factor(STAGE_RENAME[as.character(disease_stage)],
                     levels = STAGE_ORDER)]
cs[, cell_type := factor(cell_type, levels = TOP_CT)]
fwrite(cs, file.path(DATA_DIR, "figS_scdrs_gwas_panelC_data.csv"))

panel_C <- ggplot(cs, aes(x = cell_type, y = median, fill = stage)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_text(aes(label = sprintf("%.2f", median)),
            position = position_dodge(width = 0.8),
            vjust = ifelse(cs$median > 0, -0.3, 1.2),
            size = LBL_SIZE * 0.75, colour = "grey15") +
  scale_fill_manual(values = STAGE_FILL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.20, 0.15))) +
  labs(x = NULL, y = "median per-cell z-score",
       title = "C. Cell-of-origin signal lives in healthy hepatocytes  (Healthy / MASL / MASH)") +
  theme_min() +
  theme(legend.position = "right")

save_panel(panel_C, "C_disease_stage", width = 7.5, height = 3.2)

# ----------------------------------------------------------------------------
# Compose
# ----------------------------------------------------------------------------
log_msg("Composing figS_scdrs_gwas.pdf")
top <- panel_A + panel_B + plot_layout(widths = c(1.2, 1.0))
full <- top / panel_C + plot_layout(heights = c(1.3, 0.9))

out_pdf <- file.path(OUT_DIR, "figS_scdrs_gwas.pdf")
save_fig(full, out_pdf, width = 8.5, height = 7.0)
log_msg("Wrote ", out_pdf)
log_msg("Done.")
