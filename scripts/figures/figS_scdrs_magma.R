#!/usr/bin/env Rscript
# ============================================================================
# figS_scdrs_magma.R  (v2, 2026-05-14 — minimal 2-panel)
#
# KEY MESSAGE: Under canonical MAGMA-anchored scDRS (Zhang 2022 NG), most
# MASLD GWAS yield NO cell-type enrichment at FDR<0.05. The real signals
# are sparse: BBJ GGT → Hepatocytes, FinnGen NASH → Macrophages, plus a
# few Pan-UKBB ancestry-stratified hits.
#
# Two compact panels:
#   A. Per-GWAS bar of #FDR<0.05 CTs, sorted descending, colored by family.
#   B. Focused heatmap of the (CT × GWAS) hits at FDR<0.05.
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

log_msg <- function(...) cat("[figS_scdrs_magma] ", ..., "\n", sep = "")
save_panel <- function(p, slug, width, height) {
  out <- file.path(PANELS_DIR, paste0("figS_scdrs_magma_", slug, ".pdf"))
  save_fig(p, out, width = width, height = height)
  log_msg("  panel saved: ", out)
}

SIG_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures")
LONG     <- file.path(SIG_DIR, "scdrs_magma_celltype_enrichment.csv")
SUM      <- file.path(SIG_DIR, "scdrs_magma_summary.csv")

BASE_SIZE <- 7
LBL_SIZE  <- BASE_SIZE / ggplot2::.pt

theme_min <- function() {
  theme_masld(base_size = BASE_SIZE) +
    theme(
      plot.title    = element_text(size = BASE_SIZE, face = "bold"),
      plot.subtitle = element_blank(),
      panel.grid    = element_blank(),
      legend.key.size = unit(0.3, "cm"),
      legend.title    = element_text(size = BASE_SIZE - 1, face = "bold"),
      legend.text     = element_text(size = BASE_SIZE - 1)
    )
}

FAMILY_ORDER <- c(
  "Liver enzymes (EUR)", "Liver enzymes (EAS)",
  "Liver enzymes (AFR)", "Liver enzymes (CSA)",
  "NAFLD", "NASH", "HCC", "Cirrhosis", "PDFF (imaging)"
)
FAMILY_COLORS <- c(
  "Liver enzymes (EUR)" = "#1565C0",
  "Liver enzymes (EAS)" = "#7B1FA2",
  "Liver enzymes (AFR)" = "#F57C00",
  "Liver enzymes (CSA)" = "#FBC02D",
  "NAFLD"               = "#C2185B",
  "NASH"                = "#880E4F",
  "HCC"                 = "#4A148C",
  "Cirrhosis"           = "#3E2723",
  "PDFF (imaging)"      = "#00695C"
)

# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------
log_msg("loading magma enrichment + summary")
sumt <- fread(SUM)
dt   <- fread(LONG)

sumt[, family := factor(family, levels = FAMILY_ORDER)]
setorder(sumt, -n_ct_fdr05, family)
sumt[, gwas_label := factor(gwas_id, levels = unique(gwas_id))]

# ----------------------------------------------------------------------------
# Panel A — per-GWAS bar of #sig CTs (sorted, family colors)
# ----------------------------------------------------------------------------
log_msg("Panel A")
fwrite(sumt, file.path(DATA_DIR, "figS_scdrs_magma_panelA_data.csv"))

panel_A <- ggplot(sumt, aes(x = gwas_label, y = n_ct_fdr05, fill = family)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(n_ct_fdr05 > 0,
                               paste0(top_ct, "\n(z=", sprintf("%.1f", top_z), ")"),
                               "")),
            vjust = -0.25, size = LBL_SIZE * 0.7, lineheight = 0.9,
            colour = "grey15") +
  scale_fill_manual(values = FAMILY_COLORS, name = "GWAS family",
                    guide = guide_legend(ncol = 3)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8),
                     expand = expansion(mult = c(0, 0.20))) +
  labs(x = NULL, y = "# CTs at FDR < 0.05",
       title = "A. MAGMA-anchored scDRS: most MASLD GWAS yield no significant cell type") +
  theme_min() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   size = BASE_SIZE - 1),
        legend.position = "bottom")

save_panel(panel_A, "A_per_gwas_summary", width = 7.5, height = 4.0)

# ----------------------------------------------------------------------------
# Panel B — Focused heatmap of FDR<0.05 hits
# ----------------------------------------------------------------------------
# Restrict to GWAS that have ≥1 significant CT, and to the cell types
# that appear as significant in any of them. The empty cells will show
# the rest of the (CT × GWAS) grid is null.
log_msg("Panel B: focused heatmap of FDR<0.05 hits")
hit_gwas <- sumt[n_ct_fdr05 > 0]$gwas_id
hit_dt   <- dt[gwas_id %in% hit_gwas]
hit_dt[, sig := assoc_mcp < 0.05]
sig_cts  <- unique(hit_dt[sig == TRUE]$cell_type)
B_data   <- hit_dt[cell_type %in% sig_cts]

# Order GWAS by family
B_data[, family    := factor(family, levels = FAMILY_ORDER)]
B_data[, gwas_label := factor(gwas_id, levels = sumt[gwas_id %in% hit_gwas]$gwas_id)]
B_data[, cell_type := factor(cell_type, levels = rev(sig_cts))]
fwrite(B_data, file.path(DATA_DIR, "figS_scdrs_magma_panelB_data.csv"))

zmax_B <- max(abs(B_data$assoc_mcz), na.rm = TRUE)

panel_B <- ggplot(B_data, aes(x = gwas_label, y = cell_type, fill = assoc_mcz)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_point(data = B_data[sig == TRUE],
             colour = "black", size = 0.9, shape = 1, stroke = 0.5) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0, limits = c(-zmax_B, zmax_B),
                       name = "z-score") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "B. The few significant findings — CT × GWAS hits (○ = FDR<0.05)") +
  theme_min() +
  theme(axis.line  = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,
                                   size = BASE_SIZE - 1),
        legend.position = "right")

save_panel(panel_B, "B_significant_hits", width = 6.0, height = 3.4)

# ----------------------------------------------------------------------------
# Compose
# ----------------------------------------------------------------------------
log_msg("Composing figS_scdrs_magma.pdf")
full <- panel_A / panel_B + plot_layout(heights = c(1.2, 0.95))

out_pdf <- file.path(OUT_DIR, "figS_scdrs_magma.pdf")
save_fig(full, out_pdf, width = 8.0, height = 7.0)
log_msg("Wrote ", out_pdf)
log_msg("Done.")
