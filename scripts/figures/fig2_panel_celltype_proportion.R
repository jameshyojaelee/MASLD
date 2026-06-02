#!/usr/bin/env Rscript
# ============================================================================
# fig2_panel_celltype_proportion.R
# Fig 2 panel c — Cell-type proportion trajectory across F0-F4
#
# Two side-by-side compact panels sharing x = F0..F4:
#   Left:  Hepatocyte proportion loss (line + bootstrap CI ribbon)
#   Right: Non-hepatocyte renormalized stacked area (immune/stromal rise)
#
# Data: BayesPrism, 1,444 staged bulk donors across 9 cohorts.
#
# Output: figures/main/fig2_progression_sex/panels/fig2_panel_celltype_proportion.pdf
#   sized 180 x 55 mm.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PDF   <- file.path(PANEL_DIR, "fig2_panel_celltype_proportion.pdf")

PROP_FILE <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression",
  "cibersortx_celltype_expression/bayesprism_proportions.csv")
META_FILE <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
QC_FILE   <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv")

# ---------------------------------------------------------------------------
# Cell-type grouping
# ---------------------------------------------------------------------------
bp_to_umbrella <- c(
  Hepatocyte    = "Hepatocytes",
  Cholangiocyte = "Cholangiocytes",
  Endothelial   = "Endothelial",
  Stellate      = "Stellate/Fib",
  Macrophage    = "Macrophages",
  Monocyte      = "Monocytes",
  DC            = "Dendritic",
  Neutrophil    = "Granulocytes",
  T_cell        = "T cells",
  B_cell        = "B cells",
  Plasma_cell   = "Plasma cells",
  NK_cell       = "NK cells",
  Other_immune  = "Other immune"
)

umbrella_order <- c(
  "Hepatocytes", "Cholangiocytes", "Endothelial", "Stellate/Fib",
  "Macrophages", "Monocytes", "Dendritic", "Granulocytes",
  "T cells", "B cells", "Plasma cells", "NK cells", "Other immune"
)

# Palette: hepatocytes = forest green; immune = warm hues; stromal = orange/teal
umbrella_pal <- c(
  "Hepatocytes"   = "#1B5E20",
  "Cholangiocytes"= "#00BCD4",
  "Endothelial"   = "#7CB342",
  "Stellate/Fib"  = "#FF6F00",
  "Macrophages"   = "#212121",
  "Monocytes"     = "#EC407A",
  "Dendritic"     = "#6A1B9A",
  "Granulocytes"  = "#F9A825",
  "T cells"       = "#1A237E",
  "B cells"       = "#5D4037",
  "Plasma cells"  = "#8D6E63",
  "NK cells"      = "#00897B",
  "Other immune"  = "#9E9E9E"
)

# ---------------------------------------------------------------------------
# Load and prepare
# ---------------------------------------------------------------------------
prop <- fread(PROP_FILE)
meta <- fread(META_FILE)
qc   <- fread(QC_FILE)

meta <- merge(meta, qc[, .(sample_id, pass_technical)], by = "sample_id", all.x = TRUE)
meta_staged <- meta[pass_technical == TRUE & !is.na(fibrosis_stage)]
n_donors  <- nrow(meta_staged)
n_cohorts <- uniqueN(meta_staged$dataset)
cat(sprintf("[data] %d staged QC-pass donors across %d cohorts\n", n_donors, n_cohorts))

prop_long <- melt(prop, id.vars = "sample_id",
                  variable.name = "celltype_bp", value.name = "proportion")
prop_long[, celltype_bp := as.character(celltype_bp)]
prop_long[, umbrella := bp_to_umbrella[celltype_bp]]
prop_umb <- prop_long[, .(proportion = sum(proportion)), by = .(sample_id, umbrella)]
prop_umb <- merge(prop_umb,
                  meta_staged[, .(sample_id, fibrosis_stage)],
                  by = "sample_id")
prop_umb[, fibrosis_stage := factor(fibrosis_stage, levels = 0:4,
                                    labels = paste0("F", 0:4))]
prop_umb[, umbrella := factor(umbrella, levels = umbrella_order)]

boot_ci <- function(x, n_boot = 1000) {
  if (length(x) < 2) return(c(mean = mean(x), lo = NA_real_, hi = NA_real_))
  set.seed(42)
  bs <- replicate(n_boot, mean(sample(x, replace = TRUE)))
  c(mean = mean(x), lo = quantile(bs, 0.025), hi = quantile(bs, 0.975))
}

stage_summ <- prop_umb[, {
  ci <- boot_ci(proportion)
  list(mean_prop = ci[[1]], lo = ci[[2]], hi = ci[[3]])
}, by = .(umbrella, fibrosis_stage)]
stage_summ[, stage_num := as.integer(fibrosis_stage)]

# ---------------------------------------------------------------------------
# Left panel: hepatocyte loss
# ---------------------------------------------------------------------------
hep <- stage_summ[umbrella == "Hepatocytes"]

p_hep <- ggplot(hep, aes(x = stage_num, y = mean_prop)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = umbrella_pal["Hepatocytes"], alpha = 0.20) +
  geom_line(color = umbrella_pal["Hepatocytes"], linewidth = 0.9) +
  geom_point(color = umbrella_pal["Hepatocytes"], size = 1.5) +
  geom_text(aes(label = sprintf("%.0f%%", mean_prop * 100)),
            vjust = -0.8, size = 1.9, color = "grey20") +
  scale_x_continuous(breaks = 1:5, labels = paste0("F", 0:4),
                     expand = c(0.06, 0.06)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0.63, 0.93),
                     breaks = c(0.7, 0.8, 0.9)) +
  labs(title = "Hepatocyte loss",
       subtitle = sprintf("n=%d donors, %d cohorts", n_donors, n_cohorts),
       x = "Fibrosis stage", y = "Proportion") +
  theme_masld(base_size = 7) +
  theme(
    plot.title    = element_text(size = 7.5, face = "bold"),
    plot.subtitle = element_text(size = 5.5, color = "grey35"),
    panel.grid.minor = element_blank()
  )

# ---------------------------------------------------------------------------
# Right panel: non-hep renormalized stacked area
# ---------------------------------------------------------------------------
non_hep <- copy(stage_summ[umbrella != "Hepatocytes"])
non_hep[, norm_prop := mean_prop / sum(mean_prop), by = fibrosis_stage]
non_hep[, umbrella := factor(as.character(umbrella),
        levels = setdiff(umbrella_order, "Hepatocytes"))]

p_stack <- ggplot(non_hep,
                  aes(x = stage_num, y = norm_prop,
                      fill = umbrella, group = umbrella)) +
  geom_area(position = "stack", alpha = 0.90, color = "white", linewidth = 0.12) +
  scale_fill_manual(values = umbrella_pal, name = NULL,
                    guide = guide_legend(ncol = 2, keywidth = unit(0.25, "cm"),
                                        keyheight = unit(0.25, "cm"))) +
  scale_x_continuous(breaks = 1:5, labels = paste0("F", 0:4),
                     expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = c(0, 0), limits = c(0, 1)) +
  labs(title = "Non-hepatocyte remodeling",
       subtitle = "Renormalized; each stage sums to 100%",
       x = "Fibrosis stage", y = "% of non-hepatocyte cells") +
  theme_masld(base_size = 7) +
  theme(
    plot.title    = element_text(size = 7.5, face = "bold"),
    plot.subtitle = element_text(size = 5.5, color = "grey35"),
    legend.position = "right",
    legend.text   = element_text(size = 4.5),
    panel.grid    = element_blank()
  )

panel <- p_hep | p_stack

ggsave(OUT_PDF, panel,
       width  = 180 / 25.4,
       height = 55 / 25.4,
       units  = "in",
       device = cairo_pdf)
fwrite(stage_summ, file.path(DATA_DIR, "fig2_panel_celltype_proportion.csv"))
cat(sprintf("[saved] %s\n", OUT_PDF))
