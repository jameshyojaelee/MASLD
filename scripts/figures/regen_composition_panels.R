#!/usr/bin/env Rscript
# regen_composition_panels.R - regenerate celltype_composition_shift, hep_subtype_composition with F2-cellular-driver narrative.
#
# 2f: cell-type composition shift across Healthy -> Steatosis -> Steatohepatitis,
#     with Macrophages/Monocytes/Fibroblasts/Hepatocytes highlighted and other
#     umbrellas faded as a grey background. F2 transition annotated between
#     Steatosis and Steatohepatitis.
# 2g: hepatocyte meta-subtype composition across the same 3 stages. Disease-
#     Progressor wedge emphasized with the Steatosis -> Steatohepatitis step
#     annotated (the F2 switch in cellular form).
#
# After verification, mirror this block back into fig2_progression_sex.R so the
# main script stays the source of truth.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_SIZE <- 7
LBL_PT    <- 7
LBL_SIZE  <- LBL_PT / ggplot2::.pt

theme_fig2 <- function() {
  theme_masld(base_size = BASE_SIZE) +
    theme(
      plot.title    = element_text(size = BASE_SIZE),
      axis.title    = element_text(size = BASE_SIZE),
      axis.text     = element_text(size = BASE_SIZE),
      legend.title  = element_text(size = BASE_SIZE),
      legend.text   = element_text(size = BASE_SIZE),
      strip.text    = element_text(size = BASE_SIZE),
      plot.subtitle = element_blank(),
      plot.margin   = margin(3, 3, 3, 3)
    )
}

ATLAS_UMAP <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/atlas_umap_for_fig2.csv.gz")
HEP_FILE <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/f2_integration/stage_proportions_per_sample.csv")

umbrella_map <- c(
  "Hepatocytes"             = "Hepatocytes",
  "Cholangiocytes"          = "Cholangiocytes",
  "Endothelial cells"       = "Endothelial cells",
  "Fibroblasts"             = "Fibroblasts",
  "Macrophages"             = "Macrophages",
  "Mono+mono derived cells" = "Monocytes",
  "cDC1s"                   = "Dendritic cells",
  "cDC2s"                   = "Dendritic cells",
  "pDCs"                    = "Dendritic cells",
  "Mig.cDCs"                = "Dendritic cells",
  "Neutrophils"             = "Granulocytes",
  "Basophils"               = "Granulocytes",
  "T cells"                 = "T cells",
  "B cells"                 = "B cells",
  "Plasma cells"            = "Plasma cells",
  "Resident NK"             = "NK cells",
  "Circulating NK/NKT"      = "NK cells"
)
umbrella_palette <- c(
  "Hepatocytes"       = "#1B5E20",
  "Cholangiocytes"    = "#00BCD4",
  "Endothelial cells" = "#7CB342",
  "Fibroblasts"       = "#FF6F00",
  "Macrophages"       = "#212121",
  "Monocytes"         = "#EC407A",
  "Dendritic cells"   = "#6A1B9A",
  "Granulocytes"      = "#F9A825",
  "T cells"           = "#1A237E",
  "B cells"           = "#5D4037",
  "Plasma cells"      = "#8D6E63",
  "NK cells"          = "#00897B"
)
umbrella_order <- c(
  "Hepatocytes", "Cholangiocytes", "Endothelial cells", "Fibroblasts",
  "Macrophages", "Monocytes", "Dendritic cells", "Granulocytes",
  "T cells", "B cells", "Plasma cells", "NK cells"
)

# ALLOWED_PREP: unsorted only. Steatosis + Steatohepatitis exist ONLY as
# unsorted in this atlas; mixing in Healthy nuclei would create a prep-driven
# baseline bias.
ALLOWED_PREP   <- c("unsorted")
STAGE_LEVELS   <- c("Healthy", "Steatosis", "Steatohepatitis")
# F2 movers in unsorted prep (3.5x, 4x, 1.5x rises at Steatohepatitis):
#   Fibroblasts        - stellate-myofibroblast expansion (matches NMF P6)
#   Cholangiocytes     - ductular reaction
#   Endothelial cells  - LSEC capillarization
# Macrophages/Monocytes are flat in unsorted prep (cd45+ sort needed to see
# immune compartment shifts). Hepatocytes are non-monotonic. Both kept as
# grey background context lines.
HIGHLIGHT_CT   <- c("Fibroblasts", "Cholangiocytes", "Endothelial cells")
MIN_CELLS_PER_SAMPLE <- 200
MIN_HEP        <- 20

# ============================================================================
# Panel 2f - cell-type composition shift across F2 transition
# ============================================================================
message("Loading atlas UMAP...")
umap_df <- fread(ATLAS_UMAP)
umap_df[, cell_group := umbrella_map[cell_type]]

comp <- umap_df[
  preparation_method %in% ALLOWED_PREP &
  disease_stage_coarse %in% STAGE_LEVELS &
  !is.na(cell_group)
]
message(sprintf("After prep+stage filter: %s cells / %d samples",
  format(nrow(comp), big.mark = ","), uniqueN(comp$sample)))

# Per-sample cell-group counts; sample is the unit of replication.
per_sample <- comp[, .N, by = .(sample, disease_stage_coarse, cell_group)]
sample_totals <- per_sample[, .(total = sum(N)),
                            by = .(sample, disease_stage_coarse)]
per_sample <- merge(per_sample, sample_totals,
                    by = c("sample", "disease_stage_coarse"))
per_sample[, prop := N / total]

# Fill in zeros for sample x cell_group combos missing from the count table.
sample_stage <- unique(sample_totals[, .(sample, disease_stage_coarse, total)])
full_grid <- CJ(sample = sample_stage$sample, cell_group = umbrella_order)
full_grid <- merge(full_grid, sample_stage, by = "sample")
per_sample_full <- merge(full_grid,
                         per_sample[, .(sample, cell_group, N, prop)],
                         by = c("sample", "cell_group"), all.x = TRUE)
per_sample_full[is.na(prop), prop := 0]
per_sample_full[is.na(N),    N    := 0]

ok_samples <- sample_stage[total >= MIN_CELLS_PER_SAMPLE, sample]
per_sample_full <- per_sample_full[sample %in% ok_samples]
message(sprintf("Samples with >= %d cells: %d (of %d total)",
                MIN_CELLS_PER_SAMPLE, length(ok_samples), nrow(sample_stage)))

# Per-stage mean proportion (sample-as-unit; equal weight per sample).
stage_summary <- per_sample_full[, .(
  mean_prop = mean(prop, na.rm = TRUE),
  sd_prop   = sd(prop, na.rm = TRUE),
  n_samples = .N
), by = .(disease_stage_coarse, cell_group)]
stage_summary[, se_prop := sd_prop / sqrt(n_samples)]
stage_summary[, disease_stage_coarse := factor(disease_stage_coarse,
                                                levels = STAGE_LEVELS)]
stage_summary[, stage_num := as.integer(disease_stage_coarse)]
stage_summary[, cell_group := factor(cell_group, levels = umbrella_order)]
stage_summary[, highlight := cell_group %in% HIGHLIGHT_CT]

# Print sanity-check table for ALL cell groups so we can pick movers.
all_table <- dcast(stage_summary,
                  cell_group ~ disease_stage_coarse, value.var = "mean_prop")
message("All cell-type proportions per stage:")
print(all_table, digits = 3)

# Per-stage sample sizes (for the figure caption).
n_per_stage <- per_sample_full[, .(n = uniqueN(sample)),
                                by = disease_stage_coarse]
message("Samples per stage (after >=200-cell filter):")
print(n_per_stage)

p2f <- ggplot() +
  # F2 transition guide between Steatosis (x=2) and Steatohepatitis (x=3).
  geom_vline(xintercept = 2.5, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  # Background grey lines for non-highlighted umbrella types.
  geom_line(data = stage_summary[!(cell_group %in% HIGHLIGHT_CT)],
            aes(x = stage_num, y = mean_prop * 100, group = cell_group),
            color = "grey78", linewidth = 0.4) +
  geom_point(data = stage_summary[!(cell_group %in% HIGHLIGHT_CT)],
             aes(x = stage_num, y = mean_prop * 100),
             color = "grey78", size = 0.8, shape = 16) +
  # Highlighted cell types: error bars + line + points.
  geom_errorbar(data = stage_summary[cell_group %in% HIGHLIGHT_CT],
                aes(x = stage_num,
                    ymin = pmax(0, (mean_prop - se_prop) * 100),
                    ymax = (mean_prop + se_prop) * 100,
                    color = cell_group),
                width = 0.10, linewidth = 0.4) +
  geom_line(data = stage_summary[cell_group %in% HIGHLIGHT_CT],
            aes(x = stage_num, y = mean_prop * 100,
                group = cell_group, color = cell_group),
            linewidth = 0.85) +
  geom_point(data = stage_summary[cell_group %in% HIGHLIGHT_CT],
             aes(x = stage_num, y = mean_prop * 100, color = cell_group),
             size = 1.6, shape = 16) +
  scale_x_continuous(breaks = seq_along(STAGE_LEVELS), labels = STAGE_LEVELS,
                     expand = expansion(mult = 0.06)) +
  scale_color_manual(values = umbrella_palette[HIGHLIGHT_CT],
                     breaks = HIGHLIGHT_CT, name = NULL,
                     guide = guide_legend(
                       override.aes = list(linewidth = 0.85, size = 1.6),
                       keyheight = unit(7, "pt"))) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(title = "f  Stromal & epithelial expansion at F2 switch",
       x = NULL, y = "Mean % cells per sample") +
  theme_fig2() +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        panel.grid.major.x = element_blank())

save_fig(p2f, file.path(PANEL_DIR, "celltype_composition_shift.pdf"),
         width = fig_half_width * 1.05, height = 3.4)
fwrite(stage_summary, file.path(PANEL_DIR, "celltype_composition_shift_data.csv"))
message(sprintf("Saved %s", file.path(PANEL_DIR, "celltype_composition_shift.pdf")))

# ============================================================================
# Panel 2g - hepatocyte progressor expansion at F2
# ============================================================================
message("Loading hepatocyte stage proportions...")
hep <- fread(HEP_FILE)
hep[, proportion := as.numeric(proportion)]
hep <- hep[disease_stage_coarse %in% STAGE_LEVELS]

# Per-sample hepatocyte total is filled only on rows with N > 0; recover it
# and re-attach to all subtype rows for the same sample.
sample_hep_totals <- hep[!is.na(total) & total > 0,
                         .(sample_total = max(total)), by = sample]
hep <- merge(hep, sample_hep_totals, by = "sample", all.x = TRUE)
hep <- hep[!is.na(sample_total) & sample_total >= MIN_HEP]
message(sprintf("Hepatocyte samples after >=%d-hep filter: %d",
                MIN_HEP, uniqueN(hep$sample)))

meta_levels <- c("Healthy", "Neutral", "Disease-Neutral",
                 "Disease-Associated", "Disease-Progressor")
hep[, meta_subtype := factor(meta_subtype, levels = meta_levels)]
hep[, disease_stage_coarse := factor(disease_stage_coarse,
                                     levels = STAGE_LEVELS)]

hep_summary <- hep[, .(mean_prop = mean(proportion, na.rm = TRUE),
                       n_samples = .N),
                   by = .(disease_stage_coarse, meta_subtype)]
# Renormalize so each stage sums to 1 (corrects for averaging artifacts).
hep_summary[, mean_prop := mean_prop / sum(mean_prop),
            by = disease_stage_coarse]
hep_summary[, stage_num := as.integer(disease_stage_coarse)]

prog_table <- hep_summary[meta_subtype == "Disease-Progressor"]
prog_st  <- prog_table[disease_stage_coarse == "Steatosis",       mean_prop]
prog_sth <- prog_table[disease_stage_coarse == "Steatohepatitis", mean_prop]
fold_change <- prog_sth / prog_st
message(sprintf("Progressor: Steatosis %.3f -> Steatohepatitis %.3f (%.1fx)",
                prog_st, prog_sth, fold_change))

hep_pal <- c(
  "Healthy"             = "#2E7D32",
  "Neutral"             = "#A5D6A7",
  "Disease-Neutral"     = "#FFE082",
  "Disease-Associated"  = "#FB8C00",
  "Disease-Progressor"  = "#B71C1C"
)

p2g <- ggplot(hep_summary,
              aes(x = stage_num, y = mean_prop * 100, fill = meta_subtype)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 2.5, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  # Progressor expansion callout: arrow from Steatosis bar tip to Steatohep bar.
  annotate("curve", x = 2.05, xend = 2.95, y = 8, yend = 30,
           curvature = -0.25,
           arrow = arrow(length = unit(0.06, "in"), type = "closed"),
           color = "grey20", linewidth = 0.4) +
  annotate("text", x = 2.5, y = 22,
           label = sprintf("Progressor: %.1f%% → %.1f%% (%.0f×)",
                           prog_st * 100, prog_sth * 100, fold_change),
           size = LBL_SIZE * 0.9, color = "grey15",
           hjust = 0.5, vjust = 1, fontface = "italic") +
  scale_x_continuous(breaks = seq_along(STAGE_LEVELS), labels = STAGE_LEVELS,
                     expand = expansion(mult = 0.05)) +
  scale_fill_manual(values = hep_pal, breaks = meta_levels,
                    name = "Hepatocyte\nsubtype",
                    guide = guide_legend(keyheight = unit(7, "pt"))) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0, 103),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "g  Hepatocyte progressor expansion at F2",
       x = NULL, y = "% hepatocytes (sample-mean)") +
  theme_fig2() +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        panel.grid.major.x = element_blank())

save_fig(p2g, file.path(PANEL_DIR, "hep_subtype_composition.pdf"),
         width = fig_half_width * 1.15, height = 3.4)
fwrite(hep_summary, file.path(PANEL_DIR, "hep_subtype_composition_data.csv"))
message(sprintf("Saved %s", file.path(PANEL_DIR, "hep_subtype_composition.pdf")))

message("Done.")
