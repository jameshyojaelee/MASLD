#!/usr/bin/env Rscript
# gen_fig2e_embeddable.R
# Generates fig2e_embeddable.pdf — fully vectorized (NO ggrastr rasterization)
# for direct editing in Illustrator. Uses the same data / layout as fig2e.pdf
# but skips rasterize_layer() so every element is an editable vector object.
#
# Output: figures/main/fig2_progression_sex/panels/fig2e_embeddable.pdf
# Does NOT overwrite fig2e.pdf (the rasterized compositor panel).
# ============================================================================

Sys.setenv(NO_RASTERIZE = "1")   # disables ggrastr in rasterize_layer()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(MASS)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR  <- file.path(FIG2_DIR, "panels")
ATLAS_UMAP <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/atlas_umap_for_fig2.csv.gz")
OUT_PDF    <- file.path(PANEL_DIR, "fig2e_embeddable.pdf")

# Style constants (mirror archive script)
BASE_SIZE <- 7
LBL_SIZE  <- 7 / ggplot2::.pt

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

# Cell-type umbrella mapping (mirrors archive script)
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

# Load UMAP data
message("[load] ", ATLAS_UMAP)
umap_df <- fread(ATLAS_UMAP)
n_total <- nrow(umap_df)
message(sprintf("[cells] %s", comma(n_total)))

umap_df[, cell_group := umbrella_map[cell_type]]

set.seed(42)
umap_sub <- if (n_total > 150000) umap_df[sample(.N, 150000)] else copy(umap_df)
group_present       <- intersect(names(umbrella_palette), unique(umap_df$cell_group))
umbrella_palette_use <- umbrella_palette[group_present]
umap_sub[, cell_group := factor(cell_group, levels = group_present)]
umap_sub <- umap_sub[sample(.N)]

# Disease-enrichment density (MASH / Healthy log2 ratio)
ALLOWED_PREP        <- c("unsorted", "nuclei")
EXCLUDE_KDE_DATASETS <- c("GSE136103", "Liver_Atlas")
umap_pf <- umap_df[
  preparation_method %in% ALLOWED_PREP &
  !dataset %in% EXCLUDE_KDE_DATASETS
]
healthy_df <- umap_pf[disease_stage_coarse == "Healthy"]
adv_df     <- umap_pf[disease_stage_coarse == "Steatohepatitis"]

pad  <- 0.5
xlim <- range(umap_df$umap_1) + c(-pad, pad)
ylim <- range(umap_df$umap_2) + c(-pad, pad)

per_sample_avg_density <- function(df, n = 200, min_cells = 100) {
  samples <- unique(df$sample)
  z_acc <- NULL; grid_x <- NULL; grid_y <- NULL; n_used <- 0L
  for (s in samples) {
    cells <- df[sample == s]
    if (nrow(cells) < min_cells) next
    kde <- MASS::kde2d(cells$umap_1, cells$umap_2, n = n, lims = c(xlim, ylim))
    if (is.null(z_acc)) {
      z_acc <- kde$z; grid_x <- kde$x; grid_y <- kde$y
    } else {
      z_acc <- z_acc + kde$z
    }
    n_used <- n_used + 1L
  }
  list(z = z_acc / n_used, x = grid_x, y = grid_y, n_samples = n_used)
}

kde_h <- per_sample_avg_density(healthy_df)
kde_a <- per_sample_avg_density(adv_df)
message(sprintf("KDE: Healthy n=%d, MASH n=%d", kde_h$n_samples, kde_a$n_samples))
message(sprintf("[KDE] Healthy donors used: %d; MASH donors used: %d",
                kde_h$n_samples, kde_a$n_samples))

joint  <- kde_h$z + kde_a$z
floor_ <- 0.001 * max(joint)
enrich <- log2((kde_a$z + floor_) / (kde_h$z + floor_))
enrich[joint < floor_] <- NA

enrich_df <- data.table(
  expand.grid(umap_1 = kde_a$x, umap_2 = kde_a$y),
  enrich = as.vector(enrich)
)
fill_lim <- quantile(abs(enrich_df$enrich), 0.99, na.rm = TRUE)
enrich_df[, enrich := pmin(pmax(enrich, -fill_lim), fill_lim)]

group_counts <- umap_df[, .N, by = cell_group]
group_keep   <- group_counts[N >= 2000, cell_group]
ct_centroids <- umap_df[cell_group %in% group_keep,
  .(umap_1 = median(umap_1), umap_2 = median(umap_2)), by = cell_group]
ct_centroids[, cell_group := factor(cell_group, levels = group_present)]

# Build plot — rasterize_layer() is a no-op due to NO_RASTERIZE=1
p2e <- ggplot() +
  rasterize_layer(geom_raster(data = enrich_df,
                              aes(x = umap_1, y = umap_2, fill = enrich),
                              interpolate = TRUE)) +
  rasterize_layer(geom_point(data = umap_sub,
                             aes(x = umap_1, y = umap_2, color = cell_group),
                             size = 0.09, alpha = 0.65,
                             stroke = 0, shape = 16)) +
  geom_text_repel(data = ct_centroids,
            aes(x = umap_1, y = umap_2, label = cell_group),
            size = LBL_SIZE, color = "black", fontface = "bold",
            bg.color = "white", bg.r = 0.12,
            min.segment.length = 0, segment.size = 0.2,
            segment.color = "grey40",
            box.padding = 0.35, point.padding = 0.1,
            max.overlaps = Inf, force = 2, seed = 42) +
  scale_fill_gradient2(low = "#7DA0CC", mid = "white", high = "#E08A82",
                       midpoint = 0,
                       limits = c(-fill_lim, fill_lim),
                       na.value = "white",
                       name = expression(log[2] * " density ratio\n(MASH / Healthy)"),
                       guide = guide_colorbar(barwidth = 0.3, barheight = 3, order = 1)) +
  scale_color_manual(values = umbrella_palette_use, na.value = "grey80",
                     drop = FALSE,
                     guide = guide_legend(
                       override.aes = list(size = 1.5, alpha = 1, shape = 16),
                       ncol = 1, keyheight = unit(7, "pt"), order = 2)) +
  labs(title = sprintf("e  scRNA atlas (n=%s cells)", comma(n_total)),
       x = "UMAP 1", y = "UMAP 2", color = NULL) +
  theme_fig2() +
  theme(axis.text  = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.25, "cm"))

message("[save] ", OUT_PDF)
save_fig(p2e, OUT_PDF, width = fig_full_width * 0.65, height = 3.6)
message(sprintf("[done] %s (%.0f KB)", OUT_PDF,
                file.info(OUT_PDF)$size / 1024))
