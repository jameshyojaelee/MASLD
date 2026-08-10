#!/usr/bin/env Rscript
# figS05_scatac_umap_triad_integrated.R
# ==============================================================================
# Joint co-embedding UMAP of GSE244832 + GSE281367 snATAC (integrated), companion
# to the GSE244832-only figS05_scatac_umap_triad. Three panels:
#   cell_type | cohort | condition
# The COHORT panel is the point: it shows whether the two datasets occupy the same
# cell-state manifold after harmony(donor) batch correction -- the cell-level test
# of the same batch-correctability the DA framing relies on. (F_stage panel dropped:
# GSE281367 is binary MASH/NORMAL with no Kleiner F-stage.)
# Input : Analysis/ATAC/Human_External/results/coembed/umap_export_joint.tsv.gz  (from 30_coembed_umap.py)
# Output: figures/supplementary/figS05_epigenomic_spatial/figS05_scatac_umap_triad_integrated.pdf
# Env: rnaseq.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(ggrastr) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

UMAP_TSV <- file.path(BASE, "Analysis/ATAC/Human_External/results/coembed/umap_export_joint.tsv.gz")
OUT_PDF  <- file.path(FIGS05_DIR, "figS05_scatac_umap_triad_integrated.pdf")
stopifnot(file.exists(UMAP_TSV))

dat <- fread(UMAP_TSV)
# map raw SnapATAC2 marker labels -> atlas display groups (match fig1 triad aesthetic)
ct_map <- c(Hepatocyte="Hepatocytes", Stellate_Cell="Fibroblasts", Macrophage="Macrophages",
            Kupffer_Cell="Macrophages", Endothelial="Endothelial", LSEC="Endothelial",
            Cholangiocyte="Cholangiocytes", NK_T_Cell="Lymphocytes", B_Cell="Lymphocytes",
            Plasma_Cell="Lymphocytes")
dat[, ct := ifelse(cell_type %in% names(ct_map), ct_map[cell_type], "Low conf.")]
ct_colors <- c(Hepatocytes="#C9265E", Fibroblasts="#FF9800", Macrophages="#9C27B0",
               Endothelial="#1565C0", Cholangiocytes="#4CAF50", Lymphocytes="#9E9E9E",
               "Low conf."="#E0E0E0")
cohort_colors <- c(GSE244832="#2166AC", GSE281367="#B2182B")   # blue vs red — clearly distinct
cond_colors   <- c(NORMAL="#9E9E9E", MASL="#1565C0", MASH="#C9265E")

# subsample for vector size (keep ~60k), rasterize points
n <- nrow(dat); idx <- sample.int(n, min(60000, n)); d <- dat[idx]
message(sprintf("[load] %d cells -> %d plotted", n, nrow(d)))
d[, ct := factor(ct, levels = names(ct_colors))]
d[, cohort := factor(cohort, levels = names(cohort_colors))]
d[, condition := factor(condition, levels = names(cond_colors))]

base_umap <- function(mapping, cols, title) {
  ggplot(d, aes(UMAP1, UMAP2)) +
    rasterise(geom_point(mapping, size = 0.15, alpha = 0.5, stroke = 0), dpi = 400) +
    scale_color_manual(values = cols, name = NULL, drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 1.3, alpha = 1), ncol = 1)) +
    labs(x = "UMAP1", y = "UMAP2", title = title) +
    theme_masld(6) +
    theme(legend.position = "right", legend.key.size = unit(0.32, "lines"),
          axis.text = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(size = 6, face = "plain", hjust = 0))
}
p1 <- base_umap(aes(color = ct),        ct_colors,     "Cell type")
p2 <- base_umap(aes(color = cohort),    cohort_colors, "Cohort")
p3 <- base_umap(aes(color = condition), cond_colors,   "Condition")

p <- p1 + p2 + p3 + plot_layout(nrow = 1)
ggsave(OUT_PDF, p, width = 7.4, height = 2.5, useDingbats = FALSE)
message(sprintf("[wrote] %s", OUT_PDF))
message("[caption] Joint snATAC co-embedding of GSE244832 (existing) and GSE281367 (new) after ",
        "shared-feature LSI + harmony(donor). Cells colour by cell type (left), cohort (middle), ",
        "and disease condition (right). Co-mingling of the two cohorts within each cell-type region ",
        "(middle panel) indicates the cross-cohort batch is a correctable offset; per-cohort ",
        "differential accessibility is nonetheless analysed separately per cohort and treated ",
        "as supplementary corroboration (not pooled, not headline replication). This UMAP is ",
        "illustrative of shared cell states, not an inferential DA result.")
