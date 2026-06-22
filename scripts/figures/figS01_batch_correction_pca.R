#!/usr/bin/env Rscript
# KEY MESSAGE: dream's `(1|dataset)` random-intercept absorbs the cohort axis
# that dominates raw log2-CPM. PC1 captures 61.9% of variance in the raw
# matrix (cohort-driven) and collapses to 10.6% after limma::removeBatchEffect
# (the standard fixed-effect proxy for `(1|dataset)` with the disease design
# preserved). The same projection that mixes cohorts also retains the
# biological axes — disease, sex, fibrosis — confirming the mixed-model
# correction is not over-fitting.
#
# Layout: 4 metadata rows (cohort / disease / sex / fibrosis stage) × 2
# stages (raw / corrected) = 8 small panels in a compact grid.
#
# Output: figS01/panels/figS01_batch_correction_pca.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(edgeR)
  library(limma)
  library(matrixStats)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- file.path(FIGS01_DIR, "panels")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CONTROL_GRAY <- "#9E9E9E"

dge <- load_merged_dge()
stopifnot(!is.null(dge))

logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
rv     <- matrixStats::rowVars(logcpm)
top_idx <- order(rv, decreasing = TRUE)[seq_len(min(2000, length(rv)))]
logcpm_top <- logcpm[top_idx, ]

samp <- as.data.table(dge$samples, keep.rownames = "sample_id")

# Bring fibrosis_stage in from the unified metadata
meta_full <- fread(file.path(INT_META, "unified_metadata.csv"))
samp <- merge(samp,
              meta_full[, .(sample_id, fibrosis_stage)],
              by = "sample_id", all.x = TRUE, sort = FALSE)

# Mixed-effect batch removal (limma fixed-effect proxy)
design <- model.matrix(~ samp$group_binary)
logcpm_corrected <- limma::removeBatchEffect(
  logcpm_top, batch = samp$dataset, design = design
)

# PCA helper
do_pca <- function(mat, label) {
  pr  <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  pve <- 100 * (pr$sdev^2) / sum(pr$sdev^2)
  data.table(
    sample_id = colnames(mat),
    PC1 = pr$x[, 1], PC2 = pr$x[, 2],
    stage = label,
    pc1_var = pve[1], pc2_var = pve[2]
  )
}

pca_raw <- do_pca(logcpm_top,       "Raw")
pca_adj <- do_pca(logcpm_corrected, "Corrected")
pca <- rbind(pca_raw, pca_adj)
pca <- merge(pca, samp[, .(sample_id, dataset, group_binary, sex,
                           fibrosis_stage)],
             by = "sample_id")

# Short cohort labels (consistent with fig1_umap.R)
cohort_short <- c(
  GSE126848 = "GSE126848", GSE130970 = "GSE130970", GSE135251 = "GSE135251",
  GSE162694 = "GSE162694",   GSE167523 = "GSE167523", GSE174478 = "GSE174478",
  GSE193066 = "GSE193066", GSE213621 = "GSE213621", GSE240729 = "GSE240729",
  PRJNA512027 = "PRJNA512027"
)
pca[, cohort  := factor(cohort_short[dataset], levels = unname(cohort_short))]
pca[, disease := factor(group_binary, levels = c("Control", "Disease"))]
pca[sex == "" | is.na(sex), sex := NA_character_]
pca[, sex     := factor(sex, levels = c("F", "M"))]
pca[, fib     := factor(paste0("F", fibrosis_stage), levels = paste0("F", 0:4))]
pca[fibrosis_stage |> is.na(), fib := NA]
pca[, stage   := factor(stage, levels = c("Raw", "Corrected"))]

# Variance annotation (same across all metadata rows for a given stage)
ann <- pca[, .(pc1 = unique(pc1_var), pc2 = unique(pc2_var)), by = stage]
ann[, lab := sprintf("PC1=%.1f%%  PC2=%.1f%%", pc1, pc2)]

# Palettes
cohort_pal <- c(
  GSE126848="#1F77B4", GSE130970="#FF7F0E", GSE135251="#2CA02C", GSE162694="#D62728",
  GSE167523="#9467BD", GSE174478="#8C564B", GSE193066="#E377C2",
  GSE213621="#7F7F7F", GSE240729="#BCBD22", PRJNA512027="#17BECF"
)
disease_pal <- c(Control = CONTROL_GRAY, Disease = masld_colors$nash)
sex_pal     <- c(F = masld_colors$female, M = masld_colors$male)
fib_pal     <- fibrosis_stage_colors  # from publication_theme.R

# Compact base theme — small dots, no tick labels (PC coords are arbitrary).
compact_theme <- function() {
  theme_masld(base_size = 7) +
    theme(
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      axis.title      = element_text(size = 6.5),
      strip.text      = element_text(size = 7, face = "bold"),
      strip.background = element_blank(),
      plot.title      = element_blank(),
      legend.position  = "right",
      legend.title     = element_text(size = 6.5, face = "bold"),
      legend.text      = element_text(size = 5.8),
      legend.key.size  = unit(0.22, "cm"),
      legend.margin    = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0, "cm"),
      plot.margin     = margin(2, 4, 2, 4)
    )
}

# Build one row (Raw | Corrected) for a metadata variable
make_row <- function(color_var, palette, legend_title, ncol_legend = 1,
                     show_stage_strips = FALSE, add_var_annot = FALSE) {
  d <- pca[order(!is.na(get(color_var)))]  # NAs drawn first (beneath)
  p <- ggplot(d, aes(PC1, PC2, color = .data[[color_var]])) +
    rasterize_layer(geom_point(size = 0.28, alpha = 0.7, shape = 16)) +
    scale_color_manual(values = palette, name = legend_title,
                       drop = FALSE, na.value = "grey85") +
    facet_wrap(~ stage, ncol = 2, scales = "free") +
    labs(x = "PC1", y = "PC2") +
    guides(color = guide_legend(ncol = ncol_legend,
                                override.aes = list(size = 1.4, alpha = 1))) +
    compact_theme()
  if (!show_stage_strips) {
    p <- p + theme(strip.text = element_blank())
  }
  if (add_var_annot) {
    p <- p + geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab),
                       inherit.aes = FALSE, hjust = -0.08, vjust = 1.4,
                       size = 2.2, colour = "grey30")
  }
  p
}

p1 <- make_row("cohort",  cohort_pal,  "Cohort",  ncol_legend = 2,
               show_stage_strips = TRUE, add_var_annot = TRUE)
p2 <- make_row("disease", disease_pal, "Disease",   ncol_legend = 1)
p3 <- make_row("sex",     sex_pal,     "Sex",       ncol_legend = 1)
p4 <- make_row("fib",     fib_pal,     "Fibrosis",  ncol_legend = 1)

fig <- p1 / p2 / p3 / p4 +
  plot_layout(heights = c(1, 1, 1, 1)) +
  plot_annotation(
    title    = "PCA before vs after dataset batch correction",
    subtitle = sprintf(
      "Top 2,000 most-variable log2-CPM genes; %s samples · 10 cohorts · limma::removeBatchEffect proxy for dream `(1|dataset)`",
      format(ncol(logcpm_top), big.mark = ",")),
    theme    = theme(
      plot.title    = element_text(size = 9, face = "bold",
                                   family = "Helvetica"),
      plot.subtitle = element_text(size = 6.5, colour = "grey35",
                                   family = "Helvetica"),
      plot.margin   = margin(2, 2, 2, 2))
  )

out_pdf <- file.path(OUT_DIR, "figS01_batch_correction_pca.pdf")
ggsave(out_pdf, fig, width = 6.8, height = 8.4, device = cairo_pdf)

fwrite(pca[, .(sample_id, dataset, cohort, disease, sex, fib, stage,
               PC1, PC2)],
       file.path(FIGS01_DIR, "figS01_batch_correction_pca_data.csv"))

cat("=== Variance explained ===\n")
print(ann)
cat("\nWrote:\n  ", out_pdf, "\n", sep = "")
