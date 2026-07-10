#!/usr/bin/env Rscript
# fig_vae_staging.R
# VAE embedding visualizations for MASLD staging classifier
#
# Panels:
# (a) UMAP of 64-dim VAE embeddings colored by fibrosis stage (F0-F4)
# (b) UMAP colored by NAS group (4 groups)
# (c) UMAP colored by dataset (batch effect check)
# (d) VAE training curves (reconstruction + KL + contrastive losses)
# (e) Embedding vs raw feature AUROC scatter
# (f) Modality ablation lollipop (VAE dominance)
#
# Output: figures/supplementary/figS10_prediction/fig_vae_panels.pdf + individual panels

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(uwot)  # For UMAP
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts", "figures", "publication_theme.R"))
source(file.path(BASE, "scripts", "figures", "load_figure_data.R"))
SDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier")
OUTDIR <- FIGS10_DIR
dir.create(file.path(OUTDIR, "panels"), showWarnings = FALSE, recursive = TRUE)

cat("=== VAE Staging Visualizations ===\n")

# Color palettes
fibrosis_colors <- c("0" = "#2166ac", "1" = "#67a9cf", "2" = "#f7f7f7",
                      "3" = "#ef8a62", "4" = "#b2182b")
nas_colors <- c("0" = "#1b7837", "1" = "#7fbf7b", "2" = "#fc8d59", "3" = "#d73027")
# PRJNA512027 (Gerhard 2018) dropped from cohort presentation: L0/S0
# library-prep batch confounded with diagnosis. Samples filtered out of
# umap_dt below.
dataset_colors <- c(
  "GSE126848" = "#e41a1c", "GSE130970" = "#377eb8", "GSE135251" = "#4daf4a",
  "GSE162694" = "#984ea3", "GSE167523" = "#ff7f00", "GSE174478" = "#a65628",
  "GSE193066" = "#f781bf", "GSE213621" = "#999999", "GSE240729" = "#66c2a5"
)

theme_pub <- theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "plain", size = 10, margin = margin(b = 3)),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "plain"),
    legend.key.size = unit(0.35, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(3, 5, 3, 3)
  )

# ============================================================
# Load data
# ============================================================
cat("Loading embeddings...\n")
emb <- fread(file.path(SDIR, "embeddings_all_samples.csv"))
cat("  Embeddings:", nrow(emb), "samples x", ncol(emb), "columns\n")

# The first column should be sample_id or row index
# Detect which columns are embedding dimensions
emb_cols <- grep("^(dim_|V|emb_|[0-9]+$)", names(emb), value = TRUE)
if (length(emb_cols) == 0) {
  # Try: all numeric columns except sample_id
  emb_cols <- names(emb)[sapply(emb, is.numeric)]
}
cat("  Embedding dimensions:", length(emb_cols), "\n")

# Get sample IDs
id_col <- setdiff(names(emb), emb_cols)[1]
if (!is.null(id_col) && !is.na(id_col)) {
  sample_ids <- emb[[id_col]]
} else {
  sample_ids <- paste0("S", seq_len(nrow(emb)))
}

emb_mat <- as.matrix(emb[, ..emb_cols])

# Load metadata
meta <- fread(file.path(SDIR, "modeling_metadata.csv"))

# Match embeddings to metadata
meta_matched <- meta[match(sample_ids, sample_id)]
cat("  Matched:", sum(!is.na(meta_matched$sample_id)), "samples\n")

# ============================================================
# Compute UMAP
# ============================================================
cat("Computing UMAP...\n")
umap_res <- umap(emb_mat, n_neighbors = 30, min_dist = 0.3, n_components = 2,
                  metric = "euclidean", n_threads = 4, ret_model = FALSE)

umap_dt <- data.table(
  UMAP1 = umap_res[, 1],
  UMAP2 = umap_res[, 2],
  sample_id = sample_ids,
  fib_stage = meta_matched$fib_stage,
  nas_group4 = meta_matched$nas_group4,
  dataset = meta_matched$dataset,
  is_disease = meta_matched$is_disease
)
# Drop PRJNA512027 from cohort presentation (L0/S0 batch confound).
umap_dt <- umap_dt[dataset != "PRJNA512027"]

# Label fibrosis
umap_dt[, fib_label := fifelse(fib_stage >= 0, as.character(fib_stage), "NA")]
umap_dt[, nas_label := fifelse(nas_group4 >= 0, as.character(nas_group4), "NA")]

cat("  UMAP computed\n")

# ============================================================
# Panel (a): UMAP by fibrosis stage
# ============================================================
cat("Panel (a): UMAP by fibrosis stage\n")

umap_fib <- umap_dt[fib_label != "NA"]
umap_na <- umap_dt[fib_label == "NA"]

umap_theme <- theme_pub +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 7),
        plot.title = element_text(size = 9))

p_a <- ggplot() +
  geom_point(data = umap_na, aes(x = UMAP1, y = UMAP2),
             color = "gray90", size = 0.15, alpha = 0.4) +
  geom_point(data = umap_fib, aes(x = UMAP1, y = UMAP2, color = fib_label),
             size = 0.5, alpha = 0.65) +
  scale_color_manual(values = fibrosis_colors, name = "Fibrosis",
                     labels = c("F0", "F1", "F2", "F3", "F4")) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(title = "(a) Fibrosis stage", x = "UMAP 1", y = "UMAP 2") +
  umap_theme

# ============================================================
# Panel (b): UMAP by NAS group
# ============================================================
cat("Panel (b): UMAP by NAS group\n")

umap_nas <- umap_dt[nas_label != "NA"]
umap_nas_na <- umap_dt[nas_label == "NA"]

nas_labels <- c("0" = "Low (0-2)", "1" = "Moderate (3-4)",
                "2" = "High (5-6)", "3" = "Very High (7+)")

p_b <- ggplot() +
  geom_point(data = umap_nas_na, aes(x = UMAP1, y = UMAP2),
             color = "gray90", size = 0.15, alpha = 0.4) +
  geom_point(data = umap_nas, aes(x = UMAP1, y = UMAP2, color = nas_label),
             size = 0.5, alpha = 0.65) +
  scale_color_manual(values = nas_colors, name = "NAS", labels = nas_labels) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(title = "(b) NAS group", x = "UMAP 1", y = "UMAP 2") +
  umap_theme

# ============================================================
# Panel (c): UMAP by dataset (batch check)
# ============================================================
cat("Panel (c): UMAP by dataset\n")

p_c <- ggplot(umap_dt, aes(x = UMAP1, y = UMAP2, color = dataset)) +
  geom_point(size = 0.3, alpha = 0.5) +
  scale_color_manual(values = dataset_colors, name = "Cohort") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 2)) +
  labs(title = "(c) Cohort (batch check)", x = "UMAP 1", y = "UMAP 2") +
  umap_theme +
  theme(legend.text = element_text(size = 5), legend.key.size = unit(0.2, "cm"))

# ============================================================
# Panel (d): VAE training curves
# ============================================================
cat("Panel (d): Training curves\n")

curves <- fread(file.path(SDIR, "vae_training_curves.csv"))

# Melt to long format for plotting
loss_cols <- intersect(names(curves), c("recon_loss", "kl_loss", "contrastive_loss",
                                         "domain_loss", "total_loss",
                                         "reconstruction", "kl", "contrastive", "domain", "total"))
if (length(loss_cols) > 0) {
  curves_long <- melt(curves, id.vars = "epoch", measure.vars = loss_cols,
                       variable.name = "loss_type", value.name = "loss")
  curves_long <- curves_long[!is.na(loss)]

  # Clean labels
  curves_long[, loss_label := gsub("_loss$", "", loss_type)]
  curves_long[, loss_label := gsub("recon.*", "Reconstruction", loss_label)]
  curves_long[, loss_label := gsub("^kl$", "KL Divergence", loss_label)]
  curves_long[, loss_label := gsub("contrastive", "Contrastive", loss_label)]
  curves_long[, loss_label := gsub("domain", "Domain Adversarial", loss_label)]
  curves_long[, loss_label := gsub("total", "Total", loss_label)]

  loss_colors <- c("Reconstruction" = "#1f77b4", "KL Divergence" = "#ff7f0e",
                    "Contrastive" = "#2ca02c", "Domain Adversarial" = "#d62728",
                    "Total" = "#7f7f7f")

  p_d <- ggplot(curves_long, aes(x = epoch, y = loss, color = loss_label)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = loss_colors, name = "Loss") +
    scale_y_log10() +
    labs(title = "(d) VAE training convergence",
         x = "Epoch", y = "Loss (log scale)") +
    theme_pub +
    theme(legend.position = c(0.75, 0.8))
} else {
  p_d <- ggplot() + theme_void() + labs(title = "(d) Training curves [no data]")
}

# ============================================================
# Panel (e): Embedding vs Raw feature comparison
# ============================================================
cat("Panel (e): Embedding vs Raw\n")

emb_comp <- fread(file.path(SDIR, "embedding_vs_raw_comparison.csv"))

# Find the AUROC columns
auroc_cols <- names(emb_comp)
emb_auroc_col <- grep("emb|embedding", auroc_cols, value = TRUE, ignore.case = TRUE)[1]
raw_auroc_col <- grep("raw|gene", auroc_cols, value = TRUE, ignore.case = TRUE)[1]

if (!is.na(emb_auroc_col) && !is.na(raw_auroc_col)) {
  ec <- copy(emb_comp)
  setnames(ec, c(emb_auroc_col, raw_auroc_col), c("emb_auroc", "raw_auroc"),
           skip_absent = TRUE)
  ec <- ec[!is.na(emb_auroc) & !is.na(raw_auroc)]

  target_col <- intersect(names(ec), c("target", "task", "experiment"))[1]
  if (!is.na(target_col)) setnames(ec, target_col, "target", skip_absent = TRUE)

  p_e <- ggplot(ec, aes(x = raw_auroc, y = emb_auroc)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    geom_point(aes(color = target), size = 2.5, alpha = 0.8) +
    annotate("text", x = 0.55, y = 0.95, label = "Embedding\nbetter",
             size = 2.5, color = "#2ca02c", fontface = "italic") +
    annotate("text", x = 0.95, y = 0.55, label = "Raw\nbetter",
             size = 2.5, color = "#d62728", fontface = "italic") +
    coord_equal(xlim = c(0.4, 1), ylim = c(0.4, 1)) +
    labs(title = "(e) VAE embeddings vs raw features",
         x = "Raw gene AUROC", y = "VAE embedding AUROC",
         color = "Target") +
    theme_pub +
    theme(legend.position = "bottom", legend.text = element_text(size = 6))
} else {
  p_e <- ggplot() + theme_void() + labs(title = "(e) Embedding vs Raw [no data]")
}

# ============================================================
# Panel (f): Modality ablation (VAE dominance)
# ============================================================
cat("Panel (f): Modality ablation\n")

ablation <- fread(file.path(SDIR, "ablation_summary.csv"))
ablation <- ablation[order(mean_auroc)]

# Clean labels
ablation[, label := gsub("_", " + ", combination)]
ablation[, label := gsub("expr", "Expression", label)]
ablation[, label := gsub("vae", "VAE", label)]
ablation[, label := gsub("deconv", "Deconv", label)]

# Flag VAE-containing combinations
ablation[, has_vae := grepl("VAE|vae", combination, ignore.case = TRUE)]
ablation[, label := factor(label, levels = label)]

p_f <- ggplot(ablation, aes(x = mean_auroc, y = label)) +
  geom_segment(aes(xend = 0.5, yend = label), color = "gray80", linewidth = 0.4) +
  geom_point(aes(color = has_vae, size = has_vae)) +
  geom_text(aes(label = sprintf("%.3f", mean_auroc)), hjust = -0.3, size = 2.3) +
  scale_color_manual(values = c("FALSE" = "#1f77b4", "TRUE" = "#2ca02c"),
                     labels = c("Without VAE", "With VAE"), name = "") +
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 3.5), guide = "none") +
  scale_x_continuous(limits = c(0.5, 1.05)) +
  labs(title = "(f) Modality ablation — VAE dominance",
       x = "Mean QWK (Fibrosis ordinal, 6-fold LOCO)", y = "") +
  theme_pub +
  theme(legend.position = c(0.3, 0.9),
        axis.text.y = element_text(size = 6.5))

# ============================================================
# Compose 6-panel figure (3 rows x 2 cols)
# ============================================================
cat("\nComposing 6-panel VAE figure (3 UMAPs top, 3 panels bottom)...\n")

# Top row: 3 compact UMAPs side by side
# Bottom row: training curves, embedding vs raw, ablation
layout <- "
AABBCC
DDEEFF
"

composite <- p_a + p_b + p_c + p_d + p_e + p_f +
  plot_layout(design = layout, heights = c(1, 1.1)) +
  plot_annotation(
    title = "VAE Embedding Analysis for MASLD Disease Staging",
    subtitle = "64-dimensional latent space captures fibrosis biology more efficiently than 3,000 raw genes",
    theme = theme(
      plot.title = element_text(face = "plain", size = 12),
      plot.subtitle = element_text(size = 8.5, color = "gray30")
    )
  )

out_path <- file.path(OUTDIR, "fig_vae_panels.pdf")
ggsave(out_path, composite, width = 14, height = 8, device = cairo_pdf)
cat("Saved:", out_path, "\n")

# Save individual panels
ggsave(file.path(OUTDIR, "panels", "panel_vae_umap_fibrosis.pdf"), p_a, width = 4.5, height = 3.5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "panels", "panel_vae_umap_nas.pdf"), p_b, width = 4.5, height = 3.5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "panels", "panel_vae_umap_dataset.pdf"), p_c, width = 4.5, height = 3.5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "panels", "panel_vae_training.pdf"), p_d, width = 5, height = 4, device = cairo_pdf)
ggsave(file.path(OUTDIR, "panels", "panel_vae_vs_raw.pdf"), p_e, width = 5, height = 5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "panels", "panel_vae_ablation.pdf"), p_f, width = 6, height = 5, device = cairo_pdf)

cat("Individual panels saved to:", file.path(OUTDIR, "panels"), "\n")
cat("=== Done ===\n")
