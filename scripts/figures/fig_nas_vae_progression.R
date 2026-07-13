#!/usr/bin/env Rscript
# fig_nas_vae_progression.R
# Supplementary Figure: NAS-VAE Embedding & Disease Progression
#
# 6-panel layout (3 rows x 2 cols, width=11, height=13):
#   (a) NAS-VAE UMAP colored by NAS score (0-8 gradient)
#   (b) NAS-VAE UMAP colored by fibrosis stage (F0-F4)
#   (c) NAS-VAE UMAP colored by dataset (9 cohorts, batch check;
#       PRJNA512027 excluded for L0/S0 confound)
#   (d) Disease progression velocity field (arrows on UMAP)
#   (e) Bifurcation analysis (divergence score per stage)
#   (f) NAS prediction comparison (all models, horizontal bars)
#
# Output: figures/supplementary/figS10_prediction/fig_nas_vae_progression.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(uwot)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

SDIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier")
OUTDIR <- FIGS10_DIR
dir.create(file.path(OUTDIR, "panels"), showWarnings = FALSE, recursive = TRUE)

cat("=== NAS-VAE Embedding & Disease Progression Figure (6 panels) ===\n")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cat("Loading embeddings and metadata...\n")
emb <- fread(file.path(SDIR, "nas_embeddings_all_samples.csv"))
meta <- fread(file.path(SDIR, "modeling_metadata.csv"))

cat("  Embeddings:", nrow(emb), "samples x", ncol(emb) - 1, "latent dims\n")

# Extract 64-dim embedding matrix (z0..z63)
emb_cols <- grep("^z[0-9]+$", names(emb), value = TRUE)
emb_mat <- as.matrix(emb[, ..emb_cols])
sample_ids <- emb$sample_id
cat("  Detected", length(emb_cols), "embedding dimensions\n")

# Match metadata to embeddings
meta_idx <- match(sample_ids, meta$sample_id)
meta_matched <- meta[meta_idx]

# ---------------------------------------------------------------------------
# Compute UMAP from 64-dim embeddings
# ---------------------------------------------------------------------------
cat("Computing UMAP (n_neighbors=30, min_dist=0.3)...\n")
umap_res <- umap(emb_mat, n_neighbors = 30, min_dist = 0.3, n_components = 2,
                  metric = "euclidean", n_threads = 4, ret_model = FALSE)

umap_dt <- data.table(
  UMAP1     = umap_res[, 1],
  UMAP2     = umap_res[, 2],
  sample_id = sample_ids,
  nas_score = meta_matched$nas_score,
  fib_stage = meta_matched$fib_stage,
  dataset   = meta_matched$dataset
)
# Drop PRJNA512027 from cohort presentation (L0/S0 batch confound).
umap_dt <- umap_dt[dataset != "PRJNA512027"]

# Parse NAS as integer (original metadata: -1 = missing)
umap_dt[, nas_int := as.integer(nas_score)]
umap_dt[nas_int < 0 | is.na(nas_int), nas_int := NA_integer_]

# Parse fibrosis as integer
umap_dt[, fib_int := as.integer(fib_stage)]
umap_dt[fib_int < 0 | is.na(fib_int), fib_int := NA_integer_]

n_nas  <- sum(!is.na(umap_dt$nas_int))
n_fib  <- sum(!is.na(umap_dt$fib_int))
cat(sprintf("  Samples with NAS: %d; with fibrosis: %d; total: %d\n",
            n_nas, n_fib, nrow(umap_dt)))

# ---------------------------------------------------------------------------
# UMAP theme (minimal, no tick labels)
# ---------------------------------------------------------------------------
umap_theme <- theme_masld(base_size = 7) +
  theme(
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.title      = element_text(size = 6),
    legend.position = "right"
  )

# ---------------------------------------------------------------------------
# Dataset palette (9 cohorts; PRJNA512027 excluded for L0/S0 confound)
# ---------------------------------------------------------------------------
dataset_pal <- c(
  "GSE126848"   = "#e41a1c", "GSE130970"  = "#377eb8", "GSE135251" = "#4daf4a",
  "GSE162694"   = "#984ea3", "GSE167523"  = "#ff7f00", "GSE174478" = "#a65628",
  "GSE193066"   = "#f781bf", "GSE213621"  = "#999999", "GSE240729" = "#66c2a5"
)

# ============================================================
# Panel (a): UMAP colored by NAS score (0-8), blue-to-red
# ============================================================
cat("Panel (a): UMAP by NAS score (0-8 gradient)\n")

umap_nas     <- umap_dt[!is.na(nas_int)]
umap_nas_na  <- umap_dt[is.na(nas_int)]

p_a <- ggplot() +
  geom_point(data = umap_nas_na, aes(x = UMAP1, y = UMAP2),
             color = "gray85", size = 0.2, alpha = 0.35) +
  geom_point(data = umap_nas[order(nas_int)],
             aes(x = UMAP1, y = UMAP2, color = nas_int),
             size = 0.5, alpha = 0.7) +
  scale_color_gradient(low = "#2166ac", high = "#b2182b",
                       limits = c(0, 8), breaks = seq(0, 8, 2),
                       name = "NAS") +
  guides(color = guide_colorbar(barwidth = unit(0.35, "cm"),
                                 barheight = unit(2.5, "cm"))) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  umap_theme

# ============================================================
# Panel (b): UMAP colored by fibrosis stage (F0-F4)
# ============================================================
cat("Panel (b): UMAP by fibrosis stage (F0-F4)\n")

umap_fib    <- umap_dt[!is.na(fib_int)]
umap_fib_na <- umap_dt[is.na(fib_int)]

# Factor for legend ordering
umap_fib[, fib_label := factor(paste0("F", fib_int),
                                levels = c("F0", "F1", "F2", "F3", "F4"))]

fib_pal <- c("F0" = "#E3F2FD", "F1" = "#90CAF9", "F2" = "#42A5F5",
             "F3" = "#1565C0", "F4" = "#0D47A1")

p_b <- ggplot() +
  geom_point(data = umap_fib_na, aes(x = UMAP1, y = UMAP2),
             color = "gray85", size = 0.2, alpha = 0.35) +
  geom_point(data = umap_fib, aes(x = UMAP1, y = UMAP2, color = fib_label),
             size = 0.5, alpha = 0.7) +
  scale_color_manual(values = fib_pal, name = "Fibrosis") +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  umap_theme

# ============================================================
# Panel (c): UMAP colored by dataset (batch check)
# ============================================================
cat("Panel (c): UMAP by dataset (batch check)\n")

p_c <- ggplot(umap_dt, aes(x = UMAP1, y = UMAP2, color = dataset)) +
  geom_point(size = 0.3, alpha = 0.55) +
  scale_color_manual(values = dataset_pal, name = "Cohort") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1),
                               ncol = 2)) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  umap_theme +
  theme(legend.text     = element_text(size = 6),
        legend.title    = element_text(size = 6),
        legend.key.size = unit(0.22, "cm"))

# ============================================================
# Panel (d): Disease progression velocity field
# ============================================================
cat("Panel (d): Progression velocity field\n")

prog <- fread(file.path(SDIR, "progression_umap_data.csv"))

# Use NAS mode for the velocity field
prog_nas <- prog[mode == "nas"]
cat(sprintf("  Progression NAS mode: %d samples\n", nrow(prog_nas)))

# Merge NAS labels from metadata for coloring
prog_nas <- merge(prog_nas, meta[, .(sample_id, nas_score_meta = nas_score)],
                  by = "sample_id", all.x = TRUE)
prog_nas[, nas_val := as.integer(nas_score_meta)]
prog_nas[nas_val < 0 | is.na(nas_val), nas_val := NA_integer_]

# Separate samples with and without velocity
prog_vel  <- prog_nas[velocity_magnitude > 0]
prog_stat <- prog_nas[velocity_magnitude == 0]

# Scale arrows for visibility: normalize to fraction of UMAP range
umap_range <- max(diff(range(prog_nas$umap_1)), diff(range(prog_nas$umap_2)))
arrow_scale <- umap_range * 0.06 / max(abs(c(prog_vel$arrow_u, prog_vel$arrow_v)), na.rm = TRUE)

p_d <- ggplot() +
  # Static samples (terminal stages)
  geom_point(data = prog_stat, aes(x = umap_1, y = umap_2, color = nas_val),
             size = 0.4, alpha = 0.4) +
  # Mobile samples (with velocity)
  geom_point(data = prog_vel, aes(x = umap_1, y = umap_2, color = nas_val),
             size = 0.5, alpha = 0.6) +
  # Velocity arrows
  geom_segment(data = prog_vel,
               aes(x = umap_1, y = umap_2,
                   xend = umap_1 + arrow_u * arrow_scale,
                   yend = umap_2 + arrow_v * arrow_scale),
               arrow = arrow(length = unit(0.06, "cm"), type = "closed"),
               linewidth = 0.25, alpha = 0.5, color = "gray25") +
  scale_color_gradient(low = "#2166ac", high = "#b2182b",
                       limits = c(0, 8), breaks = seq(0, 8, 2),
                       na.value = "gray80", name = "NAS") +
  guides(color = guide_colorbar(barwidth = unit(0.35, "cm"),
                                 barheight = unit(2.5, "cm"))) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  umap_theme

# ============================================================
# Panel (e): Bifurcation analysis
# ============================================================
cat("Panel (e): Bifurcation analysis\n")

bif <- fread(file.path(SDIR, "progression_bifurcation.csv"))

# Select NAS and fibrosis modes
bif_plot <- bif[mode %in% c("nas", "fib")]
bif_plot[, stage_axis := fifelse(mode == "nas", "NAS", "Fibrosis")]
bif_plot[, stage_label := tp_label]

# Create combined factor for ordering: NAS stages first, then Fibrosis
bif_plot[, plot_order := seq_len(.N)]
bif_plot[, stage_label := factor(stage_label, levels = stage_label)]

# Color by bifurcation status
bif_plot[, bif_status := fifelse(is_bifurcation == TRUE | divergence_score >= 0.99,
                                  "Terminal", "Transitional")]

p_e <- ggplot(bif_plot, aes(x = stage_label, y = divergence_score, fill = bif_status)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", divergence_score)),
            vjust = -0.4, size = GEOM_TEXT_6PT, fontface = "plain") +
  # Annotate sample counts
  geom_text(aes(label = paste0("n=", n_samples), y = 0.03),
            size = GEOM_TEXT_6PT, color = "white", fontface = "plain") +
  scale_fill_manual(
    values = c("Terminal" = "#b2182b", "Transitional" = "#42A5F5"),
    name = ""
  ) +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
  facet_wrap(~ stage_axis, scales = "free_x", nrow = 1) +
  labs(x = "", y = "Divergence score") +
  theme_masld(base_size = 7) +
  theme(
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 6),
    strip.text       = element_text(size = 6, face = "plain"),
    legend.position  = "bottom",
    legend.key.size  = unit(0.3, "cm")
  )

# ============================================================
# Panel (f): NAS prediction comparison — horizontal bars
# ============================================================
cat("Panel (f): Model comparison (NAS metrics)\n")

fmc <- fread(file.path(SDIR, "foundation_model_comparison.csv"))

# Focus on NAS-related metrics: NAS 9-class QWK and NAS>=5 AUROC
# Extract NAS 9-class QWK rows
qwk_rows <- fmc[grepl("nas_9class|NAS 9-class", target, ignore.case = TRUE) &
                  !is.na(metric_qwk) & metric_qwk > 0]
# Deduplicate: keep the row with maximum QWK per model
qwk_rows <- qwk_rows[, .SD[which.max(metric_qwk)], by = model]
qwk_rows[, metric := "NAS 9-class QWK"]
qwk_rows[, value := metric_qwk]
qwk_rows[, value_sd := metric_qwk_sd]

# Extract NAS>=5 AUROC rows
auroc_rows <- fmc[grepl("nas_ge5|NAS>=5", target, ignore.case = TRUE) &
                    !is.na(metric_auroc) & metric_auroc > 0]
auroc_rows <- auroc_rows[, .SD[which.max(metric_auroc)], by = model]
auroc_rows[, metric := "NAS>=5 AUROC"]
auroc_rows[, value := metric_auroc]
auroc_rows[, value_sd := metric_auroc_sd]

comp <- rbind(
  qwk_rows[, .(model, metric, value, value_sd)],
  auroc_rows[, .(model, metric, value, value_sd)]
)

# Shorten model names
comp[, model_short := gsub(" \\+ .*", "", model)]
comp[, model_short := gsub("V3 Elastic Net", "Elastic Net (V3)", model_short)]
comp[, model_short := gsub("NAS-VAE \\+ SVM", "NAS-VAE + SVM", model)]
comp[, model_short := gsub("NAS-VAE \\+ RandomForest", "NAS-VAE + RF", model_short)]
comp[, model_short := gsub("Transfer \\(B_scratch\\)", "Transfer (scratch)", model_short)]
comp[, model_short := gsub("Transfer \\(B_finetune\\)", "Transfer (fine-tune)", model_short)]

# Flag NAS-VAE as best
comp[, is_vae := grepl("NAS-VAE", model, ignore.case = TRUE)]

# Order by value within each metric
comp[, model_short := reorder(model_short, value)]

p_f <- ggplot(comp, aes(x = value, y = model_short, fill = is_vae)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(xmin = pmax(0, value - fifelse(is.na(value_sd), 0, value_sd)),
                     xmax = pmin(1, value + fifelse(is.na(value_sd), 0, value_sd))),
                 width = 0.2, linewidth = 0.3, orientation = "y") +
  geom_text(aes(label = sprintf("%.3f", value)), hjust = -0.15, size = GEOM_TEXT_6PT,
            fontface = "plain") +
  scale_fill_manual(values = c("FALSE" = "#42A5F5", "TRUE" = "#C2185B"),
                    labels = c("Other", "NAS-VAE"),
                    name = "") +
  scale_x_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  labs(x = "Performance", y = "") +
  theme_masld(base_size = 7) +
  theme(
    strip.text       = element_text(size = 6, face = "plain"),
    axis.text.y      = element_text(size = 6),
    legend.position  = "bottom",
    legend.key.size  = unit(0.3, "cm")
  )

# ============================================================
# Compose 6-panel figure (3 rows x 2 cols)
# ============================================================
cat("\nComposing 6-panel figure (3 rows x 2 cols, 11 x 13)...\n")

layout <- "
AABB
CCDD
EEFF
"

composite <- p_a + p_b + p_c + p_d + p_e + p_f +
  plot_layout(design = layout) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(face = "plain", size = 10)
    )
  )

out_path <- file.path(OUTDIR, "fig_nas_vae_progression.pdf")
ggsave(out_path, composite, width = fig_full_width, height = fig_full_width * 13 / 11,
       device = cairo_pdf)
cat("Saved composite:", out_path, "\n")
cat("  File size:", format(file.info(out_path)$size, big.mark = ","), "bytes\n")

# --- Save individual panels ---
cat("\nSaving individual panels...\n")

panel_list <- list(
  panel_a_nas_umap       = list(p = p_a, w = 5, h = 4),
  panel_b_fib_umap       = list(p = p_b, w = 5, h = 4),
  panel_c_dataset_umap   = list(p = p_c, w = 5.5, h = 4),
  panel_d_velocity       = list(p = p_d, w = 5, h = 4),
  panel_e_bifurcation    = list(p = p_e, w = 5.5, h = 3.5),
  panel_f_model_compare  = list(p = p_f, w = 5.5, h = 4.5)
)

for (nm in names(panel_list)) {
  pinfo <- panel_list[[nm]]
  ppath <- file.path(OUTDIR, "panels", paste0(nm, ".pdf"))
  ggsave(ppath, pinfo$p, width = pinfo$w, height = pinfo$h, device = cairo_pdf)
  cat(sprintf("  %s: %s bytes\n", nm, format(file.info(ppath)$size, big.mark = ",")))
}

cat("\nAll panels saved to:", file.path(OUTDIR, "panels"), "\n")
cat("=== Done ===\n")
