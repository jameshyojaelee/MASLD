#!/usr/bin/env Rscript
# ==========================================================================
# Single-Cell Compositional Analysis + Cross-Species Concordance Figures
#
# Module 1 — CLR + LMM Compositional Analysis (panels a–c):
#   a: Forest plot — LMM beta (CLR) per cell type (MASLD vs Healthy)
#   b: CLR proportion boxplots — hepatocytes + key cell types by condition
#   c: Sample PCA of CLR-transformed proportions (colored by hepatocyte CLR)
#
# Module 2 — Cross-Species Effect Size Normalization (panels d–e):
#   d: Human sc hepatocyte vs Mouse bulk logFC scatter (SD-normalized)
#   e: Human bulk vs Mouse bulk logFC scatter (reference; SD-normalized)
#
# Output: figures/fig_sc_compositional_crossspecies.pdf (combined)
#         figures/panels/fig_sc_{a,b,c,d,e}.pdf (individual)
# ==========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(ggrepel)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

SC_DIR     <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2")
COMP_DIR   <- file.path(SC_DIR, "compositional")
MOUSE_DIR  <- file.path(SC_DIR, "mouse_sc")
PANEL_DIR  <- file.path(FIG_OUT, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

OUT <- file.path(FIG_OUT, "fig_sc_compositional_crossspecies.pdf")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
lmm   <- fread(file.path(COMP_DIR, "lmm_binary_results.csv"))
clr   <- fread(file.path(COMP_DIR, "clr_values.csv"))
norm  <- fread(file.path(MOUSE_DIR, "species_effect_norm_hepatocytes.csv"))

# Prettify cell type labels
prettify <- function(x) {
  x <- gsub("\\.", " ", x)
  x <- gsub("Mono mono", "Mono+mono", x)
  x
}
lmm[, ct_label := prettify(cell_type)]

# Sort by beta for forest plot
lmm <- lmm[order(beta)]
lmm[, ct_label := factor(ct_label, levels = ct_label)]

# Color scheme: up (magenta) / down (blue) / NS (gray)
lmm[, color_group := "NS"]
lmm[significant == TRUE & beta > 0, color_group := "Up in MASLD"]
lmm[significant == TRUE & beta < 0, color_group := "Down in MASLD"]

ct_colors <- c(
  "Up in MASLD"   = masld_colors$up,
  "Down in MASLD" = masld_colors$down,
  "NS"            = masld_colors$ns
)

# ===========================================================================
# Panel A — Forest plot
# ===========================================================================
panel_a <- ggplot(lmm, aes(x = beta, y = ct_label, color = color_group)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3, linetype = "dashed") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                orientation = "y", width = 0.35, linewidth = 0.45) +
  geom_point(aes(shape = significant), size = 1.8) +
  scale_color_manual(values = ct_colors, name = NULL) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), guide = "none") +
  annotate("text",
           x = max(lmm$ci_hi, na.rm = TRUE) * 0.95,
           y = 1.5,
           label = paste0("Hepatocytes\nbeta = ", round(lmm[cell_type == "Hepatocytes", beta], 3),
                          "\npadj = ", round(lmm[cell_type == "Hepatocytes", padj], 2)),
           size = 2, hjust = 1, color = "grey30", fontface = "italic") +
  labs(
    title    = "Cell type proportions: MASLD vs Healthy",
    subtitle = paste0("CLR-LMM | dataset random effect | n=",
                      lmm$n_samples[1], " samples"),
    x = "beta (CLR units; positive = increased in MASLD)",
    y = NULL
  ) +
  theme_masld(base_size = 7) +
  theme(
    legend.position  = c(0.98, 0.12),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.key.size  = unit(0.25, "cm")
  )

ggsave(file.path(PANEL_DIR, "fig_sc_a_forest.pdf"),
       panel_a, width = 6.5, height = 4.5)
message("Saved panel a")

# ===========================================================================
# Panel B — CLR boxplots for key cell types
# ===========================================================================
key_cts  <- c("Hepatocytes", "Cholangiocytes", "Endothelial cells",
              "Macrophages", "Fibroblasts", "T cells")
ct_cols  <- setdiff(colnames(clr), c("sample", "dataset", "condition_harmonized"))
key_cols <- intersect(key_cts, ct_cols)

box_dt <- melt(
  clr[condition_harmonized %in% c("Healthy", "MASLD"),
       c("sample", "condition_harmonized", key_cols), with = FALSE],
  id.vars      = c("sample", "condition_harmonized"),
  variable.name = "cell_type",
  value.name    = "clr_value"
)
box_dt[, cell_type := factor(cell_type, levels = key_cols)]

# Add significance annotation from LMM
sig_dt <- lmm[prettify(cell_type) %in% key_cols,
              .(cell_type = prettify(cell_type), padj, significant)]
box_dt <- merge(box_dt, sig_dt, by = "cell_type", all.x = TRUE)
box_dt[, sig_label := ifelse(padj < 0.001, "***",
                      ifelse(padj < 0.01,  "**",
                      ifelse(padj < 0.05,  "*", "ns")))]

# Per-cell-type y-max for annotation placement
y_max_dt <- box_dt[, .(y_max = max(clr_value, na.rm = TRUE) + 0.1 * diff(range(clr_value, na.rm = TRUE))),
                   by = cell_type]
box_dt <- merge(box_dt, y_max_dt, by = "cell_type")

panel_b <- ggplot(box_dt,
                  aes(x = condition_harmonized, y = clr_value,
                      fill = condition_harmonized)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.35, width = 0.55) +
  geom_jitter(width = 0.12, alpha = 0.35, size = 0.3, color = "grey30") +
  geom_text(
    data = unique(box_dt[, .(cell_type, sig_label, y_max)]),
    aes(x = 1.5, y = y_max, label = sig_label),
    inherit.aes = FALSE, size = 2.5, vjust = 0
  ) +
  scale_fill_manual(
    values = c("Healthy" = masld_colors$control, "MASLD" = masld_colors$masld),
    name   = NULL
  ) +
  scale_x_discrete(labels = c("Healthy" = "Ctrl", "MASLD" = "MASLD")) +
  facet_wrap(~ cell_type, scales = "free_y", nrow = 1) +
  labs(
    title    = "CLR-transformed cell type proportions",
    subtitle = "* padj<0.05  ** padj<0.01  *** padj<0.001  ns = not significant",
    x = NULL, y = "CLR score"
  ) +
  theme_masld(base_size = 7) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(size = 6),
    axis.text.x      = element_text(size = 5.5)
  )

ggsave(file.path(PANEL_DIR, "fig_sc_b_boxplot.pdf"),
       panel_b, width = 8, height = 3)
message("Saved panel b")

# ===========================================================================
# Panel C — Sample PCA of CLR proportions
# ===========================================================================
clr_mat  <- as.matrix(clr[condition_harmonized %in% c("Healthy", "MASLD"), ..ct_cols])
clr_meta <- clr[condition_harmonized %in% c("Healthy", "MASLD")]

pca_res  <- prcomp(clr_mat, scale. = FALSE, center = TRUE)
pvar     <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- data.table(
  PC1       = pca_res$x[, 1],
  PC2       = pca_res$x[, 2],
  condition = clr_meta$condition_harmonized,
  dataset   = clr_meta$dataset,
  hep_clr   = clr_mat[, "Hepatocytes"]
)

# Centroids per condition
centroids <- pca_df[, .(PC1 = mean(PC1), PC2 = mean(PC2)), by = condition]

panel_c <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = hep_clr, shape = condition),
             size = 1.8, stroke = 0.2, color = "grey50", alpha = 0.85) +
  geom_point(data = centroids,
             aes(x = PC1, y = PC2, color = condition),
             size = 3.5, shape = 4, stroke = 1.2) +
  scale_fill_gradient2(
    low = masld_colors$down, mid = "white", high = masld_colors$up,
    midpoint = 0, name = "Hepatocyte\nCLR score",
    guide = guide_colorbar(barheight = 2.5, barwidth = 0.4)
  ) +
  scale_shape_manual(
    values = c("Healthy" = 21, "MASLD" = 24), name = "Condition"
  ) +
  scale_color_manual(
    values = c("Healthy" = masld_colors$control, "MASLD" = masld_colors$masld),
    guide = "none"
  ) +
  labs(
    title    = "Sample composition PCA",
    subtitle = "CLR-transformed cell type proportions | × = centroid",
    x        = paste0("PC1 (", pvar[1], "%)"),
    y        = paste0("PC2 (", pvar[2], "%)")
  ) +
  theme_masld(base_size = 7) +
  theme(legend.position = "right")

ggsave(file.path(PANEL_DIR, "fig_sc_c_pca.pdf"),
       panel_c, width = 4.5, height = 3.5)
message("Saved panel c")

# ===========================================================================
# Panel D — Human sc hepatocyte vs Mouse bulk (SD-normalized)
# ===========================================================================
norm_ab <- norm[!is.na(human_sc_logFC_norm) & !is.na(mouse_bulk_logFC_norm)]
rho_d   <- round(cor(norm_ab$human_sc_logFC_norm,
                     norm_ab$mouse_bulk_logFC_norm,
                     method = "spearman"), 3)
n_d     <- nrow(norm_ab)
conc_d  <- round(100 * mean(sign(norm_ab$human_sc_logFC_norm) ==
                            sign(norm_ab$mouse_bulk_logFC_norm)), 1)

# Highlight Conserved
norm_ab[, point_class := ifelse(is_conserved_core == TRUE, "Conserved", "Other")]

panel_d <- ggplot(norm_ab, aes(x = human_sc_logFC_norm,
                               y = mouse_bulk_logFC_norm)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_point(data = norm_ab[point_class == "Other"],
             color = masld_colors$ns, size = 0.3, alpha = 0.35, shape = 16) +
  geom_point(data = norm_ab[point_class == "Conserved"],
             color = masld_colors$up, size = 0.7, alpha = 0.75, shape = 16) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              color = "grey20", fill = "grey85", linewidth = 0.6) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
           label = sprintf("Spearman rho = %.3f\nn = %s genes\n%s%% concordant",
                           rho_d, format(n_d, big.mark = ","), conc_d),
           size = 2, family = "mono") +
  labs(
    title    = "Human sc vs Mouse bulk (Hepatocytes)",
    subtitle = "SD-normalized logFC | Red = Conserved (723 genes)",
    x = "Human sc hepatocyte logFC / SD",
    y = "Mouse bulk hepatocyte logFC / SD"
  ) +
  theme_masld(base_size = 7)

ggsave(file.path(PANEL_DIR, "fig_sc_d_sc_vs_bulk.pdf"),
       panel_d, width = 3.5, height = 3.5)
message("Saved panel d")

# ===========================================================================
# Panel E — Human bulk vs Mouse bulk (SD-normalized; reference comparison)
# ===========================================================================
norm_bb <- norm[!is.na(human_bulk_logFC_norm) & !is.na(mouse_bulk_logFC_norm)]
rho_e   <- round(cor(norm_bb$human_bulk_logFC_norm,
                     norm_bb$mouse_bulk_logFC_norm,
                     method = "spearman"), 3)
n_e     <- nrow(norm_bb)
conc_e  <- round(100 * mean(sign(norm_bb$human_bulk_logFC_norm) ==
                            sign(norm_bb$mouse_bulk_logFC_norm)), 1)

norm_bb[, point_class := ifelse(is_conserved_core == TRUE, "Conserved", "Other")]

panel_e <- ggplot(norm_bb, aes(x = human_bulk_logFC_norm,
                               y = mouse_bulk_logFC_norm)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_point(data = norm_bb[point_class == "Other"],
             color = masld_colors$ns, size = 0.3, alpha = 0.35, shape = 16) +
  geom_point(data = norm_bb[point_class == "Conserved"],
             color = masld_colors$up, size = 0.7, alpha = 0.75, shape = 16) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              color = "grey20", fill = "grey85", linewidth = 0.6) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
           label = sprintf("Spearman rho = %.3f\nn = %s genes\n%s%% concordant",
                           rho_e, format(n_e, big.mark = ","), conc_e),
           size = 2, family = "mono") +
  labs(
    title    = "Human bulk vs Mouse bulk (reference)",
    subtitle = "SD-normalized integrated logFC | Red = Conserved",
    x = "Human bulk integrated logFC / SD",
    y = "Mouse bulk integrated logFC / SD"
  ) +
  theme_masld(base_size = 7)

ggsave(file.path(PANEL_DIR, "fig_sc_e_bulk_vs_bulk.pdf"),
       panel_e, width = 3.5, height = 3.5)
message("Saved panel e")

# ===========================================================================
# Panel F — Hepatocyte CLR score: Healthy vs MASLD (focused boxplot)
# ===========================================================================
hep_box <- data.table(
  condition = clr_meta$condition_harmonized,
  dataset   = clr_meta$dataset,
  hep_clr   = clr_mat[, "Hepatocytes"]
)

hep_t    <- t.test(hep_clr ~ condition, data = hep_box)
hep_p    <- format(hep_t$p.value, digits = 2, scientific = TRUE)
hep_mean <- hep_box[, .(mean_clr = round(mean(hep_clr), 3)), by = condition]

panel_f <- ggplot(hep_box, aes(x = condition, y = hep_clr, fill = condition)) +
  geom_boxplot(width = 0.5, outlier.size = 0.5, linewidth = 0.4,
               outlier.alpha = 0.4) +
  geom_jitter(width = 0.12, alpha = 0.35, size = 0.4, color = "grey30") +
  scale_fill_manual(
    values = c("Healthy" = masld_colors$control, "MASLD" = masld_colors$masld),
    guide  = "none"
  ) +
  scale_x_discrete(labels = c("Healthy" = "Healthy\n(n=34)",
                               "MASLD"   = "MASLD\n(n=188)")) +
  annotate("text", x = 1.5, y = max(hep_box$hep_clr) + 0.1,
           label = paste0("p = ", hep_p, "\n(t-test)"),
           size = 2, hjust = 0.5) +
  annotate("segment", x = 1, xend = 2,
           y = max(hep_box$hep_clr) + 0.05,
           yend = max(hep_box$hep_clr) + 0.05,
           linewidth = 0.3) +
  labs(
    title    = "Hepatocyte proportion stable in MASLD",
    subtitle = paste0("CLR score (beta = ",
                      round(lmm[cell_type == "Hepatocytes", beta], 3),
                      ", padj = ",
                      round(lmm[cell_type == "Hepatocytes", padj], 2), ")"),
    x = NULL, y = "Hepatocyte CLR score"
  ) +
  theme_masld(base_size = 7)

ggsave(file.path(PANEL_DIR, "fig_sc_f_hepatocyte_clr.pdf"),
       panel_f, width = 2.5, height = 3.5)
message("Saved panel f")

# ===========================================================================
# Combined figure: (a) forest | (b) boxplots / (c) PCA | (d+e) scatter
# ===========================================================================
combined <- (
  panel_a /
  (panel_b) /
  ((panel_c | panel_f) / (panel_d | panel_e))
) +
  plot_annotation(
    title = "Single-Cell Compositional Analysis + Cross-Species Concordance",
    tag_levels = "a",
    theme = theme_masld(base_size = 8)
  ) &
  theme(plot.tag = element_text(face = "bold", size = 8))

ggsave(OUT, combined, width = 10, height = 14)
message("\nCombined figure saved: ", OUT)

# ===========================================================================
# Summary table printed to console
# ===========================================================================
message("\n=== Key numbers for figure legends ===")
message("Module 1 — CLR+LMM:")
message("  n samples: ", nrow(clr_meta), " (unsorted/nuclei, Healthy+MASLD)")
message("  n datasets: ", length(unique(clr_meta$dataset)))
message("  Hepatocytes: beta=", round(lmm[cell_type=="Hepatocytes", beta], 3),
        " padj=", round(lmm[cell_type=="Hepatocytes", padj], 2))
message("  Sig cell types: ",
        paste(lmm[significant==TRUE, cell_type], collapse=", "))
message("")
message("Module 2 — Cross-species:")
message("  Human sc vs Mouse bulk: rho=", rho_d, " n=", format(n_d, big.mark=","),
        " concordant=", conc_d, "%")
message("  Human bulk vs Mouse bulk: rho=", rho_e, " n=", format(n_e, big.mark=","),
        " concordant=", conc_e, "%")
message("  Conserved genes used: ", sum(norm_bb$is_conserved_core, na.rm=TRUE))
