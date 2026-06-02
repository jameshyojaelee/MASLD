#!/usr/bin/env Rscript
#' 307c: Pseudotime Figures v3 — Clean 4-panel supplementary figure.
#'
#' a: Per-cell-type UMAP colored by pseudotime
#' b: Pseudotime by disease condition (Healthy vs MASLD, the two powered groups)
#' c: Top gene trajectories per cell type (data-driven, not hardcoded)
#' d: Bulk fibrosis signature concordance heatmap

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PT_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudotime")
FIG_DIR <- file.path(BASE, "figures/supplementary/figS02_progression")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

theme_src <- file.path(BASE, "scripts/figures/publication_theme.R")
if (file.exists(theme_src)) source(theme_src)

CELL_TYPES <- c("Hepatocytes", "Macrophages", "Fibroblasts",
                "Endothelial_cells", "Cholangiocytes")
CT_LABELS <- c(Hepatocytes = "Hepatocytes", Macrophages = "Macrophages",
               Fibroblasts = "Fibroblasts", Endothelial_cells = "Endothelial",
               Cholangiocytes = "Cholangiocytes")

cond_colors <- c(Healthy = "#2196F3", MASLD = "#E91E63")

cat("=== 307c: Pseudotime Figures v3 ===\n")

# =========================================================================
# Load
# =========================================================================
cpt <- fread(file.path(PT_DIR, "consensus_pseudotime_all.csv"))
names(cpt)[1] <- "cell"

meta_list <- lapply(CELL_TYPES, function(ct) {
  f <- file.path(PT_DIR, paste0(ct, "_metadata.csv"))
  if (!file.exists(f)) return(NULL)
  m <- fread(f); names(m)[1] <- "cell"; m$ct <- ct; m
})
meta_all <- rbindlist(meta_list, fill = TRUE)

dat <- merge(meta_all, cpt[, .(cell, consensus_pseudotime, cell_type)],
             by = "cell", all.x = TRUE)

conc <- fread(file.path(PT_DIR, "bulk_sc_concordance.csv"))

# =========================================================================
# Panel A: UMAPs (subsampled, rasterized)
# =========================================================================
cat("Panel A...\n")
set.seed(42)
dp <- dat[!is.na(consensus_pseudotime),
          .SD[sample(.N, min(.N, 15000))], by = ct]
# Shuffle to avoid overplotting bias
dp <- dp[sample(.N)]

pa <- ggplot(dp, aes(x = UMAP_1, y = UMAP_2, color = consensus_pseudotime)) +
  geom_point(size = 0.05, alpha = 0.6, shape = 16) +
  scale_color_viridis_c(option = "magma", name = "Pseudotime",
                         breaks = c(0, 0.5, 1)) +
  facet_wrap(~factor(ct, levels = CELL_TYPES, labels = CT_LABELS),
             nrow = 1, scales = "free") +
  theme_masld(base_size = 7) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(0.8, "cm"), legend.key.height = unit(0.15, "cm"),
        panel.spacing = unit(0.1, "cm"),
        plot.margin = margin(2, 2, 2, 2)) +
  labs(x = NULL, y = NULL) +
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))

# =========================================================================
# Panel B: Pseudotime by Healthy vs MASLD (the two powered conditions)
# =========================================================================
cat("Panel B...\n")

# Only Healthy and MASLD have enough cells for meaningful comparison
dat_hm <- dat[condition %in% c("Healthy", "MASLD") & !is.na(consensus_pseudotime)]
dat_hm[, condition := factor(condition, levels = c("Healthy", "MASLD"))]

# Count cells per condition per type
n_tab <- dat_hm[, .N, by = .(ct, condition)]
n_tab[, label := paste0("n=", format(N, big.mark = ","))]

pb <- ggplot(dat_hm, aes(x = condition, y = consensus_pseudotime, fill = condition)) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
               alpha = 0.9, linewidth = 0.3) +
  geom_text(data = n_tab, aes(x = condition, label = label),
            y = -0.08, size = 1.8, inherit.aes = FALSE) +
  scale_fill_manual(values = cond_colors, guide = "none") +
  facet_wrap(~factor(ct, levels = CELL_TYPES, labels = CT_LABELS),
             nrow = 1) +
  theme_masld(base_size = 7) +
  theme(strip.text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 7),
        panel.spacing = unit(0.3, "cm")) +
  coord_cartesian(ylim = c(-0.12, 1.05), clip = "off") +
  labs(x = NULL, y = "Consensus\npseudotime")

# Add Wilcoxon p-values
for (ct_name in CELL_TYPES) {
  h_vals <- dat_hm[ct == ct_name & condition == "Healthy", consensus_pseudotime]
  m_vals <- dat_hm[ct == ct_name & condition == "MASLD", consensus_pseudotime]
  if (length(h_vals) > 10 & length(m_vals) > 10) {
    wt <- wilcox.test(h_vals, m_vals)
    es <- median(m_vals) - median(h_vals)
    cat(sprintf("  %s: Wilcoxon p=%.2e, Δmedian=%.3f\n", ct_name, wt$p.value, es))
  }
}

# =========================================================================
# Panel C: Top gene trajectories (data-driven from pseudotime_corr)
# =========================================================================
cat("Panel C...\n")

# Load gene dynamics + correlations to pick the best genes
gene_panels <- list()

for (ct in CELL_TYPES) {
  dyn_f <- file.path(PT_DIR, paste0("gene_dynamics_", ct, ".csv"))
  corr_f <- file.path(PT_DIR, paste0("pseudotime_corr_", ct, ".csv"))
  if (!file.exists(dyn_f) || !file.exists(corr_f)) next

  dyn <- fread(dyn_f); names(dyn)[1] <- "gene"
  corr <- fread(corr_f)

  # Pick top 2 positive + top 2 negative correlated genes
  # Exclude ribosomal (RPL/RPS), mitochondrial (MT-), and ENSG IDs
  corr_clean <- corr[!grepl("^RPL|^RPS|^MT-|^ENSG", gene) & padj < 0.05]
  top_pos <- head(corr_clean[spearman_rho > 0][order(-spearman_rho)], 2)
  top_neg <- head(corr_clean[spearman_rho < 0][order(spearman_rho)], 2)
  sel_genes <- c(top_pos$gene, top_neg$gene)

  for (g in sel_genes) {
    row <- dyn[gene == g]
    if (nrow(row) == 0) next
    vals <- as.numeric(row[1, -1, with = FALSE])
    bins <- seq(0, 1, length.out = length(vals))
    # Smooth
    fit <- tryCatch(loess(vals ~ bins, span = 0.25), error = function(e) NULL)
    if (is.null(fit)) next
    pred <- predict(fit, bins)
    rho_val <- corr[gene == g, spearman_rho]
    gene_panels[[paste0(ct, "_", g)]] <- data.table(
      ct = ct, gene = g, pseudotime = bins, expression = pred,
      direction = ifelse(rho_val > 0, "Up with disease", "Down with disease")
    )
  }
}
traj_df <- rbindlist(gene_panels)

if (nrow(traj_df) > 0) {
  # Use direct labeling instead of shared legend
  # Get endpoint positions for labels
  label_pos <- traj_df[, .SD[which.max(pseudotime)], by = .(ct, gene)]

  pc <- ggplot(traj_df, aes(x = pseudotime, y = expression,
                              color = direction, group = gene)) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    geom_text(data = label_pos, aes(label = gene), size = 1.8,
              hjust = 0, nudge_x = 0.02, show.legend = FALSE) +
    scale_color_manual(values = c("Up with disease" = "#C2185B",
                                   "Down with disease" = "#1565C0"),
                       name = NULL) +
    facet_wrap(~factor(ct, levels = CELL_TYPES, labels = CT_LABELS),
               nrow = 1, scales = "free_y") +
    theme_masld(base_size = 7) +
    theme(strip.text = element_text(size = 8, face = "bold"),
          legend.position = "bottom",
          legend.key.size = unit(0.3, "cm"),
          panel.spacing = unit(0.3, "cm")) +
    coord_cartesian(clip = "off") +
    labs(x = "Pseudotime →", y = "Expression\n(smoothed)")
} else {
  pc <- ggplot() + theme_void()
}

# =========================================================================
# Panel D: Bulk-SC concordance heatmap (larger, cleaner)
# =========================================================================
cat("Panel D...\n")

conc_pal <- conc[pseudotime_method == "palantir"]
if (nrow(conc_pal) == 0) conc_pal <- conc[pseudotime_method == "dpt"]

# Clean names
conc_pal[, sig_clean := gsub("transition_", "", signature)]
conc_pal[, sig_clean := gsub("_DOWN", " ↓", sig_clean)]
conc_pal[, sig_clean := gsub("_UP", " ↑", sig_clean)]
conc_pal[, sig_clean := gsub("_to_", " → ", sig_clean)]

# Keep fibrosis transitions only (NAS are redundant / noisy)
fib_sigs <- conc_pal[grepl("^F[0-4]", sig_clean)]
sig_order <- c("F0 → F1 ↑", "F0 → F1 ↓", "F1 → F2 ↑", "F1 → F2 ↓",
               "F2 → F3 ↑", "F2 → F3 ↓", "F3 → F4 ↑", "F3 → F4 ↓")
fib_sigs <- fib_sigs[sig_clean %in% sig_order]
fib_sigs[, sig_clean := factor(sig_clean, levels = rev(sig_order))]
fib_sigs[, ct_label := CT_LABELS[cell_type]]
fib_sigs[, ct_label := factor(ct_label,
                                levels = c("Hepatocytes", "Macrophages", "Fibroblasts",
                                           "Endothelial", "Cholangiocytes"))]

pd_heatmap <- ggplot(fib_sigs,
                      aes(x = ct_label, y = sig_clean, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 2.2,
            color = ifelse(abs(fib_sigs$spearman_rho) > 0.3, "white", "gray30")) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
                        midpoint = 0, limits = c(-0.65, 0.65),
                        name = "Spearman ρ") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
        axis.text.y = element_text(size = 6.5),
        legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.6, "cm"),
        plot.title = element_text(size = 8, face = "bold")) +
  labs(x = NULL, y = NULL,
       title = "Bulk fibrosis signatures vs\nsc pseudotime")

# =========================================================================
# Assemble: A / B / C / D (4 rows)
# =========================================================================
cat("Assembling...\n")

fig <- (pa / pb / pc / pd_heatmap) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

# Adjust heights: UMAPs need more, heatmap less
fig <- fig + plot_layout(heights = c(1.2, 1, 1, 1.3))

out <- file.path(FIG_DIR, "figS_pseudotime_v3.pdf")
save_fig(fig, out, width = 7.09, height = 9)
cat("Saved:", out, "\n")

# Individual panels
save_fig(pa, file.path(FIG_DIR, "v3_panel_a_umap.pdf"), width = 7.09, height = 1.8)
save_fig(pb, file.path(FIG_DIR, "v3_panel_b_condition.pdf"), width = 7.09, height = 2)
save_fig(pc, file.path(FIG_DIR, "v3_panel_c_genes.pdf"), width = 7.09, height = 2)
save_fig(pd_heatmap, file.path(FIG_DIR, "v3_panel_d_heatmap.pdf"), width = 4.5, height = 3.5)

cat("=== 307c COMPLETE ===\n")
