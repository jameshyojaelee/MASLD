#!/usr/bin/env Rscript
#' 307b: Pseudotime Figures v2 — Informative panels only.
#'
#' Panel A: Per-cell-type UMAP colored by consensus pseudotime (5 types)
#' Panel B: Pseudotime vs disease condition (validates trajectory = disease axis)
#' Panel C: Key gene trajectories along pseudotime (marker genes per cell type)
#' Panel D: Palantir entropy (differentiation potential) by condition
#' Panel E: Bulk-SC concordance heatmap (fibrosis/NAS signatures x cell types)

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

# Publication theme
theme_src <- file.path(BASE, "scripts/figures/publication_theme.R")
if (file.exists(theme_src)) {
  source(theme_src)
} else {
  theme_masld <- function(base_size = 7) theme_classic(base_size = base_size)
  save_fig <- function(p, f, w = 7, h = 5, ...) ggsave(f, p, width = w, height = h, device = cairo_pdf)
}

CELL_TYPES <- c("Hepatocytes", "Macrophages", "Fibroblasts",
                "Endothelial_cells", "Cholangiocytes")

CT_LABELS <- c("Hepatocytes" = "Hepatocytes", "Macrophages" = "Macrophages",
               "Fibroblasts" = "Fibroblasts", "Endothelial_cells" = "Endothelial",
               "Cholangiocytes" = "Cholangiocytes")

# Color palette for conditions (disease progression order)
cond_colors <- c("Healthy" = "#2196F3", "NAFLD" = "#81D4FA",
                 "MASLD" = "#F48FB1", "NASH" = "#C2185B", "Cirrhotic" = "#4A148C")
cond_order <- c("Healthy", "NAFLD", "MASLD", "NASH", "Cirrhotic")

cat("=== 307b: Pseudotime Figures v2 ===\n")

# =========================================================================
# Load data
# =========================================================================

# Consensus pseudotime
cpt <- fread(file.path(PT_DIR, "consensus_pseudotime_all.csv"))
names(cpt)[1] <- "cell"

# Per-cell-type metadata (UMAP coords)
meta_list <- list()
for (ct in CELL_TYPES) {
  f <- file.path(PT_DIR, paste0(ct, "_metadata.csv"))
  if (file.exists(f)) {
    m <- fread(f)
    names(m)[1] <- "cell"
    m$cell_type_clean <- ct
    meta_list[[ct]] <- m
  }
}
meta_all <- rbindlist(meta_list, fill = TRUE)

# Merge UMAP + pseudotime
dat <- merge(meta_all, cpt[, .(cell, consensus_pseudotime, cell_type)],
             by = "cell", all.x = TRUE)

# Palantir results (for entropy)
palantir_list <- list()
for (ct in CELL_TYPES) {
  f <- file.path(PT_DIR, paste0("palantir_pseudotime_", ct, ".csv"))
  if (file.exists(f)) {
    p <- fread(f)
    names(p)[1] <- "cell"
    p$cell_type_clean <- ct
    palantir_list[[ct]] <- p
  }
}
palantir_all <- rbindlist(palantir_list, fill = TRUE)

# Gene dynamics
dynamics_list <- list()
for (ct in CELL_TYPES) {
  f <- file.path(PT_DIR, paste0("gene_dynamics_", ct, ".csv"))
  if (file.exists(f)) {
    d <- fread(f)
    names(d)[1] <- "gene"
    dynamics_list[[ct]] <- d
  }
}

# Bulk-SC concordance
conc <- fread(file.path(PT_DIR, "bulk_sc_concordance.csv"))

# Gene correlations (for selecting marker genes)
corr_list <- list()
for (ct in CELL_TYPES) {
  f <- file.path(PT_DIR, paste0("pseudotime_corr_", ct, ".csv"))
  if (file.exists(f)) {
    c <- fread(f)
    c$cell_type <- ct
    corr_list[[ct]] <- c
  }
}

cat("Data loaded.\n")

# =========================================================================
# Panel A: UMAP colored by pseudotime (5 cell types)
# =========================================================================
cat("--- Panel A: UMAPs ---\n")

# Subsample for plotting speed
set.seed(42)
dat_plot <- dat[, .SD[sample(.N, min(.N, 20000))], by = cell_type_clean]

pa <- ggplot(dat_plot[!is.na(consensus_pseudotime)],
             aes(x = UMAP_1, y = UMAP_2, color = consensus_pseudotime)) +
  geom_point(size = 0.1, alpha = 0.5) +
  scale_color_viridis_c(option = "magma", name = "Pseudotime") +
  facet_wrap(~factor(cell_type_clean, levels = CELL_TYPES, labels = CT_LABELS),
             nrow = 1, scales = "free") +
  theme_masld(base_size = 7) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        legend.key.width = unit(0.4, "cm"), legend.key.height = unit(0.2, "cm"),
        legend.position = "bottom") +
  labs(x = NULL, y = NULL)

# =========================================================================
# Panel B: Pseudotime by disease condition (boxplots)
# =========================================================================
cat("--- Panel B: Pseudotime vs condition ---\n")

# Filter to conditions with enough cells
dat_cond <- dat[condition %in% cond_order & !is.na(consensus_pseudotime)]
dat_cond[, condition := factor(condition, levels = cond_order)]

# Aggregate to median per condition per cell type for cleaner display
pb <- ggplot(dat_cond,
             aes(x = condition, y = consensus_pseudotime, fill = condition)) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8, linewidth = 0.3) +
  scale_fill_manual(values = cond_colors, guide = "none") +
  facet_wrap(~factor(cell_type_clean, levels = CELL_TYPES, labels = CT_LABELS),
             nrow = 1, scales = "free_y") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text = element_text(size = 8, face = "bold")) +
  labs(x = NULL, y = "Consensus pseudotime")

# =========================================================================
# Panel C: Key gene trajectories (ribbon plots)
# =========================================================================
cat("--- Panel C: Gene trajectories ---\n")

# Select biologically meaningful marker genes per cell type
marker_genes <- list(
  Hepatocytes = c("CRP", "SAA2", "TTR", "CYP3A4", "APOB", "FASN"),
  Macrophages = c("HLA-DRA", "HLA-DPA1", "CD163", "TREM2", "SPP1", "ALB"),
  Fibroblasts = c("COL1A1", "ACTA2", "DCN", "VIM", "PDGFRA", "COL3A1"),
  Endothelial_cells = c("PECAM1", "VWF", "CLEC4G", "ACKR1", "PLVAP", "CD34"),
  Cholangiocytes = c("KRT19", "KRT7", "EPCAM", "SOX9", "HNF1B", "SPP1")
)

# Build trajectory curves from gene_dynamics
traj_data <- list()
for (ct in names(dynamics_list)) {
  d <- dynamics_list[[ct]]
  genes_avail <- intersect(marker_genes[[ct]], d$gene)
  if (length(genes_avail) == 0) next

  for (g in genes_avail) {
    row <- d[gene == g]
    if (nrow(row) == 0) next
    vals <- as.numeric(row[1, -1, with = FALSE])
    # Smooth with loess
    bins <- seq(0, 1, length.out = length(vals))
    fit <- tryCatch(loess(vals ~ bins, span = 0.2), error = function(e) NULL)
    if (is.null(fit)) next
    pred <- predict(fit, bins)
    traj_data[[paste0(ct, "_", g)]] <- data.table(
      cell_type = ct, gene = g,
      pseudotime = bins, expression = pred
    )
  }
}
traj_df <- rbindlist(traj_data)

# Pick top 3 genes per cell type (most dynamic range)
if (nrow(traj_df) > 0) {
  gene_range <- traj_df[, .(dynamic_range = max(expression) - min(expression)), by = .(cell_type, gene)]
  top_genes <- gene_range[, head(.SD[order(-dynamic_range)], 3), by = cell_type]
  traj_filt <- traj_df[paste(cell_type, gene) %in% paste(top_genes$cell_type, top_genes$gene)]

  pc <- ggplot(traj_filt,
               aes(x = pseudotime, y = expression, color = gene)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~factor(cell_type, levels = CELL_TYPES, labels = CT_LABELS),
               nrow = 1, scales = "free_y") +
    theme_masld(base_size = 7) +
    theme(strip.text = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 5),
          legend.key.height = unit(0.3, "cm")) +
    labs(x = "Pseudotime", y = "Expression (smoothed)", color = NULL) +
    scale_color_brewer(palette = "Set1")
} else {
  pc <- ggplot() + theme_void() + labs(title = "No gene dynamics available")
}

# =========================================================================
# Panel D: Palantir entropy by condition
# =========================================================================
cat("--- Panel D: Differentiation potential ---\n")

if (nrow(palantir_all) > 0 && "palantir_entropy" %in% names(palantir_all)) {
  pal_cond <- palantir_all[condition %in% cond_order]
  pal_cond[, condition := factor(condition, levels = cond_order)]

  pd <- ggplot(pal_cond[!is.na(palantir_entropy)],
               aes(x = condition, y = palantir_entropy, fill = condition)) +
    geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8, linewidth = 0.3) +
    scale_fill_manual(values = cond_colors, guide = "none") +
    facet_wrap(~factor(cell_type_clean, levels = CELL_TYPES, labels = CT_LABELS),
               nrow = 1, scales = "free_y") +
    theme_masld(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text = element_text(size = 8, face = "bold")) +
    labs(x = NULL, y = "Differentiation potential\n(Palantir entropy)")
} else {
  pd <- ggplot() + theme_void() + labs(title = "No Palantir entropy data")
}

# =========================================================================
# Panel E: Bulk-SC concordance heatmap
# =========================================================================
cat("--- Panel E: Bulk-SC heatmap ---\n")

if (nrow(conc) > 0) {
  # Use Palantir method, reshape to matrix
  conc_pal <- conc[pseudotime_method == "palantir"]
  if (nrow(conc_pal) == 0) conc_pal <- conc[pseudotime_method == "dpt"]

  # Clean signature names
  conc_pal[, sig_clean := gsub("transition_", "", signature)]
  conc_pal[, sig_clean := gsub("_", " → ", sig_clean)]
  conc_pal[, sig_clean := gsub(" → UP", " ↑", sig_clean)]
  conc_pal[, sig_clean := gsub(" → DOWN", " ↓", sig_clean)]

  # Order by fibrosis stage
  sig_order <- c("F0 → F1 ↑", "F0 → F1 ↓", "F1 → F2 ↑", "F1 → F2 ↓",
                 "F2 → F3 ↑", "F2 → F3 ↓", "F3 → F4 ↑", "F3 → F4 ↓")
  conc_pal <- conc_pal[sig_clean %in% sig_order]
  conc_pal[, sig_clean := factor(sig_clean, levels = rev(sig_order))]

  pe <- ggplot(conc_pal,
               aes(x = factor(cell_type, levels = CELL_TYPES, labels = CT_LABELS),
                   y = sig_clean, fill = spearman_rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 2) +
    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
                         midpoint = 0, limits = c(-0.6, 0.6),
                         name = "Spearman ρ") +
    theme_masld(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.5, "cm")) +
    labs(x = NULL, y = NULL, title = "Bulk fibrosis signature vs sc pseudotime")
} else {
  pe <- ggplot() + theme_void() + labs(title = "No concordance data")
}

# =========================================================================
# Combine and save
# =========================================================================
cat("--- Assembling figure ---\n")

# Full figure: A over B over C over (D | E)
fig <- (pa / pb / pc / (pd | pe)) +
  plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(1, 1, 1, 1.2))

out_path <- file.path(FIG_DIR, "figS_pseudotime_v2.pdf")
save_fig(fig, out_path, width = 7.09, height = 10)
cat("Saved:", out_path, "\n")

# Also save panels individually for flexibility
save_fig(pa, file.path(FIG_DIR, "panel_umap_pseudotime.pdf"), width = 7.09, height = 2)
save_fig(pb, file.path(FIG_DIR, "panel_pseudotime_by_condition.pdf"), width = 7.09, height = 2)
save_fig(pc, file.path(FIG_DIR, "panel_gene_trajectories.pdf"), width = 7.09, height = 2)
save_fig(pd, file.path(FIG_DIR, "panel_differentiation_potential.pdf"), width = 7.09, height = 2)
save_fig(pe, file.path(FIG_DIR, "panel_bulk_sc_heatmap.pdf"), width = 4, height = 3)

cat("\n=== 307b: Figures v2 COMPLETE ===\n")
