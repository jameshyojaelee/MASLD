#!/usr/bin/env Rscript
# ==========================================================================
# Supplementary Figure: Single-Cell Validation of Bulk Findings
# Shows that hepatocyte pseudobulk DE recapitulates bulk RNA-seq results
#
# Panels:
#   a: Hepatocyte pseudobulk logFC vs bulk dream logFC (scatter)
#   b: Per-cell-type DEG overlap with bulk dream (horizontal bar)
#   c: Cell-type-specific unique DEGs (genes only found in sc, not bulk)
#   d: Hepatocyte fraction vs deconvolution estimate (MuSiC validation)
# ==========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

dir.create(file.path(FIGS03_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(FIGS03_DIR, "panels", "figS_sc_validation.pdf")

# ==========================================================================
# Load data
# ==========================================================================
pb_de_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
dream <- load_dream_results()

# Cell-type palette
ct_palette <- c(
  "Hepatocytes"          = "#0D47A1",
  "Cholangiocytes"       = "#1565C0",
  "Endothelial cells"    = "#2E7D32",
  "Fibroblasts"          = "#F57F17",
  "Macrophages"          = "#C2185B",
  "Mono+mono derived cells" = "#E91E63",
  "T cells"              = "#7B1FA2",
  "B cells"              = "#9C27B0",
  "Resident NK"          = "#00695C",
  "Circulating NK/NKT"   = "#00897B",
  "Plasma cells"         = "#5D4037",
  "Neutrophils"          = "#FF6F00",
  "cDC1s"                = "#AD1457",
  "cDC2s"                = "#D81B60",
  "pDCs"                 = "#6A1B9A",
  "Basophils"            = "#78909C"
)

# ==========================================================================
# Panel a: Hepatocyte pseudobulk vs bulk dream scatter
# ==========================================================================
hep_de_file <- file.path(pb_de_dir, "Hepatocytes_de.csv")

if (file.exists(hep_de_file) && !is.null(dream)) {
  hep_de <- fread(hep_de_file)
  hep_de <- add_symbols(hep_de, "gene")

  # Keep lfsr + shrunk_logFC so is_dream_deg() can apply the canonical ashr gate.
  dream_merge_cols <- intersect(c("symbol", "dream_logFC", "dream_padj",
                                  "dream_shrunk_logFC", "dream_lfsr"), names(dream))
  merged <- merge(
    dream[, ..dream_merge_cols],
    hep_de[, .(symbol, hep_logFC = logFC, hep_padj = padj)],
    by = "symbol", allow.cartesian = FALSE
  )

  merged[, bulk_sig := is_dream_deg(merged)]
  merged[, sig_class := fcase(
    bulk_sig & hep_padj < 0.05, "Both significant",
    bulk_sig,                    "Bulk only",
    hep_padj < 0.05,             "Hepatocyte only",
    default = "Neither"
  )]
  merged[, bulk_sig := NULL]

  rho <- cor(merged$dream_logFC, merged$hep_logFC,
             use = "complete.obs", method = "spearman")
  r <- cor(merged$dream_logFC, merged$hep_logFC,
           use = "complete.obs", method = "pearson")
  n_both <- sum(merged$sig_class == "Both significant")
  n_total <- nrow(merged)

  # Direction concordance for significant genes
  sig_both <- merged[sig_class == "Both significant"]
  dir_conc <- if (nrow(sig_both) > 0) {
    mean(sign(sig_both$dream_logFC) == sign(sig_both$hep_logFC)) * 100
  } else NA

  merged[, sig_class := factor(sig_class,
    levels = c("Neither", "Bulk only", "Hepatocyte only", "Both significant"))]
  setorder(merged, sig_class)

  sig_colors <- c(
    "Both significant"  = masld_colors$up,
    "Bulk only"         = "#42A5F5",
    "Hepatocyte only"   = "#7B1FA2",
    "Neither"           = masld_colors$ns
  )

  # Label top discordant/concordant genes
  top_genes <- merged[sig_class == "Both significant"][order(-abs(dream_logFC))][1:min(15, .N)]

  p_a <- ggplot(merged, aes(x = dream_logFC, y = hep_logFC, color = sig_class)) +
    rasterize_layer(geom_point(size = 0.15, alpha = 0.3, stroke = 0, shape = 16)) +
    scale_color_manual(values = sig_colors) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.3, color = "grey40") +
    geom_smooth(data = merged, aes(x = dream_logFC, y = hep_logFC),
                method = "lm", linewidth = 0.4, color = "black",
                se = FALSE, inherit.aes = FALSE) +
    geom_text_repel(data = top_genes, aes(label = symbol),
                    size = 1.8, max.overlaps = 20, color = "black",
                    segment.size = 0.2) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("rho = %.3f\nr = %.3f\nn = %s shared\nDir. conc. = %.0f%%",
                             rho, r, format(n_both, big.mark = ","),
                             ifelse(is.na(dir_conc), 0, dir_conc)),
             hjust = -0.1, vjust = 1.3, size = 2, color = "gray20") +
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
    theme_masld() +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.2, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(x = expression("Bulk integrated log"[2]*"FC"),
         y = expression("Hepatocyte pseudobulk log"[2]*"FC"),
         color = NULL, title = "Hepatocyte sc-DE vs bulk concordance")
} else {
  p_a <- placeholder("Hepatocyte pseudobulk DE not available")
}

# ==========================================================================
# Panel b: Per-cell-type overlap with bulk dream DEGs
# ==========================================================================
de_files <- list.files(pb_de_dir, pattern = "_de\\.csv$", full.names = TRUE)

if (length(de_files) > 0 && !is.null(dream)) {
  bulk_degs <- dream[is_dream_deg(dream), symbol]

  overlap_dt <- rbindlist(lapply(de_files, function(f) {
    ct <- gsub("_de\\.csv$", "", basename(f))
    ct_clean <- gsub("_", " ", ct)
    dt <- fread(f)
    # Skip the all-cell pseudobulk matrix (wide, no per-celltype logFC/padj pair).
    if (!all(c("logFC", "padj") %in% names(dt))) return(NULL)
    dt <- add_symbols(dt, "gene")
    ct_degs <- dt[padj < 0.05, symbol]
    n_ct <- length(ct_degs)
    n_overlap <- length(intersect(ct_degs, bulk_degs))
    n_ct_only <- length(setdiff(ct_degs, bulk_degs))
    jaccard <- length(intersect(ct_degs, bulk_degs)) /
               length(union(ct_degs, bulk_degs))
    data.table(cell_type = ct_clean, n_total = n_ct,
               n_overlap = n_overlap, n_unique = n_ct_only,
               jaccard = jaccard)
  }))

  overlap_dt <- overlap_dt[n_total > 0]
  overlap_dt[, cell_type := factor(cell_type,
    levels = overlap_dt[order(n_overlap)]$cell_type)]

  overlap_long <- melt(overlap_dt,
    id.vars = c("cell_type", "n_total", "jaccard"),
    measure.vars = c("n_overlap", "n_unique"),
    variable.name = "type", value.name = "n")
  overlap_long[, type := fifelse(type == "n_overlap",
                                  "Shared with bulk", "Cell-type specific")]

  p_b <- ggplot(overlap_long, aes(x = n, y = cell_type, fill = type)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c("Shared with bulk" = masld_colors$down,
                                  "Cell-type specific" = masld_colors$up)) +
    theme_masld() +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.2, "cm")) +
    labs(x = "Number of DEGs (padj < 0.05)", y = NULL, fill = NULL,
         title = "Pseudobulk DEG overlap with bulk integrated")
} else {
  p_b <- placeholder("Pseudobulk DE files not found")
}

# ==========================================================================
# Panel c: Direction concordance per cell type
# ==========================================================================
if (length(de_files) > 0 && !is.null(dream)) {
  conc_dt <- rbindlist(lapply(de_files, function(f) {
    ct <- gsub("_de\\.csv$", "", basename(f))
    ct_clean <- gsub("_", " ", ct)
    dt <- fread(f)
    # Skip non-per-celltype DE tables that lack standard logFC/padj columns
    # (e.g. allcell_pseudobulk_de.csv uses `lfc` + a wide cell-type matrix).
    if (!all(c("logFC", "padj") %in% names(dt))) return(NULL)
    dt <- add_symbols(dt, "gene")
    dream_cols <- intersect(c("symbol", "dream_logFC", "dream_padj",
                              "dream_shrunk_logFC", "dream_lfsr"), names(dream))
    m <- merge(dt[, .(symbol, ct_logFC = logFC, ct_padj = padj)],
               dream[, ..dream_cols],
               by = "symbol")
    shared <- m[ct_padj < 0.05 & is_dream_deg(m)]
    if (nrow(shared) == 0) return(NULL)
    rho <- cor(shared$ct_logFC, shared$dream_logFC, method = "spearman")
    dir_agree <- mean(sign(shared$ct_logFC) == sign(shared$dream_logFC)) * 100
    data.table(cell_type = ct_clean, rho = rho, dir_concordance = dir_agree,
               n_shared = nrow(shared))
  }))

  if (nrow(conc_dt) > 0) {
    conc_dt <- conc_dt[!is.na(rho)]
    conc_dt[, cell_type := factor(cell_type, levels = conc_dt[order(rho)]$cell_type)]

    p_c <- ggplot(conc_dt, aes(x = rho, y = cell_type)) +
      geom_segment(aes(x = 0, xend = rho, y = cell_type, yend = cell_type),
                   linewidth = 0.4, color = "gray60") +
      geom_point(aes(size = n_shared), color = masld_colors$hep_intrinsic, shape = 16) +
      geom_text(aes(label = sprintf("%.0f%%", dir_concordance)),
                hjust = -0.3, size = 2, color = "gray30") +
      scale_size_continuous(range = c(1, 4), name = "Shared DEGs") +
      scale_x_continuous(limits = c(-0.2, 1.1)) +
      theme_masld() +
      theme(legend.position = "right",
            legend.key.size = unit(0.25, "cm")) +
      labs(x = expression("Spearman " * rho * " (vs bulk integrated)"),
           y = NULL,
           title = "LFC concordance per cell type")
  } else {
    p_c <- placeholder("No shared significant genes")
  }
} else {
  p_c <- placeholder("Data not available")
}

# ==========================================================================
# Panel d: sc-derived hepatocyte fraction vs MuSiC estimates
# ==========================================================================
sc_props_file <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/cell_type_proportions.csv")
deconv_dir <- file.path(BASE, "Analysis/Deconvolution/results")

if (file.exists(sc_props_file)) {
  sc_props <- fread(sc_props_file)

  # Compute per-dataset mean hepatocyte fraction from sc
  ct_cols <- setdiff(names(sc_props), c("sample", "dataset"))
  if ("Hepatocytes" %in% ct_cols) {
    sc_hep <- sc_props[, .(sc_hep_mean = mean(Hepatocytes, na.rm = TRUE),
                           sc_hep_sd = sd(Hepatocytes, na.rm = TRUE),
                           n_samples = .N), by = dataset]

    # Load MuSiC estimates per dataset
    music_hep <- rbindlist(lapply(list.dirs(deconv_dir, recursive = FALSE), function(d) {
      ds <- basename(d)
      f <- file.path(d, paste0(ds, "_music_prop_weighted.tsv"))
      if (!file.exists(f)) return(NULL)
      props <- fread(f)
      if ("Hepatocytes" %in% names(props)) {
        data.table(dataset = ds,
                   music_hep_mean = mean(props$Hepatocytes, na.rm = TRUE),
                   music_hep_sd = sd(props$Hepatocytes, na.rm = TRUE))
      }
    }))

    if (nrow(music_hep) > 0) {
      comp <- merge(sc_hep, music_hep, by = "dataset")
      r_val <- cor(comp$sc_hep_mean, comp$music_hep_mean, method = "pearson")

      p_d <- ggplot(comp, aes(x = sc_hep_mean * 100, y = music_hep_mean * 100)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    linewidth = 0.3, color = "gray50") +
        geom_errorbar(aes(ymin = (music_hep_mean - music_hep_sd) * 100,
                          ymax = (music_hep_mean + music_hep_sd) * 100),
                      width = 1, linewidth = 0.3, color = "gray60") +
        geom_errorbarh(aes(xmin = (sc_hep_mean - sc_hep_sd) * 100,
                           xmax = (sc_hep_mean + sc_hep_sd) * 100),
                       height = 1, linewidth = 0.3, color = "gray60") +
        geom_point(size = 2.5, color = masld_colors$hep_intrinsic, shape = 16) +
        geom_text_repel(aes(label = dataset), size = 2, max.overlaps = 20) +
        annotate("text", x = -Inf, y = Inf,
                 label = sprintf("r = %.3f", r_val),
                 hjust = -0.1, vjust = 1.5, size = 2.5) +
        coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
        theme_masld() +
        labs(x = "sc-derived hepatocyte fraction (%)",
             y = "MuSiC-estimated hepatocyte fraction (%)",
             title = "sc vs deconvolution hepatocyte estimate")
    } else {
      p_d <- placeholder("MuSiC results not found")
    }
  } else {
    p_d <- placeholder("Hepatocytes not in sc proportions")
  }
} else {
  p_d <- placeholder("sc proportions file not found")
}

# ==========================================================================
# Assemble
# ==========================================================================
fig <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(fig, OUT, width = fig_full_width, height = 8)
message("Single-cell validation figure saved to ", OUT)
