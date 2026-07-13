##############################################################################
# Supplementary Figure: COLOC Threshold Sensitivity Analysis
#
# 6 panels (3 rows x 2 cols):
#   (a) p12 prior sensitivity heatmap — gene count stability
#   (b) PP.H4 rank stability scatter — across p12 priors
#   (c) Min SNP threshold sensitivity — line plot per PP.H4 threshold
#   (d) Palindromic SNP removal — distribution histogram
#   (e) PP.H4 distribution — with MHC and LD cluster annotations
#   (f) SNP count distribution — in COLOC windows vs hits
#
# Demonstrates robustness of COLOC results to analytic parameter choices.
#
# Input files:
#   GWAS/finemapping/results/coloc_threshold_audit/p12_sensitivity.csv
#   GWAS/finemapping/results/coloc_threshold_audit/maf_palindromic_audit.csv
#   GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv
#   GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#
# Output:
#   figures/supplementary/figS04_coloc/figS_coloc_threshold_sensitivity.pdf
#   figures/supplementary/figS04_coloc/panelX_*.pdf  (individual panels)
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
  library(viridis)
  library(RColorBrewer)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# publication_theme.R defines placeholder() — but define locally as fallback
if (!exists("placeholder")) {
  placeholder <- function(label) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = label, size = GEOM_TEXT_6PT,
               color = "black") +
      theme_void()
  }
}

FINEMAPPING   <- file.path(BASE, "GWAS/finemapping/results")
AUDIT_DIR     <- file.path(FINEMAPPING, "coloc_threshold_audit")
SUSIE_DIR     <- file.path(FINEMAPPING, "susie_coloc")

OUT       <- file.path(FIGS04_DIR, "figS_coloc_threshold_sensitivity.pdf")
PANEL_DIR <- file.path(FIGS04_DIR, "panels")
dir.create(FIGS04_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Publication color constants
# ---------------------------------------------------------------------------
COL_CURRENT   <- "#C2185B"   # highlight current setting (magenta)
COL_CONSERV   <- "#1565C0"   # conservative (blue)
COL_PERMISSIVE<- "#E91E63"   # permissive (pink)
COL_BOTH      <- "#9E9E9E"   # both / shared (grey)
COL_MHC       <- "#C62828"   # MHC region (red)
COL_LD        <- "#F57F17"   # LD cluster (amber)
COL_CLEAN     <- "#BDBDBD"   # clean / unflagged (grey)

# ---------------------------------------------------------------------------
# Shared data loading (done once, shared across panels)
# ---------------------------------------------------------------------------
cat("Loading data files...\n")

p12_f    <- file.path(AUDIT_DIR, "p12_sensitivity.csv")
pal_f    <- file.path(AUDIT_DIR, "maf_palindromic_audit.csv")
all_f    <- file.path(SUSIE_DIR, "susie_coloc_all_gwas.csv")
gene_f   <- file.path(SUSIE_DIR, "gene_level_coloc.csv")

p12_dat  <- if (file.exists(p12_f))  fread(p12_f)  else NULL
pal_dat  <- if (file.exists(pal_f))  fread(pal_f)  else NULL
all_dat  <- if (file.exists(all_f))  fread(all_f)  else NULL
gene_dat <- if (file.exists(gene_f)) fread(gene_f) else NULL

if (!is.null(p12_dat))  cat(sprintf("  p12_sensitivity.csv:       %d rows\n",  nrow(p12_dat)))
if (!is.null(pal_dat))  cat(sprintf("  maf_palindromic_audit.csv: %d rows\n",  nrow(pal_dat)))
if (!is.null(all_dat))  cat(sprintf("  susie_coloc_all_gwas.csv:  %d rows\n",  nrow(all_dat)))
if (!is.null(gene_dat)) cat(sprintf("  gene_level_coloc.csv:      %d rows\n",  nrow(gene_dat)))

# Pre-initialize all panels as placeholders
p_a <- placeholder("(a) p12 prior sensitivity\n[data not found]")
p_b <- placeholder("(b) PP.H4 rank stability\n[data not found]")
p_c <- placeholder("(c) Min SNP threshold sensitivity\n[data not found]")
p_d <- placeholder("(d) Palindromic SNP burden\n[data not found]")
p_e <- placeholder("(e) PP.H4 distribution by region\n[data not found]")
p_f <- placeholder("(f) SNP count distribution\n[data not found]")

# ═══════════════════════════════════════════════════════════════════════════
# Panel (a): p12 prior sensitivity — gene count heatmap
# ═══════════════════════════════════════════════════════════════════════════
cat("Panel (a): p12 prior sensitivity heatmap...\n")
tryCatch({
  if (is.null(p12_dat)) stop("p12_sensitivity.csv not found")

  # PP.H4 thresholds of interest
  pp4_thresholds <- c(0.5, 0.8, 0.9)
  # p12 columns available in the file
  p12_cols <- c("pp4_1e7", "pp4_1e6", "pp4_5e6", "pp4_1e5", "pp4_5e5")
  p12_vals <- c(1e-7, 1e-6, 5e-6, 1e-5, 5e-5)
  p12_labels <- c("1×10⁻⁷", "1×10⁻⁶", "5×10⁻⁶\n(current)", "1×10⁻⁵", "5×10⁻⁵")

  # Build heatmap data: for each (pp4_threshold, p12_col), count unique genes
  heat_rows <- list()
  for (thr in pp4_thresholds) {
    for (i in seq_along(p12_cols)) {
      col <- p12_cols[i]
      if (!col %in% names(p12_dat)) next
      n_genes <- p12_dat[get(col) >= thr, uniqueN(ensembl)]
      heat_rows[[length(heat_rows) + 1]] <- data.table(
        pp4_thr = thr,
        p12_val = p12_vals[i],
        p12_lab = p12_labels[i],
        p12_col = col,
        n_genes = n_genes,
        is_current = (p12_vals[i] == 5e-6)
      )
    }
  }
  heat <- rbindlist(heat_rows)

  # Ordered factors for display
  heat[, pp4_thr_f := factor(paste0("PP.H4 ≥ ", pp4_thr),
                             levels = paste0("PP.H4 ≥ ", rev(pp4_thresholds)))]
  heat[, p12_f := factor(p12_lab, levels = p12_labels)]

  # Highlight current column (p12 = 5e-6)
  current_col_lab <- p12_labels[p12_vals == 5e-6]

  p_a <- ggplot(heat, aes(x = p12_f, y = pp4_thr_f, fill = n_genes)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_tile(data = heat[is_current == TRUE],
              color = COL_CURRENT, fill = NA, linewidth = 1.2) +
    geom_text(aes(label = n_genes), size = GEOM_TEXT_6PT, fontface = "plain",
              color = ifelse(heat$n_genes > max(heat$n_genes) * 0.6, "white", "black")) +
    scale_fill_viridis_c(option = "mako", direction = -1,
                         name = "Unique\ngenes", labels = scales::comma) +
    labs(
      x = "p12 prior value",
      y = NULL
    ) +
    theme_masld(base_size = 7) +
    theme(
      axis.text.x  = element_text(size = 6, lineheight = 0.85),
      panel.border = element_blank(),
      axis.line    = element_blank(),
      axis.ticks   = element_blank(),
      legend.position = "right"
    )

  message("[caption] Panel (a): p12 prior sensitivity, gene counts. Box = current setting (p12 = 5x10^-6).")
  cat("  Panel (a) built successfully.\n")
}, error = function(e) {
  cat(sprintf("  Panel (a) FAILED: %s\n", conditionMessage(e)))
  p_a <<- placeholder(sprintf("(a) p12 prior sensitivity\nError: %s", conditionMessage(e)))
})

# ═══════════════════════════════════════════════════════════════════════════
# Panel (b): PP.H4 rank stability scatter across p12 values
# ═══════════════════════════════════════════════════════════════════════════
cat("Panel (b): PP.H4 rank stability scatter...\n")
tryCatch({
  if (is.null(p12_dat)) stop("p12_sensitivity.csv not found")

  # Use columns: pp4_5e6 (current), pp4_1e6 (conservative), pp4_1e5 (permissive)
  needed <- c("ensembl", "gene", "pp4_5e6", "pp4_1e6", "pp4_1e5")
  missing_cols <- setdiff(needed, names(p12_dat))
  if (length(missing_cols) > 0) stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))

  # One row per unique gene-GWAS window: average if multiple windows per gene
  scatter <- p12_dat[, .(
    pp4_current    = mean(pp4_5e6, na.rm = TRUE),
    pp4_conserv    = mean(pp4_1e6, na.rm = TRUE),
    pp4_permissive = mean(pp4_1e5, na.rm = TRUE)
  ), by = .(ensembl, gene)]

  # Classify concordance relative to 0.5 threshold
  scatter[, category_conserv := fcase(
    pp4_current >= 0.5 & pp4_conserv >= 0.5,   "Both ≥0.5",
    pp4_current >= 0.5 & pp4_conserv < 0.5,    "Only current",
    pp4_current < 0.5  & pp4_conserv >= 0.5,   "Only conservative",
    default = "Neither"
  )]
  scatter[, category_permissive := fcase(
    pp4_current >= 0.5 & pp4_permissive >= 0.5, "Both ≥0.5",
    pp4_current >= 0.5 & pp4_permissive < 0.5,  "Only current",
    pp4_current < 0.5  & pp4_permissive >= 0.5, "Only permissive",
    default = "Neither"
  )]

  # Compute Spearman rho
  rho_conserv <- cor(scatter$pp4_current, scatter$pp4_conserv,
                     method = "spearman", use = "complete.obs")
  rho_permissive <- cor(scatter$pp4_current, scatter$pp4_permissive,
                        method = "spearman", use = "complete.obs")

  cat(sprintf("    Spearman rho (current vs conservative p12): %.4f\n", rho_conserv))
  cat(sprintf("    Spearman rho (current vs permissive p12):   %.4f\n", rho_permissive))

  # Color map
  color_map <- c(
    "Both ≥0.5"        = COL_BOTH,
    "Only current"      = COL_CURRENT,
    "Only conservative" = COL_CONSERV,
    "Only permissive"   = COL_PERMISSIVE,
    "Neither"           = "#E0E0E0"
  )

  # Long format for overlay (two comparison layers)
  scatter_long <- rbind(
    scatter[, .(gene, ensembl, pp4_current, pp4_other = pp4_conserv,
                category = category_conserv, comparison = "p12 = 1×10⁻⁶ (conservative)")],
    scatter[, .(gene, ensembl, pp4_current, pp4_other = pp4_permissive,
                category = category_permissive, comparison = "p12 = 1×10⁻⁵ (permissive)")]
  )

  # Annotation text
  annot_df <- data.table(
    comparison = c("p12 = 1×10⁻⁶ (conservative)", "p12 = 1×10⁻⁵ (permissive)"),
    rho        = c(rho_conserv, rho_permissive),
    label      = sprintf("ρ = %.3f", c(rho_conserv, rho_permissive))
  )

  p_b <- ggplot(scatter_long[category != "Neither"],
                aes(x = pp4_current, y = pp4_other, color = category)) +
    geom_point(data = scatter_long[category == "Neither"],
               aes(x = pp4_current, y = pp4_other),
               color = "#E0E0E0", size = 0.5, alpha = 0.4, inherit.aes = FALSE) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "solid",
                color = "black", linewidth = 0.4) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = 0.5, linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    geom_text(data = annot_df,
              aes(x = 0.08, y = 0.91, label = label),
              size = GEOM_TEXT_6PT, color = "black", hjust = 0, inherit.aes = FALSE) +
    facet_wrap(~ comparison, nrow = 1) +
    scale_color_manual(values = color_map, name = NULL,
                       guide = guide_legend(override.aes = list(size = 2))) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1")) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1")) +
    labs(
      x = "PP.H4 (current p12 = 5×10⁻⁶)",
      y = "PP.H4 (alternative p12)"
    ) +
    theme_masld(base_size = 7) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.25, "cm"))

  message("[caption] Panel (b): PP.H4 rank stability across p12 priors. x-axis = current (p12 = 5x10^-6); dashed lines at PP.H4 = 0.5.")
  cat("  Panel (b) built successfully.\n")
}, error = function(e) {
  cat(sprintf("  Panel (b) FAILED: %s\n", conditionMessage(e)))
  p_b <<- placeholder(sprintf("(b) PP.H4 rank stability\nError: %s", conditionMessage(e)))
})

# ═══════════════════════════════════════════════════════════════════════════
# Panel (c): Minimum SNP threshold sensitivity — line plot
# ═══════════════════════════════════════════════════════════════════════════
cat("Panel (c): Min SNP threshold sensitivity...\n")
tryCatch({
  if (is.null(all_dat)) stop("susie_coloc_all_gwas.csv not found")

  needed_cols <- c("gene", "ensembl", "PP.H4.abf", "n_snps")
  missing_cols <- setdiff(needed_cols, names(all_dat))
  if (length(missing_cols) > 0) stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))

  pp4_thresholds <- c(0.5, 0.8, 0.9)
  min_snp_grid   <- c(10, 50, 100, 200, 500)
  current_snp    <- 100

  # For each (min_snps, pp4_threshold), count unique genes passing both
  snp_rows <- list()
  for (ms in min_snp_grid) {
    sub <- all_dat[n_snps >= ms]
    for (thr in pp4_thresholds) {
      n_genes <- sub[PP.H4.abf >= thr, uniqueN(ensembl)]
      snp_rows[[length(snp_rows) + 1]] <- data.table(
        min_snps  = ms,
        pp4_thr   = thr,
        n_genes   = n_genes,
        is_current = (ms == current_snp)
      )
    }
  }
  snp_dt <- rbindlist(snp_rows)
  snp_dt[, pp4_label := paste0("PP.H4 ≥ ", pp4_thr)]

  # Palette for the three PP.H4 thresholds
  thr_pal <- c(
    "PP.H4 ≥ 0.5" = "#1565C0",
    "PP.H4 ≥ 0.8" = "#E91E63",
    "PP.H4 ≥ 0.9" = "#880E4F"
  )

  # Annotation: gene counts at current threshold
  annot_snp <- snp_dt[is_current == TRUE]

  p_c <- ggplot(snp_dt, aes(x = min_snps, y = n_genes,
                             color = pp4_label, group = pp4_label)) +
    geom_vline(xintercept = current_snp, linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    geom_text(data = annot_snp,
              aes(x = min_snps, y = n_genes, label = n_genes, color = pp4_label),
              size = GEOM_TEXT_6PT, vjust = -0.9, fontface = "plain", show.legend = FALSE) +
    annotate("text", x = current_snp, y = max(snp_dt$n_genes) * 0.97,
             label = "current\n(n=100)", size = GEOM_TEXT_6PT, color = "black",
             hjust = -0.05, lineheight = 0.85) +
    scale_color_manual(values = thr_pal, name = NULL) +
    scale_x_continuous(breaks = min_snp_grid,
                       labels = as.character(min_snp_grid)) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.05, 0.12))) +
    labs(
      x = "Minimum SNPs per COLOC window",
      y = "Unique genes (any GWAS)"
    ) +
    theme_masld(base_size = 7) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.25, "cm"))

  message("[caption] Panel (c): sensitivity to minimum SNP count threshold per COLOC window.")
  cat("  Panel (c) built successfully.\n")
}, error = function(e) {
  cat(sprintf("  Panel (c) FAILED: %s\n", conditionMessage(e)))
  p_c <<- placeholder(sprintf("(c) Min SNP threshold\nError: %s", conditionMessage(e)))
})

# ═══════════════════════════════════════════════════════════════════════════
# Panel (d): Palindromic SNP removal — histogram of pct_palindromic
# ═══════════════════════════════════════════════════════════════════════════
cat("Panel (d): Palindromic SNP burden histogram...\n")
tryCatch({
  if (is.null(pal_dat)) stop("maf_palindromic_audit.csv not found")
  if (!"pct_palindromic" %in% names(pal_dat)) stop("pct_palindromic column missing")

  mean_pal <- mean(pal_dat$pct_palindromic, na.rm = TRUE)
  med_pal  <- median(pal_dat$pct_palindromic, na.rm = TRUE)

  # Total palindromic removed across all windows
  total_palindromic <- if ("n_palindromic" %in% names(pal_dat)) {
    sum(pal_dat$n_palindromic, na.rm = TRUE)
  } else NA_integer_

  subtitle_text <- if (!is.na(total_palindromic)) {
    sprintf("%s A/T and C/G variants removed across all windows",
            scales::comma(total_palindromic))
  } else {
    "A/T and C/G variants removed across all windows"
  }

  p_d <- ggplot(pal_dat[!is.na(pct_palindromic)],
                aes(x = pct_palindromic)) +
    geom_histogram(bins = 40, fill = "#1565C0", color = "white",
                   alpha = 0.85, linewidth = 0.2) +
    geom_vline(xintercept = mean_pal, linetype = "solid",
               color = COL_CURRENT, linewidth = 0.7) +
    geom_vline(xintercept = med_pal, linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    annotate("text",
             x = mean_pal + 0.3, y = Inf,
             label = sprintf("Mean: %.1f%%", mean_pal),
             size = GEOM_TEXT_6PT, color = COL_CURRENT, vjust = 1.4, hjust = 0,
             fontface = "plain") +
    annotate("text",
             x = med_pal - 0.3, y = Inf,
             label = sprintf("Median: %.1f%%", med_pal),
             size = GEOM_TEXT_6PT, color = "black", vjust = 2.8, hjust = 1) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
    labs(
      x = "Palindromic SNPs (%)",
      y = "COLOC windows"
    ) +
    theme_masld(base_size = 7)

  message(sprintf("[caption] Panel (d): palindromic SNP burden per COLOC window. %s", subtitle_text))
  cat("  Panel (d) built successfully.\n")
}, error = function(e) {
  cat(sprintf("  Panel (d) FAILED: %s\n", conditionMessage(e)))
  p_d <<- placeholder(sprintf("(d) Palindromic SNP burden\nError: %s", conditionMessage(e)))
})

# ═══════════════════════════════════════════════════════════════════════════
# Panel (e): PP.H4 distribution with MHC and LD cluster annotations
# ═══════════════════════════════════════════════════════════════════════════
cat("Panel (e): PP.H4 distribution by region type...\n")
tryCatch({
  if (is.null(gene_dat)) stop("gene_level_coloc.csv not found")

  needed_cols <- c("coloc_best_pp4", "is_mhc", "ld_cluster_flag")
  missing_cols <- setdiff(needed_cols, names(gene_dat))
  if (length(missing_cols) > 0) stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))

  # Filter to PP.H4 > 0.1 for visibility; classify by region type
  plot_e <- gene_dat[!is.na(coloc_best_pp4) & coloc_best_pp4 > 0.1]
  plot_e <- copy(plot_e)

  plot_e[, region_type := fcase(
    is_mhc == TRUE,                                          "MHC region",
    !is.na(ld_cluster_flag) & is_mhc == FALSE,              "LD cluster",
    default = "Unflaged"
  )]

  # Count by region type for legend label
  counts <- plot_e[, .N, by = region_type]
  count_map <- setNames(counts$N, counts$region_type)

  region_pal <- c(
    "Unflaged"   = COL_CLEAN,
    "LD cluster" = COL_LD,
    "MHC region" = COL_MHC
  )

  # Threshold lines
  thr_lines <- data.table(
    xint  = c(0.5, 0.8, 0.9),
    label = c("0.5", "0.8", "0.9")
  )

  p_e <- ggplot(plot_e, aes(x = coloc_best_pp4, fill = region_type)) +
    geom_histogram(bins = 50, alpha = 0.8, color = "white",
                   linewidth = 0.1, position = "stack") +
    geom_vline(data = thr_lines, aes(xintercept = xint),
               linetype = "dashed", color = "black",
               linewidth = 0.4, inherit.aes = FALSE) +
    geom_text(data = thr_lines,
              aes(x = xint, y = Inf, label = xint),
              vjust = 1.4, size = GEOM_TEXT_6PT, color = "black",
              inherit.aes = FALSE) +
    scale_fill_manual(
      values = region_pal,
      name   = NULL,
      labels = function(x) {
        n <- count_map[x]
        ifelse(!is.na(n), sprintf("%s (n=%s)", x, scales::comma(n)), x)
      }
    ) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    labs(
      x = "Best PP.H4 (any GWAS)",
      y = "Genes"
    ) +
    theme_masld(base_size = 7) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.25, "cm"))

  message("[caption] Panel (e): PP.H4 distribution by region type. Genes with best PP.H4 > 0.1; dashed lines at 0.5, 0.8, 0.9.")
  cat("  Panel (e) built successfully.\n")
}, error = function(e) {
  cat(sprintf("  Panel (e) FAILED: %s\n", conditionMessage(e)))
  p_e <<- placeholder(sprintf("(e) PP.H4 distribution\nError: %s", conditionMessage(e)))
})

# ═══════════════════════════════════════════════════════════════════════════
# Panel (f): SNP count distribution in COLOC windows vs hits
# ═══════════════════════════════════════════════════════════════════════════
cat("Panel (f): SNP count distribution in COLOC windows...\n")
tryCatch({
  if (is.null(all_dat)) stop("susie_coloc_all_gwas.csv not found")

  needed_cols <- c("PP.H4.abf", "n_snps")
  missing_cols <- setdiff(needed_cols, names(all_dat))
  if (length(missing_cols) > 0) stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))

  # Build three groups: all windows, PP.H4>0.5 hits, PP.H4>0.8 hits
  # De-duplicate by gene+gwas_name for per-window SNP counts
  dedup_key <- intersect(c("ensembl", "gwas_name"), names(all_dat))
  if (length(dedup_key) < 2) dedup_key <- c("gene", "gwas_name")

  windows_all  <- unique(all_dat[, c(dedup_key, "PP.H4.abf", "n_snps"), with = FALSE])

  snp_long <- rbind(
    windows_all[, .(n_snps, group = "All windows")],
    windows_all[PP.H4.abf >= 0.5, .(n_snps, group = "PP.H4 ≥ 0.5")],
    windows_all[PP.H4.abf >= 0.8, .(n_snps, group = "PP.H4 ≥ 0.8")]
  )
  snp_long[, group := factor(group,
                              levels = c("All windows", "PP.H4 ≥ 0.5", "PP.H4 ≥ 0.8"))]

  group_pal <- c(
    "All windows" = "#BDBDBD",
    "PP.H4 ≥ 0.5" = "#1565C0",
    "PP.H4 ≥ 0.8" = COL_CURRENT
  )

  # Trim at 99th percentile for display
  snp_p99 <- quantile(windows_all$n_snps, 0.99, na.rm = TRUE)

  # Median per group for annotation
  group_medians <- snp_long[, .(median_snps = median(n_snps, na.rm = TRUE),
                                 n = .N), by = group]
  cat(sprintf("    Group medians:\n"))
  group_medians[, cat(sprintf("      %s: median=%d (n=%s)\n",
                               group, round(median_snps), scales::comma(n))),
                by = group]

  p_f <- ggplot(snp_long[!is.na(n_snps) & n_snps <= snp_p99],
                aes(x = n_snps, fill = group, color = group)) +
    geom_histogram(bins = 50, alpha = 0.55, linewidth = 0.2,
                   position = "identity") +
    geom_vline(xintercept = 10,  linetype = "dotted",
               color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 100, linetype = "dashed",
               color = "black",  linewidth = 0.5) +
    annotate("text", x = 12, y = Inf,
             label = "n=10\n(old)", size = GEOM_TEXT_6PT, color = "black",
             vjust = 1.4, hjust = 0, lineheight = 0.85) +
    annotate("text", x = 103, y = Inf,
             label = "n=100\n(current)", size = GEOM_TEXT_6PT, color = "black",
             vjust = 1.4, hjust = 0, lineheight = 0.85) +
    scale_fill_manual(values  = group_pal, name = NULL) +
    scale_color_manual(values = group_pal, name = NULL, guide = "none") +
    scale_x_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.12))) +
    labs(
      x = "SNPs per COLOC window (trimmed at 99th pct.)",
      y = "Windows"
    ) +
    theme_masld(base_size = 7) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.25, "cm"))

  message("[caption] Panel (f): SNP count distribution in COLOC windows. Dashed = current threshold (n=100); dotted = old threshold (n=10).")
  cat("  Panel (f) built successfully.\n")
}, error = function(e) {
  cat(sprintf("  Panel (f) FAILED: %s\n", conditionMessage(e)))
  p_f <<- placeholder(sprintf("(f) SNP count distribution\nError: %s", conditionMessage(e)))
})

# ═══════════════════════════════════════════════════════════════════════════
# Save individual panels
# ═══════════════════════════════════════════════════════════════════════════
cat("Saving individual panels...\n")

panel_specs <- list(
  list(plot = p_a, name = "panelA_p12_heatmap",          w = 3.5, h = 2.8),
  list(plot = p_b, name = "panelB_rank_stability",         w = 5.0, h = 2.8),
  list(plot = p_c, name = "panelC_minsnp_sensitivity",     w = 3.5, h = 2.8),
  list(plot = p_d, name = "panelD_palindromic_histogram",  w = 3.5, h = 2.8),
  list(plot = p_e, name = "panelE_pp4_distribution",       w = 3.5, h = 2.8),
  list(plot = p_f, name = "panelF_snp_count_distribution", w = 3.5, h = 2.8)
)

for (spec in panel_specs) {
  tryCatch({
    out_path <- file.path(PANEL_DIR, paste0(spec$name, ".pdf"))
    pdf_dev  <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
    ggplot2::ggsave(out_path, spec$plot,
                    width  = spec$w,
                    height = spec$h,
                    dpi    = 300,
                    device = pdf_dev)
    cat(sprintf("  Saved: %s\n", basename(out_path)))
  }, error = function(e) {
    cat(sprintf("  WARNING: Failed to save %s: %s\n", spec$name, conditionMessage(e)))
  })
}

# ═══════════════════════════════════════════════════════════════════════════
# Assemble combined figure with patchwork
# ═══════════════════════════════════════════════════════════════════════════
cat("Assembling composite figure...\n")
message("[caption] COLOC threshold sensitivity analysis (panels a-f).")

combined <- tryCatch({
  (p_a | p_b) /
  (p_c | p_d) /
  (p_e | p_f) +
    patchwork::plot_annotation(
      tag_levels = "a",
      theme = theme(
        plot.tag   = element_text(size = 6, face = "plain", family = "Helvetica")
      )
    ) &
    theme(plot.margin = margin(3, 4, 3, 4))
}, error = function(e) {
  cat(sprintf("  Patchwork assembly FAILED: %s\n", conditionMessage(e)))
  NULL
})

# ═══════════════════════════════════════════════════════════════════════════
# Save composite figure
# ═══════════════════════════════════════════════════════════════════════════
cat(sprintf("Saving composite figure to:\n  %s\n", OUT))

tryCatch({
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
  if (!is.null(combined)) {
    ggplot2::ggsave(OUT, combined,
                    width  = fig_full_width,
                    height = 8.27,
                    dpi    = 300,
                    device = pdf_dev)
    cat(sprintf("SUCCESS: Composite figure saved (%.2f x 8.27 in).\n", fig_full_width))
  } else {
    # Fallback: save panels manually with pdf()
    pdf_dev(OUT, width = fig_full_width, height = 8.27)
    gridExtra::grid.arrange(
      ggplotGrob(p_a), ggplotGrob(p_b),
      ggplotGrob(p_c), ggplotGrob(p_d),
      ggplotGrob(p_e), ggplotGrob(p_f),
      ncol = 2
    )
    grDevices::dev.off()
    cat(sprintf("SUCCESS: Fallback composite figure saved.\n"))
  }
}, error = function(e) {
  cat(sprintf("ERROR: Failed to save composite figure: %s\n", conditionMessage(e)))
})

cat("\nfigS_coloc_threshold_sensitivity.R complete.\n")
