#!/usr/bin/env Rscript
# 29_compare_three_methods.R
# ---------------------------------------------------------------------------
# Three-way deconvolution method comparison across the 5 canonical MASLD
# "pooled (cohort-adjusted)" cohorts:  MuSiC vs BayesPrism vs Rectangle.
#
# Extends the load_props() pattern of 16_compare_methods.R from two methods to
# three, and adds quantitative cross-method agreement + composition diagnostics.
#
# READS existing per-cohort proportion TSVs (never modifies them) from
#   Analysis/Deconvolution/results/{DATASET}/
#     {DATASET}_music_prop_weighted.tsv     (MuSiC,      closed simplex, 16 CT)
#     {DATASET}_bayesprism_proportions.tsv  (BayesPrism, closed simplex, 16 CT)
#     {DATASET}_rectangle_proportions.tsv   (Rectangle,  open simplex,   16 CT)  <- produced later
#     {DATASET}_rectangle_unknown.tsv       (Rectangle Unknown fraction)         <- produced later
# and WRITES all outputs (CSV + PDF) to a NEW dir:
#   Analysis/Deconvolution/rectangle_comparison/
#
# Missing files are guarded (warn + skip) exactly like 16_compare_methods.R, so
# this runs today on MuSiC+BayesPrism and gains the Rectangle arm once its
# proportions exist.
#
# Outputs:
#   three_method_concordance.csv          (tidy: per cell type x method-pair:
#                                          Pearson r, Spearman rho, Lin's CCC, RMSE, n)
#   three_method_concordance_heatmap.pdf  (cell type x pair, one page per metric)
#   scatter_rectangle_vs_methods.pdf      (Rectangle vs BayesPrism / MuSiC, faceted by CT)
#   bland_altman_hep_mac.pdf              (Hepatocyte + Macrophage fractions)
#   mean_composition_stacked.pdf          (mean composition per cohort, per method)
#   rectangle_unknown_by_cohort.pdf       (Rectangle Unknown-fraction distribution)
#   rectangle_unknown_summary.csv         (per-cohort Unknown-mass summary)
#
# Figure rules: PDF only (useDingbats=FALSE); base 6 pt Helvetica, no bold; all
# text black; captions via message() (never subtitle); control/reference = #9E9E9E;
# no lollipops; no 3D.  Env: rnaseq (ggplot2 / dplyr / tidyr / data.table).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

pdf.options(useDingbats = FALSE)

# ---------------------------------------------------------------------------
# Paths + constants
# ---------------------------------------------------------------------------
project_root <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RESULTS_BASE <- file.path(project_root, "Analysis/Deconvolution/results")
# Output dir is overridable (for a preview/dry run) but defaults to the
# deliverable dir. NEW dir only — never writes over existing results.
OUTPUT_DIR <- Sys.getenv("RECT_CMP_OUTDIR",
  unset = file.path(project_root, "Analysis/Deconvolution/rectangle_comparison"))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 5 canonical pooled cohorts (control-bearing; contribute to the Disease-vs-Control analysis)
MEGA_COHORTS <- c("GSE126848", "GSE135251", "GSE130970", "GSE213621", "GSE162694")

# 16 canonical cell-type strings (shared contract; identical column order across methods)
CANON_CT <- c("Endothelial cells", "Hepatocytes", "Plasma cells", "T cells",
              "Cholangiocytes", "Fibroblasts", "Macrophages",
              "Circulating NK/NKT", "Resident NK", "Mono+mono derived cells",
              "Basophils", "B cells", "cDC1s", "cDC2s", "pDCs", "Neutrophils")

METHODS      <- c("music", "bayesprism", "rectangle")
METHOD_LABEL <- c(music = "MuSiC", bayesprism = "BayesPrism", rectangle = "Rectangle")
CONTROL_GREY <- "#9E9E9E"   # reference / identity lines + Unknown slice

# Per-cell-type fill palette (data encoding, not text -> color is allowed here).
# Deterministic 16-colour categorical set; Hepatocytes anchored to a stable hue.
CT_PALETTE <- setNames(
  grDevices::colorRampPalette(
    c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
      "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
  )(length(CANON_CT)),
  CANON_CT)

# ---------------------------------------------------------------------------
# Loaders (guard missing files: warn + skip, per 16_compare_methods.R)
# ---------------------------------------------------------------------------
prop_filename <- function(dataset, method) {
  switch(method,
    music      = paste0(dataset, "_music_prop_weighted.tsv"),
    bayesprism = paste0(dataset, "_bayesprism_proportions.tsv"),
    rectangle  = paste0(dataset, "_rectangle_proportions.tsv"),
    stop("unknown method: ", method))
}

# Long-format proportions for one (dataset, method); NULL if file missing.
load_props <- function(dataset, method) {
  fp <- file.path(RESULTS_BASE, dataset, prop_filename(dataset, method))
  if (!file.exists(fp)) {
    warning(sprintf("Missing %s proportions: %s", METHOD_LABEL[[method]], fp))
    return(NULL)
  }
  props <- read.delim(fp, row.names = 1, check.names = FALSE)
  keep <- intersect(CANON_CT, colnames(props))
  if (length(keep) == 0L) {
    warning(sprintf("No canonical cell-type columns in %s", fp))
    return(NULL)
  }
  props <- props[, keep, drop = FALSE]
  props$Sample  <- rownames(props)
  props$Method  <- method
  props$Dataset <- dataset
  tidyr::pivot_longer(props, cols = dplyr::all_of(keep),
                      names_to = "CellType", values_to = "Proportion")
}

# Rectangle Unknown fraction for one dataset; NULL if missing. Robust to either a
# rownames-indexed single column or a (sample_id, Unknown) two-column layout.
load_unknown <- function(dataset) {
  fp <- file.path(RESULTS_BASE, dataset, paste0(dataset, "_rectangle_unknown.tsv"))
  if (!file.exists(fp)) {
    warning(sprintf("Missing Rectangle Unknown fractions: %s", fp))
    return(NULL)
  }
  u <- read.delim(fp, check.names = FALSE, stringsAsFactors = FALSE)
  if ("Unknown" %in% colnames(u)) {
    sid <- if ("sample_id" %in% colnames(u)) u[["sample_id"]] else u[[1]]
    return(data.frame(Sample = as.character(sid),
                      Unknown = as.numeric(u[["Unknown"]]),
                      Dataset = dataset, stringsAsFactors = FALSE))
  }
  # Fallback: rownames-indexed single value column
  u2 <- read.delim(fp, row.names = 1, check.names = FALSE)
  data.frame(Sample = rownames(u2), Unknown = as.numeric(u2[[1]]),
             Dataset = dataset, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Lin's Concordance Correlation Coefficient (DescTools::CCC unavailable in rnaseq)
#   CCC = 2 * cov(x,y) / (var(x) + var(y) + (mean(x) - mean(y))^2)   [population moments]
# ---------------------------------------------------------------------------
lin_ccc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 3L) return(NA_real_)
  mx <- mean(x); my <- mean(y)
  vx  <- mean((x - mx)^2)
  vy  <- mean((y - my)^2)
  sxy <- mean((x - mx) * (y - my))
  denom <- vx + vy + (mx - my)^2
  if (denom == 0) return(NA_real_)
  2 * sxy / denom
}

rmse <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 1L) return(NA_real_)
  sqrt(mean((x[ok] - y[ok])^2))
}

safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  # guard zero-variance vectors (e.g. a cell type never detected by a method)
  if (stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

# ---------------------------------------------------------------------------
# Load all methods x cohorts, assemble long + wide tables
# ---------------------------------------------------------------------------
message("Loading proportions for ", length(MEGA_COHORTS), " cohorts x ",
        length(METHODS), " methods ...")
long_list <- list()
for (ds in MEGA_COHORTS) {
  for (m in METHODS) {
    long_list[[paste(ds, m, sep = "__")]] <- load_props(ds, m)
  }
}
long <- data.table::rbindlist(long_list, use.names = TRUE, fill = TRUE)

methods_present <- if (nrow(long) > 0) intersect(METHODS, unique(long$Method)) else character(0)
if (length(methods_present) < 2L) {
  stop("Need >=2 methods with data to compare; found: ",
       paste(methods_present, collapse = ", "),
       ". (Rectangle proportions are produced later; MuSiC+BayesPrism should exist.)")
}
message("Methods with data: ",
        paste(METHOD_LABEL[methods_present], collapse = ", "))
if (!"rectangle" %in% methods_present) {
  message("NOTE: Rectangle proportions not found yet -> Rectangle panels/pairs ",
          "are skipped; MuSiC-vs-BayesPrism reference comparison still produced.")
}

# Wide: one row per (Dataset, Sample, CellType); one column per method.
wide <- tidyr::pivot_wider(long, id_cols = c("Dataset", "Sample", "CellType"),
                           names_from = "Method", values_from = "Proportion")
wide <- as.data.table(wide)
for (m in METHODS) if (!m %in% colnames(wide)) wide[[m]] <- NA_real_

message("Matched (sample x cell-type) rows: ", nrow(wide),
        " across ", length(unique(wide$Sample)), " samples.")

# ---------------------------------------------------------------------------
# 1. Cross-method concordance (per cell type + pooled) for each method-pair
# ---------------------------------------------------------------------------
PAIRS <- list(
  c("rectangle", "bayesprism"),
  c("rectangle", "music"),
  c("bayesprism", "music")
)

concordance_rows <- list()
for (pr in PAIRS) {
  a <- pr[1]; b <- pr[2]
  pair_label <- sprintf("%s vs %s", METHOD_LABEL[[a]], METHOD_LABEL[[b]])
  if (!(a %in% methods_present && b %in% methods_present)) {
    message("Skip concordance pair (missing method): ", pair_label)
    next
  }
  for (ct in c("ALL cell types", CANON_CT)) {
    sub <- if (identical(ct, "ALL cell types")) wide else wide[CellType == ct]
    xa <- sub[[a]]; xb <- sub[[b]]
    ok <- is.finite(xa) & is.finite(xb)
    concordance_rows[[paste(pair_label, ct, sep = "||")]] <- data.table(
      pair         = pair_label,
      method_a     = METHOD_LABEL[[a]],
      method_b     = METHOD_LABEL[[b]],
      cell_type    = ct,
      n            = sum(ok),
      pearson_r    = safe_cor(xa, xb, "pearson"),
      spearman_rho = safe_cor(xa, xb, "spearman"),
      lin_ccc      = lin_ccc(xa, xb),
      rmse         = rmse(xa, xb)
    )
  }
}
concordance <- data.table::rbindlist(concordance_rows, use.names = TRUE, fill = TRUE)

if (nrow(concordance) > 0) {
  data.table::fwrite(concordance, file.path(OUTPUT_DIR, "three_method_concordance.csv"))
  message("Wrote three_method_concordance.csv (", nrow(concordance), " rows).")
} else {
  message("No concordance rows produced (insufficient overlapping methods).")
}

# ---------------------------------------------------------------------------
# Heatmap: cell type x pair, one page per metric
# ---------------------------------------------------------------------------
if (nrow(concordance) > 0) {
  metrics <- list(
    list(col = "pearson_r",    name = "Pearson r",    lo = 0, hi = 1, low = "white", high = "#1E88E5"),
    list(col = "spearman_rho", name = "Spearman rho", lo = 0, hi = 1, low = "white", high = "#1E88E5"),
    list(col = "lin_ccc",      name = "Lin's CCC",    lo = 0, hi = 1, low = "white", high = "#2E7D32"),
    list(col = "rmse",         name = "RMSE",         lo = NA, hi = NA, low = "white", high = "#C62828")
  )
  ct_levels <- rev(c("ALL cell types", CANON_CT))
  hm <- copy(concordance)
  hm[, cell_type := factor(cell_type, levels = ct_levels)]
  hm[, pair := factor(pair, levels = unique(concordance$pair))]

  pdf(file.path(OUTPUT_DIR, "three_method_concordance_heatmap.pdf"), width = 6.5, height = 5.5)
  for (mt in metrics) {
    dd <- hm[, .(pair, cell_type, val = get(mt$col))]
    lims <- if (is.na(mt$lo)) range(dd$val, na.rm = TRUE) else c(mt$lo, mt$hi)
    p <- ggplot(dd, aes(x = pair, y = cell_type, fill = val)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = ifelse(is.na(val), "", sprintf("%.2f", val))),
                size = 1.9, color = "black") +
      scale_fill_gradient(low = mt$low, high = mt$high, limits = lims,
                          na.value = CONTROL_GREY, name = mt$name) +
      labs(title = paste0("Cross-method agreement: ", mt$name), x = NULL, y = NULL) +
      theme_minimal(base_family = "Helvetica", base_size = 6) +
      theme(
        plot.title  = element_text(size = 7, face = "plain", color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 6, color = "black"),
        legend.text  = element_text(size = 6, color = "black"),
        panel.grid   = element_blank()
      )
    print(p)
  }
  invisible(dev.off())
  message("Wrote three_method_concordance_heatmap.pdf (one page per metric).")
}

# ---------------------------------------------------------------------------
# 2. Scatter panels: Rectangle vs each other method (faceted + colored by cell type)
#    plus BayesPrism-vs-MuSiC reference. Identity line in control grey.
# ---------------------------------------------------------------------------
scatter_pairs <- Filter(function(pr) pr[1] %in% methods_present && pr[2] %in% methods_present, PAIRS)
if (length(scatter_pairs) > 0) {
  pdf(file.path(OUTPUT_DIR, "scatter_rectangle_vs_methods.pdf"), width = 7.2, height = 6.8)
  for (pr in scatter_pairs) {
    a <- pr[1]; b <- pr[2]
    dd <- wide[is.finite(get(a)) & is.finite(get(b))]
    if (nrow(dd) == 0) next
    r_all   <- safe_cor(dd[[a]], dd[[b]], "pearson")
    ccc_all <- lin_ccc(dd[[a]], dd[[b]])
    p <- ggplot(dd, aes(x = .data[[b]], y = .data[[a]], color = CellType)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = CONTROL_GREY, linewidth = 0.3) +
      geom_point(alpha = 0.5, size = 0.4) +
      facet_wrap(~ CellType, scales = "free", ncol = 4) +
      scale_color_manual(values = CT_PALETTE, guide = "none") +
      labs(
        title = sprintf("%s (y) vs %s (x): per-cell-type fractions",
                        METHOD_LABEL[[a]], METHOD_LABEL[[b]]),
        x = paste(METHOD_LABEL[[b]], "proportion"),
        y = paste(METHOD_LABEL[[a]], "proportion")
      ) +
      theme_bw(base_family = "Helvetica", base_size = 6) +
      theme(
        plot.title  = element_text(size = 7, face = "plain", color = "black"),
        strip.text  = element_text(size = 5.5, color = "black"),
        strip.background = element_rect(fill = "grey92", color = NA),
        axis.text   = element_text(size = 5, color = "black"),
        axis.title  = element_text(size = 6, color = "black"),
        panel.grid.minor = element_blank()
      )
    print(p)
    message(sprintf("  scatter %s vs %s: n=%d, overall Pearson r=%.3f, CCC=%.3f",
                    METHOD_LABEL[[a]], METHOD_LABEL[[b]], nrow(dd),
                    ifelse(is.na(r_all), NA_real_, r_all),
                    ifelse(is.na(ccc_all), NA_real_, ccc_all)))
  }
  invisible(dev.off())
  message("Wrote scatter_rectangle_vs_methods.pdf.")
}

# ---------------------------------------------------------------------------
# 3. Bland-Altman (mean vs difference) for Hepatocyte + Macrophage fractions
#    for each Rectangle pair (and the BP-vs-MuSiC reference).
# ---------------------------------------------------------------------------
ba_pairs <- Filter(function(pr) pr[1] %in% methods_present && pr[2] %in% methods_present, PAIRS)
ba_cts   <- intersect(c("Hepatocytes", "Macrophages"), unique(wide$CellType))
ba_rows  <- list()
for (pr in ba_pairs) {
  a <- pr[1]; b <- pr[2]
  for (ct in ba_cts) {
    sub <- wide[CellType == ct & is.finite(get(a)) & is.finite(get(b))]
    if (nrow(sub) == 0) next
    ba_rows[[paste(a, b, ct, sep = "__")]] <- data.table(
      pair      = sprintf("%s - %s", METHOD_LABEL[[a]], METHOD_LABEL[[b]]),
      cell_type = ct,
      mean_frac = (sub[[a]] + sub[[b]]) / 2,
      diff_frac = sub[[a]] - sub[[b]]
    )
  }
}
if (length(ba_rows) > 0) {
  ba <- data.table::rbindlist(ba_rows, use.names = TRUE, fill = TRUE)
  ba_stats <- ba[, .(md = mean(diff_frac), sd = stats::sd(diff_frac)), by = .(pair, cell_type)]
  ba_stats[, `:=`(loa_lo = md - 1.96 * sd, loa_hi = md + 1.96 * sd)]
  p_ba <- ggplot(ba, aes(x = mean_frac, y = diff_frac)) +
    geom_hline(yintercept = 0, color = CONTROL_GREY, linewidth = 0.3) +
    geom_point(alpha = 0.5, size = 0.5, color = "black") +
    geom_hline(data = ba_stats, aes(yintercept = md), linewidth = 0.3, color = "#C62828") +
    geom_hline(data = ba_stats, aes(yintercept = loa_lo), linetype = "dashed",
               linewidth = 0.3, color = CONTROL_GREY) +
    geom_hline(data = ba_stats, aes(yintercept = loa_hi), linetype = "dashed",
               linewidth = 0.3, color = CONTROL_GREY) +
    facet_grid(cell_type ~ pair) +
    labs(title = "Bland-Altman: Hepatocyte & Macrophage fractions",
         x = "Mean of the two methods' proportion",
         y = "Difference (method A - method B)") +
    theme_bw(base_family = "Helvetica", base_size = 6) +
    theme(
      plot.title = element_text(size = 7, face = "plain", color = "black"),
      strip.text = element_text(size = 6, color = "black"),
      strip.background = element_rect(fill = "grey92", color = NA),
      axis.text  = element_text(size = 5, color = "black"),
      axis.title = element_text(size = 6, color = "black"),
      panel.grid.minor = element_blank()
    )
  pdf(file.path(OUTPUT_DIR, "bland_altman_hep_mac.pdf"),
      width = 1.8 + 1.8 * length(unique(ba$pair)), height = 4.2)
  print(p_ba)
  invisible(dev.off())
  message("Wrote bland_altman_hep_mac.pdf.")
} else {
  message("Bland-Altman skipped (no Hepatocyte/Macrophage overlap across a method pair).")
}

# ---------------------------------------------------------------------------
# 4. Stacked composition bars: mean composition per cohort, per method.
#    Rectangle gets an explicit Unknown slice (control grey) so open- vs closed-
#    simplex methods are visually comparable (all sum to ~1).
# ---------------------------------------------------------------------------
comp <- long[, .(mean_prop = mean(Proportion, na.rm = TRUE)),
             by = .(Dataset, Method, CellType)]

# Append Rectangle Unknown mean as an extra slice, if available.
if ("rectangle" %in% methods_present) {
  unk_list <- lapply(MEGA_COHORTS, load_unknown)
  unk <- data.table::rbindlist(Filter(Negate(is.null), unk_list), use.names = TRUE, fill = TRUE)
  if (nrow(unk) > 0) {
    unk_mean <- as.data.table(unk)[, .(mean_prop = mean(Unknown, na.rm = TRUE)), by = .(Dataset)]
    unk_mean[, `:=`(Method = "rectangle", CellType = "Unknown")]
    comp <- rbind(comp, unk_mean[, .(Dataset, Method, CellType, mean_prop)], use.names = TRUE)
  }
}

ct_fill <- c(CT_PALETTE, Unknown = CONTROL_GREY)
comp[, CellType := factor(CellType, levels = c(CANON_CT, "Unknown"))]
comp[, MethodLab := factor(METHOD_LABEL[Method], levels = METHOD_LABEL[methods_present])]

p_comp <- ggplot(comp, aes(x = Dataset, y = mean_prop, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  facet_wrap(~ MethodLab, nrow = 1) +
  scale_fill_manual(values = ct_fill, name = "Cell type") +
  labs(title = "Mean cell-type composition per cohort, by deconvolution method",
       x = NULL, y = "Mean proportion") +
  theme_bw(base_family = "Helvetica", base_size = 6) +
  theme(
    plot.title  = element_text(size = 7, face = "plain", color = "black"),
    strip.text  = element_text(size = 6, color = "black"),
    strip.background = element_rect(fill = "grey92", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    axis.title  = element_text(size = 6, color = "black"),
    legend.title = element_text(size = 6, color = "black"),
    legend.text  = element_text(size = 5, color = "black"),
    legend.key.size = unit(0.25, "cm"),
    panel.grid.minor = element_blank()
  )
pdf(file.path(OUTPUT_DIR, "mean_composition_stacked.pdf"),
    width = 2.0 + 1.5 * length(methods_present), height = 4.2)
print(p_comp)
invisible(dev.off())
message("Wrote mean_composition_stacked.pdf.")

# ---------------------------------------------------------------------------
# 5. Rectangle Unknown-content summary (distribution per cohort + note)
# ---------------------------------------------------------------------------
if ("rectangle" %in% methods_present) {
  unk_list <- lapply(MEGA_COHORTS, load_unknown)
  unk <- data.table::rbindlist(Filter(Negate(is.null), unk_list), use.names = TRUE, fill = TRUE)
  if (nrow(unk) > 0) {
    unk <- as.data.table(unk)
    unk_summary <- unk[, .(
      n      = .N,
      mean   = mean(Unknown, na.rm = TRUE),
      median = stats::median(Unknown, na.rm = TRUE),
      sd     = stats::sd(Unknown, na.rm = TRUE),
      min    = min(Unknown, na.rm = TRUE),
      max    = max(Unknown, na.rm = TRUE)
    ), by = .(Dataset)]
    data.table::fwrite(unk_summary, file.path(OUTPUT_DIR, "rectangle_unknown_summary.csv"))
    message("Wrote rectangle_unknown_summary.csv.")

    p_unk <- ggplot(unk, aes(x = Dataset, y = Unknown)) +
      geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 0.4, color = CONTROL_GREY) +
      stat_summary(fun = mean, geom = "point", size = 1.2, color = "black") +
      stat_summary(fun = mean, geom = "errorbar",
                   fun.min = function(z) mean(z) - stats::sd(z),
                   fun.max = function(z) mean(z) + stats::sd(z),
                   width = 0.25, linewidth = 0.3, color = "black") +
      labs(title = "Rectangle Unknown fraction per cohort (mean +/- SD)",
           x = NULL, y = "Unknown fraction") +
      theme_bw(base_family = "Helvetica", base_size = 6) +
      theme(
        plot.title  = element_text(size = 7, face = "plain", color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color = "black"),
        axis.text.y = element_text(size = 5, color = "black"),
        axis.title  = element_text(size = 6, color = "black"),
        panel.grid.minor = element_blank()
      )
    pdf(file.path(OUTPUT_DIR, "rectangle_unknown_by_cohort.pdf"), width = 4.0, height = 3.4)
    print(p_unk)
    invisible(dev.off())
    message("Wrote rectangle_unknown_by_cohort.pdf.")

    grand_mean <- mean(unk$Unknown, na.rm = TRUE)
    message(sprintf(
      "Rectangle assigns on average %.1f%% of mass to Unknown (range %.1f-%.1f%% across cohorts). MuSiC and BayesPrism are closed-simplex (Unknown == 0 by construction), so the 16 cell-type fractions absorb this mass and are not directly comparable to Rectangle's re-normalised open simplex.",
      100 * grand_mean,
      100 * min(unk_summary$mean, na.rm = TRUE),
      100 * max(unk_summary$mean, na.rm = TRUE)))
  } else {
    message("Rectangle Unknown files not found -> Unknown summary skipped.")
  }
} else {
  message("Rectangle absent -> Unknown-content summary skipped (produced once Rectangle runs).")
}

message("\n29_compare_three_methods.R complete. Outputs in: ", OUTPUT_DIR)
