#!/usr/bin/env Rscript
# plot_stage_unique_degs_heatmap.R
# Heatmaps showing ONLY stage-unique DEGs, split by which stage they are unique to.
# Handles extreme count imbalance (NAS0 ~5k vs NAS1 ~2) by capping per-stage
# representation at TOP_N genes (ranked by tau specificity) while showing all
# genes for stages below that threshold.
#
# Produces 2 heatmaps:
#   1) NAS axis: columns = NAS0_vs_rest .. NAS7_vs_rest, rows = unique DEGs
#   2) Fibrosis axis: columns = F0_vs_rest .. F4_vs_rest, rows = unique DEGs
#
# Row-split by unique_stage. Stages with 0 unique genes are omitted.
# Row titles show: "NAS0 (4,868 — top 25 shown)" or "NAS1 (2)".

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
})

cat("=== Stage-Unique DEG Heatmaps ===\n")

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- file.path(FIG2_DIR, "panels")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_N  <- 25          # max genes per stage
PADJ_THRESH <- 0.1    # significance threshold used in stage_unique_signatures

# ── Load data ─────────────────────────────────────────────────────────────────
staging_dir <- file.path(INTEGRATION, "results/staging_classifier")

sigs <- fread(file.path(staging_dir, "stage_unique_signatures.csv"))
ovr_nas <- fread(file.path(staging_dir, "one_vs_rest_nas_dream.csv"))
ovr_fib <- fread(file.path(staging_dir, "one_vs_rest_fibrosis_dream.csv"))

cat("Loaded:", nrow(sigs), "stage-unique genes\n")
cat("OVR NAS:", nrow(ovr_nas), "rows |  OVR Fibrosis:", nrow(ovr_fib), "rows\n")

# ── Gene symbol mapping (uses multi-evidence atlas via load_figure_data.R) ────
sigs <- add_symbols(sigs, "gene")
ovr_nas <- add_symbols(ovr_nas, "gene")
ovr_fib <- add_symbols(ovr_fib, "gene")

# Shorten any remaining bare Ensembl IDs for display
shorten_ensg <- function(x) {
  is_ensg <- grepl("^ENSG", x)
  x[is_ensg] <- sub("^ENSG0+", "ENSG..", x[is_ensg])
  x
}

# ── Known MASLD markers for annotation ────────────────────────────────────────
target_markers <- c(
  "PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7", "GCKR",
  "CXCL10", "CCL2", "TREM2", "SPP1", "IL1B",
  "COL1A1", "TIMP1", "TGFB1", "ACTA2", "LUM",
  "FASN", "SCD", "PPARA", "CYP7A1", "AKR1B10"
)

# ── Color scheme (matches publication_theme.R) ────────────────────────────────
col_fun <- colorRamp2(
  c(-2.5, 0, 2.5),
  c(masld_colors$down, "white", masld_colors$up)
)

# ── Helper: select top-N per stage by tau, build matrix ───────────────────────
build_unique_matrix <- function(sigs_axis, ovr_dt, axis_name, stage_col) {
  # axis_name = "NAS" or "Fibrosis"
  # stage_col = "stage_label" in OVR data

  # Count total per stage for row titles
  stage_counts <- sigs_axis[, .N, by = unique_stage]
  setnames(stage_counts, "N", "total")

  # Select top-N per stage by tau (descending)
  selected <- sigs_axis[order(-tau),
                         head(.SD, TOP_N),
                         by = unique_stage]
  cat(axis_name, "- selected", nrow(selected), "genes across",
      uniqueN(selected$unique_stage), "stages\n")

  # Build logFC matrix from OVR dream: gene x stage_label
  ovr_sub <- ovr_dt[gene %in% selected$gene]
  mat_dt <- dcast(ovr_sub, gene + symbol ~ stage_label, value.var = "logFC",
                  fun.aggregate = function(x) x[1])

  # Deduplicate symbols: keep highest-tau version
  tau_map <- selected[, .(gene, tau)]
  mat_dt <- merge(mat_dt, tau_map, by = "gene", all.x = TRUE)
  mat_dt <- mat_dt[order(-tau)]
  mat_dt <- mat_dt[!duplicated(symbol)]
  mat_dt[, tau := NULL]

  # Add unique_stage mapping for row-split
  stage_map <- selected[, .(gene, unique_stage)]
  mat_dt <- merge(mat_dt, stage_map, by = "gene", all.x = TRUE)

  # Order columns naturally; drop reference stage (NAS0/F0)
  stage_labels <- sort(unique(ovr_dt$stage_label))
  stage_labels <- stage_labels[!stage_labels %in% c("NAS0_vs_rest", "F0_vs_rest")]
  stage_labels <- stage_labels[stage_labels %in% names(mat_dt)]

  # Build numeric matrix
  mat <- as.matrix(mat_dt[, ..stage_labels])
  rownames(mat) <- shorten_ensg(mat_dt$symbol)
  mat[is.na(mat)] <- 0

  # Prettify column names: "NAS0_vs_rest" -> "NAS 0"
  pretty_cols <- gsub("_vs_rest", "", colnames(mat))
  pretty_cols <- gsub("NAS", "NAS ", pretty_cols)
  pretty_cols <- gsub("^F", "F", pretty_cols)
  colnames(mat) <- pretty_cols

  # Build row-split factor with informative titles
  row_stage <- mat_dt$unique_stage
  # Create labels with counts
  stage_labels_pretty <- sapply(row_stage, function(s) {
    total <- stage_counts[unique_stage == s, total]
    shown <- min(total, TOP_N)
    short <- gsub("_vs_rest", "", s)
    if (total > TOP_N) {
      sprintf("%s  (n=%s; top %d shown)", short, format(total, big.mark = ","), TOP_N)
    } else {
      sprintf("%s  (n=%d)", short, total)
    }
  })

  # Ordered factor: stages in natural numeric order
  stage_order <- gsub("_vs_rest", "", sort(unique(row_stage)))
  label_order <- unique(stage_labels_pretty[match(sort(unique(row_stage)), row_stage)])
  row_split <- factor(stage_labels_pretty, levels = label_order)

  list(mat = mat, row_split = row_split, symbols = rownames(mat))
}

# ── Helper: draw heatmap ─────────────────────────────────────────────────────
draw_stage_unique_heatmap <- function(mat, row_split, out_file, title) {
  if (nrow(mat) == 0) { cat("SKIP (empty):", out_file, "\n"); return() }

  mat_cap <- pmax(pmin(mat, 2.5), -2.5)

  # Mark known MASLD genes on right side
  found_markers <- intersect(target_markers, rownames(mat_cap))
  right_ha <- NULL
  if (length(found_markers) > 0) {
    idx <- match(found_markers, rownames(mat_cap))
    right_ha <- rowAnnotation(
      genes = anno_mark(at = idx, labels = found_markers,
                        labels_gp = gpar(fontsize = 7, fontface = "italic"),
                        padding = unit(1, "mm"))
    )
  }

  # Direction annotation (left): up/down at peak stage
  peak_col <- apply(abs(mat_cap), 1, which.max)
  peak_val <- mat_cap[cbind(seq_len(nrow(mat_cap)), peak_col)]
  direction <- ifelse(peak_val > 0, "Up", "Down")
  dir_colors <- c("Up" = masld_colors$up, "Down" = masld_colors$down)
  ha_left <- rowAnnotation(
    Direction = direction,
    col = list(Direction = dir_colors),
    border = TRUE,
    simple_anno_size = unit(3, "mm"),
    annotation_name_gp = gpar(fontsize = 7)
  )

  # Dynamic height: ~2.5mm per gene + gaps + title space
  n_genes <- nrow(mat_cap)
  n_stages <- nlevels(row_split)
  ht_height <- max(6, n_genes * 0.1 + n_stages * 0.3 + 2)  # inches
  ht_width <- ncol(mat_cap) * 0.55 + 5  # inches

  pdf(out_file, width = ht_width, height = ht_height)
  ht <- Heatmap(
    mat_cap,
    name = "logFC",
    col = col_fun,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 5),
    show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 9),
    left_annotation = ha_left,
    right_annotation = right_ha,
    row_split = row_split,
    row_title_gp = gpar(fontsize = 7),
    row_title_rot = 0,
    row_gap = unit(1.5, "mm"),
    border = TRUE,
    column_title = title,
    column_title_gp = gpar(fontsize = 11, fontface = "plain"),
    use_raster = nrow(mat_cap) > 50,
    raster_quality = 4,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 7),
      legend_height = unit(3, "cm")
    )
  )
  draw(ht, padding = unit(c(2, 12, 2, 2), "mm"))
  dev.off()
  cat("Saved:", out_file, "(", nrow(mat_cap), "genes )\n")
}

# ── NAS axis (exclude NAS0 — reference group) ────────────────────────────────
nas_sigs <- sigs[staging_axis == "NAS" & is_stage_unique == TRUE & unique_stage != "NAS0_vs_rest"]
cat("\nNAS unique DEGs by stage:\n")
print(nas_sigs[, .N, by = unique_stage][order(unique_stage)])

nas_data <- build_unique_matrix(nas_sigs, ovr_nas, "NAS", "stage_label")
draw_stage_unique_heatmap(
  nas_data$mat, nas_data$row_split,
  file.path(OUT_DIR, "nas_stage_unique_degs_heatmap.pdf"),
  "NAS Stage-Unique DEGs (One-vs-Rest logFC)"
)

# ── Fibrosis axis (exclude F0 — reference group) ─────────────────────────────
fib_sigs <- sigs[staging_axis == "Fibrosis" & is_stage_unique == TRUE & unique_stage != "F0_vs_rest"]
cat("\nFibrosis unique DEGs by stage:\n")
print(fib_sigs[, .N, by = unique_stage][order(unique_stage)])

fib_data <- build_unique_matrix(fib_sigs, ovr_fib, "Fibrosis", "stage_label")
draw_stage_unique_heatmap(
  fib_data$mat, fib_data$row_split,
  file.path(OUT_DIR, "fib_stage_unique_degs_heatmap.pdf"),
  "Fibrosis Stage-Unique DEGs (One-vs-Rest logFC)"
)

cat("\n=== Done ===\n")
