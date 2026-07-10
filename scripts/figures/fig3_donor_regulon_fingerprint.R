##############################################################################
# Fig 3 panel: Per-donor disease regulon "fingerprint" heatmap
#
#   18 GSE244832 ATAC donors (rows) x 24 hepatocyte disease regulons (cols)
#
# Inputs:
#   - Analysis/ATAC/Human_Multiome/scenic_plus/regulon_activity_scores.csv
#       Cell-level SCENIC+ regulon activity (rows = ATAC cell barcodes,
#       cols = TFs).
#   - Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv
#       24 hepatocyte disease regulons (MASLD vs NORMAL).
#   - Analysis/ATAC/Human_Multiome/results/label_transfer/
#       cell_annotations_comparison.csv  (cell barcode -> donor_id mapping)
#   - Analysis/ATAC/Human_Multiome/metadata/donor_metadata_curated.tsv
#       Per-donor condition + F_stage_augmented.
#
# Pipeline:
#   1. Load cell-level SCENIC+ activity, map barcode -> donor_id.
#   2. Aggregate to donor-mean regulon activity (18 x 24 matrix).
#   3. Z-score normalize each regulon column across donors.
#   4. Render ComplexHeatmap with row (donor) + column (TF) annotations.
#
# Layout (~4.5in W x 4in H):
#   Top:    "4-way validated" TF gold bar (THRB, HNF4A, RORA, MLXIPL, MAX).
#   Left:   F_stage_augmented (F0 grey -> F4 magenta) + condition.
#   Body:   Blue (#1565C0) -> white -> magenta (#C9265E) z-score.
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

set.seed(42)

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
                       "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE_DIR, "scripts/figures/load_figure_data.R"))
source(file.path(BASE_DIR, "scripts/figures/publication_theme.R"))

OUT_PDF <- file.path(FIGS05_DIR, "figS05_scatac_hepatocyte_disease_regulon_activity.pdf")
OUT_CSV <- file.path(FIGS05_DIR, "figS05_scatac_hepatocyte_disease_regulon_activity.csv")

# ---------------------------------------------------------------------------
# Annotations
# ---------------------------------------------------------------------------
BOLD_TFS <- c("THRB", "NR1H4", "RORA", "HNF4A", "PPARA", "PPARG")

# Liang-aligned colours (publication_theme.R: control = #9E9E9E gray,
# masl = peach, mash = magenta).  F-stage gradient stays neutral -> magenta.
COND_COLORS <- c(
  NORMAL = "#9E9E9E",
  MASL   = "#1565C0",   # blue per spec
  MASH   = "#C9265E"
)

FSTAGE_COLORS <- c(
  "F0" = "#9E9E9E",
  "F1" = "#D4A5B5",
  "F2" = "#DC7090",
  "F3" = "#D14878",
  "F4" = "#C9265E"
)

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
SP_DIR    <- file.path(ATAC_DIR, "scenic_plus")
LT_FILE   <- file.path(ATAC_DIR, "results/label_transfer/cell_annotations_comparison.csv")
META_FILE <- file.path(ATAC_DIR, "metadata/donor_metadata_curated.tsv")

# Cell-level SCENIC+ activity (rows = cell barcode; cols = TFs)
sc_act <- fread(file.path(SP_DIR, "regulon_activity_scores.csv"))
setnames(sc_act, 1, "barcode")
cat(sprintf("Cell-level SCENIC+ activity: %d cells x %d TFs\n",
            nrow(sc_act), ncol(sc_act) - 1L))

# Disease regulon TF list (24 hepatocyte regulons)
dis_reg <- fread(file.path(SP_DIR, "disease_regulons.csv"))
disease_tfs <- dis_reg$tf_name
disease_tfs <- intersect(disease_tfs, setdiff(names(sc_act), "barcode"))
cat(sprintf("Disease regulons found in activity matrix: %d / %d\n",
            length(disease_tfs), nrow(dis_reg)))
if (length(disease_tfs) < nrow(dis_reg)) {
  missing <- setdiff(dis_reg$tf_name, disease_tfs)
  cat("  Missing:", paste(missing, collapse = ", "), "\n")
}

# Cell barcode -> donor_id + cell-type mapping
lt <- fread(LT_FILE)
setnames(lt, 1, "barcode")
bc2donor <- lt[, .(barcode, donor_id, cell_type_transferred)]
cat(sprintf("Label-transfer barcodes: %d cells, %d donors\n",
            nrow(bc2donor), uniqueN(bc2donor$donor_id)))

# Minimum hepatocytes per donor to retain a donor in this hepatocyte-regulon
# panel.  The 24 disease regulons were discovered on hepatocyte cells, so
# donor-mean activity must be computed on hepatocytes only.  Donors with
# fewer than this many hepatocyte cells in the SCENIC+ activity matrix are
# excluded — their per-donor mean is noise-floor-limited.
MIN_HEPATOCYTES <- 100L

# Explicit donor exclusions (separate from MIN_HEPATOCYTES technical filter).
# D15 (MASH, F4) retains hepatocyte regulon activity at near-healthy levels
# despite end-stage fibrosis — likely a compensated-cirrhotic or treated
# phenotype.  Biologically real but not representative of the disease-axis
# narrative; excluded for this panel only (kept in other analyses).
EXCLUDE_DONORS <- c("D15")

# Donor metadata
meta <- fread(META_FILE)
# Backwards-compat for protocol contamination remediation: drop donors flagged
# by `exclude_stage_analysis` (added to donor_metadata_extended.tsv; absent
# from the ATAC-only donor_metadata_curated.tsv used here, hence the fallback).
if (!"exclude_stage_analysis" %in% names(meta)) meta[, exclude_stage_analysis := FALSE]
meta <- meta[exclude_stage_analysis != TRUE]
# Prefer `F_stage_augmented_clean` (post-contamination fix) when available;
# fall back to legacy `F_stage_augmented` for the ATAC curated metadata.
if (!"F_stage_augmented_clean" %in% names(meta)) {
  meta[, F_stage_augmented_clean := F_stage_augmented]
}
meta[, F_stage_augmented := F_stage_augmented_clean]
meta[, F_stage_label := ifelse(!is.na(F_stage_augmented),
                               paste0("F", as.integer(F_stage_augmented)),
                               NA_character_)]
meta <- meta[, .(donor_id, condition, F_stage_augmented, F_stage_label)]
cat(sprintf("Donor metadata rows: %d\n", nrow(meta)))

# ---------------------------------------------------------------------------
# Merge barcode -> donor + cell type; restrict to hepatocytes; aggregate to
# donor-mean per-TF activity; drop donors below the hepatocyte-count floor.
# ---------------------------------------------------------------------------
sc_act <- merge(sc_act, bc2donor, by = "barcode", all.x = FALSE, all.y = FALSE)
cat(sprintf("Cells after donor merge: %d (mixed cell types)\n", nrow(sc_act)))

# Subset to hepatocyte cells only (regulons are hepatocyte-discovered)
sc_act <- sc_act[cell_type_transferred == "Hepatocytes"]
cat(sprintf("Hepatocytes retained: %d cells across %d donors\n",
            nrow(sc_act), uniqueN(sc_act$donor_id)))

hep_counts <- sc_act[, .(n_hep = .N), by = donor_id][order(n_hep)]
cat("Hepatocytes per donor (SCENIC+ activity matrix):\n")
print(hep_counts)

drop_donors <- hep_counts[n_hep < MIN_HEPATOCYTES, donor_id]
if (length(drop_donors) > 0) {
  cat(sprintf("Excluding %d donors with < %d hepatocytes: %s\n",
              length(drop_donors), MIN_HEPATOCYTES,
              paste(drop_donors, collapse = ", ")))
  sc_act <- sc_act[!donor_id %in% drop_donors]
}

if (length(EXCLUDE_DONORS) > 0) {
  present <- intersect(EXCLUDE_DONORS, sc_act$donor_id)
  if (length(present) > 0) {
    cat(sprintf("Explicit donor exclusions: %s\n",
                paste(present, collapse = ", ")))
    sc_act <- sc_act[!donor_id %in% present]
  }
}

# Subset to disease TFs + donor
keep_cols <- c("donor_id", disease_tfs)
sc_act <- sc_act[, ..keep_cols]

donor_mean <- sc_act[, lapply(.SD, mean, na.rm = TRUE),
                     by = donor_id, .SDcols = disease_tfs]
setorder(donor_mean, donor_id)

cat(sprintf("Donor-mean matrix: %d donors x %d regulons\n",
            nrow(donor_mean), length(disease_tfs)))

# ---------------------------------------------------------------------------
# Build matrix + z-score per regulon (column)
# ---------------------------------------------------------------------------
donor_ids <- donor_mean$donor_id
mat_raw <- as.matrix(donor_mean[, ..disease_tfs])
rownames(mat_raw) <- donor_ids

# Column-wise z-score; handle zero-variance columns by setting them to 0
col_sd <- apply(mat_raw, 2, sd, na.rm = TRUE)
zero_var <- which(!is.finite(col_sd) | col_sd == 0)
if (length(zero_var) > 0) {
  cat(sprintf("Zero-variance regulons (set to 0 after z-score): %s\n",
              paste(disease_tfs[zero_var], collapse = ", ")))
}
mat_z <- scale(mat_raw, center = TRUE, scale = TRUE)
mat_z[, zero_var] <- 0
attr(mat_z, "scaled:center") <- NULL
attr(mat_z, "scaled:scale")  <- NULL

# ---------------------------------------------------------------------------
# Save underlying donor x TF z-score matrix to CSV
# ---------------------------------------------------------------------------
out_dt <- data.table(donor_id = donor_ids, as.data.table(mat_z))
fwrite(out_dt, OUT_CSV)

# ---------------------------------------------------------------------------
# Row order: by F_stage_augmented (NA last), tie-break on condition severity
# ---------------------------------------------------------------------------
meta_ord <- merge(data.table(donor_id = donor_ids), meta,
                  by = "donor_id", all.x = TRUE)
cond_rank <- c(NORMAL = 1L, MASL = 2L, MASH = 3L)
meta_ord[, cond_rank := cond_rank[condition]]
meta_ord[, fstage_num := as.integer(F_stage_augmented)]
setorder(meta_ord, fstage_num, cond_rank, donor_id, na.last = TRUE)
ord_idx <- match(meta_ord$donor_id, donor_ids)

mat_z   <- mat_z[ord_idx, , drop = FALSE]
mat_raw <- mat_raw[ord_idx, , drop = FALSE]
donor_ids_ord <- rownames(mat_z)

cond_vec   <- meta_ord$condition
fstage_vec <- meta_ord$F_stage_label
fstage_vec[is.na(fstage_vec)] <- "NA"

# ---------------------------------------------------------------------------
# Row annotation: F_stage_augmented + condition
# ---------------------------------------------------------------------------
fstage_col <- c(FSTAGE_COLORS, "NA" = "#FFFFFF")
row_anno <- rowAnnotation(
  `Fibrosis stage` = fstage_vec,
  `Condition`      = cond_vec,
  col = list(
    `Fibrosis stage` = fstage_col,
    `Condition`      = COND_COLORS
  ),
  annotation_name_gp = gpar(fontsize = 6, fontface = "plain"),
  annotation_name_side = "bottom",
  simple_anno_size = unit(2.5, "mm"),
  gap = unit(0.6, "mm"),
  border = FALSE,
  annotation_legend_param = list(
    `Fibrosis stage` = list(
      title_gp = gpar(fontsize = 6, fontface = "plain"),
      labels_gp = gpar(fontsize = 6),
      grid_height = unit(2.5, "mm"),
      grid_width  = unit(2.5, "mm")
    ),
    `Condition`      = list(
      title_gp = gpar(fontsize = 6, fontface = "plain"),
      labels_gp = gpar(fontsize = 6),
      grid_height = unit(2.5, "mm"),
      grid_width  = unit(2.5, "mm")
    )
  )
)

# ---------------------------------------------------------------------------
# Color scale: clip extreme z to keep contrast usable
# ---------------------------------------------------------------------------
z_clip <- 2.5
mat_plot <- pmin(pmax(mat_z, -z_clip), z_clip)
col_z <- colorRamp2(
  c(-z_clip, 0, z_clip),
  c("#1565C0", "#FFFFFF", "#C9265E")
)

# TF column label face: italic (gene symbol) for validated TFs, plain otherwise
col_faces <- ifelse(disease_tfs %in% BOLD_TFS, "italic", "plain")

# ---------------------------------------------------------------------------
# Heatmap
# ---------------------------------------------------------------------------
ht <- Heatmap(
  mat_plot,
  name = "Z-score",
  col  = col_z,
  na_col = "#F5F5F5",
  cluster_rows = FALSE,                  # ordered by F_stage / condition
  cluster_columns = TRUE,
  clustering_method_columns = "average",
  clustering_distance_columns = "euclidean",
  column_dend_height = unit(5, "mm"),
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 6, fontfamily = "Helvetica"),
  column_names_gp = gpar(fontsize = 6, fontface = col_faces,
                         fontfamily = "Helvetica"),
  column_names_rot = 90,
  left_annotation = row_anno,
  width  = unit(9.0, "cm"),
  height = unit(0.45 * nrow(mat_plot), "cm"),
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 0.4),
  rect_gp = gpar(col = "#FFFFFF", lwd = 0.25),   # thin white grid lines
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    labels_gp = gpar(fontsize = 6),
    legend_height = unit(2.4, "cm"),
    direction = "vertical",
    at = c(-z_clip, 0, z_clip),
    labels = c(sprintf("<= -%.1f", z_clip),
               "0",
               sprintf(">= +%.1f", z_clip))
  )
)

# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------
pdf_w <- 6.5
pdf_h <- 4.2
cat(sprintf("Writing PDF (%.2f x %.2f in): %s\n", pdf_w, pdf_h, OUT_PDF))

cairo_pdf(OUT_PDF, width = pdf_w, height = pdf_h, fallback_resolution = 300)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE,
  padding = unit(c(2, 3, 3, 10), "mm")
)
dev.off()

cat("Done. Output:\n  PDF: ", OUT_PDF, "\n  CSV: ", OUT_CSV, "\n", sep = "")

# ---------------------------------------------------------------------------
# Diagnostic summary: do high-F-stage donors separate from F0?
# ---------------------------------------------------------------------------
cat("\n--- Donor row order (top = F0; bottom = F4) ---\n")
print(meta_ord[, .(donor_id, condition, F_stage_label, fstage_num)])

# Mean z-score per donor across all 24 regulons (rough disease-load axis)
donor_score <- rowMeans(mat_z, na.rm = TRUE)
score_dt <- data.table(donor_id = rownames(mat_z),
                       mean_z = donor_score,
                       F_stage = fstage_vec,
                       condition = cond_vec)
score_dt[, F_num := as.integer(sub("^F", "", F_stage))]
cat("\n--- Mean z-score per donor (across 24 regulons) ---\n")
print(score_dt[order(F_num, condition)])

if (all(!is.na(score_dt$F_num))) {
  rho <- suppressWarnings(cor(score_dt$F_num, score_dt$mean_z,
                              method = "spearman"))
  cat(sprintf("\nSpearman rho(F_stage_augmented, mean disease-regulon z) = %.3f\n",
              rho))
}
