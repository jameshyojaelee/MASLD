##############################################################################
# disease_master_regulators_celltype: Cross-modality disease master-regulator
#   TFs (rows) x cell types (cols) heatmap, colored by regulon-activity EFFECT
#   SIZE.
#
#   This panel shows the DIRECTION and MAGNITUDE of SCENIC+ regulon-activity
#   change (MASLD vs normal) for cross-modality-defined disease master
#   regulators across 5 cell types. It is an effect-size display, NOT a
#   significance gate: single-cohort donor-level differential significance is
#   underpowered (n=18 donors). The previously used FDR-gated SCENIC+
#   disease_regulons.csv files are EMPTY (the honest donor-level n=18 DE test
#   finds 0; the old non-zero set was cell-level pseudoreplication).
#
#   TF panel definition (cross-modality, defensible):
#     scenic_plus/disease_master_regulators.csv -- hepatocyte SCENIC+ regulon
#     TFs that are bulk MASLD DEGs (n=846) OR COLOC hits (HNF4A/RORA/THRB
#     included). `tf_name` defines the kept set.
#
#   Per-cell-type regulon-activity effect sizes (heatmap fill) come from the
#   populated (NOT disease-filtered) per-CT regulon tables:
#     - hepatocyte_regulons.csv     (Hepatocyte)
#     - macrophage_regulons.csv     (Macrophage)
#     - stellate_regulons.csv       (Stellate)
#     - endothelial_regulons.csv    (Endothelial)
#     - cholangiocyte_regulons.csv  (Cholangiocyte)
#   Each carries tf_name, regulon_activity_diff, activity_padj per TF per CT
#   (long format: one row per regulon target gene; collapsed to one row/TF).
#
#   Aesthetic (Liang Cell 2024 refinement, 2026-05-18):
#     - Row + column dendrograms on the effect-size matrix (guarded for n>=2)
#     - Column-side colored annotation bar (one color per CT)
#     - White-centered diverging scale anchored at #1565C0 / #C9265E
#     - Drug-target row-side annotation bar (gold rectangles, NOT row borders)
#     - Hepatocyte headline TFs (THRB/HNF4A/RORA/MLXIPL/MAX) in bold italic
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

OUT_DIR <- file.path(FIG3_DIR, "panels")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PDF <- file.path(OUT_DIR, "disease_master_regulators_celltype.pdf")
OUT_CSV <- file.path(OUT_DIR, "disease_master_regulators_celltype.csv")

# ---------------------------------------------------------------------------
# Annotations
# ---------------------------------------------------------------------------
DRUG_TARGET_TFS  <- c("THRB", "NR1H4", "PPARA", "PPARG")  # FDA / clinical-stage MASLD
# Hepatocyte headline master regulators (bold-italic row labels). Whichever of
# these are present in the master-regulator panel are emphasized; HNF4A/RORA/THRB
# are guaranteed members of the hepatocyte master-regulator set.
HEADLINE_TFS     <- c("THRB", "HNF4A", "RORA", "MLXIPL", "MAX")

# Cross-modality master-regulator TF set (defines the kept rows). Not a
# significance gate; single-cohort donor-level regulon-DE is underpowered (n=18).
MASTER_REG_F <- file.path(ATAC_DIR, "scenic_plus", "disease_master_regulators.csv")

# Liang-style warm/cool CT palette (warm = positive findings / hepatic-myeloid,
# cool = baseline / vascular-biliary).
CT_COLORS <- c(
  "Hepatocyte"    = "#C9265E",  # warm magenta
  "Macrophage"    = "#F28E2B",  # warm orange
  "Stellate"      = "#A01753",  # warm deep red
  "Endothelial"   = "#2E8B8B",  # cool teal
  "Cholangiocyte" = "#1565C0"   # cool blue
)

# ---------------------------------------------------------------------------
# Load per-CT regulons
# ---------------------------------------------------------------------------
SP_DIR <- file.path(ATAC_DIR, "scenic_plus")
# Populated (NOT disease-filtered) per-CT regulon tables. The old
# {ct}_disease_regulons.csv files are empty (FDR-gated donor-level n=18 DE).
ct_files <- list(
  Hepatocyte    = file.path(SP_DIR, "hepatocyte_regulons.csv"),
  Macrophage    = file.path(SP_DIR, "macrophage_regulons.csv"),
  Stellate      = file.path(SP_DIR, "stellate_regulons.csv"),
  Endothelial   = file.path(SP_DIR, "endothelial_regulons.csv"),
  Cholangiocyte = file.path(SP_DIR, "cholangiocyte_regulons.csv")
)

load_ct <- function(ct, f) {
  if (!file.exists(f)) {
    warning("Missing: ", f); return(NULL)
  }
  dt <- fread(f)
  needed <- c("tf_name", "regulon_activity_diff", "activity_padj")
  if (!all(needed %in% names(dt))) {
    warning("Missing columns in ", f); return(NULL)
  }
  dt[, cell_type := ct]
  dt[, .(cell_type, tf_name, regulon_activity_diff, activity_padj)]
}

reg_list <- mapply(load_ct, names(ct_files), ct_files,
                   SIMPLIFY = FALSE, USE.NAMES = TRUE)
reg_all <- rbindlist(reg_list, use.names = TRUE, fill = TRUE)
cat("Loaded regulon rows per CT:\n"); print(table(reg_all$cell_type))

# Collapse duplicate TF rows within a CT (keep most significant)
reg_all <- reg_all[order(cell_type, tf_name, activity_padj)]
reg_all <- reg_all[, .SD[1], by = .(cell_type, tf_name)]

# ---------------------------------------------------------------------------
# Restrict to the cross-modality disease master-regulator TF set
#   (defensible "disease-associated" definition; NOT a per-CT significance gate)
# ---------------------------------------------------------------------------
if (!file.exists(MASTER_REG_F))
  stop("Master-regulator file not found: ", MASTER_REG_F)
master_tfs <- unique(fread(MASTER_REG_F)$tf_name)
keep_tfs   <- intersect(master_tfs, unique(reg_all$tf_name))
cat(sprintf("Cross-modality master-regulator TFs: %d (present in per-CT regulons: %d)\n",
            length(master_tfs), length(keep_tfs)))
cat("Kept TFs: ", paste(sort(keep_tfs), collapse = ", "), "\n")
hl_present <- intersect(HEADLINE_TFS, keep_tfs)
cat("Headline TFs present: ", paste(hl_present, collapse = ", "), "\n")

reg_keep <- reg_all[tf_name %in% keep_tfs]

# ---------------------------------------------------------------------------
# Build matrices (TF x CT)
# ---------------------------------------------------------------------------
ct_all <- c("Hepatocyte", "Macrophage", "Stellate", "Endothelial", "Cholangiocyte")
ct_present <- intersect(ct_all, unique(reg_keep$cell_type))

mat_delta <- dcast(reg_keep, tf_name ~ cell_type,
                   value.var = "regulon_activity_diff")
mat_padj  <- dcast(reg_keep, tf_name ~ cell_type,
                   value.var = "activity_padj")

tf_vec <- mat_delta$tf_name
mat_delta <- as.matrix(mat_delta[, ..ct_present]); rownames(mat_delta) <- tf_vec
mat_padj  <- as.matrix(mat_padj[,  ..ct_present]);  rownames(mat_padj)  <- tf_vec

# Save the underlying data (delta and padj side-by-side)
delta_df <- data.table(tf_name = rownames(mat_delta))
for (cc in ct_present) delta_df[[paste0("delta_", cc)]] <- mat_delta[, cc]
for (cc in ct_present) delta_df[[paste0("padj_",  cc)]] <- mat_padj[,  cc]
fwrite(delta_df, OUT_CSV)

# ---------------------------------------------------------------------------
# Clustering distances (handle NA by 0-imputation for distance only;
# matrix values displayed remain unimputed). Clustering on the effect-size
# matrix. Guard hclust for n>=2 rows/cols (single-row matrices break dist()).
# ---------------------------------------------------------------------------
clust_mat <- mat_delta
clust_mat[is.na(clust_mat)] <- 0

row_clust <- if (nrow(clust_mat) >= 2)
  hclust(dist(clust_mat, method = "euclidean"), method = "average") else FALSE
col_clust <- if (ncol(clust_mat) >= 2)
  hclust(dist(t(clust_mat), method = "euclidean"), method = "average") else FALSE

# ---------------------------------------------------------------------------
# Color scale: white-centered diverging with Liang endpoints
# ---------------------------------------------------------------------------
clamp_range <- 0.20
mat_plot <- pmin(pmax(mat_delta, -clamp_range), clamp_range)

col_delta <- colorRamp2(
  c(-clamp_range, 0, clamp_range),
  c("#1565C0", "white", "#C9265E")
)

# ---------------------------------------------------------------------------
# Cell function: gray out absent (NA) TF x CT combinations only.
# No significance stars: this is an EFFECT-SIZE display. Single-cohort
# donor-level regulon-activity DE is underpowered (n=18), so per-CT activity_padj
# is uninformative and intentionally NOT overlaid as significance markers.
# ---------------------------------------------------------------------------
star_fun <- function(j, i, x, y, width, height, fill) {
  d <- mat_delta[i, j]
  if (is.na(d)) {
    grid.rect(x, y, width, height, gp = gpar(fill = "#F5F5F5", col = NA))
    return(invisible())
  }
}

# ---------------------------------------------------------------------------
# Annotations: CT identity bar (top) + Drug-target bar (left)
# ---------------------------------------------------------------------------
col_anno <- HeatmapAnnotation(
  `Cell type` = ct_present,
  col = list(`Cell type` = CT_COLORS[ct_present]),
  show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = 5, fontface = "bold"),
  annotation_name_side = "right",
  simple_anno_size = unit(2.8, "mm"),
  border = FALSE
)

is_drug <- rownames(mat_delta) %in% DRUG_TARGET_TFS
drug_anno <- rowAnnotation(
  `Drug target` = anno_simple(
    ifelse(is_drug, "yes", "no"),
    col = c("yes" = "#E6A817", "no" = "#FFFFFF00"),  # transparent white for "no"
    border = FALSE,
    width = unit(2.8, "mm")
  ),
  annotation_name_gp = gpar(fontsize = 5, fontface = "bold"),
  annotation_name_side = "top",
  show_legend = FALSE
)

# Row label faces: bold italic for hepatocyte headline master regulators
# (THRB/HNF4A/RORA/MLXIPL/MAX), regular italic otherwise.
row_faces <- ifelse(rownames(mat_delta) %in% HEADLINE_TFS,
                    "bold.italic", "italic")

# ---------------------------------------------------------------------------
# Build heatmap
# ---------------------------------------------------------------------------
ht <- Heatmap(
  mat_plot,
  name = "Regulon activity effect size",
  col = col_delta,
  na_col = "#F5F5F5",
  cluster_rows = row_clust,
  cluster_columns = col_clust,
  row_dend_width  = unit(8, "mm"),
  column_dend_height = unit(8, "mm"),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7, fontface = row_faces,
                      fontfamily = "Helvetica"),
  column_names_gp = gpar(fontsize = 8, fontface = "bold",
                         fontfamily = "Helvetica"),
  column_names_rot = 45,
  top_annotation = col_anno,
  left_annotation = drug_anno,
  cell_fun = star_fun,
  rect_gp = gpar(col = NA),                    # NO cell borders
  width  = unit(3.2, "cm"),
  height = unit(0.42 * nrow(mat_plot), "cm"),
  border = FALSE,
  heatmap_legend_param = list(
    title = expression(Delta * " activity\n(MASLD - normal)"),
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5),
    legend_height = unit(2.2, "cm"),
    at = c(-clamp_range, -clamp_range/2, 0, clamp_range/2, clamp_range),
    labels = c(sprintf("<= -%.2f", clamp_range), sprintf("-%.2f", clamp_range/2),
               "0", sprintf("+%.2f", clamp_range/2),
               sprintf(">= +%.2f", clamp_range))
  )
)

# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------
n_rows <- nrow(mat_plot)
pdf_h  <- max(4.6, 0.20 * n_rows + 3.0)
pdf_w  <- 5.0

cat(sprintf("Writing PDF (%d rows, %.2f x %.2f in): %s\n",
            n_rows, pdf_w, pdf_h, OUT_PDF))

pdf(OUT_PDF, width = pdf_w, height = pdf_h, useDingbats = FALSE)

# Push a 2-row viewport: heatmap (top) + footnote (bottom)
grid.newpage()
pushViewport(viewport(layout = grid.layout(
  nrow = 2, ncol = 1,
  heights = unit.c(unit(1, "null"), unit(10, "mm"))
)))

# --- Heatmap viewport ---
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht,
     column_title = "Regulon-activity effect size of cross-modality disease master regulators",
     column_title_gp = gpar(fontsize = 7.5, fontface = "bold"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(2, 2, 6, 2), "mm"),
     newpage = FALSE)
upViewport()

# --- Footnote viewport ---
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.text(
  paste0("Fill = SCENIC+ regulon-activity effect size (MASLD - normal); cross-modality master ",
         "regulators (hepatocyte regulon TFs that are bulk DEGs or COLOC hits)."),
  x = 0.02, y = 0.82, just = c("left", "top"),
  gp = gpar(fontsize = 5.5, col = "grey25", fontface = "plain")
)
grid.text(
  paste0("Effect-size display only; single-cohort donor-level regulon DE is underpowered (n=18). ",
         "Bold italic = hepatocyte headline TFs (THRB, HNF4A, RORA, MLXIPL, MAX). Gray = TF absent in cell type."),
  x = 0.02, y = 0.40, just = c("left", "top"),
  gp = gpar(fontsize = 5.5, col = "grey45", fontface = "italic")
)
upViewport()
upViewport()

dev.off()

cat("Done. Output:\n  PDF: ", OUT_PDF, "\n  CSV: ", OUT_CSV, "\n", sep = "")

# Print quick summary of which CTs show the strongest regulon-activity effect per TF
sig_summary <- reg_keep[tf_name %in% rownames(mat_plot),
  .(n_ct_present  = .N,
    max_abs_delta = max(abs(regulon_activity_diff), na.rm = TRUE)),
  by = tf_name][order(-max_abs_delta, -n_ct_present)]
cat("\nTop TFs by max |regulon-activity effect size|:\n")
print(sig_summary, nrows = 30)
