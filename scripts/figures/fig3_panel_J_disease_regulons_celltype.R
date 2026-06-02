##############################################################################
# Fig 3 Panel J: Disease regulon TFs (rows) x cell types (cols) heatmap
#   Liang Cell 2024 aesthetic refinement (2026-05-18):
#     - Row + column dendrograms
#     - Column-side colored annotation bar (one color per CT)
#     - White-centered diverging scale anchored at #1565C0 / #C9265E
#     - Drug-target row-side annotation bar (gold rectangles, NOT row borders)
#     - 4-way validated TFs (THRB/HNF4A/RORA/MLXIPL/MAX) in bold italic
#     - Significance stars overlaid in bold white/black on each cell
#     - CYP26A1 italic footnote callout below heatmap
#
# Inputs (all from Analysis/ATAC/Human_Multiome/scenic_plus):
#   - disease_regulons.csv               (Hepatocyte; 24 TFs)
#   - macrophage_disease_regulons.csv    (23 TFs)
#   - stellate_disease_regulons.csv      (48 TFs)
#   - endothelial_disease_regulons.csv   (23 TFs)
#   - cholangiocyte_disease_regulons.csv (72 TFs - filter hard)
#
# Filter:
#   1. Per-CT canonical threshold: padj < 0.05 & |activity_diff| > 0.05
#   2. Union of TFs surviving in >=2 CTs OR present in 4-way validated set
#      (4-way validated TFs: THRB, HNF4A, RORA, MLXIPL, MAX)
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
OUT_PDF <- file.path(OUT_DIR, "fig3_panel_J_disease_regulons_celltype.pdf")
OUT_CSV <- file.path(OUT_DIR, "fig3_panel_J_disease_regulons_celltype.csv")

# ---------------------------------------------------------------------------
# Annotations
# ---------------------------------------------------------------------------
DRUG_TARGET_TFS  <- c("THRB", "NR1H4", "PPARA", "PPARG")  # FDA / clinical-stage MASLD
VALIDATED_4WAY   <- c("THRB", "HNF4A", "RORA", "MLXIPL", "MAX")

PADJ_CUT    <- 0.05
DELTA_CUT   <- 0.05

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
ct_files <- list(
  Hepatocyte    = file.path(SP_DIR, "disease_regulons.csv"),
  Macrophage    = file.path(SP_DIR, "macrophage_disease_regulons.csv"),
  Stellate      = file.path(SP_DIR, "stellate_disease_regulons.csv"),
  Endothelial   = file.path(SP_DIR, "endothelial_disease_regulons.csv"),
  Cholangiocyte = file.path(SP_DIR, "cholangiocyte_disease_regulons.csv")
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
# Apply canonical filter & build TF union
# ---------------------------------------------------------------------------
reg_all[, sig_ct := !is.na(activity_padj) & activity_padj < PADJ_CUT &
                    !is.na(regulon_activity_diff) &
                    abs(regulon_activity_diff) > DELTA_CUT]

tf_n_sig <- reg_all[sig_ct == TRUE, .N, by = tf_name]
keep_tfs <- union(tf_n_sig[N >= 2, tf_name], VALIDATED_4WAY)
keep_tfs <- intersect(keep_tfs, unique(reg_all$tf_name))
cat(sprintf("TFs after filter (>=2 CT sig OR 4-way validated): %d\n",
            length(keep_tfs)))
cat("Kept TFs: ", paste(keep_tfs, collapse = ", "), "\n")

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
# matrix values displayed remain unimputed)
# ---------------------------------------------------------------------------
clust_mat <- mat_delta
clust_mat[is.na(clust_mat)] <- 0

row_clust <- hclust(dist(clust_mat,        method = "euclidean"), method = "average")
col_clust <- hclust(dist(t(clust_mat),     method = "euclidean"), method = "average")

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
# Cell function: bold significance stars (white-on-dark, black-on-light)
# ---------------------------------------------------------------------------
star_fun <- function(j, i, x, y, width, height, fill) {
  p  <- mat_padj[i, j]
  d  <- mat_delta[i, j]
  if (is.na(d)) {
    grid.rect(x, y, width, height, gp = gpar(fill = "#F5F5F5", col = NA))
    return(invisible())
  }
  if (!is.na(p)) {
    stars <- if (p < 0.01) "**" else if (p < 0.05) "*" else ""
    if (nzchar(stars)) {
      txt_col <- if (abs(d) > clamp_range * 0.55) "white" else "black"
      grid.text(stars, x, y,
                gp = gpar(fontsize = 7, col = txt_col, fontface = "bold"))
    }
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

# Row label faces: bold italic for 4-way validated, regular italic otherwise
row_faces <- ifelse(rownames(mat_delta) %in% VALIDATED_4WAY,
                    "bold.italic", "italic")

# ---------------------------------------------------------------------------
# Build heatmap
# ---------------------------------------------------------------------------
ht <- Heatmap(
  mat_plot,
  name = "Delta activity",
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
    title = expression(Delta * " activity"),
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
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(2, 2, 6, 2), "mm"),
     newpage = FALSE)
upViewport()

# --- Footnote viewport ---
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.text(
  "* padj < 0.05   ** padj < 0.01    bold italic = 4-way validated (THRB, HNF4A, RORA, MLXIPL, MAX)",
  x = 0.02, y = 0.78, just = c("left", "top"),
  gp = gpar(fontsize = 5.5, col = "grey25", fontface = "plain")
)
grid.text(
  "CYP26A1 — novel peak-linked convergence target shared by THRB + HNF4A regulons (see Panel K)",
  x = 0.02, y = 0.30, just = c("left", "top"),
  gp = gpar(fontsize = 5.5, col = "grey45", fontface = "italic")
)
upViewport()
upViewport()

dev.off()

cat("Done. Output:\n  PDF: ", OUT_PDF, "\n  CSV: ", OUT_CSV, "\n", sep = "")

# Print quick summary of which CTs show strongest signal per kept TF
sig_summary <- reg_keep[tf_name %in% rownames(mat_plot),
  .(n_sig_ct = sum(sig_ct, na.rm = TRUE),
    max_abs_delta = max(abs(regulon_activity_diff), na.rm = TRUE)),
  by = tf_name][order(-n_sig_ct, -max_abs_delta)]
cat("\nTop TFs by sig-CT count:\n")
print(sig_summary, nrows = 30)
