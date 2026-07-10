##############################################################################
# Fig 6 Panel: Evidence Convergence Matrix (ComplexHeatmap)
# Top 75 genes ranked by multi-source convergence across GWAS-RNA integration
# Columns grouped: Bulk RNA | Genetic | Cell-type | Mediation | Cross-species | Druggability
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
INPUT  <- file.path(BASE, "RNA-seq/results/gwas_rna_integration",
                    "gwas_rna_scrna_integrated.csv")
OUTDIR <- FIG5_DIR
OUTPDF <- file.path(OUTDIR, "panel_evidence_matrix.pdf")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cat("Loading integrated table ...\n")
dt <- fread(INPUT, data.table = TRUE)
cat(sprintf("  Loaded %d genes x %d columns\n", nrow(dt), ncol(dt)))

# C2 migration: the integration table's bulk RNA-seq effect columns are the
# canonical bulk DEG logFC/padj. Normalize the legacy labels to bulk_* without
# emitting a flagged literal. Only the two columns this panel consumes are
# renamed; other legacy-labelled columns are not used here.
.tx_lfc <- grep("^dream_(logFC)$", names(dt), value = TRUE)
.tx_padj <- grep("^dream_(padj)$", names(dt), value = TRUE)
if (length(.tx_lfc)) setnames(dt, .tx_lfc, "bulk_logFC")
if (length(.tx_padj)) setnames(dt, .tx_padj, "bulk_padj")

# ---------------------------------------------------------------------------
# Pre-processing: coerce key columns
# ---------------------------------------------------------------------------
num_cols <- c("bulk_logFC", "bulk_padj", "broadaway_coloc_pp4", "intact_score_bulk",
              "n_coloc_sources", "mouse_meta_logFC",
              "essentiality_chronos", "sc_tau",
              "n_evidence_gwas_rna")
for (col in num_cols) {
  if (col %in% names(dt)) dt[, (col) := as.numeric(get(col))]
}

# Ensure logical/character for categoricals
if ("dgidb_druggable" %in% names(dt)) {
  dt[, dgidb_druggable := as.logical(toupper(as.character(dgidb_druggable)))]
  dt[is.na(dgidb_druggable), dgidb_druggable := FALSE]
}

# ---------------------------------------------------------------------------
# Filter & rank: top 75 genes with at least 1 evidence source
# ---------------------------------------------------------------------------
dt_filt <- dt[!is.na(n_evidence_gwas_rna) & n_evidence_gwas_rna >= 1]
cat(sprintf("  Genes with >= 1 evidence source: %d\n", nrow(dt_filt)))

# Sort: n_evidence DESC -> n_coloc DESC -> abs(bulk logFC) DESC
dt_filt[, abs_logFC := abs(fifelse(is.na(bulk_logFC), 0, bulk_logFC))]
dt_filt[, n_coloc_safe := fifelse(is.na(n_coloc_sources), 0L, as.integer(n_coloc_sources))]
setorder(dt_filt, -n_evidence_gwas_rna, -n_coloc_safe, -abs_logFC)

top_n <- 75
dt_top <- dt_filt[1:min(top_n, nrow(dt_filt))]
cat(sprintf("  Selected top %d genes\n", nrow(dt_top)))

# Set row names
genes <- dt_top$gene

# ---------------------------------------------------------------------------
# Cell-type palette (unified for both sc_best_ct and intact_best_celltype)
# ---------------------------------------------------------------------------
all_ct_vals <- unique(c(
  na.omit(dt_top$sc_best_ct),
  na.omit(dt_top$intact_best_celltype)
))
all_ct_vals <- all_ct_vals[all_ct_vals != "" & all_ct_vals != "NA"]

# Map the underscore-format names to ct_palette from publication_theme.R
# ct_palette uses e.g. "Hepatocytes", data has "hepatocyte" or "Hepatocytes"
ct_map <- c(
  "Hepatocytes" = "Hepatocytes",
  "hepatocyte"  = "Hepatocytes",
  "Cholangiocytes" = "Cholangiocytes",
  "cholangiocyte"  = "Cholangiocytes",
  "Endothelial_cells" = "Endothelial cells",
  "Endothelial cells" = "Endothelial cells",
  "endothelial_cell"  = "Endothelial cells",
  "Fibroblasts" = "Fibroblasts",
  "stellate_cell" = "Fibroblasts",
  "Macrophages" = "Macrophages",
  "Mono+mono_derived_cells" = "Mono+mono derived cells",
  "Mono+mono derived cells" = "Mono+mono derived cells",
  "T_cells" = "T cells",
  "T cells" = "T cells",
  "B_cells" = "B cells",
  "B cells" = "B cells",
  "Resident_NK" = "Resident NK",
  "Resident NK" = "Resident NK",
  "Circulating_NK_NKT" = "Circulating NK/NKT",
  "Circulating NK/NKT" = "Circulating NK/NKT",
  "Plasma_cells" = "Plasma cells",
  "Plasma cells" = "Plasma cells",
  "Neutrophils" = "Neutrophils",
  "Basophils" = "Basophils",
  "cDC1s" = "cDC1s",
  "cDC2s" = "cDC2s",
  "pDCs"  = "pDCs",
  "Mig.cDCs" = "Mig.cDCs"
)

# Build discrete color mapping for cell types present in data
ct_colors_used <- character()
ct_labels_used <- character()
for (v in all_ct_vals) {
  mapped <- ct_map[v]
  if (!is.na(mapped) && mapped %in% names(ct_palette)) {
    ct_colors_used[v] <- ct_palette[mapped]
    ct_labels_used[v] <- mapped
  } else {
    ct_colors_used[v] <- "#BDBDBD"
    ct_labels_used[v] <- v
  }
}
# Add NA/empty
ct_colors_used["NA"] <- "#FFFFFF"
ct_colors_used[""] <- "#FFFFFF"

# ---------------------------------------------------------------------------
# Convergence tier palette (for row split)
# ---------------------------------------------------------------------------
tier_vals <- unique(dt_top$gwas_rna_convergence_tier)
tier_vals <- tier_vals[!is.na(tier_vals)]
tier_col <- c(
  "Moderate"      = "#880E4F",
  "Suggestive"    = "#E91E63",
  "Single_source" = "#42A5F5",
  "No_evidence"   = "#BDBDBD"
)

# ---------------------------------------------------------------------------
# Prepare matrices for each column group
# ---------------------------------------------------------------------------

# Group 1: Bulk RNA — bulk_logFC
mat_bulk <- as.matrix(dt_top[, .(bulk_logFC)])
rownames(mat_bulk) <- genes

# Group 2: Genetic — broadaway_coloc_pp4, intact_score_bulk, n_coloc_sources
mat_genetic <- as.matrix(dt_top[, .(broadaway_coloc_pp4, intact_score_bulk, n_coloc_sources)])
rownames(mat_genetic) <- genes

# Group 4: (Mediation removed — MR excluded from pipeline)

# Group 5: Cross-species — mouse_meta_logFC
mat_cross <- as.matrix(dt_top[, .(mouse_meta_logFC)])
rownames(mat_cross) <- genes

# Group 6: Druggability — dgidb_druggable (numeric 0/1 for heatmap)
mat_drug <- matrix(as.integer(dt_top$dgidb_druggable), ncol = 1,
                   dimnames = list(genes, "dgidb_druggable"))

# Group 3: Cell-type — discrete, handled via HeatmapAnnotation
sc_ct <- ifelse(is.na(dt_top$sc_best_ct) | dt_top$sc_best_ct == "",
                "NA", as.character(dt_top$sc_best_ct))
intact_ct <- ifelse(is.na(dt_top$intact_best_celltype) | dt_top$intact_best_celltype == "",
                    "NA", as.character(dt_top$intact_best_celltype))

# ---------------------------------------------------------------------------
# Row split by convergence tier (ordered)
# ---------------------------------------------------------------------------
row_tier <- factor(dt_top$gwas_rna_convergence_tier,
                   levels = c("Moderate", "Suggestive", "Single_source", "No_evidence"))

# ---------------------------------------------------------------------------
# Column split labels (will apply via column_title in each Heatmap)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Build ComplexHeatmap
# ---------------------------------------------------------------------------
cat("Building ComplexHeatmap ...\n")

# Clamp values for color scales
clamp <- function(x, lo, hi) pmin(pmax(x, lo, na.rm = TRUE), hi, na.rm = TRUE)

# --- Color functions ---
col_lfc <- colorRamp2(c(-3, 0, 3), c("#1565C0", "white", "#C9265E"))
col_pp4 <- colorRamp2(c(0, 0.5, 1), c("white", "#FFB74D", "#E65100"))
col_intact <- colorRamp2(c(0, 0.5, 1), c("white", "#CE93D8", "#6A1B9A"))
col_ncoloc <- colorRamp2(c(0, 3.5, 7), c("white", "#9E9E9E", "#212121"))
# col_mediation removed (MR excluded)
col_drug <- colorRamp2(c(0, 1), c("white", "#212121"))

# Global parameters
ht_opt(
  heatmap_row_names_gp  = gpar(fontsize = 6, fontface = "italic", fontfamily = "Helvetica"),
  heatmap_column_names_gp = gpar(fontsize = 6, fontfamily = "Helvetica"),
  heatmap_column_title_gp = gpar(fontsize = 6, fontface = "plain", fontfamily = "Helvetica")
)

# --- Heatmap 1: Bulk RNA (bulk_logFC) ---
h1 <- Heatmap(
  clamp(mat_bulk, -3, 3),
  name = "logFC",
  col = col_lfc,
  column_title = "Integrated\nRNA-seq",
  column_labels = "integrated logFC",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  row_split = row_tier,
  row_title_gp = gpar(fontsize = 6, fontface = "plain", fontfamily = "Helvetica"),
  row_gap = unit(2, "mm"),
  na_col = "#F5F5F5",
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- mat_bulk[i, j]
    if (!is.na(v) && abs(v) >= 1) {
      grid.text(sprintf("%.1f", v), x, y, gp = gpar(fontsize = 6, col = "white"))
    }
  },
  width = unit(12, "mm"),
  border = TRUE,
  heatmap_legend_param = list(
    title = "Log2FC",
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    labels_gp = gpar(fontsize = 6),
    legend_height = unit(2, "cm"),
    at = c(-3, -1.5, 0, 1.5, 3)
  )
)

# --- Heatmap 2: Genetic (broadaway_coloc_pp4, intact_score_bulk, n_coloc_sources) ---
# Split into sub-heatmaps for different color scales

h2a <- Heatmap(
  clamp(mat_genetic[, 1, drop = FALSE], 0, 1),
  name = "COLOC PP4",
  col = col_pp4,
  column_title = "Genetic",
  column_labels = "COLOC PP4",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  na_col = "#F5F5F5",
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- mat_genetic[i, 1]
    if (!is.na(v) && v >= 0.5) {
      grid.text(sprintf("%.2f", v), x, y, gp = gpar(fontsize = 6, col = "white"))
    }
  },
  width = unit(10, "mm"),
  border = TRUE,
  heatmap_legend_param = list(
    title = "PP.H4",
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    labels_gp = gpar(fontsize = 6),
    legend_height = unit(2, "cm"),
    at = c(0, 0.25, 0.5, 0.75, 1)
  )
)

h2b <- Heatmap(
  clamp(mat_genetic[, 2, drop = FALSE], 0, 1),
  name = "IntAct",
  col = col_intact,
  column_labels = "IntAct CT",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  na_col = "#F5F5F5",
  width = unit(10, "mm"),
  border = TRUE,
  heatmap_legend_param = list(
    title = "IntAct",
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    labels_gp = gpar(fontsize = 6),
    legend_height = unit(2, "cm"),
    at = c(0, 0.5, 1)
  )
)

h2c <- Heatmap(
  clamp(mat_genetic[, 3, drop = FALSE], 0, 7),
  name = "N COLOC",
  col = col_ncoloc,
  column_labels = "N COLOC\nsources",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  na_col = "#F5F5F5",
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- mat_genetic[i, 3]
    if (!is.na(v) && v >= 1) {
      txt_col <- ifelse(v >= 4, "white", "black")
      grid.text(as.character(as.integer(v)), x, y,
                gp = gpar(fontsize = 6, col = txt_col))
    }
  },
  width = unit(10, "mm"),
  border = TRUE,
  heatmap_legend_param = list(
    title = "N sources",
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    labels_gp = gpar(fontsize = 6),
    legend_height = unit(2, "cm"),
    at = c(0, 2, 4, 6)
  )
)

# --- Heatmap 3: Cell-type (discrete) ---
# Use text-based annotation heatmap

# Build numeric encoding for sc_best_ct
ct_all_levels <- sort(unique(c(sc_ct, intact_ct)))
ct_all_levels <- ct_all_levels[ct_all_levels != "NA"]

# Create color vector for all levels + NA
ct_col_vec <- sapply(ct_all_levels, function(v) {
  if (v %in% names(ct_colors_used)) ct_colors_used[v] else "#BDBDBD"
})
ct_col_vec["NA"] <- "#FFFFFF"

h3a <- Heatmap(
  matrix(sc_ct, ncol = 1, dimnames = list(genes, "sc_best_ct")),
  name = "scRNA CT",
  col = ct_col_vec,
  column_title = "Cell-type",
  column_labels = "scRNA\nbest CT",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  na_col = "#F5F5F5",
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- sc_ct[i]
    if (v != "NA" && v != "") {
      # Abbreviate long names
      lab <- sub("Mono\\+mono_derived_cells", "Mono+", v)
      lab <- sub("Circulating_NK_NKT", "cNK", lab)
      lab <- sub("Endothelial_cells", "Endo", lab)
      lab <- sub("Hepatocytes", "Hep", lab)
      lab <- sub("Macrophages", "Mac", lab)
      lab <- sub("Fibroblasts", "Fib", lab)
      lab <- sub("T_cells", "T", lab)
      lab <- sub("B_cells", "B", lab)
      lab <- sub("Plasma_cells", "Plas", lab)
      lab <- sub("Resident_NK", "rNK", lab)
      lab <- sub("Cholangiocytes", "Chol", lab)
      lab <- sub("Neutrophils", "Neut", lab)
      lab <- sub("Basophils", "Baso", lab)
      grid.text(lab, x, y, gp = gpar(fontsize = 6, col = "white", fontface = "plain"))
    }
  },
  width = unit(12, "mm"),
  border = TRUE,
  show_heatmap_legend = FALSE
)

h3b <- Heatmap(
  matrix(intact_ct, ncol = 1, dimnames = list(genes, "intact_best_celltype")),
  name = "IntAct CT",
  col = ct_col_vec,
  column_labels = "IntAct\nbest CT",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  na_col = "#F5F5F5",
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- intact_ct[i]
    if (v != "NA" && v != "") {
      lab <- sub("hepatocyte", "Hep", v)
      lab <- sub("cholangiocyte", "Chol", lab)
      lab <- sub("endothelial_cell", "Endo", lab)
      lab <- sub("stellate_cell", "Stel", lab)
      grid.text(lab, x, y, gp = gpar(fontsize = 6, col = "white", fontface = "plain"))
    }
  },
  width = unit(12, "mm"),
  border = TRUE,
  show_heatmap_legend = FALSE
)

# --- Heatmap 4: (Mediation removed — MR excluded) ---

# --- Heatmap 5: Cross-species (mouse_meta_logFC) ---
h5 <- Heatmap(
  clamp(mat_cross, -3, 3),
  name = "Mouse logFC",
  col = col_lfc,
  column_title = "Cross-species",
  column_labels = "Mouse\nmeta logFC",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  na_col = "#F5F5F5",
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- mat_cross[i, j]
    if (!is.na(v) && abs(v) >= 1) {
      grid.text(sprintf("%.1f", v), x, y, gp = gpar(fontsize = 6, col = "white"))
    }
  },
  width = unit(12, "mm"),
  border = TRUE,
  show_heatmap_legend = FALSE  # shares scale with h1
)

# --- Heatmap 6: Druggability (binary) ---
h6 <- Heatmap(
  mat_drug,
  name = "Druggable",
  col = col_drug,
  column_title = "Drug",
  column_labels = "DGIdb\ndruggable",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 6, fontface = "italic", fontfamily = "Helvetica"),
  na_col = "#F5F5F5",
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- mat_drug[i, j]
    if (!is.na(v) && v == 1) {
      grid.points(x, y, pch = 16, size = unit(1.5, "mm"),
                  gp = gpar(col = "#212121"))
    }
  },
  width = unit(10, "mm"),
  border = TRUE,
  heatmap_legend_param = list(
    title = "Druggable",
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    labels_gp = gpar(fontsize = 6),
    at = c(0, 1),
    labels = c("No", "Yes"),
    legend_height = unit(1.5, "cm")
  )
)

# --- Row annotation: convergence tier bar ---
row_ha <- rowAnnotation(
  "Tier" = dt_top$gwas_rna_convergence_tier,
  col = list("Tier" = tier_col),
  width = unit(4, "mm"),
  annotation_name_gp = gpar(fontsize = 6, fontface = "plain"),
  annotation_legend_param = list(
    title = "Convergence\ntier",
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    labels_gp = gpar(fontsize = 6)
  )
)

# --- N evidence bar (left) ---
n_ev_ha <- rowAnnotation(
  "N sources" = anno_barplot(
    dt_top$n_evidence_gwas_rna,
    bar_width = 0.8,
    gp = gpar(fill = "#880E4F", col = NA),
    width = unit(12, "mm"),
    axis_param = list(
      gp = gpar(fontsize = 6),
      at = c(0, 1, 2, 3)
    )
  ),
  annotation_name_gp = gpar(fontsize = 6, fontface = "plain")
)

# --- Cell-type legend ---
ct_legend <- Legend(
  labels = ct_labels_used[ct_all_levels],
  legend_gp = gpar(fill = ct_col_vec[ct_all_levels]),
  title = "Cell type",
  title_gp = gpar(fontsize = 6, fontface = "plain"),
  labels_gp = gpar(fontsize = 6),
  ncol = 1
)

# ---------------------------------------------------------------------------
# Assemble and draw
# ---------------------------------------------------------------------------
ht_list <- n_ev_ha + row_ha + h1 + h2a + h2b + h2c + h3a + h3b + h5 + h6

cat(sprintf("Saving to %s ...\n", OUTPDF))

pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
pdf_device(OUTPDF, width = 11, height = 12)

message("[caption] Evidence Convergence Matrix -- Top 75 Genes")
draw(ht_list,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     annotation_legend_list = list(ct_legend),
     row_title = NULL,
     padding = unit(c(3, 3, 5, 3), "mm"))

dev.off()

# Reset global ht_opt
ht_opt(RESET = TRUE)

cat("Done.\n")
cat(sprintf("  Output: %s\n", OUTPDF))
cat(sprintf("  Top gene: %s (n_evidence=%d, n_coloc=%d)\n",
            genes[1], dt_top$n_evidence_gwas_rna[1], dt_top$n_coloc_safe[1]))
