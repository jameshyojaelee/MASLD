# Fig 3 panel — chromVAR TF × (cell type × condition) cascade heatmap
# 18 curated TFs (rows) × 4 MASLD-core cell types × 3 conditions
# (NORMAL / MASL / MASH) = 12 columns.  Reveals the cellular regulatory
# cascade: which TF programs activate / deactivate in which cell type
# at which stage of progression.

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
PER_DONOR_TSV <- file.path(ATAC_DIR, "results/chromvar_v2/chromvar_per_donor.tsv.gz")
META_FILE     <- file.path(ATAC_DIR, "metadata/donor_metadata_curated.tsv")
OUT_PDF       <- file.path(FIGS05_DIR, "figS05_scatac_chromvar_celltype_progression.pdf")
OUT_CSV       <- sub("\\.pdf$", ".csv", OUT_PDF)

# -----------------------------------------------------------------------------
# Curated TF list — matches the cell-type and progression dotplot panels
# -----------------------------------------------------------------------------
TF_GROUPS <- list(
  `Hepatocyte-related` =
    c(HNF4A  = "HNF4A",  HNF1A = "HNF1A",
      FOXA1  = "FOXA1",  FOXA2 = "FOXA2"),
  `NR drug targets` =
    c(THRB   = "THRB",   NR1H4 = "Nr1H4",
      PPARA  = "Ppara",  PPARG = "PPARG",
      RORA   = "RORA",   RXRA  = "Rxra"),
  `Metabolic` =
    c(MLXIPL = "MLXIPL", `NR1H3 (LXRα)` = "Nr1h3", KLF15 = "KLF15"),
  `Stress / inflammation` =
    c(XBP1   = "XBP1",   ETS1  = "ETS1"),
  `Fibrogenic / EMT` =
    c(SNAI1  = "SNAI1",  ZEB1  = "ZEB1",  TCF4  = "TCF4")
)

BOLD_TFS     <- c("THRB", "NR1H4", "RORA", "HNF4A", "PPARA", "PPARG")
DROP_LOW_HEP <- c("D05", "D07", "D18")
EXCLUDE_OUTL <- c("D15")
CT_KEEP      <- c("Hepatocytes", "Cholangiocytes", "Fibroblasts", "Macrophages")
COND_LEVELS  <- c("NORMAL", "MASL", "MASH")

CT_COLORS <- c(
  Hepatocytes    = "#E07A5F",   # warm orange
  Cholangiocytes = "#81B29A",   # green
  Fibroblasts    = "#3D5A80",   # deep blue
  Macrophages    = "#C9A227"    # ochre
)

# Condition colors — kept distinct from the blue↔magenta deviation scale
COND_COLORS <- c(
  NORMAL = "#B8B8B8",
  MASL   = "#F4A261",
  MASH   = "#8C5E58"
)

# -----------------------------------------------------------------------------
# Load + filter
# -----------------------------------------------------------------------------
tf_lookup <- rbindlist(lapply(names(TF_GROUPS), function(g) {
  v <- TF_GROUPS[[g]]
  data.table(display = names(v), motif = unname(v), group = g)
}))

pd <- fread(cmd = sprintf("zcat %s", PER_DONOR_TSV))
setnames(pd, "donor_id_atac", "donor_id")
pd <- pd[cell_type %in% CT_KEEP]

meta <- fread(META_FILE)
# Backwards-compat for protocol contamination remediation: drop donors flagged
# by `exclude_stage_analysis` (added to donor_metadata_extended.tsv; absent
# from the ATAC-only donor_metadata_curated.tsv used here, hence the fallback).
if (!"exclude_stage_analysis" %in% names(meta)) meta[, exclude_stage_analysis := FALSE]
meta <- meta[exclude_stage_analysis != TRUE]
pd <- merge(pd, meta[, .(donor_id, condition)], by = "donor_id", all.x = FALSE)

drop_set <- c(DROP_LOW_HEP, EXCLUDE_OUTL)
pd <- pd[!donor_id %in% drop_set]
cat(sprintf("Donors retained: %d (%s excluded)\n",
            uniqueN(pd$donor_id), paste(drop_set, collapse = ",")))

# Restrict to curated TFs + collapse motif variants by TF_symbol
pd <- pd[TF %in% tf_lookup$motif | TF_symbol %in% tf_lookup$motif]
pd <- pd[, .(mean_deviation = mean(mean_deviation, na.rm = TRUE)),
         by = .(donor_id, cell_type, condition, TF_symbol)]
pd <- merge(pd, tf_lookup, by.x = "TF_symbol", by.y = "motif", all.x = TRUE)
pd <- pd[!is.na(display)]
cat(sprintf("TFs retained: %d\n", uniqueN(pd$display)))

# -----------------------------------------------------------------------------
# Aggregate to cell_type × condition × TF mean deviation, plus N donors
# -----------------------------------------------------------------------------
agg <- pd[, .(mean_dev = mean(mean_deviation, na.rm = TRUE),
              n_donors = uniqueN(donor_id)),
          by = .(display, group, cell_type, condition)]
agg[, col_key := paste(cell_type, condition, sep = "|")]

# Row + column ordering
y_order <- unlist(lapply(TF_GROUPS, names), use.names = FALSE)
y_order <- intersect(y_order, agg$display)
col_keys <- as.vector(outer(CT_KEEP, COND_LEVELS,
                            function(ct, cond) paste(ct, cond, sep = "|")))

mat <- dcast(agg, display ~ col_key, value.var = "mean_dev")
mat <- mat[match(y_order, display)]
mat_m <- as.matrix(mat[, ..col_keys])
rownames(mat_m) <- mat$display

# Clip values for display
z_clip <- 0.6
mat_plot <- pmin(pmax(mat_m, -z_clip), z_clip)

# -----------------------------------------------------------------------------
# Annotations
# -----------------------------------------------------------------------------
col_ct   <- sub("\\|.*", "", colnames(mat_plot))
col_cond <- sub(".*\\|", "", colnames(mat_plot))

col_anno <- HeatmapAnnotation(
  `Cell type` = col_ct,
  `Condition` = factor(col_cond, levels = COND_LEVELS),
  col = list(
    `Cell type` = CT_COLORS[CT_KEEP],
    `Condition` = COND_COLORS
  ),
  annotation_name_gp   = gpar(fontsize = 6, fontface = "bold"),
  annotation_name_side = "left",
  simple_anno_size     = unit(2.5, "mm"),
  gap                  = unit(0.6, "mm"),
  border               = FALSE,
  show_legend          = TRUE,
  annotation_legend_param = list(
    `Cell type` = list(title_gp = gpar(fontsize = 6, fontface = "bold"),
                       labels_gp = gpar(fontsize = 5),
                       grid_height = unit(2.5, "mm"),
                       grid_width  = unit(2.5, "mm")),
    `Condition` = list(title_gp = gpar(fontsize = 6, fontface = "bold"),
                       labels_gp = gpar(fontsize = 5),
                       grid_height = unit(2.5, "mm"),
                       grid_width  = unit(2.5, "mm"))
  )
)

# Column-gap pattern: 0 within a cell-type block, big gap between blocks
n_per_block <- length(COND_LEVELS)
col_split <- factor(col_ct, levels = CT_KEEP)

# Y-axis label aesthetics: bold for user-specified TFs
y_levels <- rownames(mat_plot)
y_stem   <- sub(" \\(.*", "", y_levels)
y_face   <- ifelse(y_stem %in% BOLD_TFS, "bold", "plain")

# Row group separators (functional TF categories)
row_group <- tf_lookup$group[match(y_levels, tf_lookup$display)]
row_split <- factor(row_group,
                    levels = unique(tf_lookup$group))

# -----------------------------------------------------------------------------
# Heatmap
# -----------------------------------------------------------------------------
col_z <- colorRamp2(c(-z_clip, 0, z_clip),
                    c("#1565C0", "#FFFFFF", "#C9265E"))

ht <- Heatmap(
  mat_plot,
  name = "Motif deviation z",
  col  = col_z,
  na_col = "#F5F5F5",
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  column_split = col_split,
  row_split    = row_split,
  row_gap      = unit(1.2, "mm"),
  column_gap   = unit(2.0, "mm"),
  row_title_gp = gpar(fontsize = 6, fontface = "bold"),
  row_title_side = "left",
  row_title_rot = 0,
  column_title  = NULL,
  row_names_side = "left",
  row_names_gp   = gpar(fontsize = 6, fontface = y_face),
  show_column_names = FALSE,
  top_annotation = col_anno,
  width  = unit(7.0, "cm"),
  height = unit(0.42 * nrow(mat_plot), "cm"),
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 0.4),
  rect_gp   = gpar(col = "#FFFFFF", lwd = 0.25),
  heatmap_legend_param = list(
    title = "Motif deviation z",
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5),
    legend_height = unit(2.0, "cm"),
    direction = "vertical",
    at = c(-z_clip, 0, z_clip),
    labels = c(sprintf("<= -%.1f", z_clip),
               "0",
               sprintf(">= +%.1f", z_clip))
  )
)

# -----------------------------------------------------------------------------
# Render
# -----------------------------------------------------------------------------
pdf_w <- 5.5
pdf_h <- 4.5

cairo_pdf(OUT_PDF, width = pdf_w, height = pdf_h, fallback_resolution = 300)
draw(
  ht,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE,
  padding = unit(c(2, 3, 3, 8), "mm")
)
dev.off()

# -----------------------------------------------------------------------------
# Source CSV
# -----------------------------------------------------------------------------
fwrite(agg[, .(display, group, cell_type, condition, mean_dev, n_donors)],
       OUT_CSV)

message(sprintf("[fig3_chromvar_celltype_progression] wrote %s", OUT_PDF))
message(sprintf("  TFs              : %d", nrow(mat_plot)))
message(sprintf("  Columns          : %d (%s × %s)",
                ncol(mat_plot), paste(CT_KEEP, collapse=","),
                paste(COND_LEVELS, collapse=",")))
message(sprintf("  Source CSV       : %s", OUT_CSV))
