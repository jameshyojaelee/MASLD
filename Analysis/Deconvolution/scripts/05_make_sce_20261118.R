# Wrap stratified 100,000-cell NUCLEUS-FIX extract into a SingleCellExperiment RDS.
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(data.table)
})

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CAND <- file.path(ROOT, "Analysis/Deconvolution/candidates/loco_20261118")
TMP <- file.path(CAND, "tmp_sce")

mat <- Matrix::readMM(file.path(TMP, "counts.mtx"))
mat <- as(mat, "dgCMatrix")

obs <- fread(file.path(TMP, "obs.tsv"), data.table = FALSE)
var <- fread(file.path(TMP, "var.tsv"), data.table = FALSE)

stopifnot(nrow(mat) == nrow(var), ncol(mat) == nrow(obs), ncol(mat) == 100000L)
stopifnot(!anyDuplicated(var$symbol), !anyDuplicated(obs$cell_id))

rownames(mat) <- var$symbol
colnames(mat) <- obs$cell_id
rownames(var) <- var$symbol
rownames(obs) <- obs$cell_id
obs$celltype <- obs$cell_type

sce <- SingleCellExperiment(
  assays = list(counts = mat),
  colData = obs,
  rowData = var
)

out_rds <- file.path(CAND, "reference_nucleus_fix_100k_sce.rds")
saveRDS(sce, out_rds, compress = TRUE)
cat("Saved SCE:", out_rds, "dim:", paste(dim(sce), collapse = "x"), "\n")
print(table(colData(sce)$dataset))
print(table(colData(sce)$cell_type))
