#!/usr/bin/env Rscript
# pseudobulk_de_coarse_lineage.R
# Coarse 7-LINEAGE scRNA pseudobulk DE (MASLD vs Healthy), for the redesigned
# Fig 3E genetics<->expression disconnect panel.
#
# Mirrors pseudobulk_de.R exactly (same h5ad metadata, same limma-voom
# ~dataset+condition, same filters) but aggregates the FINE per-cell-type
# pseudobulk RAW COUNTS into 7 lineages first (raw counts are additive).
#
# Lineages: Hepatocytes / Cholangiocytes / Endothelial / Mesenchymal /
#           Myeloid / T_NK / B_Plasma.
#
# Input:  results_gpu_v2/pseudobulk/{FineType}_pseudobulk.csv (genes x samples, raw counts)
# Output: results_gpu_v2/pseudobulk_de_coarse/{Lineage}_de.csv (limma-voom DE)

suppressPackageStartupMessages({
  library(rhdf5); library(limma); library(edgeR); library(data.table)
})

BASE     <- Sys.getenv("MASLD_PROJECT_ROOT",
                       "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
scvi_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2")
h5ad_file<- file.path(scvi_dir, "integrated_atlas.h5ad")
pb_dir   <- file.path(scvi_dir, "pseudobulk")
out_dir  <- file.path(scvi_dir, "pseudobulk_de_coarse")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 7-lineage map: canonical FINE cell-type stem -> lineage --------------
LINEAGE_MAP <- c(
  "Hepatocytes"             = "Hepatocytes",
  "Cholangiocytes"          = "Cholangiocytes",
  "Endothelial_cells"       = "Endothelial",
  "Fibroblasts"             = "Mesenchymal",
  "Macrophages"             = "Myeloid",
  "Mono+mono_derived_cells" = "Myeloid",
  "cDC1s" = "Myeloid", "cDC2s" = "Myeloid", "pDCs" = "Myeloid",
  "Mig.cDCs" = "Myeloid", "Basophils" = "Myeloid", "Neutrophils" = "Myeloid",
  "T_cells"            = "T_NK",
  "Circulating_NK_NKT" = "T_NK",
  "Resident_NK"        = "T_NK",
  "B_cells"      = "B_Plasma",
  "Plasma_cells" = "B_Plasma")

# ---- 1. sample metadata from h5ad (identical to pseudobulk_de.R) ----------
message("Extracting sample metadata from h5ad...")
stopifnot(file.exists(h5ad_file))
items <- h5ls(h5ad_file)
read_categorical <- function(h5file, key) {
  gp <- paste0("/obs/", key)
  obs_groups <- subset(items, group == "/obs" & otype == "H5I_GROUP")$name
  if (key %in% obs_groups) {
    codes <- h5read(h5file, paste0(gp, "/codes"))
    cats  <- h5read(h5file, paste0(gp, "/categories"))
    return(cats[codes + 1L])
  }
  vals <- h5read(h5file, gp)
  if (is.list(vals) && "codes" %in% names(vals)) return(vals$categories[vals$codes + 1L])
  as.vector(vals)
}
sample_vec  <- read_categorical(h5ad_file, "sample")
dataset_vec <- read_categorical(h5ad_file, "dataset")
cond_vec <- tryCatch(read_categorical(h5ad_file, "condition_harmonized"),
                     error = function(e) read_categorical(h5ad_file, "condition"))
prep_vec <- tryCatch(read_categorical(h5ad_file, "preparation_method"),
                     error = function(e) rep("unsorted", length(sample_vec)))
cell_meta <- data.table(sample = sample_vec, dataset = dataset_vec,
                        condition = cond_vec, preparation_method = prep_vec)
sample_meta <- cell_meta[, .(dataset = dataset[1],
  condition = names(sort(table(condition), decreasing = TRUE))[1],
  preparation_method = preparation_method[1]), by = sample]
sample_meta <- sample_meta[preparation_method %in% c("nuclei","cd45_negative","unsorted")]
sample_meta <- sample_meta[condition %in% c("Healthy", "MASLD")]
message("  DE-eligible samples: ", nrow(sample_meta), " (",
        paste(names(table(sample_meta$condition)), table(sample_meta$condition),
              sep="=", collapse=", "), ")")

# ---- 2. aggregate fine pseudobulk raw counts -> lineage count matrix ------
pb_files <- list.files(pb_dir, pattern = "_pseudobulk\\.csv$", full.names = TRUE)
# canonical stem: strip suffix; keep only files whose stem is in the map
stem_of <- function(f) sub("_pseudobulk\\.csv$", "", basename(f))
pb_files <- pb_files[stem_of(pb_files) %in% names(LINEAGE_MAP)]

aggregate_lineage <- function(files) {
  # union genes x union samples, sum counts (missing -> 0)
  mats <- lapply(files, function(f) {
    dt <- fread(f); g <- dt[[1]]
    m <- as.matrix(dt[, -1, with = FALSE]); storage.mode(m) <- "numeric"
    rownames(m) <- g; m
  })
  all_genes <- sort(unique(unlist(lapply(mats, rownames))))
  all_samps <- sort(unique(unlist(lapply(mats, colnames))))
  acc <- matrix(0, length(all_genes), length(all_samps),
                dimnames = list(all_genes, all_samps))
  for (m in mats) acc[rownames(m), colnames(m)] <- acc[rownames(m), colnames(m)] + m
  acc
}

for (ln in unique(LINEAGE_MAP)) {
  fine_files <- pb_files[LINEAGE_MAP[stem_of(pb_files)] == ln]
  if (!length(fine_files)) { message("skip ", ln, " (no files)"); next }
  message("\n--- ", ln, " <- ", paste(stem_of(fine_files), collapse=", "), " ---")
  counts_mat <- aggregate_lineage(fine_files)

  # match samples to DE-eligible metadata
  midx <- match(colnames(counts_mat), sample_meta$sample)
  keep <- !is.na(midx)
  counts_mat <- counts_mat[, keep, drop = FALSE]
  matched <- sample_meta[midx[keep]]
  cond_tab <- table(matched$condition)
  valid <- names(cond_tab[cond_tab >= 2])
  if (length(valid) < 2) { message("  skip: <2 conditions w/ >=2 samples"); next }
  kc <- matched$condition %in% valid
  counts_mat <- counts_mat[, kc, drop = FALSE]; matched <- matched[kc]

  # gene filter: >=5 counts in >=3 samples
  counts_mat <- counts_mat[rowSums(counts_mat >= 5) >= 3, , drop = FALSE]
  message("  samples=", ncol(counts_mat), " (",
          paste(names(table(matched$condition)), table(matched$condition), sep="=", collapse=", "),
          "); genes=", nrow(counts_mat))
  if (nrow(counts_mat) < 100) { message("  skip: too few genes"); next }

  condition <- factor(matched$condition, levels = c("Healthy","MASLD"))
  dataset   <- factor(matched$dataset)
  design <- if (length(unique(matched$dataset)) > 1)
              model.matrix(~ dataset + condition) else model.matrix(~ condition)
  dge <- calcNormFactors(DGEList(counts = counts_mat))
  v   <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  ci  <- which(colnames(design) == "conditionMASLD")
  res <- setDT(topTable(fit, coef = ci, number = Inf, sort.by = "none"))
  res[, gene := rownames(v$E)]
  setnames(res, c("adj.P.Val","P.Value","t"), c("padj","pvalue","t_stat"), skip_absent = TRUE)
  res[, `:=`(cell_type = ln, contrast = "MASLD_vs_Healthy")]
  fwrite(res, file.path(out_dir, paste0(ln, "_de.csv")))
  message("  saved ", ln, "_de.csv (", nrow(res), " genes, ",
          sum(res$padj < 0.05, na.rm = TRUE), " sig padj<0.05; mean|logFC|=",
          round(mean(abs(res$logFC), na.rm = TRUE), 3), ")")
}
message("\nCoarse-lineage pseudobulk DE complete -> ", out_dir)
