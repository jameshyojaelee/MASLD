#!/usr/bin/env Rscript
# 310b: Cross-modal — Spatial ssGSEA + GWAS-ATAC TF activity
#
# Spatial: Compute ssGSEA scores for each hepatocyte subtype's top 50
# markers across Visium spots from GSE192741.
#
# GWAS-ATAC: Score 12 disrupted disease regulon TFs across hepatocyte
# subtypes using regulon target gene expression.
#
# Input:
#   hepatocyte_subtypes/subtype_markers.csv                (from 309)
#   hepatocyte_subtypes/hepatocyte_subtype_metadata.csv    (from 309)
#   hepatocyte_subtypes/subtype_disease_enrichment.csv     (from 309)
#   hepatocyte_subtypes/hepatocyte_atlas.h5ad              (from 308)
#   Analysis/Spatial/results/preprocessed/merged_spatial.h5ad  (GSE192741)
#   GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv
#   Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv
#
# Output:
#   crossmodal/spatial/subtype_spatial_scores.csv
#   crossmodal/spatial/subtype_spatial_zonation_correlation.csv
#   crossmodal/gwas_atac/subtype_tf_activity.csv
#   crossmodal/gwas_atac/subtype_gwas_atac_enrichment.csv
#
# Environment: rnaseq (CPU)

suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(rhdf5)
  library(Matrix)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SC_DIR      <- file.path(BASE, "Analysis/SingleCell")
RESULTS     <- file.path(SC_DIR, "results_gpu_v2")
SUB_DIR     <- file.path(RESULTS, "hepatocyte_subtypes")
SPATIAL_DIR <- file.path(BASE, "Analysis/Spatial/results")
GWAS_ATAC_DIR <- file.path(BASE, "GWAS/finemapping/results/gwas_atac")
SCENIC_DIR  <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus")

OUT_SPATIAL <- file.path(SUB_DIR, "crossmodal", "spatial")
OUT_GWAS    <- file.path(SUB_DIR, "crossmodal", "gwas_atac")
dir.create(OUT_SPATIAL, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_GWAS, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helper: read h5ad categorical column (categories + codes)
# ---------------------------------------------------------------------------
read_h5ad_categorical <- function(h5file, col_path) {
  cats  <- as.character(h5read(h5file, paste0(col_path, "/categories")))
  codes <- h5read(h5file, paste0(col_path, "/codes"))
  cats[codes + 1L]
}

# ---------------------------------------------------------------------------
# Helper: read a sparse matrix from h5ad group path (CSR/CSC aware)
# Matches 310a pattern for consistency.
# ---------------------------------------------------------------------------
read_sparse_h5 <- function(h5file, group_path, n_obs, n_var) {
  data_vec <- h5read(h5file, paste0(group_path, "/data"))
  indices  <- h5read(h5file, paste0(group_path, "/indices"))
  indptr   <- h5read(h5file, paste0(group_path, "/indptr"))

  enc_type <- tryCatch(
    h5readAttributes(h5file, group_path)[["encoding-type"]],
    error = function(e) NULL
  )

  if (!is.null(enc_type) && grepl("csr", enc_type, ignore.case = TRUE)) {
    mat <- sparseMatrix(
      j = as.integer(indices) + 1L,
      p = as.integer(indptr),
      x = as.numeric(data_vec),
      dims = c(n_obs, n_var),
      repr = "C"
    )
  } else if (!is.null(enc_type) && grepl("csc", enc_type, ignore.case = TRUE)) {
    mat <- sparseMatrix(
      i = as.integer(indices) + 1L,
      p = as.integer(indptr),
      x = as.numeric(data_vec),
      dims = c(n_obs, n_var)
    )
  } else {
    # Default: infer from indptr length
    if (length(indptr) == n_obs + 1L) {
      mat <- sparseMatrix(
        j = as.integer(indices) + 1L,
        p = as.integer(indptr),
        x = as.numeric(data_vec),
        dims = c(n_obs, n_var),
        repr = "C"
      )
    } else {
      mat <- sparseMatrix(
        i = as.integer(indices) + 1L,
        p = as.integer(indptr),
        x = as.numeric(data_vec),
        dims = c(n_obs, n_var)
      )
    }
  }
  return(mat)
}


# ═══════════════════════════════════════════════════════════════════════════
# PART A: Spatial ssGSEA
# ═══════════════════════════════════════════════════════════════════════════
message("=== PART A: Spatial ssGSEA ===\n")

# Load subtype markers
markers <- fread(file.path(SUB_DIR, "subtype_markers.csv"))
subtypes <- sort(unique(markers$subtype))
message("Subtypes: ", paste(subtypes, collapse = ", "))

# Build gene sets: top 50 markers per subtype (already ordered by significance in 309)
gene_sets <- list()
for (st in subtypes) {
  genes <- markers[subtype == st, names]
  if (length(genes) > 50) genes <- genes[1:50]
  gene_sets[[paste0("Hep_", st)]] <- genes
}
message("Gene sets built: ", length(gene_sets), " subtypes, ",
        paste(sapply(gene_sets, length), collapse = "/"), " genes each")

# Load spatial expression matrix from h5ad via rhdf5
# merged_spatial.h5ad is ~24K spots x 10K genes (CSR)
message("Loading spatial data...")
sp_h5    <- file.path(SPATIAL_DIR, "preprocessed/merged_spatial.h5ad")
stopifnot(file.exists(sp_h5))

sp_genes    <- as.character(h5read(sp_h5, "var/_index"))
sp_barcodes <- as.character(h5read(sp_h5, "obs/_index"))
n_spots     <- length(sp_barcodes)
n_sp_genes  <- length(sp_genes)

message("Reading sparse expression matrix...")
sp_sparse <- read_sparse_h5(sp_h5, "X", n_obs = n_spots, n_var = n_sp_genes)
rownames(sp_sparse) <- sp_barcodes
colnames(sp_sparse) <- sp_genes
message("Spatial matrix: ", n_spots, " spots x ", n_sp_genes, " genes")

# Read sample/condition info (categorical columns)
sp_sample  <- read_h5ad_categorical(sp_h5, "obs/sample_id")
sp_dataset <- read_h5ad_categorical(sp_h5, "obs/dataset")

# Transpose to genes x spots for GSVA
sp_mat <- t(sp_sparse)  # Now genes x spots (dense required by GSVA)
rm(sp_sparse); gc()

# Convert to dense for GSVA (24K spots x 10K genes is manageable)
sp_mat <- as.matrix(sp_mat)
message("Dense spatial matrix: ", nrow(sp_mat), " genes x ", ncol(sp_mat), " spots")

# Filter gene sets to genes present in spatial data
for (nm in names(gene_sets)) {
  gene_sets[[nm]] <- gene_sets[[nm]][gene_sets[[nm]] %in% sp_genes]
}
gene_sets <- gene_sets[sapply(gene_sets, length) >= 5]
message("Gene sets with >=5 spatial genes: ", length(gene_sets))

if (length(gene_sets) > 0) {
  message("Running ssGSEA on spatial data...")
  # GSVA >= 1.50 (Bioconductor 3.18+) uses ssgseaParam objects
  gsva_version <- packageVersion("GSVA")
  message("  GSVA version: ", as.character(gsva_version))

  if (gsva_version >= "1.50") {
    param <- ssgseaParam(exprData = sp_mat, geneSets = gene_sets)
    gsva_res <- gsva(param, verbose = FALSE)
  } else {
    gsva_res <- gsva(sp_mat, gene_sets, method = "ssgsea", verbose = FALSE)
  }
  # gsva_res: gene_sets x spots

  # Build output data.table (spots as rows)
  spatial_scores <- as.data.table(t(gsva_res))
  spatial_scores[, spot := sp_barcodes]
  spatial_scores[, dataset := sp_dataset]

  fwrite(spatial_scores, file.path(OUT_SPATIAL, "subtype_spatial_scores.csv"))
  message("Saved spatial scores: ", nrow(spatial_scores), " spots x ",
          length(gene_sets), " subtypes")

  # Zonation correlation: periportal / pericentral markers as reference axis
  zonation_pp <- c("ASS1", "HAL", "CPS1", "ALDOB", "PCK1")
  zonation_pc <- c("GLUL", "CYP2E1", "CYP1A2", "CYP3A4")
  pp_present <- zonation_pp[zonation_pp %in% sp_genes]
  pc_present <- zonation_pc[zonation_pc %in% sp_genes]

  if (length(pp_present) >= 2 && length(pc_present) >= 2) {
    pp_score <- colMeans(sp_mat[pp_present, , drop = FALSE])
    pc_score <- colMeans(sp_mat[pc_present, , drop = FALSE])

    zon_corr_rows <- list()
    for (nm in rownames(gsva_res)) {
      subtype_score <- gsva_res[nm, ]
      rho_pp <- cor(subtype_score, pp_score, method = "spearman",
                    use = "complete.obs")
      rho_pc <- cor(subtype_score, pc_score, method = "spearman",
                    use = "complete.obs")
      zon_corr_rows[[nm]] <- data.table(
        subtype = nm,
        rho_periportal = rho_pp,
        rho_pericentral = rho_pc
      )
    }
    zon_corr <- rbindlist(zon_corr_rows)
    fwrite(zon_corr, file.path(OUT_SPATIAL,
                               "subtype_spatial_zonation_correlation.csv"))
    message("Zonation correlation:")
    print(zon_corr)
  } else {
    message("WARNING: Insufficient zonation markers in spatial data ",
            "(PP: ", length(pp_present), ", PC: ", length(pc_present), ")")
  }
} else {
  message("WARNING: No gene sets had >=5 genes in spatial data")
}

# Free spatial memory
rm(sp_mat); if (exists("gsva_res")) rm(gsva_res); gc()


# ═══════════════════════════════════════════════════════════════════════════
# PART B: GWAS-ATAC TF Activity
# ═══════════════════════════════════════════════════════════════════════════
message("\n=== PART B: GWAS-ATAC TF Activity ===\n")

# Load disrupted TFs (only those in disease regulons)
motif <- fread(file.path(GWAS_ATAC_DIR, "motif_disruption_scores.csv"))
disrupted_tfs <- unique(motif[motif_in_disease_regulon == TRUE, tf_name])
message("Disrupted disease-regulon TFs: ", paste(disrupted_tfs, collapse = ", "))

# Load hepatocyte regulons (target gene lists from SCENIC+)
regulons <- fread(file.path(SCENIC_DIR, "hepatocyte_regulons.csv"))
message("Hepatocyte regulons: ", length(unique(regulons$tf_name)), " TFs")

# Build regulon gene sets for disrupted TFs (require >= 3 target genes)
regulon_sets <- list()
for (tf in disrupted_tfs) {
  targets <- unique(regulons[tf_name == tf, target_gene])
  if (length(targets) >= 3) {
    regulon_sets[[tf]] <- targets
  }
}
message("Regulon gene sets for disrupted TFs: ", length(regulon_sets),
        " (of ", length(disrupted_tfs), " disrupted TFs)")

# Load per-cell metadata (need subtype labels + cell barcodes)
message("Loading hepatocyte subtype metadata...")
cell_meta <- fread(file.path(SUB_DIR, "hepatocyte_subtype_metadata.csv"))
# Script 309 uses pandas to_csv(index=True) -> first column is cell barcode
if ("V1" %in% names(cell_meta)) {
  setnames(cell_meta, "V1", "cell_barcode")
} else if (!("cell_barcode" %in% names(cell_meta))) {
  setnames(cell_meta, names(cell_meta)[1], "cell_barcode")
}
cell_subtypes <- sort(unique(cell_meta$hepatocyte_subtype))
message("  ", nrow(cell_meta), " cells, ", length(cell_subtypes), " subtypes")

# Load normalized hepatocyte expression from h5ad for TF scoring
# This is the full 657K x 37K matrix — requires ~128GB RAM
message("Loading hepatocyte expression matrix for TF scoring...")
hep_h5 <- file.path(SUB_DIR, "hepatocyte_atlas.h5ad")
stopifnot(file.exists(hep_h5))

hep_genes    <- as.character(h5read(hep_h5, "var/_index"))
hep_barcodes <- as.character(h5read(hep_h5, "obs/_index"))
n_hep_cells  <- length(hep_barcodes)
n_hep_genes  <- length(hep_genes)
message("  Dimensions: ", n_hep_cells, " cells x ", n_hep_genes, " genes")

message("Reading sparse expression matrix (this may take a few minutes)...")
hep_sparse <- read_sparse_h5(hep_h5, "X",
                              n_obs = n_hep_cells, n_var = n_hep_genes)
colnames(hep_sparse) <- hep_genes
rownames(hep_sparse) <- hep_barcodes
message("  Loaded: ", format(object.size(hep_sparse), units = "GB"))

# Map barcodes to subtypes
barcode_to_subtype <- setNames(cell_meta$hepatocyte_subtype,
                               cell_meta$cell_barcode)

# Compute regulon activity scores per subtype
message("Computing regulon activity scores per subtype...")
tf_activity_full <- list()

for (tf in names(regulon_sets)) {
  targets <- regulon_sets[[tf]]
  targets_present <- targets[targets %in% hep_genes]
  if (length(targets_present) < 3) {
    message("  Skipping ", tf, " — only ", length(targets_present),
            " target genes found")
    next
  }

  # Extract target gene columns and compute mean expression per cell
  target_expr <- hep_sparse[, targets_present, drop = FALSE]
  cell_scores <- Matrix::rowMeans(target_expr)

  # Mean per subtype
  for (st in cell_subtypes) {
    st_cells <- names(barcode_to_subtype)[barcode_to_subtype == st]
    st_cells <- st_cells[st_cells %in% names(cell_scores)]
    if (length(st_cells) == 0) next

    tf_activity_full[[paste0(tf, "_", st)]] <- data.table(
      tf = tf,
      subtype = st,
      mean_activity = mean(cell_scores[st_cells]),
      n_targets = length(targets_present),
      n_cells = length(st_cells)
    )
  }
  message("  ", tf, ": ", length(targets_present), " targets scored")
}

rm(hep_sparse); gc()

if (length(tf_activity_full) > 0) {
  tf_dt <- rbindlist(tf_activity_full)
  fwrite(tf_dt, file.path(OUT_GWAS, "subtype_tf_activity.csv"))
  message("Saved TF activity: ", nrow(tf_dt), " rows")

  # Pivot for display: TFs x subtypes
  tf_wide <- dcast(tf_dt, tf ~ subtype, value.var = "mean_activity")
  message("TF activity matrix:")
  print(tf_wide)

  # Statistical tests: Kruskal-Wallis across subtypes for each TF
  # + Wilcoxon MASLD-enriched vs rest
  enrich_dt <- fread(file.path(SUB_DIR, "subtype_disease_enrichment.csv"))
  masld_subtypes <- enrich_dt[enrichment_class == "MASLD_enriched", subtype]
  message("MASLD-enriched subtypes: ", paste(masld_subtypes, collapse = ", "))

  test_rows <- list()
  for (tf_name_i in unique(tf_dt$tf)) {
    tf_sub <- tf_dt[tf == tf_name_i]
    if (nrow(tf_sub) < 3) next

    # Kruskal-Wallis across all subtypes
    kw <- kruskal.test(mean_activity ~ subtype, data = tf_sub)

    # Wilcoxon: MASLD-enriched subtypes vs rest
    if (length(masld_subtypes) > 0 &&
        length(masld_subtypes) < nrow(tf_sub)) {
      masld_vals <- tf_sub[subtype %in% masld_subtypes, mean_activity]
      other_vals <- tf_sub[!subtype %in% masld_subtypes, mean_activity]

      if (length(masld_vals) >= 1 && length(other_vals) >= 1) {
        wt <- wilcox.test(masld_vals, other_vals, alternative = "greater",
                          exact = FALSE)
        test_rows[[tf_name_i]] <- data.table(
          tf = tf_name_i,
          kw_stat = kw$statistic,
          kw_pval = kw$p.value,
          wilcox_stat = wt$statistic,
          wilcox_pval = wt$p.value,
          masld_mean = mean(masld_vals),
          other_mean = mean(other_vals)
        )
      }
    }
  }

  if (length(test_rows) > 0) {
    test_dt <- rbindlist(test_rows)
    fwrite(test_dt, file.path(OUT_GWAS, "subtype_gwas_atac_enrichment.csv"))
    message("GWAS-ATAC enrichment tests:")
    print(test_dt)
  } else {
    message("WARNING: No TFs had enough subtypes for statistical testing")
  }
} else {
  message("WARNING: No TFs could be scored (check regulon target overlap)")
}

message("\n310b complete.")
