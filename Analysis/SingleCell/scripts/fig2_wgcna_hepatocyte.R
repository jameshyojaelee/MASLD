#!/usr/bin/env Rscript
# fig3_wgcna_hepatocyte.R
# WGCNA co-expression modules on hepatocyte pseudobulk counts.
# Shows hepatocyte disease programs are organized into modules
# with quantitative correlation to MASLD severity (NAS, fibrosis).
#
# Gene selection: top 5000 most variable genes by IQR (avoids Ensembl-symbol
# mapping issues between DE results and pseudobulk row names).
# Disease correlation: matched to unified_metadata.csv via sample SRR IDs.

# Install WGCNA if needed (idempotent)
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  message("Installing WGCNA...")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("WGCNA", update = FALSE, ask = FALSE, quiet = TRUE)
}

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(edgeR)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 8)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SC_DIR <- file.path(BASE, "Analysis/SingleCell")
RESULTS <- file.path(SC_DIR, "results_gpu_v2")
OUT_DIR <- file.path(RESULTS, "fig2_data")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PB_FILE   <- file.path(RESULTS, "pseudobulk", "Hepatocytes_pseudobulk.csv")
# cell_type_proportions.csv has sc SRR IDs matching pseudobulk sample columns
# (condition, condition_harmonized, dataset columns present)
META_FILE <- file.path(RESULTS, "cell_type_proportions.csv")

# ---------------------------------------------------------------------------
# 1. Load pseudobulk counts
# ---------------------------------------------------------------------------
message("Loading pseudobulk counts...")
pb_raw <- fread(PB_FILE)

# First column is gene symbols (no header), rest are SRR sample IDs
gene_ids <- pb_raw[[1]]
pb_mat   <- as.matrix(pb_raw[, -1, with = FALSE])
rownames(pb_mat) <- gene_ids
colnames(pb_mat) <- names(pb_raw)[-1]
message(sprintf("Pseudobulk matrix: %d genes x %d samples", nrow(pb_mat), ncol(pb_mat)))

# Remove genes with all-zero counts
nonzero_genes <- rowSums(pb_mat > 0) >= 3
pb_mat <- pb_mat[nonzero_genes, ]
message(sprintf("After zero-filter: %d genes", nrow(pb_mat)))

# ---------------------------------------------------------------------------
# 2. TMM normalize + log2(CPM + 1)
# ---------------------------------------------------------------------------
message("Normalizing (TMM + logCPM)...")
dge <- DGEList(counts = pb_mat)
dge <- calcNormFactors(dge, method = "TMM")
lcpm <- cpm(dge, log = TRUE, prior.count = 1)

# ---------------------------------------------------------------------------
# 3. Select top 5000 variable genes by IQR (avoids Ensembl-symbol mapping)
# ---------------------------------------------------------------------------
message("Selecting top variable genes by IQR...")
gene_iqr <- apply(lcpm, 1, IQR)
top_genes <- names(sort(gene_iqr, decreasing = TRUE))[1:min(5000, nrow(lcpm))]
datExpr   <- t(lcpm[top_genes, ])   # WGCNA expects samples x genes
message(sprintf("Working matrix: %d samples x %d genes", nrow(datExpr), ncol(datExpr)))

# ---------------------------------------------------------------------------
# 4. Sample / gene QC
# ---------------------------------------------------------------------------
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  message(sprintf("After QC: %d samples x %d genes", nrow(datExpr), ncol(datExpr)))
}

# ---------------------------------------------------------------------------
# 5. Load metadata and build trait matrix
# ---------------------------------------------------------------------------
message("Loading metadata...")
meta <- tryCatch(fread(META_FILE), error = function(e) { message("Metadata not found: ", e$message); NULL })

sample_ids <- rownames(datExpr)  # SRR IDs

build_trait_matrix <- function(meta, sample_ids) {
  n <- length(sample_ids)
  traits <- data.frame(row.names = sample_ids)

  if (is.null(meta)) {
    message("No metadata available — skipping trait matrix")
    return(traits)
  }

  # Match by sample column (sc cell_type_proportions.csv uses "sample" = SRR ID)
  id_col <- NULL
  for (col in c("sample", "sample_id", "srr_id", "Run", "SRR")) {
    if (col %in% names(meta)) { id_col <- col; break }
  }
  if (is.null(id_col)) {
    message("No ID column found in metadata")
    return(traits)
  }

  meta_match <- meta[get(id_col) %chin% sample_ids]
  message(sprintf("Metadata matched: %d / %d samples", nrow(meta_match), n))

  if (nrow(meta_match) < 10) {
    message("Too few matched samples for trait correlation")
    return(traits)
  }

  idx <- match(sample_ids, meta_match[[id_col]])

  # Condition: 1 = MASLD/NASH/NAFLD, 0 = Control/Healthy
  if ("diagnosis_harmonized" %in% names(meta_match)) {
    diag_vals <- meta_match$diagnosis_harmonized[idx]
    traits$condition <- as.numeric(diag_vals %in% c("NASH", "MASLD", "NAFL", "Borderline",
                                                      "NAFLD", "Fibrosis"))
  } else if ("condition" %in% names(meta_match)) {
    cond_vals <- toupper(meta_match$condition[idx])
    traits$condition <- as.numeric(grepl("MASLD|NASH|NAFLD|DISEASE", cond_vals))
  }

  # NAS score
  for (col in c("NAS", "nas_score", "nas")) {
    if (col %in% names(meta_match)) {
      traits$nas_score <- suppressWarnings(as.numeric(meta_match[[col]][idx]))
      break
    }
  }

  # Fibrosis stage
  for (col in c("fibrosis_stage", "Fibrosis_stage", "fibrosis")) {
    if (col %in% names(meta_match)) {
      traits$fibrosis_stage <- suppressWarnings(as.numeric(meta_match[[col]][idx]))
      break
    }
  }

  # Remove all-NA columns
  traits <- traits[, colSums(!is.na(traits)) >= 10, drop = FALSE]
  message(sprintf("Trait columns with sufficient data: %s",
                  paste(names(traits), collapse = ", ")))
  traits
}

datTraits <- build_trait_matrix(meta, sample_ids)

# ---------------------------------------------------------------------------
# 6. Soft threshold selection
# ---------------------------------------------------------------------------
message("Selecting soft threshold power...")
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0,
                          networkType = "signed hybrid")

soft_power <- sft$powerEstimate
if (is.na(soft_power) || soft_power < 1) {
  good_fit <- sft$fitIndices[!is.na(sft$fitIndices$SFT.R.sq) &
                               sft$fitIndices$SFT.R.sq > 0.80, "Power"]
  soft_power <- if (length(good_fit) > 0) min(good_fit) else 6
}
message(sprintf("Soft threshold power: %d", soft_power))

# ---------------------------------------------------------------------------
# 7. Blockwise WGCNA
# ---------------------------------------------------------------------------
message("Running blockwise WGCNA...")
net <- blockwiseModules(
  datExpr,
  power             = soft_power,
  networkType       = "signed hybrid",
  TOMType           = "signed",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = FALSE,
  pamRespectsDendro = FALSE,
  maxBlockSize      = 8000,
  nThreads          = 8,
  verbose           = 0,
  saveTOMs          = FALSE
)
n_modules <- length(setdiff(unique(net$colors), "grey"))
message(sprintf("Modules found (excl. grey): %d", n_modules))

# ---------------------------------------------------------------------------
# 8. Module membership (kME) and gene table
# ---------------------------------------------------------------------------
MEs <- net$MEs
kME <- cor(datExpr, MEs, use = "pairwise.complete.obs")

module_genes <- data.table(
  gene         = colnames(datExpr),
  module_color = net$colors
)
module_genes <- module_genes[module_color != "grey"]

module_genes[, kME := mapply(function(g, m) {
  me_col <- paste0("ME", m)
  if (me_col %in% colnames(kME) && g %in% rownames(kME)) kME[g, me_col] else NA_real_
}, gene, module_color)]

out_genes <- file.path(OUT_DIR, "wgcna_module_genes.csv")
fwrite(module_genes, out_genes)
message(sprintf("Saved: %s (%d genes)", out_genes, nrow(module_genes)))

# ---------------------------------------------------------------------------
# 9. Module-trait correlation (only if datTraits has columns)
# ---------------------------------------------------------------------------
module_summary <- data.table(
  module_color = setdiff(unique(net$colors), "grey"),
  n_genes      = tabulate(factor(net$colors[net$colors != "grey"],
                                  levels = setdiff(unique(net$colors), "grey")))
)

if (ncol(datTraits) > 0) {
  message("Computing module-trait correlations...")
  module_trait_cor  <- cor(MEs, datTraits, use = "pairwise.complete.obs")
  module_trait_pmat <- corPvalueStudent(module_trait_cor, nrow(datExpr))

  pvals_flat <- as.vector(module_trait_pmat)
  padj_flat  <- p.adjust(pvals_flat, method = "BH")
  padj_mat   <- matrix(padj_flat, nrow = nrow(module_trait_pmat),
                        dimnames = dimnames(module_trait_pmat))

  for (trait in colnames(datTraits)) {
    me_rows <- rownames(module_trait_cor)
    module_colors_me <- sub("^ME", "", me_rows)

    r_col    <- paste0("r_",    trait)
    padj_col <- paste0("padj_", trait)

    module_summary[[r_col]]    <- NA_real_
    module_summary[[padj_col]] <- NA_real_

    for (i in seq_len(nrow(module_summary))) {
      mc     <- module_summary$module_color[i]
      me_key <- paste0("ME", mc)
      if (me_key %in% me_rows) {
        module_summary[[r_col]][i]    <- module_trait_cor[me_key, trait]
        module_summary[[padj_col]][i] <- padj_mat[me_key, trait]
      }
    }
  }
} else {
  message("No trait data available — skipping correlation")
}

# Add hub genes per module
hub_genes_per_mod <- chooseTopHubInEachModule(datExpr, net$colors, omitColors = "grey")
module_summary[, top_hub_gene := hub_genes_per_mod[module_color]]

# Top 5 hub genes (by kME)
module_summary[, top5_hubs := sapply(module_color, function(mc) {
  sub_dt <- module_genes[module_color == mc][!is.na(kME)][order(-kME)]
  paste(head(sub_dt$gene, 5), collapse = ";")
})]

out_trait <- file.path(OUT_DIR, "wgcna_module_trait_correlation.csv")
fwrite(module_summary, out_trait)
message(sprintf("Saved: %s (%d modules, %d trait columns)",
                out_trait, nrow(module_summary), ncol(datTraits)))

# Hub gene table (all modules, top 500 by kME)
hub_tbl <- module_genes[!is.na(kME)][order(module_color, -kME)]
out_hub <- file.path(OUT_DIR, "wgcna_hub_genes.csv")
fwrite(hub_tbl[, .(gene, module_color, kME)][1:min(500, .N)], out_hub)
message(sprintf("Saved: %s", out_hub))

message("=== WGCNA complete ===")
