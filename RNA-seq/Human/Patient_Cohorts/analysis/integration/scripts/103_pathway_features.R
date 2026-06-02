#!/usr/bin/env Rscript
# 103_pathway_features.R
# ---------------------------------------------------------------------------
# Compute ssGSEA and GSVA pathway activity scores for all 1,444 samples.
# Uses Hallmark (50), KEGG (~186), and Reactome (~674) gene sets from MSigDB.
#
# Input:  results/integration/merged_dge.rds (34,453 genes x 1,444 samples)
# Output: results/staging_classifier/pathway_scores_ssgsea.rds (matrix)
#         results/staging_classifier/pathway_scores_gsva.rds   (matrix)
#         results/staging_classifier/pathway_feature_names.csv
#
# Usage: Rscript 103_pathway_features.R
# SLURM: cpu, 8 CPUs, 64G, 48h
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(msigdbr)
})

# Install GSVA if not present
if (!requireNamespace("GSVA", quietly = TRUE)) {
  cat("GSVA not found. Installing from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("GSVA", ask = FALSE, update = FALSE)
}
library(GSVA)

set.seed(42)

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results/integration")
ME     <- file.path(BASE, "RNA-seq/results/multi_evidence")
OUTDIR <- file.path(INT, "results/staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 103: Pathway Features (ssGSEA / GSVA) ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ============================================================
# STEP 1: Load expression data and compute TMM logCPM
# ============================================================
cat("--- Step 1: Load DGE and compute logCPM ---\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
cat("DGE loaded:", nrow(dge), "genes x", ncol(dge), "samples\n")

dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
cat("logCPM matrix:", nrow(logcpm), "x", ncol(logcpm), "\n\n")

# ============================================================
# STEP 2: Map Ensembl IDs to gene symbols
# ============================================================
cat("--- Step 2: Map Ensembl IDs to gene symbols ---\n")

# Load symbol mapping from multi-evidence atlas
atlas <- fread(file.path(ME, "multi_evidence_atlas.csv"),
               select = c("ensembl_id", "human_symbol"))
atlas <- atlas[!is.na(human_symbol) & human_symbol != ""]
atlas[, ensembl_clean := sub("\\..*", "", ensembl_id)]
atlas <- atlas[!duplicated(ensembl_clean)]

# Clean rownames of logCPM
ensembl_clean <- sub("\\..*", "", rownames(logcpm))
idx <- match(ensembl_clean, atlas$ensembl_clean)
mapped <- !is.na(idx)

cat("Genes mapped to symbols:", sum(mapped), "of", nrow(logcpm),
    sprintf("(%.1f%%)\n", 100 * sum(mapped) / nrow(logcpm)))

# Subset to mapped genes
logcpm_sym <- logcpm[mapped, ]
symbols <- atlas$human_symbol[idx[mapped]]

# Handle duplicate symbols: keep the one with highest mean expression
mean_expr <- rowMeans(logcpm_sym)
dup_syms <- symbols[duplicated(symbols)]
if (length(dup_syms) > 0) {
  cat("Resolving", length(unique(dup_syms)), "duplicate symbols (keeping highest mean expression)\n")
  keep <- rep(TRUE, length(symbols))
  for (s in unique(dup_syms)) {
    which_dup <- which(symbols == s)
    best <- which_dup[which.max(mean_expr[which_dup])]
    keep[setdiff(which_dup, best)] <- FALSE
  }
  logcpm_sym <- logcpm_sym[keep, ]
  symbols <- symbols[keep]
}

rownames(logcpm_sym) <- symbols
cat("Final symbol-mapped matrix:", nrow(logcpm_sym), "genes x", ncol(logcpm_sym), "samples\n\n")

# ============================================================
# STEP 3: Load MSigDB gene sets
# ============================================================
cat("--- Step 3: Load MSigDB gene sets ---\n")

# Hallmark
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
cat("Hallmark:", length(hallmark_sets), "gene sets\n")

# KEGG (try KEGG_MEDICUS first, fall back to KEGG)
kegg_df <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
  error = function(e) {
    cat("  CP:KEGG_MEDICUS not found, falling back to CP:KEGG\n")
    msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG")
  }
)
kegg_sets <- if (!is.null(kegg_df) && nrow(kegg_df) > 0) {
  split(kegg_df$gene_symbol, kegg_df$gs_name)
} else {
  cat("WARNING: No KEGG gene sets found.\n")
  list()
}
cat("KEGG:", length(kegg_sets), "gene sets\n")

# Reactome
reactome_df <- msigdbr(species = "Homo sapiens", collection = "C2",
                        subcollection = "CP:REACTOME")
reactome_sets <- split(reactome_df$gene_symbol, reactome_df$gs_name)
cat("Reactome:", length(reactome_sets), "gene sets\n")

# Combine all
all_sets <- c(hallmark_sets, kegg_sets, reactome_sets)
cat("\nTotal gene sets:", length(all_sets), "\n")

# Filter to sets with at least 10 and at most 500 genes present in our data
set_sizes <- sapply(all_sets, function(gs) sum(gs %in% rownames(logcpm_sym)))
all_sets <- all_sets[set_sizes >= 10 & set_sizes <= 500]
cat("Gene sets after size filter (10-500 mapped genes):", length(all_sets), "\n\n")

# ============================================================
# STEP 4: Run ssGSEA
# ============================================================
cat("--- Step 4: Run ssGSEA ---\n")
cat("This may take 30-60 minutes for", ncol(logcpm_sym), "samples...\n")
t0 <- Sys.time()

# GSVA >= 1.50 uses gsvaParam/ssgseaParam objects; detect API version
gsva_version <- packageVersion("GSVA")
cat("GSVA version:", as.character(gsva_version), "\n")

if (gsva_version >= "1.50") {
  # New API (Bioconductor 3.18+)
  ssgsea_param <- ssgseaParam(exprData = logcpm_sym, geneSets = all_sets,
                               normalize = TRUE)
  ssgsea_scores <- gsva(ssgsea_param, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(8))
} else {
  # Legacy API
  ssgsea_scores <- gsva(logcpm_sym, all_sets, method = "ssgsea",
                         kcdf = "Gaussian", parallel.sz = 8, verbose = TRUE)
}

t1 <- Sys.time()
cat("\nssGSEA completed in", round(difftime(t1, t0, units = "mins"), 1), "minutes\n")
cat("ssGSEA score matrix:", nrow(ssgsea_scores), "pathways x", ncol(ssgsea_scores), "samples\n")
cat("Score range:", round(range(ssgsea_scores), 3), "\n\n")

# ============================================================
# STEP 5: Run GSVA
# ============================================================
cat("--- Step 5: Run GSVA ---\n")
cat("This may take 30-60 minutes for", ncol(logcpm_sym), "samples...\n")
t0 <- Sys.time()

if (gsva_version >= "1.50") {
  gsva_param <- gsvaParam(exprData = logcpm_sym, geneSets = all_sets,
                           kcdf = "Gaussian")
  gsva_scores <- gsva(gsva_param, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(8))
} else {
  gsva_scores <- gsva(logcpm_sym, all_sets, method = "gsva",
                       kcdf = "Gaussian", parallel.sz = 8, verbose = TRUE)
}

t1 <- Sys.time()
cat("\nGSVA completed in", round(difftime(t1, t0, units = "mins"), 1), "minutes\n")
cat("GSVA score matrix:", nrow(gsva_scores), "pathways x", ncol(gsva_scores), "samples\n")
cat("Score range:", round(range(gsva_scores), 3), "\n\n")

# ============================================================
# STEP 6: Save outputs
# ============================================================
cat("--- Step 6: Save outputs ---\n")

saveRDS(ssgsea_scores, file.path(OUTDIR, "pathway_scores_ssgsea.rds"))
cat("Saved:", file.path(OUTDIR, "pathway_scores_ssgsea.rds"), "\n")

saveRDS(gsva_scores, file.path(OUTDIR, "pathway_scores_gsva.rds"))
cat("Saved:", file.path(OUTDIR, "pathway_scores_gsva.rds"), "\n")

# Save pathway feature names
pathway_names <- data.table(
  pathway = rownames(ssgsea_scores),
  collection = fifelse(grepl("^HALLMARK_", rownames(ssgsea_scores)), "Hallmark",
                fifelse(grepl("^KEGG_|^KEGG_MEDICUS_", rownames(ssgsea_scores)), "KEGG",
                fifelse(grepl("^REACTOME_", rownames(ssgsea_scores)), "Reactome", "Other")))
)
fwrite(pathway_names, file.path(OUTDIR, "pathway_feature_names.csv"))
cat("Saved:", file.path(OUTDIR, "pathway_feature_names.csv"), "\n")

cat("\n=== Summary ===\n")
cat("ssGSEA:", nrow(ssgsea_scores), "pathways x", ncol(ssgsea_scores), "samples\n")
cat("GSVA:", nrow(gsva_scores), "pathways x", ncol(gsva_scores), "samples\n")
cat("Collections: Hallmark =", pathway_names[collection == "Hallmark", .N],
    ", KEGG =", pathway_names[collection == "KEGG", .N],
    ", Reactome =", pathway_names[collection == "Reactome", .N], "\n")
cat("\nFinished:", as.character(Sys.time()), "\n")
