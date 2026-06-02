#!/usr/bin/env Rscript
# 191_pseudotime_signatures.R
# Transfer scRNA pseudotime transition signatures to bulk RNA-seq via ssGSEA.
#
# Strategy: Extract top genes from each pseudotime transition program
# (gene_dynamics files: 2000 genes x 100 bins per cell type), define
# UP/DOWN gene sets per transition, score on bulk expression via ssGSEA.
#
# SLURM: sbatch --partition=cpu --cpus-per-task=8 --mem=64G --time=48:00:00
# Env: micromamba activate rnaseq

library(data.table)
library(GSVA)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUTDIR <- file.path(INT, "results/multiprogram")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PSEUDOTIME_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudotime")

cat("=== 191: Pseudotime Transition Signature Scoring ===\n")

# ---------------------------------------------------------------------------
# 1. Load bulk expression matrix + symbol-to-Ensembl mapping
# ---------------------------------------------------------------------------
cat("Loading bulk expression matrix...\n")
dge <- readRDS(file.path(INT, "results/integration/merged_dge.rds"))
# DGEList: extract logCPM
expr <- edgeR::cpm(dge, log = TRUE)
cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(expr), ncol(expr)))

# Build symbol-to-Ensembl mapping from dream results
# Bulk uses Ensembl IDs (ENSG...) as rownames; pseudotime uses gene symbols
cat("Building symbol-to-Ensembl mapping...\n")
dream <- fread(file.path(INT, "results/integration/dream_results_ashr.csv"),
               select = c("gene", "symbol"))
dream <- unique(dream[symbol != "" & !is.na(symbol)])
symbol_to_ensembl <- setNames(dream$gene, dream$symbol)
cat(sprintf("  Mapping: %d symbols -> Ensembl IDs\n", length(symbol_to_ensembl)))

# Also create Ensembl-to-symbol for labeling
# Strip version suffix for more flexible matching
ensembl_base <- sub("\\.[0-9]+$", "", rownames(expr))
names(ensembl_base) <- rownames(expr)

# ---------------------------------------------------------------------------
# 2. Build transition gene sets from scRNA pseudotime gene dynamics
# ---------------------------------------------------------------------------
cat("Building transition gene sets from pseudotime gene dynamics...\n")

CELL_TYPES <- c("Hepatocytes", "Macrophages", "Fibroblasts",
                "Endothelial_cells", "Cholangiocytes")

# For each cell type, gene_dynamics has 2000 genes x 100 pseudotime bins.
# Early bins (0-30) = early disease; late bins (70-99) = advanced disease.
# Transition genes: early UP = high in bins 0-30 & low in 70-99 (and vice versa).
N_TOP <- 100  # top genes per transition direction

gene_sets <- list()

for (ct in CELL_TYPES) {
  dyn_path <- file.path(PSEUDOTIME_DIR, sprintf("gene_dynamics_%s.csv", ct))
  if (!file.exists(dyn_path)) {
    cat(sprintf("  WARNING: %s not found, skipping\n", dyn_path))
    next
  }
  dyn <- fread(dyn_path)
  genes <- dyn$gene
  dyn_mat <- as.matrix(dyn[, -1, with = FALSE])  # 2000 x 100 bins

  # Early program: mean of bins 0-30
  early_score <- rowMeans(dyn_mat[, 1:31], na.rm = TRUE)
  # Late program: mean of bins 70-100
  late_score <- rowMeans(dyn_mat[, 71:100], na.rm = TRUE)
  # Transition delta
  delta <- late_score - early_score

  # Late-UP genes (increase with pseudotime = progression-associated)
  top_up <- genes[order(delta, decreasing = TRUE)[1:N_TOP]]
  # Late-DOWN genes (decrease with pseudotime = early/protective)
  top_down <- genes[order(delta, decreasing = FALSE)[1:N_TOP]]

  ct_short <- gsub("_cells$", "", ct)
  gene_sets[[paste0(ct_short, "_progression_UP")]] <- top_up
  gene_sets[[paste0(ct_short, "_progression_DOWN")]] <- top_down

  cat(sprintf("  %s: %d UP, %d DOWN transition genes\n",
    ct_short, length(top_up), length(top_down)))
}

cat(sprintf("Total gene sets: %d\n", length(gene_sets)))

# ---------------------------------------------------------------------------
# 3. Score bulk samples via ssGSEA
# ---------------------------------------------------------------------------
cat("Running ssGSEA on bulk expression...\n")

# Map gene symbols to Ensembl IDs and filter to those in bulk expression
gene_sets_filtered <- lapply(gene_sets, function(gs) {
  # Map symbols to Ensembl IDs
  ensembl_ids <- symbol_to_ensembl[gs]
  ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]
  # Keep only those present in bulk expression matrix
  intersect(ensembl_ids, rownames(expr))
})
keep <- sapply(gene_sets_filtered, length) >= 10
gene_sets_filtered <- gene_sets_filtered[keep]
cat(sprintf("  Gene sets with >=10 genes in bulk: %d / %d\n",
  sum(keep), length(gene_sets)))

# Run ssGSEA
gsva_param <- ssgseaParam(expr, gene_sets_filtered, normalize = TRUE)
scores <- gsva(gsva_param, verbose = FALSE)

# Transpose to samples x gene_sets
scores_df <- as.data.table(t(scores), keep.rownames = "sample_id")
cat(sprintf("  Transition signature scores: %d samples x %d signatures\n",
  nrow(scores_df), ncol(scores_df) - 1))

# ---------------------------------------------------------------------------
# 4. Save
# ---------------------------------------------------------------------------
out_path <- file.path(OUTDIR, "pseudotime_transition_scores.csv")
fwrite(scores_df, out_path)
cat(sprintf("Written to: %s\n", out_path))

# Also save the gene sets for reference
gene_set_df <- rbindlist(lapply(names(gene_sets_filtered), function(nm) {
  data.table(gene_set = nm, gene = gene_sets_filtered[[nm]])
}))
fwrite(gene_set_df, file.path(OUTDIR, "pseudotime_transition_gene_sets.csv"))
cat(sprintf("Gene sets written to: %s\n",
  file.path(OUTDIR, "pseudotime_transition_gene_sets.csv")))

cat("Done.\n")
