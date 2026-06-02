#!/usr/bin/env Rscript
# build_kallisto_dge.R
# ---------------------------------------------------------------------------
# Build merged_dge.rds from kallisto gene counts so that pillar worker scripts
# (CPSS, bootstrap, k-fold, permutation, PVCA, SVA) produce results with
# unversioned gene IDs matching the canonical dream_results.csv.
#
# Mirrors the DGEList construction in 04b_dream_kallisto_enhanced.R but
# includes ALL cohorts (not just mega) and saves the DGEList as merged_dge.rds.
# Backs up the existing STAR merged_dge.rds first.
#
# Output: results/integration/merged_dge.rds (overwrite with kallisto version)
#         results/integration/merged_dge_star_backup.rds (STAR backup)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(yaml)
})

PROJECT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
BASE <- file.path(PROJECT, "RNA-seq/Human/Patient_Cohorts")
RDIR <- file.path(BASE, "analysis/integration/results/integration")
KALL <- file.path(PROJECT, "RNA-seq/results/kallisto")

# --- Backup existing STAR merged_dge.rds ---
existing <- file.path(RDIR, "merged_dge.rds")
backup   <- file.path(RDIR, "merged_dge_star_backup.rds")
if (file.exists(existing) && !file.exists(backup)) {
  file.copy(existing, backup)
  cat("Backed up STAR merged_dge.rds ->", backup, "\n")
} else if (file.exists(backup)) {
  cat("STAR backup already exists:", backup, "\n")
}

# --- Load kallisto counts ---
cat("Loading kallisto counts...\n")
kc <- fread(file.path(KALL, "all_cohorts_gene_counts.tsv.gz"))
gene_ids <- kc$gene_id
kc[, gene_id := NULL]
cnts <- as.matrix(kc)
rownames(cnts) <- gene_ids
storage.mode(cnts) <- "double"
cnts[is.na(cnts)] <- 0
cnts <- round(cnts)
cat("  Raw kallisto matrix:", nrow(cnts), "genes x", ncol(cnts), "samples\n")

# Strip ENSG version suffix to match canonical dream_results.csv gene IDs
rownames(cnts) <- sub("\\..*$", "", rownames(cnts))

# --- Load canonical sample metadata from existing DGE ---
cat("Loading canonical metadata...\n")
star_dge <- readRDS(backup)
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))

# QC filter: the star_dge already has pass_technical applied (via 03_integrate_counts.R)
# Use ALL cohorts (not just mega) so pillar B and other scripts can filter as needed
ds_samples <- star_dge$samples
ds_samples$sample_id <- rownames(ds_samples)
common_ids <- intersect(ds_samples$sample_id, colnames(cnts))
cat("  STAR samples:", nrow(ds_samples),
    " | kallisto samples:", ncol(cnts),
    " | common:", length(common_ids), "\n")

cnts <- cnts[, common_ids]
ds_samples <- ds_samples[match(common_ids, ds_samples$sample_id), ]

# --- Build DGEList ---
group <- factor(ds_samples$group_binary, levels = c("Control", "Disease"))
y <- DGEList(counts = cnts, samples = data.frame(
  sample_id    = common_ids,
  dataset      = ds_samples$dataset,
  group_binary = ds_samples$group_binary,
  condition    = ds_samples$condition,
  sex          = ds_samples$sex,
  stringsAsFactors = FALSE))
rownames(y$samples) <- common_ids

# Apply filterByExpr (same design as 03_integrate_counts.R)
design_filter <- model.matrix(~ 0 + group_binary, data = y$samples)
keep <- filterByExpr(y, design = design_filter)
cat("  Genes passing filterByExpr:", sum(keep), "/", nrow(y), "\n")
y <- y[keep, , keep.lib.sizes = FALSE]

# TMM normalization
y <- calcNormFactors(y, method = "TMM")
cat("  Final DGE:", nrow(y), "genes x", ncol(y), "samples\n")

# --- Verify gene overlap with dream_results.csv ---
dream_f <- file.path(RDIR, "dream_results.csv")
if (file.exists(dream_f)) {
  dream <- fread(dream_f, select = "gene")
  overlap <- length(intersect(rownames(y), dream$gene))
  dream_only <- length(setdiff(dream$gene, rownames(y)))
  dge_only <- length(setdiff(rownames(y), dream$gene))
  cat(sprintf("\n  Validation vs dream_results.csv:\n"))
  cat(sprintf("    DGE genes: %d, dream genes: %d\n", nrow(y), nrow(dream)))
  cat(sprintf("    Overlap: %d, dream-only: %d, DGE-only: %d\n",
              overlap, dream_only, dge_only))
  # The DGE has all cohorts; dream was run on mega-only. DGE should be a superset.
  if (dream_only > 0) {
    cat("  WARNING: dream_results.csv has", dream_only,
        "genes NOT in the new DGE (may be from mega-only filterByExpr)\n")
  }
}

# --- Save ---
saveRDS(y, existing)
cat("\nSaved kallisto merged_dge.rds:", existing, "\n")
cat("  genes:", nrow(y), " samples:", ncol(y), "\n")
cat("Done.\n")
