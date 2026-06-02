#!/usr/bin/env Rscript
# 09_build_kallisto_dge_9cohort.R
# ---------------------------------------------------------------------------
# Build merged_dge.rds and merged_counts_raw.rds from the 9-cohort kallisto
# counts (wave1 5 mega + wave2 4 staging cohorts; PRJNA512027 excluded).
#
# Mirrors 03_integrate_counts.R logic:
#   1. Load kallisto gene counts (all_cohorts_gene_counts.tsv.gz)
#   2. Match to canonical metadata (meta_matched.rds)
#   3. Apply pass_technical QC filter
#   4. Build DGEList with filterByExpr + RLE normalization
#
# Outputs (overwrites existing with kallisto versions):
#   results/integration/merged_counts_raw.rds  (all 9 cohorts, pre-QC)
#   results/integration/merged_dge.rds         (QC-filtered + normalized)
#   results/integration/merged_dge_star_backup.rds  (backup of STAR version)
#   results/integration/merged_counts_raw_star_backup.rds
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

PROJECT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
BASE <- file.path(PROJECT, "RNA-seq/Human/Patient_Cohorts")
RDIR <- file.path(BASE, "analysis/integration/results/integration")
KALL <- file.path(PROJECT, "RNA-seq/results/kallisto")
INT  <- file.path(BASE, "analysis/integration")

# ── Backup existing STAR files ──────────────────────────────────────────────
for (f in c("merged_dge.rds", "merged_counts_raw.rds")) {
  src <- file.path(RDIR, f)
  bak <- file.path(RDIR, sub("\\.rds$", "_star_backup.rds", f))
  if (file.exists(src) && !file.exists(bak)) {
    file.copy(src, bak)
    cat("Backed up", f, "->", basename(bak), "\n")
  } else if (file.exists(bak)) {
    cat("Backup already exists:", basename(bak), "\n")
  }
}

# ── Load kallisto counts ────────────────────────────────────────────────────
cat("\n[1/5] Loading kallisto counts...\n")
kc <- fread(file.path(KALL, "all_cohorts_gene_counts.tsv.gz"))
gene_ids_versioned <- kc$gene_id
kc[, gene_id := NULL]
cnts <- as.matrix(kc)
rownames(cnts) <- gene_ids_versioned
storage.mode(cnts) <- "double"
cnts[is.na(cnts)] <- 0
cnts <- round(cnts)
cat("  Raw kallisto matrix:", nrow(cnts), "genes x", ncol(cnts), "samples\n")

# Strip ENSG version suffix for downstream compatibility
rownames(cnts) <- sub("\\..*$", "", rownames(cnts))
# Handle duplicates (rare: different versions of same gene -> sum)
if (any(duplicated(rownames(cnts)))) {
  dups <- unique(rownames(cnts)[duplicated(rownames(cnts))])
  cat("  Summing", length(dups), "duplicated gene IDs (multi-version)\n")
  cnts_dt <- as.data.table(cnts, keep.rownames = "gene_id")
  cnts_dt <- cnts_dt[, lapply(.SD, sum), by = gene_id]
  gene_ids_clean <- cnts_dt$gene_id
  cnts_dt[, gene_id := NULL]
  cnts <- as.matrix(cnts_dt)
  rownames(cnts) <- gene_ids_clean
}

# ── Load canonical metadata ─────────────────────────────────────────────────
cat("\n[2/5] Loading metadata...\n")
meta <- readRDS(file.path(RDIR, "meta_matched.rds"))
# Exclude PRJNA512027 (dropped from paper 2026-05-15)
meta <- meta[meta$dataset != "PRJNA512027", ]
cat("  Metadata:", nrow(meta), "samples across",
    length(unique(meta$dataset)), "datasets\n")

# Match kallisto samples to metadata (inner join)
common_ids <- intersect(colnames(cnts), meta$sample_id)
cat("  Kallisto samples:", ncol(cnts),
    "| Metadata samples:", nrow(meta),
    "| Common:", length(common_ids), "\n")

if (length(common_ids) == 0) {
  stop("No samples in common between kallisto counts and metadata!")
}

cnts <- cnts[, common_ids]
meta <- meta[match(common_ids, meta$sample_id), ]
stopifnot(all(meta$sample_id == colnames(cnts)))

cat("  Per-dataset sample counts (pre-QC):\n")
print(table(meta$dataset))

# ── Save merged_counts_raw.rds ──────────────────────────────────────────────
cat("\n[3/5] Saving merged_counts_raw.rds...\n")
saveRDS(cnts, file.path(RDIR, "merged_counts_raw.rds"))
cat("  Saved:", nrow(cnts), "genes x", ncol(cnts), "samples\n")

# ── Apply QC filter (pass_technical) ────────────────────────────────────────
cat("\n[4/5] Applying QC filter...\n")
qc_file <- file.path(INT, "qc/sample_qc_report.csv")
if (!file.exists(qc_file)) {
  stop("QC report not found: ", qc_file)
}
qc <- fread(qc_file)
pass_ids <- qc[pass_technical == TRUE, sample_id]
pass_ids <- intersect(pass_ids, common_ids)
cat("  QC-passing samples:", length(pass_ids), "/", length(common_ids), "\n")

cnts_qc <- cnts[, pass_ids]
meta_qc <- meta[match(pass_ids, meta$sample_id), ]
stopifnot(all(meta_qc$sample_id == colnames(cnts_qc)))

cat("  Per-dataset sample counts (post-QC):\n")
print(table(meta_qc$dataset))

# ── Build DGEList ───────────────────────────────────────────────────────────
cat("\n[5/5] Building DGEList...\n")

# For disease-only cohorts (GSE167523, GSE174478, GSE193066, GSE240729),
# group_binary may be all "Disease". Set a safe default.
grp <- meta_qc$group_binary
if (is.null(grp) || all(is.na(grp))) grp <- rep("Disease", nrow(meta_qc))
grp[is.na(grp)] <- "Disease"  # disease-only cohorts

dge <- DGEList(counts = cnts_qc)
dge$samples$dataset      <- meta_qc$dataset
dge$samples$condition    <- factor(meta_qc$condition)
dge$samples$group_binary <- factor(grp, levels = c("Control", "Disease"))
dge$samples$sex          <- meta_qc$sex

# filterByExpr
design_filter <- model.matrix(~ 0 + group_binary, data = dge$samples)
keep <- filterByExpr(dge, design = design_filter)
cat("  Genes passing filterByExpr:", sum(keep), "/", nrow(dge), "\n")
dge <- dge[keep, , keep.lib.sizes = FALSE]

# RLE normalization (same as 03_integrate_counts.R)
dge <- calcNormFactors(dge, method = "RLE")

cat("\n===== MERGED DGE SUMMARY (KALLISTO, 9 COHORTS) =====\n")
cat("Samples:", ncol(dge), "\n")
cat("Genes:", nrow(dge), "\n")
cat("Datasets:", paste(sort(unique(dge$samples$dataset)), collapse = ", "), "\n")
cat("Group distribution:\n")
print(table(dge$samples$group_binary, dge$samples$dataset))

# ── Save ────────────────────────────────────────────────────────────────────
saveRDS(dge, file.path(RDIR, "merged_dge.rds"))
cat("\nSaved:", file.path(RDIR, "merged_dge.rds"), "\n")
cat("  genes:", nrow(dge), " samples:", ncol(dge), "\n")

# ── Validate against dream results ─────────────────────────────────────────
dream_f <- file.path(RDIR, "dream_results.csv")
if (file.exists(dream_f)) {
  dream <- fread(dream_f, select = "gene")
  overlap <- length(intersect(rownames(dge), dream$gene))
  dream_only <- length(setdiff(dream$gene, rownames(dge)))
  cat(sprintf("\nValidation vs dream_results.csv:\n"))
  cat(sprintf("  DGE genes: %d, dream genes: %d, overlap: %d, dream-only: %d\n",
              nrow(dge), nrow(dream), overlap, dream_only))
}

cat("\nDone.\n")
