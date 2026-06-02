#!/usr/bin/env Rscript
# 220_process_pxd052937.R
# Process DIA-MS plasma proteomics data from Spectronaut output (PXD052937).
#
# Input:  data/PXD052937/.../Souria 72 Plasma DirectDIA_Report_Protein Quant (Pivot).xls
#         data/PXD052937/.../Souria 72 Plasma DirectDIA_ConditionSetup.tsv
# Output: results/multiprogram/ms_plasma_matrix.csv   (72 samples x N proteins)
#         results/multiprogram/ms_plasma_metadata.csv  (sample_id, condition, run_label)
#
# Steps:
#   1. Read protein quant pivot table (TSV despite .xls extension)
#   2. Parse sample metadata from condition setup
#   3. Log2-transform quantities (pseudocount for zeros/NAs)
#   4. Map UniProt IDs to HGNC gene symbols via biomaRt
#   5. Median-impute remaining missing values per protein
#   6. Write outputs
#
# SLURM: --partition=cpu --cpus-per-task=4 --mem=32G --time=48:00:00
# Env: micromamba activate rnaseq

suppressPackageStartupMessages({
  library(data.table)
  library(biomaRt)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

DATA_DIR <- file.path(BASE,
  "data/PXD052937/20220808_140452_Souria 72 Plasma DirectDIA")
QUANT_FILE <- file.path(DATA_DIR,
  "Souria 72 Plasma DirectDIA_Report_Protein Quant (Pivot).xls")
COND_FILE <- file.path(DATA_DIR,
  "Souria 72 Plasma DirectDIA_ConditionSetup.tsv")

RDIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram")
dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 220: Process PXD052937 DIA-MS Plasma Proteomics ===\n")

# ── 1. Read protein quant pivot table ─────────────────────────────────────
cat("Reading protein quant pivot table...\n")
quant_raw <- fread(QUANT_FILE, sep = "\t", header = TRUE,
                   na.strings = c("", "NA", "NaN", "Filtered"))

cat(sprintf("  Dimensions: %d proteins x %d columns\n",
            nrow(quant_raw), ncol(quant_raw)))

# First column is PG.ProteinGroups; remaining columns are per-sample quantities
protein_col <- names(quant_raw)[1]
sample_cols  <- names(quant_raw)[-1]

# Extract run label from column header: "[N] RUNLABEL.PG.Quantity" -> RUNLABEL
# e.g. "[1] FL07202022_160_P116.htrms.PG.Quantity" -> "FL07202022_160_P116.htrms"
parse_run_label <- function(col_name) {
  # Strip leading "[N] " prefix and trailing ".PG.Quantity"
  s <- sub("^\\[\\d+\\]\\s+", "", col_name)
  sub("\\.PG\\.Quantity$", "", s)
}

run_labels <- vapply(sample_cols, parse_run_label, character(1))
cat(sprintf("  %d sample columns parsed\n", length(run_labels)))

# Extract UniProt IDs (first accession in semicolon-delimited group)
uniprot_ids <- sub(";.*", "", quant_raw[[protein_col]])

# Build numeric matrix: proteins x samples
quant_mat <- as.matrix(quant_raw[, .SD, .SDcols = sample_cols])
rownames(quant_mat) <- uniprot_ids
colnames(quant_mat) <- run_labels
storage.mode(quant_mat) <- "numeric"

cat(sprintf("  Non-missing values: %d / %d (%.1f%%)\n",
            sum(!is.na(quant_mat)), length(quant_mat),
            100 * mean(!is.na(quant_mat))))

# ── 2. Parse sample metadata ───────────────────────────────────────────────
cat("Reading condition setup...\n")
cond_raw <- fread(COND_FILE, sep = "\t", header = TRUE,
                  skip = 0, fill = TRUE)

# Column names: "#", "Is Reference", "Run Label", "Condition", ...
# Guard against comment-prefixed header
if (names(cond_raw)[1] == "#") {
  setnames(cond_raw, "#", "row_num")
} else {
  # fread may have read the # as part of the first col name
  setnames(cond_raw, names(cond_raw)[1], "row_num")
}

# Ensure expected columns exist
needed <- c("Run Label", "Condition")
missing_cols <- setdiff(needed, names(cond_raw))
if (length(missing_cols) > 0) {
  stop("ConditionSetup missing columns: ", paste(missing_cols, collapse = ", "))
}

metadata <- data.table(
  run_label = cond_raw[["Run Label"]],
  condition = cond_raw[["Condition"]]
)
metadata <- metadata[!is.na(run_label) & run_label != "Run Label"]

cat(sprintf("  %d samples in condition setup\n", nrow(metadata)))
cat("  Condition distribution:\n")
print(metadata[, .N, by = condition][order(condition)])

# ── 3. Log2-transform (pseudocount = min non-zero / 2) ────────────────────
cat("Log2-transforming quantities...\n")
min_nonzero <- min(quant_mat[quant_mat > 0], na.rm = TRUE)
pseudocount <- min_nonzero / 2
cat(sprintf("  Pseudocount: %.4f (half of minimum non-zero value)\n", pseudocount))

# Treat zeros as missing before log transform
quant_mat[quant_mat == 0] <- NA
log2_mat <- log2(quant_mat + pseudocount)

# ── 4. Map UniProt IDs to HGNC gene symbols via biomaRt ───────────────────
cat("Mapping UniProt IDs to gene symbols via biomaRt...\n")

# Use only the primary (first) accession per protein group
primary_uniprot <- rownames(log2_mat)
# Remove isoform suffixes (e.g., P12345-2 -> P12345) for biomaRt lookup
primary_uniprot_clean <- sub("-\\d+$", "", primary_uniprot)

symbol_map <- tryCatch({
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                     host = "https://www.ensembl.org")
  bm <- getBM(
    attributes = c("uniprotswissprot", "hgnc_symbol"),
    filters    = "uniprotswissprot",
    values     = unique(primary_uniprot_clean),
    mart       = ensembl
  )
  setDT(bm)
  # Keep only rows with a valid gene symbol
  bm <- bm[hgnc_symbol != "" & !is.na(hgnc_symbol)]
  # One symbol per UniProt: prefer shorter (canonical) symbol on ties
  bm <- bm[order(nchar(hgnc_symbol))][!duplicated(uniprotswissprot)]
  cat(sprintf("  biomaRt: mapped %d / %d UniProt IDs\n",
              nrow(bm), length(unique(primary_uniprot_clean))))
  bm
}, error = function(e) {
  cat(sprintf("  WARNING: biomaRt unavailable (%s). Falling back to UniProt IDs as symbols.\n",
              conditionMessage(e)))
  data.table(uniprotswissprot = character(0), hgnc_symbol = character(0))
})

# Build a lookup vector: UniProt -> symbol
if (nrow(symbol_map) > 0) {
  lookup <- setNames(symbol_map$hgnc_symbol, symbol_map$uniprotswissprot)
  gene_symbols <- lookup[primary_uniprot_clean]
  # Fall back to UniProt ID where mapping failed
  unmapped <- is.na(gene_symbols)
  gene_symbols[unmapped] <- primary_uniprot_clean[unmapped]
  names(gene_symbols) <- NULL
} else {
  gene_symbols <- primary_uniprot_clean
}

cat(sprintf("  %d proteins with gene symbols, %d using UniProt fallback\n",
            sum(gene_symbols != primary_uniprot_clean),
            sum(gene_symbols == primary_uniprot_clean)))

# Handle duplicate gene symbols: keep the protein with fewest NAs per symbol
log2_dt <- as.data.table(log2_mat)
log2_dt[, gene_symbol := gene_symbols]
log2_dt[, n_missing := rowSums(is.na(log2_mat))]
log2_dt <- log2_dt[order(n_missing)][!duplicated(gene_symbol)]
log2_dt[, n_missing := NULL]
cat(sprintf("  %d unique gene symbols retained (from %d proteins)\n",
            nrow(log2_dt), nrow(log2_mat)))

# ── 5. Median imputation per protein ──────────────────────────────────────
cat("Applying per-protein median imputation...\n")
n_before <- sum(is.na(as.matrix(log2_dt[, .SD, .SDcols = run_labels])))

sample_cols_present <- intersect(run_labels, names(log2_dt))
for (i in seq_len(nrow(log2_dt))) {
  row_vals <- as.numeric(log2_dt[i, ..sample_cols_present])
  med_val  <- median(row_vals, na.rm = TRUE)
  if (!is.na(med_val)) {
    for (col in sample_cols_present) {
      if (is.na(log2_dt[[col]][i])) set(log2_dt, i, col, med_val)
    }
  }
}

n_after <- sum(is.na(as.matrix(log2_dt[, .SD, .SDcols = sample_cols_present])))
cat(sprintf("  Imputed %d missing values (%.1f%% of matrix)\n",
            n_before - n_after,
            100 * (n_before - n_after) / (nrow(log2_dt) * length(sample_cols_present))))

# ── 6. Build output: samples x proteins matrix ────────────────────────────
cat("Building output matrices...\n")

# Transpose: samples (rows) x proteins (columns)
protein_matrix <- t(as.matrix(log2_dt[, .SD, .SDcols = sample_cols_present]))
colnames(protein_matrix) <- log2_dt$gene_symbol
rownames(protein_matrix) <- sample_cols_present

# Convert to data.table with sample_id column
out_mat <- as.data.table(protein_matrix, keep.rownames = "sample_id")

cat(sprintf("  Output matrix: %d samples x %d proteins\n",
            nrow(out_mat), ncol(out_mat) - 1L))

# Build metadata output: align to actual sample columns
# run_label in metadata uses the raw run label from condition file; match by stripping
# the ".htrms" suffix variation to a common key
normalize_run <- function(x) sub("\\.htrms$", "", x)

metadata[, run_key := normalize_run(run_label)]
sample_key <- data.table(
  sample_id = sample_cols_present,
  run_key   = normalize_run(sample_cols_present)
)

out_meta <- merge(sample_key, metadata[, .(run_key, condition, run_label)],
                  by = "run_key", all.x = TRUE)
out_meta[, run_key := NULL]
setcolorder(out_meta, c("sample_id", "condition", "run_label"))
setorder(out_meta, sample_id)

cat(sprintf("  Metadata: %d samples mapped to conditions (%d unmapped)\n",
            sum(!is.na(out_meta$condition)),
            sum(is.na(out_meta$condition))))

# ── 7. Write outputs ───────────────────────────────────────────────────────
mat_path  <- file.path(RDIR, "ms_plasma_matrix.csv")
meta_path <- file.path(RDIR, "ms_plasma_metadata.csv")

fwrite(out_mat,  mat_path)
fwrite(out_meta, meta_path)

cat(sprintf("Wrote: %s\n", mat_path))
cat(sprintf("Wrote: %s\n", meta_path))
cat("=== 220: DONE ===\n")
