#!/usr/bin/env Rscript
# 17_annotate_gene_symbols.R
# ---------------------------------------------------------------------------
# 1. Parse GENCODE v49 GTF → ENSG ↔ gene_name ↔ gene_type mapping
# 2. Download & cache human↔mouse 1:1 ortholog table via biomaRt
# 3. Add gene_symbol, gene_type, mouse columns to all disease signature CSVs
# 4. Re-run sanity checks from script 16 using symbols
# 5. Re-compute library overlap with mapped orthologs
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(biomaRt)
  library(ggplot2)
})

# ============================================================
#  Paths
# ============================================================
BASE   <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT    <- file.path(BASE, "analysis/integration")
DS_DIR <- file.path(INT, "results/disease_signatures")
INT_DIR <- file.path(INT, "results/integration")
MOUSE_UI <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration"

HUMAN_GTF <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
MOUSE_GTF <- file.path(MOUSE_UI, "../Public_Diet_Models/reference/raw/gencode.vM38.annotation.gtf")

CACHE_DIR <- file.path(INT, "results/gene_annotation")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# Helper: parse GENCODE GTF using awk (no rtracklayer needed)
parse_gtf_genes <- function(gtf_path, cache_path) {
  if (file.exists(cache_path)) {
    cat("  Loading cached:", basename(cache_path), "\n")
    return(fread(cache_path))
  }

  cat("  Parsing GTF with awk:", basename(gtf_path), "\n")
  # awk extracts gene_id, gene_name, gene_type from gene-level records
  is_gz <- grepl("\\.gz$", gtf_path)
  reader <- if (is_gz) paste("zcat", shQuote(gtf_path)) else paste("cat", shQuote(gtf_path))

  cmd <- paste0(reader, " | awk -F'\\t' '$3 == \"gene\" {",
    "match($9, /gene_id \"([^\"]+)\"/, gid); ",
    "match($9, /gene_name \"([^\"]+)\"/, gn); ",
    "match($9, /gene_type \"([^\"]+)\"/, gt); ",
    "print gid[1] \"\\t\" gn[1] \"\\t\" gt[1]",
    "}'")

  raw <- fread(cmd = cmd, header = FALSE, col.names = c("gene_id", "symbol", "gene_type"))
  raw[, gene_base := gsub("\\..*", "", gene_id)]
  raw <- unique(raw)
  fwrite(raw, cache_path, sep = "\t")
  cat("  Saved:", cache_path, "(", nrow(raw), "genes)\n")
  return(raw)
}

# ============================================================
#  Step 1: Parse GENCODE GTFs for gene ID → symbol
# ============================================================
cat("=== STEP 1: Parsing GENCODE GTFs ===\n")

h_annot <- parse_gtf_genes(HUMAN_GTF, file.path(CACHE_DIR, "human_ensg_to_symbol.tsv"))
cat("  Human genes in GTF:", nrow(h_annot), "\n")

cat("\n")
m_annot <- parse_gtf_genes(MOUSE_GTF, file.path(CACHE_DIR, "mouse_ensmusg_to_symbol.tsv"))
# Rename columns for mouse
if (!"mouse_gene_id" %in% names(m_annot)) {
  setnames(m_annot, c("gene_id", "symbol", "gene_type", "gene_base"),
           c("mouse_gene_id", "mouse_symbol", "mouse_gene_type", "mouse_gene_base"))
}
cat("  Mouse genes in GTF:", nrow(m_annot), "\n")

# ============================================================
#  Step 2: Download & cache human↔mouse ortholog mapping
# ============================================================
cat("\n=== STEP 2: Ortholog mapping ===\n")

cache_ortho <- file.path(CACHE_DIR, "ortholog_mapping.tsv")
if (file.exists(cache_ortho)) {
  cat("  Loading cached ortholog mapping...\n")
  orthologs <- fread(cache_ortho)
} else {
  cat("  Querying Ensembl for human↔mouse orthologs...\n")
  human_mart <- tryCatch({
    useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
               mirror = "useast")
  }, error = function(e) {
    cat("  useast mirror failed, trying www...\n")
    useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  })

  ortho_raw <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name",
                    "mmusculus_homolog_ensembl_gene",
                    "mmusculus_homolog_associated_gene_name",
                    "mmusculus_homolog_orthology_type",
                    "mmusculus_homolog_perc_id",
                    "mmusculus_homolog_perc_id_r1"),
    mart = human_mart
  )
  orthologs <- as.data.table(ortho_raw)
  setnames(orthologs, c("human_gene_id", "human_symbol",
                         "mouse_gene_id", "mouse_symbol",
                         "orthology_type", "perc_id_human_to_mouse",
                         "perc_id_mouse_to_human"))

  # Keep only one2one with valid mouse IDs
  orthologs <- orthologs[mouse_gene_id != "" & orthology_type == "ortholog_one2one"]
  orthologs <- unique(orthologs)

  fwrite(orthologs, cache_ortho, sep = "\t")
  cat("  Saved:", cache_ortho, "\n")
}
cat("  One-to-one orthologs:", nrow(orthologs), "\n")

# ============================================================
#  Step 3: Annotate all disease signature CSVs
# ============================================================
cat("\n=== STEP 3: Annotating CSVs ===\n")

# Helper: add symbol/type/mouse columns to a CSV (idempotent)
annotate_csv <- function(filepath, gene_col = "gene", save = TRUE) {
  if (!file.exists(filepath)) {
    cat("  SKIP (not found):", basename(filepath), "\n")
    return(NULL)
  }
  dt <- fread(filepath)
  if (!gene_col %in% names(dt)) {
    cat("  SKIP (no gene col):", basename(filepath), "\n")
    return(dt)
  }

  # Remove existing annotation columns (idempotent re-run)
  drop <- intersect(names(dt), c("symbol", "gene_type", "mouse_gene_id",
                                   "mouse_symbol", "gene_base"))
  if (length(drop) > 0) dt[, (drop) := NULL]

  # Strip version for matching
  dt[, gene_base := gsub("\\..*", "", get(gene_col))]

  # Merge human annotation
  dt <- merge(dt, h_annot[, .(gene_base, symbol, gene_type)],
              by = "gene_base", all.x = TRUE)

  # Merge ortholog info
  dt <- merge(dt, orthologs[, .(human_gene_id, mouse_gene_id, mouse_symbol)],
              by.x = "gene_base", by.y = "human_gene_id", all.x = TRUE)

  # Clean up: move symbol/type to front
  front_cols <- c(gene_col, "symbol", "gene_type", "mouse_gene_id", "mouse_symbol")
  other_cols <- setdiff(names(dt), c(front_cols, "gene_base"))
  setcolorder(dt, c(front_cols, other_cols))
  dt[, gene_base := NULL]

  # Report
  mapped <- sum(!is.na(dt$symbol))
  ortho  <- sum(!is.na(dt$mouse_gene_id))
  cat(sprintf("  %s: %d genes, %d mapped (%.1f%%), %d with mouse ortholog\n",
    basename(filepath), nrow(dt), mapped, mapped/nrow(dt)*100, ortho))

  if (save) fwrite(dt, filepath)
  return(dt)
}

# Disease signatures directory
ds_files <- c(
  "nafl_vs_nash_dream.csv", "nafl_vs_nash_meta.csv", "nafl_vs_nash_consensus.csv",
  "fibrosis_dream.csv", "fibrosis_slopes_meta.csv", "fibrosis_early_late.csv",
  "nas_components_all.csv", "unified_disease_signatures.csv"
)
cat("\nDisease signature files:\n")
for (f in ds_files) {
  annotate_csv(file.path(DS_DIR, f))
}

# Integration directory
int_files <- c("consensus_degs.csv", "dream_results.csv")
cat("\nIntegration files:\n")
for (f in int_files) {
  annotate_csv(file.path(INT_DIR, f))
}

# ============================================================
#  Step 4: Sanity checks with gene symbols
# ============================================================
cat("\n=== STEP 4: Sanity checks with gene symbols ===\n")

# Reload annotated files
nn_dream <- fread(file.path(DS_DIR, "nafl_vs_nash_dream.csv"))
fib_dream <- fread(file.path(DS_DIR, "fibrosis_dream.csv"))
unified <- fread(file.path(DS_DIR, "unified_disease_signatures.csv"))

# Known NASH markers
nash_markers <- c("CYP7A1", "CYP2E1", "PNPLA3", "FASN", "SCD", "SREBF1",
                   "CCL2", "TNF", "IL6", "IL1B", "TGFB1")
cat("\nKnown NASH markers in NAFL-vs-NASH dream:\n")
for (marker in nash_markers) {
  hit <- nn_dream[symbol == marker]
  if (nrow(hit) > 0) {
    cat(sprintf("  %s: LFC=%.3f, padj=%.2e, gene_type=%s\n",
      hit$symbol[1], hit$logFC[1], hit$adj.P.Val[1], hit$gene_type[1]))
  } else {
    cat(sprintf("  %s: NOT FOUND\n", marker))
  }
}

# Known fibrosis markers
fib_markers <- c("COL1A1", "COL3A1", "ACTA2", "TGFB1", "TIMP1", "LOX",
                  "COL1A2", "MMP2", "SERPINE1")
cat("\nKnown fibrosis markers in fibrosis dream:\n")
for (marker in fib_markers) {
  hit <- fib_dream[symbol == marker]
  if (nrow(hit) > 0) {
    cat(sprintf("  %s: slope=%.3f, padj=%.2e\n",
      hit$symbol[1], hit$logFC[1], hit$adj.P.Val[1]))
  } else {
    cat(sprintf("  %s: NOT FOUND\n", marker))
  }
}

# Stage classification by known markers
cat("\nStage classification of known markers:\n")
all_known <- c(nash_markers, fib_markers)
known_in_unified <- unified[symbol %in% all_known]
if (nrow(known_in_unified) > 0) {
  print(known_in_unified[, .(symbol, stage_category, nn_lfc, nn_padj, fib_slope, fib_padj)])
}

# ============================================================
#  Step 5: Library overlap with ortholog mapping
# ============================================================
cat("\n=== STEP 5: Library overlap via ortholog mapping ===\n")

lib_file <- file.path(INT_DIR, "library_intersection.csv")
if (file.exists(lib_file)) {
  lib <- fread(lib_file)
  cat("Library genes:", nrow(lib), "\n")

  # Library uses mouse gene IDs — map to human orthologs
  if ("mouse_gene_id" %in% names(lib)) {
    lib[, mouse_base := gsub("\\..*", "", mouse_gene_id)]
  } else {
    # First column might be the gene ID
    lib[, mouse_base := gsub("\\..*", "", lib[[1]])]
  }

  # Map library mouse genes → human orthologs
  lib_mapped <- merge(lib,
    orthologs[, .(mouse_gene_id, human_gene_id, human_symbol)],
    by.x = "mouse_base", by.y = "mouse_gene_id", all.x = TRUE)

  n_mapped <- sum(!is.na(lib_mapped$human_gene_id))
  cat(sprintf("  Library genes with human ortholog: %d / %d (%.1f%%)\n",
    n_mapped, nrow(lib), n_mapped/nrow(lib)*100))

  # Now check overlap with human disease DEGs
  human_nn_sig <- nn_dream[adj.P.Val < 0.05, gsub("\\..*", "", gene)]
  human_fib_sig <- fib_dream[adj.P.Val < 0.05, gsub("\\..*", "", gene)]

  lib_human_ids <- lib_mapped[!is.na(human_gene_id), unique(human_gene_id)]

  nn_overlap <- sum(lib_human_ids %in% human_nn_sig)
  fib_overlap <- sum(lib_human_ids %in% human_fib_sig)

  cat(sprintf("  Library genes in NAFL-vs-NASH DEGs: %d / %d (%.1f%%)\n",
    nn_overlap, length(lib_human_ids), nn_overlap/length(lib_human_ids)*100))
  cat(sprintf("  Library genes in fibrosis DEGs: %d / %d (%.1f%%)\n",
    fib_overlap, length(lib_human_ids), fib_overlap/length(lib_human_ids)*100))

  # Stage category of library genes
  lib_in_unified <- unified[gsub("\\..*", "", gene) %in% lib_human_ids]
  cat("\nLibrary genes by stage category:\n")
  if (nrow(lib_in_unified) > 0) {
    print(lib_in_unified[, .N, by = stage_category][order(-N)])
  } else {
    cat("  No library genes found in unified table\n")
  }
}

# ============================================================
#  Summary stats
# ============================================================
cat("\n=== ANNOTATION SUMMARY ===\n")
cat("  Human genes in GENCODE v49 GTF:", nrow(h_annot), "\n")
cat("  Mouse genes in GENCODE vM38 GTF:", nrow(m_annot), "\n")
cat("  One-to-one orthologs:", nrow(orthologs), "\n")

# Gene type breakdown for unified disease signatures
cat("\nGene type breakdown (unified disease signatures):\n")
type_counts <- unified[, .N, by = gene_type][order(-N)]
print(head(type_counts, 10))

cat("\n=== Script 17 complete ===\n")
