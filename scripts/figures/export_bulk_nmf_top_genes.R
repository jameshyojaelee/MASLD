#!/usr/bin/env Rscript
# Export top-50 genes per bulk NMF k=6 program from the cached NMF fit.
# Output: figures/misc/nmf_umap_exploration/nmf_k6_top50_genes.tsv
# Columns: program, rank, ensembl_or_symbol, symbol
#
# Idempotent: re-running overwrites the TSV but is fast (~30s).

suppressPackageStartupMessages({
  library(data.table)
  library(NMF)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

CACHE  <- file.path(BASE, "RNA-seq/results/subtypes/nmf_results_cache_clean.rds")
META   <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
OUTDIR <- file.path(BASE, "figures/misc/nmf_umap_exploration")
OUT    <- file.path(OUTDIR, "nmf_k6_top50_genes.tsv")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

K_TARGET <- 6
TOPN     <- 50

cat(sprintf("Reading cache: %s\n", CACHE))
cache <- readRDS(CACHE)

# Cache layout: list(nmf_results=list("3"=fit, "4"=fit, ...), metrics=..., mat_nn=..., mat=...)
fits <- cache$nmf_results
if (!as.character(K_TARGET) %in% names(fits)) {
  stop(sprintf("k=%d not in cache$nmf_results; available: %s",
               K_TARGET, paste(names(fits), collapse=",")))
}
fit <- fits[[as.character(K_TARGET)]]
W <- basis(fit)  # gene x program
cat(sprintf("W shape: %d genes x %d programs\n", nrow(W), ncol(W)))

# Gene IDs come from row names; could be ENSEMBL or symbols. Detect.
gene_ids <- rownames(W)
is_ensembl <- mean(grepl("^ENSG", gene_ids)) > 0.5
cat(sprintf("Gene-id format: %s (frac ENSG = %.2f)\n",
            ifelse(is_ensembl, "ENSEMBL", "symbol"),
            mean(grepl("^ENSG", gene_ids))))

# Build symbol mapping if needed.
sym_lookup <- NULL
if (file.exists(META)) {
  meta <- fread(META)
  # gencode metadata likely has columns: gene_id, gene_name (symbol)
  cols <- intersect(c("gene_id", "ensembl_id", "ensembl_gene_id"), names(meta))
  scol <- intersect(c("gene_name", "symbol", "gene_symbol"), names(meta))
  if (length(cols) > 0 && length(scol) > 0) {
    sym_lookup <- setNames(meta[[scol[1]]], sub("\\..*$","", meta[[cols[1]]]))
  }
}

# Build top-N table.
out <- rbindlist(lapply(seq_len(ncol(W)), function(p) {
  loadings <- W[, p]
  ord <- order(loadings, decreasing = TRUE)[seq_len(min(TOPN, nrow(W)))]
  ids <- gene_ids[ord]
  ids_clean <- sub("\\..*$","", ids)  # strip ENSEMBL version suffix if present

  if (is_ensembl && !is.null(sym_lookup)) {
    syms <- sym_lookup[ids_clean]
    syms[is.na(syms) | syms == ""] <- ids_clean[is.na(syms) | syms == ""]
  } else {
    syms <- ids_clean
  }

  data.table(
    program = sprintf("P%d", p),
    rank = seq_along(ord),
    gene_id = ids_clean,
    symbol = syms,
    loading = loadings[ord]
  )
}))

cat(sprintf("Writing %s (%d rows)\n", OUT, nrow(out)))
fwrite(out, OUT, sep = "\t")

# Quick sanity: print top-5 of each program
cat("\n--- Top-5 per program ---\n")
for (p in unique(out$program)) {
  syms <- out[program == p][1:5, symbol]
  cat(sprintf("%s: %s\n", p, paste(syms, collapse=", ")))
}
