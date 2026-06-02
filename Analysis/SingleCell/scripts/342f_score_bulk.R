#!/usr/bin/env Rscript
# 342f_score_bulk.R - Score polyploid signatures on bulk RNA-seq across 1,469 samples
# via GSVA / ssGSEA on ComBat-corrected log-CPM.
#
# Inputs:
#   RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/corrected_logcpm.rds
#   RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv
#   data/ploidy_signatures/{richter,katsuda,yin,consensus}_*.tsv
# Outputs:
#   Analysis/SingleCell/results_gpu_v2/ploidy/signature_scores_bulk_persample.csv
#   Analysis/SingleCell/results_gpu_v2/ploidy/signature_scores_bulk_summary.csv

suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
setwd(ROOT)
SIG_DIR <- "data/ploidy_signatures"
OUT_DIR <- "Analysis/SingleCell/results_gpu_v2/ploidy"

cat("[", format(Sys.time()), "] Loading bulk log-CPM\n", sep="")
expr <- readRDS("RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/corrected_logcpm.rds")
cat("  bulk matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("  first 5 gene IDs:", paste(rownames(expr)[1:5], collapse=", "), "\n")
cat("  expr range:", round(min(expr),2), "to", round(max(expr),2), "\n")

meta <- fread("RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
cat("  metadata rows:", nrow(meta), "\n")

# Match samples
common <- intersect(colnames(expr), meta$sample_id)
cat("  matched samples:", length(common), "\n")
expr <- expr[, common]
meta <- meta[match(common, sample_id)]

# Bulk matrix uses Ensembl IDs (versioned); signatures use gene symbols.
# Build Ensembl -> gene_symbol mapping via gencode v49 metadata
cat("\n[", format(Sys.time()), "] Mapping bulk Ensembl IDs to gene symbols\n", sep="")
gm <- fread("data/gencode_v49_gene_metadata.tsv.gz")
ens2sym <- setNames(gm$gene_name, gm$gene_id)
expr_genes_ens <- rownames(expr)
expr_genes_sym <- ens2sym[expr_genes_ens]
miss <- sum(is.na(expr_genes_sym))
cat("  mapped:", sum(!is.na(expr_genes_sym)), "/", length(expr_genes_ens), " (missing", miss, ")\n")
# Collapse to gene symbol by taking max expression per symbol
sym_unique <- unique(expr_genes_sym[!is.na(expr_genes_sym)])
cat("  unique gene symbols:", length(sym_unique), "\n")
# Build collapsed matrix: for each symbol, take row with highest mean expression
sym_groups <- split(seq_along(expr_genes_ens), expr_genes_sym)
expr_sym <- matrix(0, nrow = length(sym_unique), ncol = ncol(expr),
                   dimnames = list(sym_unique, colnames(expr)))
for (s in sym_unique) {
  idx <- sym_groups[[s]]
  if (length(idx) == 1) {
    expr_sym[s, ] <- expr[idx, ]
  } else {
    # take the row with highest mean expression
    means <- rowMeans(expr[idx, , drop = FALSE])
    expr_sym[s, ] <- expr[idx[which.max(means)], ]
  }
}
expr <- expr_sym
cat("  collapsed matrix:", dim(expr), "\n")

# Load signatures (human gene symbols)
load_genes <- function(f, col = "human_symbol") {
  d <- fread(f)
  if (!col %in% names(d)) col <- intersect(c("gene_symbol", "Gene", "human_symbol"), names(d))[1]
  unique(toupper(as.character(d[[col]])))
}
sigs <- list(
  richter_up   = load_genes(file.path(SIG_DIR, "richter2021_polyploid_up.tsv")),
  richter_down = load_genes(file.path(SIG_DIR, "richter2021_polyploid_down.tsv")),
  katsuda_up   = load_genes(file.path(SIG_DIR, "katsuda2019_polyploid_up.tsv")),
  katsuda_down = load_genes(file.path(SIG_DIR, "katsuda2019_polyploid_down.tsv")),
  yin_all      = load_genes(file.path(SIG_DIR, "yin2024_polyploid_candidates.tsv")),
  yin_conserved = load_genes(file.path(SIG_DIR, "yin2024_conserved_polyploid.tsv")),
  consensus_up = load_genes(file.path(SIG_DIR, "consensus_polyploid_up.tsv")),
  consensus_down = load_genes(file.path(SIG_DIR, "consensus_polyploid_down.tsv"))
)
cat("\nSignatures (post-filter to expression matrix gene names):\n")
expr_genes <- rownames(expr)
expr_genes_upper <- toupper(expr_genes)
for (k in names(sigs)) {
  sigs[[k]] <- intersect(sigs[[k]], expr_genes_upper)
  cat("  ", k, ":", length(sigs[[k]]), "\n", sep="")
}

# Convert gene name case to match expr rownames
sym_map <- setNames(expr_genes, toupper(expr_genes))
sigs_native <- lapply(sigs, function(g) sym_map[g][!is.na(sym_map[g])])
print(sapply(sigs_native, length))

# GSVA / ssGSEA
cat("\n[", format(Sys.time()), "] Running ssGSEA via GSVA\n", sep="")
# Need >=2 genes per set; drop small ones
sigs_use <- sigs_native[sapply(sigs_native, length) >= 2]
cat("  signatures to score:", length(sigs_use), "\n")

# Use the GSVA v2 API with ssgseaParam
ssgsea_param <- ssgseaParam(exprData = expr, geneSets = sigs_use, minSize = 1, maxSize = 1000)
scores <- gsva(ssgsea_param, verbose = FALSE)
cat("  scored:", dim(scores), "(sig × samples)\n")

# Composite indices
ssc <- as.data.frame(t(scores))
ssc$sample_id <- rownames(ssc)
ssc$polyploid_index_richter   <- ssc[["richter_up"]] - ssc[["richter_down"]]
ssc$polyploid_index_katsuda   <- ssc[["katsuda_up"]] - ssc[["katsuda_down"]]
ssc$polyploid_index_consensus <- ifelse(
  "consensus_up" %in% colnames(ssc) & "consensus_down" %in% colnames(ssc),
  ssc[["consensus_up"]] - ssc[["consensus_down"]], NA_real_
)

# Merge with metadata
out_dt <- merge(meta, as.data.table(ssc), by = "sample_id")
fwrite(out_dt, file.path(OUT_DIR, "signature_scores_bulk_persample.csv"))
cat("\nWrote", file.path(OUT_DIR, "signature_scores_bulk_persample.csv"),
    "(", nrow(out_dt), "samples ×", ncol(out_dt), "cols )\n")

# Sample summary by fibrosis stage
cat("\n=== Polyploid index by fibrosis stage ===\n")
sig_cols <- grep("^polyploid_index|^richter|^katsuda|^yin|^consensus", names(out_dt), value=TRUE)
sum_stage <- out_dt[!is.na(fibrosis_stage), lapply(.SD, mean, na.rm=TRUE),
                    by = fibrosis_stage, .SDcols = sig_cols][order(fibrosis_stage)]
print(sum_stage)
fwrite(sum_stage, file.path(OUT_DIR, "signature_scores_bulk_by_fibrosis.csv"))

# Sample summary by NAS
out_dt[, nas_bin := cut(nas_score, breaks=c(-Inf,1,3,5,Inf), labels=c("0-1","2-3","4-5","6+"))]
sum_nas <- out_dt[!is.na(nas_bin), lapply(.SD, mean, na.rm=TRUE),
                  by = nas_bin, .SDcols = sig_cols][order(nas_bin)]
cat("\n=== Polyploid index by NAS bin ===\n")
print(sum_nas)
fwrite(sum_nas, file.path(OUT_DIR, "signature_scores_bulk_by_nas.csv"))

cat("\n[", format(Sys.time()), "] DONE\n", sep="")
