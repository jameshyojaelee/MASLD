#!/usr/bin/env Rscript
# 342a_richter_DE.R - Richter 2021 polyploid DE via DESeq2 with paired design.
# Pseudobulk per (individual × ploidy), 4n vs 2n in hepatocytes.

suppressPackageStartupMessages({
  library(DESeq2)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(data.table)
})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
setwd(ROOT)

H5AD    <- "data/external/richter2021/richter2021_hepatocytes.h5ad"
OUT_DIR <- "data/ploidy_signatures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("[", format(Sys.time()), "] Loading h5ad via zellkonverter\n", sep="")
sce <- readH5AD(H5AD, X_name = "counts", reader = "python")  # X_name guesses; .X may be only matrix
# inspect
cat("  cells:", ncol(sce), "  genes:", nrow(sce), "\n")
cat("  colData cols:", paste(colnames(colData(sce)), collapse=", "), "\n")

# Filter
keep <- (colData(sce)$passing_QC == "pass" &
         colData(sce)$inferred_cell_label == "Hepatocytes" &
         colData(sce)$Ploidy %in% c("2n","4n"))
sce <- sce[, keep]
cat("  after filter:", ncol(sce), "\n")

# Pseudobulk by (Individual × Ploidy)
ind <- gsub("[^A-Za-z0-9]", "_", as.character(colData(sce)$ID.Individual))
ploidy <- as.character(colData(sce)$Ploidy)
pb_key <- paste0(ind, "@@", ploidy)
counts_mat <- assay(sce, 1L)  # first assay
counts_mat <- round(as.matrix(counts_mat))  # snRNA-seq2 floats -> int counts

groups <- unique(pb_key)
pb_counts <- matrix(0L, nrow = nrow(counts_mat), ncol = length(groups),
                    dimnames = list(rownames(counts_mat), groups))
pb_meta <- data.table(sample = groups, individual = NA_character_, ploidy = NA_character_,
                      n_cells = NA_integer_)
for (g in groups) {
  cells <- pb_key == g
  pb_counts[, g] <- rowSums(counts_mat[, cells, drop=FALSE])
  parts <- strsplit(g, "@@", fixed = TRUE)[[1]]
  pb_meta[sample == g, `:=`(individual = parts[1], ploidy = parts[2], n_cells = sum(cells))]
}

# Drop tiny groups (<5 cells)
keep_samples <- pb_meta$n_cells >= 5
pb_counts <- pb_counts[, keep_samples]
pb_meta   <- pb_meta[keep_samples]
cat("  pseudobulk samples:", ncol(pb_counts), "(>=5 cells each)\n")
print(pb_meta)

# Gene filter
ok_genes <- rowSums(pb_counts > 5) >= 3
pb_counts <- pb_counts[ok_genes, ]
cat("  genes kept (>5 in >=3 samples):", nrow(pb_counts), "\n")

# DESeq2 with individual + ploidy
pb_meta[, individual := factor(individual)]
pb_meta[, ploidy := factor(ploidy, levels = c("2n","4n"))]
dds <- DESeqDataSetFromMatrix(countData = pb_counts,
                              colData   = data.frame(pb_meta, row.names = pb_meta$sample),
                              design    = ~ individual + ploidy)
dds <- DESeq(dds, fitType = "local", quiet = TRUE)
res <- results(dds, name = "ploidy_4n_vs_2n")
cat("\nDE summary:\n"); print(summary(res))

# Map Ensembl -> gene_name via rowData
gene_names <- rowData(sce)$gene_name[match(rownames(res), rownames(sce))]
biotype    <- rowData(sce)$biotype[match(rownames(res), rownames(sce))]
res_dt <- data.table(
  gene_ensembl   = rownames(res),
  gene_symbol    = as.character(gene_names),
  biotype        = as.character(biotype),
  log2FoldChange = res$log2FoldChange,
  lfcSE          = res$lfcSE,
  stat           = res$stat,   # Wald statistic = LFC/SE (rank metric for GSEA)
  pvalue         = res$pvalue,
  padj           = res$padj,
  baseMean       = res$baseMean
)
res_dt <- res_dt[order(padj)]
fwrite(res_dt, file.path(OUT_DIR, "richter2021_de_full.tsv"), sep = "\t")
cat("\nWrote", file.path(OUT_DIR, "richter2021_de_full.tsv"), "(", nrow(res_dt), "genes)\n")

# Signatures - require human ortholog via mouse-symbol uppercase rule (standard convention)
# Mouse gene_symbol -> human via uppercase. Mouse-specific genes (Gm, Rik, lowercase-only) filter out.
res_dt[, human_symbol := toupper(gene_symbol)]
# Filter mouse-specific transcripts that won't transfer:
# - empty/NA symbols
# - Gm-prefixed (mouse gene model predictions)
# - Rik-suffixed (RIKEN cDNAs)
# - "BC" prefixed (BC clone mouse names)
bad <- is.na(res_dt$gene_symbol) | res_dt$gene_symbol == "" |
       grepl("^Gm[0-9]+$", res_dt$gene_symbol) |
       grepl("Rik$", res_dt$gene_symbol) |
       grepl("^BC[0-9]+$", res_dt$gene_symbol) |
       grepl("^[0-9]", res_dt$gene_symbol)  # numeric-prefix names
res_dt[bad, human_symbol := NA_character_]
fwrite(res_dt, file.path(OUT_DIR, "richter2021_de_full.tsv"), sep = "\t")  # re-save with human_symbol

# Filter to genes with human ortholog candidate + protein_coding + decent expression
res_clean <- res_dt[!is.na(stat) & !is.na(log2FoldChange) & !is.na(human_symbol) &
                    biotype == "protein_coding" & baseMean > 10]
cat("genes after filter (protein-coding, named-transferable, baseMean>10):", nrow(res_clean), "\n")
# Rank by Wald statistic (= LFC/SE) — standard pre-ranking metric for GSEA-style gene sets
# Top N up + top N down without strict FDR filter; FDR enrichment characterized via top-K cumulative dist
N_SIG <- 100  # signature size per direction
sig_up <- res_clean[stat > 0][order(-stat)][1:N_SIG]
sig_dn <- res_clean[stat < 0][order(stat)][1:N_SIG]
sig_up <- sig_up[!is.na(gene_ensembl)]
sig_dn <- sig_dn[!is.na(gene_ensembl)]
cat("\nN sig up (top by Wald):", nrow(sig_up), "  N sig down:", nrow(sig_dn), "\n")
cat("up: Wald stat range", range(sig_up$stat), "padj range", range(sig_up$padj, na.rm=TRUE), "\n")
cat("dn: Wald stat range", range(sig_dn$stat), "padj range", range(sig_dn$padj, na.rm=TRUE), "\n")
cat("\nSignature: up_in_4n =", nrow(sig_up), "  down_in_4n =", nrow(sig_dn), "\n")
cat("Top 10 up genes:", paste(head(sig_up$gene_symbol, 10), collapse=", "), "\n")
cat("Top 10 down genes:", paste(head(sig_dn$gene_symbol, 10), collapse=", "), "\n")

fwrite(sig_up[, .(gene_ensembl, gene_symbol, human_symbol, log2FoldChange, lfcSE, stat, pvalue, padj, biotype, baseMean)],
       file.path(OUT_DIR, "richter2021_polyploid_up.tsv"), sep = "\t")
fwrite(sig_dn[, .(gene_ensembl, gene_symbol, human_symbol, log2FoldChange, lfcSE, stat, pvalue, padj, biotype, baseMean)],
       file.path(OUT_DIR, "richter2021_polyploid_down.tsv"), sep = "\t")

# Mechanism marker overlap - use HUMAN ortholog symbols
mm <- fread("data/ploidy_signatures/mechanism_markers.tsv")
mm_human <- toupper(mm$gene)
up_syms <- toupper(sig_up$human_symbol)
dn_syms <- toupper(sig_dn$human_symbol)
ovl_up <- intersect(mm_human, up_syms)
ovl_dn <- intersect(mm_human, dn_syms)
cat("\nMechanism marker overlap:\n")
cat("  up_in_4n:  ", paste(ovl_up, collapse=", "), "\n")
cat("  down_in_4n:", paste(ovl_dn, collapse=", "), "\n")
cat("  total:", length(ovl_up) + length(ovl_dn), "of", length(mm_human), "markers\n")
cat("[", format(Sys.time()), "] DONE\n", sep="")
