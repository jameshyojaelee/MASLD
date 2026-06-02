#!/usr/bin/env Rscript
# 483_coloc_enrichment.R
# Per-program hypergeometric enrichment of top loading genes in COLOC PP.H4.susie > 0.5 set.
# Output per-program enrichment across 28 GWAS (if broken down by GWAS).

suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
coloc_f <- file.path(root, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
out_dir <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/validation/coloc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
name <- args[1]; k <- as.integer(args[2])
top_f <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/cnmf_annot", name,
                   sprintf("program_topgenes.k%d.tsv", k))
stopifnot(file.exists(top_f), file.exists(coloc_f))

top <- fread(top_f)
coloc <- fread(coloc_f)
# Detect column names — prefer coloc_best_susie_pp4
susie_col <- intersect(c("coloc_best_susie_pp4", "coloc_best_pp4", "PP.H4.susie"), names(coloc))[1]
if (is.na(susie_col)) stop(sprintf("No SuSiE PP4 col found in %s", coloc_f))
# gene name column
gene_col <- intersect(c("gene", "gene_name", "gene_symbol", "symbol"), names(coloc))[1]
if (is.na(gene_col)) stop("No gene_name in coloc file")

universe <- unique(coloc[[gene_col]])
coloc05 <- unique(coloc[[gene_col]][coloc[[susie_col]] > 0.5])
coloc08 <- unique(coloc[[gene_col]][coloc[[susie_col]] > 0.8])

cat(sprintf("[483] universe=%d, coloc>0.5=%d, coloc>0.8=%d\n",
            length(universe), length(coloc05), length(coloc08)))

enrich <- function(top_genes, set_name, set_genes) {
  overlap <- length(intersect(top_genes, set_genes))
  set_a <- length(top_genes); set_b <- length(set_genes); univ <- length(universe)
  p <- phyper(overlap - 1, set_b, univ - set_b, set_a, lower.tail = FALSE)
  data.table(set = set_name, overlap = overlap, top_n = set_a, set_size = set_b, universe = univ, pval = p)
}

results <- list()
for (p in unique(top$program)) {
  tg <- top[program == p]$gene_name
  tg <- tg[tg %in% universe]
  r5 <- enrich(tg, "coloc_pp4_gt_0.5", coloc05)
  r8 <- enrich(tg, "coloc_pp4_gt_0.8", coloc08)
  r5[, program := p]; r8[, program := p]
  results[[length(results) + 1]] <- rbind(r5, r8)
}
out <- rbindlist(results)
out[, qval := p.adjust(pval, method = "BH"), by = set]
fwrite(out, file.path(out_dir, sprintf("coloc_enrichment_%s_k%d.tsv", name, k)), sep = "\t")
cat(sprintf("[483] wrote enrichment: %d rows. Programs with q<0.05 in pp4>0.5: %d\n",
            nrow(out), out[set == "coloc_pp4_gt_0.5" & qval < 0.05, .N]))
