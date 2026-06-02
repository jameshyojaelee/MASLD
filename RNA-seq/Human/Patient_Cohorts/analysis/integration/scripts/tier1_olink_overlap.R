#!/usr/bin/env Rscript
# tier1_olink_overlap.R
# Re-tabulate Tier 1 (1,905) DEG × Olink panel overlap.

suppressPackageStartupMessages({ library(data.table) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity")

dream <- fread(file.path(RDIR, "dream_results.csv"))
gene_meta <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
olink <- fread(file.path(BASE, "Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt"),
               select = c("Assay"))
olink_proteins <- unique(olink$Assay)
cat(sprintf("Olink panel: %d unique assays\n", length(olink_proteins)))

# Strip Ensembl version suffix and join symbol
strip_version <- function(x) sub("[.][0-9]+$", "", x)

dream[, base_id := strip_version(gene)]
dream <- merge(dream, gene_meta[, .(ensembl_base, gene_name)],
               by.x = "base_id", by.y = "ensembl_base", all.x = TRUE)

for (lfc in c(0.3, 0.5)) {
  deg <- dream[padj < 0.05 & abs(logFC) > lfc]
  syms <- deg[!is.na(gene_name), gene_name]
  ov   <- intersect(syms, olink_proteins)
  cat(sprintf("|LFC|>%.1f: %d DEGs, %d with symbol; overlap with Olink panel = %d (%.1f%% of DEGs)\n",
              lfc, nrow(deg), length(syms), length(ov),
              100 * length(ov) / length(syms)))
  fwrite(data.table(symbol = sort(ov)),
         file.path(OUT_DIR, sprintf("dream_olink_overlap_lfc%s.csv", sub("[.]", "", as.character(lfc)))))
}

# Tier 1 specifically — save the per-gene merge of dream + Olink-membership flag
deg_t1 <- dream[padj < 0.05 & abs(logFC) > 0.5]
deg_t1[, in_olink := !is.na(gene_name) & gene_name %in% olink_proteins]
fwrite(deg_t1[, .(gene, base_id, gene_name, logFC, padj, in_olink)],
       file.path(OUT_DIR, "tier1_olink_overlap_per_gene.csv"))
cat(sprintf("\nTier 1 in Olink: %d / %d\n", sum(deg_t1$in_olink), nrow(deg_t1)))
