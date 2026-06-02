#!/usr/bin/env Rscript
# tier1_plasma_bridge.R
# Reproduce tissue_plasma_bridge with Tier 1 (1,551 at padj<0.05, |logFC|>0.5, kallisto canonical)
# AND with the broader |LFC|>0.3 set, against the SAME 1,267-protein bridge
# Olink panel that produced the original "812".

suppressPackageStartupMessages({ library(data.table) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
BRIDGE <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier/tissue_plasma_bridge.csv"))
cat("Bridge rows:", nrow(BRIDGE), "\n")
cat("Bridge in_plasma=TRUE:", sum(BRIDGE$in_plasma), "\n")

dream <- fread(file.path(RDIR, "dream_results.csv"))
gene_meta <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
strip_v <- function(x) sub("[.][0-9]+$", "", x)
dream[, base_id := strip_v(gene)]
dream <- merge(dream, gene_meta[, .(ensembl_base, gene_name)],
               by.x = "base_id", by.y = "ensembl_base", all.x = TRUE)

# Tier 1 (canonical)
deg_t1  <- dream[padj < 0.05 & abs(logFC) > 0.5, gene_name]
deg_t1  <- deg_t1[!is.na(deg_t1)]
ov_t1   <- intersect(deg_t1,  BRIDGE[in_plasma == TRUE, human_symbol])

# Broader (|LFC|>0.3) for comparison
deg_03  <- dream[padj < 0.05 & abs(logFC) > 0.3, gene_name]
deg_03  <- deg_03[!is.na(deg_03)]
ov_03   <- intersect(deg_03, BRIDGE[in_plasma == TRUE, human_symbol])

# padj<0.1 no LFC (matches the original 16,258 used by 812 bridge)
deg_p1  <- dream[padj < 0.1, gene_name]
deg_p1  <- deg_p1[!is.na(deg_p1)]
ov_p1   <- intersect(deg_p1, BRIDGE[in_plasma == TRUE, human_symbol])

cat(sprintf("\nUsing the SAME bridge Olink panel (1,267 proteins):\n"))
cat(sprintf("  Tier 1 (1,551 @ padj<0.05, |LFC|>0.5, kallisto): %d in plasma (%.1f%% of Tier 1)\n",
            length(ov_t1),  100 * length(ov_t1)  / length(deg_t1)))
cat(sprintf("  |LFC|>0.3 (4,370, kallisto canonical):       %d in plasma (%.1f%% of |LFC|>0.3)\n",
            length(ov_03),  100 * length(ov_03)  / length(deg_03)))
cat(sprintf("  padj<0.1 (15,732, original method):         %d in plasma (matches original 812 ballpark)\n",
            length(ov_p1)))

# Save Tier 1 list
fwrite(data.table(symbol = sort(ov_t1)), file.path(OUT_DIR, "tier1_plasma_bridge_overlap.csv"))
cat(sprintf("\nSaved: tier1_plasma_bridge_overlap.csv (%d genes)\n", length(ov_t1)))
