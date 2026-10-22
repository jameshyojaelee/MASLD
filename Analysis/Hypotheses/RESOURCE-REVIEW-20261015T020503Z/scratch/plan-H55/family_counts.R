# H55 planning: fix the gene family sizes from expression levels only.
# Reads counts and TPM; never uses any exposure. Seed not needed (no random step).
suppressMessages({library(edgeR); library(data.table)})
B <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
S <- paste0(B, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/plan-H55/")
meta <- fread(paste0(S, "design_meta_counts_only.tsv"))
meta <- meta[!is.na(SMCENTER) & SMCENTER != ""]
cat("model donors (DTHHRDY and SMCENTER present):", nrow(meta), "\n")
cnt <- fread(paste0(S, "inputs/gene_reads_2017-06-05_v8_liver.gct.gz"), skip = 2)
tpm <- fread(paste0(B, "data/external/allelic_refs/gtex_v8/gene_tpm_2017-06-05_v8_liver.gct.gz"), skip = 2)
cat("count rows", nrow(cnt), "tpm rows", nrow(tpm), "same Name order", identical(cnt$Name, tpm$Name), "\n")
cnt <- cnt[!grepl("_PAR_Y$", Name)]; tpm <- tpm[!grepl("_PAR_Y$", Name)]
ens <- sub("\\..*$", "", cnt$Name)
cat("after dropping _PAR_Y rows:", nrow(cnt), "dup unversioned", sum(duplicated(ens)), "\n")
M <- as.matrix(cnt[, meta$SAMPID, with = FALSE]); rownames(M) <- ens
Tm <- as.matrix(tpm[, meta$SAMPID, with = FALSE]); rownames(Tm) <- sub("\\..*$", "", tpm$Name)
cat("integer counts:", all(M == round(M)), " library size range (M):", round(range(colSums(M))/1e6, 1), "\n")
y <- DGEList(M)
keep <- filterByExpr(y)   # group = NULL, design = NULL: min sample size from all donors
cat("filterByExpr pass (all genes):", sum(keep), "\n")
deg <- fread(paste0(B, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"))
deg[, ens := sub("\\..*$", "", gene)]
can <- deg[padj < 0.05 & abs(logFC) > 0.5]
cat("canonical", nrow(can), "up", sum(can$logFC > 0), "down", sum(can$logFC < 0), "dup ens", sum(duplicated(can$ens)), "\n")
inG <- can$ens %in% ens
cat("canonical in GTEx v26 IDs:", sum(inG), " lost:", sum(!inG), "\n")
fam <- can[inG & can$ens %in% names(keep)[keep]]
cat("PRIMARY FAMILY canonical & filterByExpr:", nrow(fam), "up", sum(fam$logFC > 0), "down", sum(fam$logFC < 0), "\n")
med <- apply(Tm[can$ens[inG], ], 1, median)
cat("canonical in GTEx with median TPM >= 1:", sum(med >= 1), " >=0.1:", sum(med >= 0.1), " ==0:", sum(med == 0), "\n")
cat("overlap filterByExpr family vs median TPM>=1:", length(intersect(fam$ens, names(med)[med >= 1])), "\n")
cat("bulk |logFC| median in family:", round(median(abs(fam$logFC)), 3), " in lost-by-filter:", round(median(abs(can[inG & !(ens %in% fam$ens)]$logFC)), 3), "\n")
# composition marker overlap
mk <- fread(paste0(B, "data/external/liver_celltype_markers_combined.tsv"))
mk <- mk[source == "Guilliams2022_consensus_ttest"]
immune <- c("B cells","T cells","Macrophages","Mono+mono derived cells","Neutrophils","Resident NK","Circulating NK/NKT","Plasma cells","cDC1s","cDC2s","pDCs","Basophils")
sets <- list(hepatocyte = mk[cell_type == "Hepatocytes"]$gene, immune = unique(mk[cell_type %in% immune]$gene),
             stromal = unique(mk[cell_type %in% c("Fibroblasts","Endothelial cells")]$gene))
sym <- cnt$Description; names(sym) <- ens
for (n in names(sets)) {
  g <- ens[sym %in% sets[[n]]]
  cat("marker set", n, "symbols", length(sets[[n]]), "in GTEx", length(g), "pass filter", sum(keep[g]), "in primary family", sum(g %in% fam$ens), "\n")
}
# positive-control composite: H8 steatosis q<0.05, not canonical
h8 <- fread(paste0(B, "Analysis/Hypotheses/H8-20261009T015529Z/run-20261009T020920Z/meta_steatosis.tsv.gz"))
q <- h8[q_meta < 0.05]
cat("H8 steatosis meta genes", nrow(h8), "q<0.05", nrow(q), "up", sum(q$mu > 0), "down", sum(q$mu < 0), "\n")
pc <- q[!(gene_id %in% can$ens) & gene_id %in% ens]
cat("H8 q<0.05 not canonical, in GTEx:", nrow(pc), " and pass filterByExpr:", sum(keep[pc$gene_id]), " up", sum(pc$mu > 0 & keep[pc$gene_id]), "down", sum(pc$mu < 0 & keep[pc$gene_id]), "\n")
cat("PLIN2 canonical:", "ENSG00000147872" %in% can$ens, " PLIN2 passes filter:", keep["ENSG00000147872"], " PLIN2 in H8 q<0.05:", "ENSG00000147872" %in% q$gene_id, "\n")
cat("H8 q<0.05 among primary family:", sum(fam$ens %in% q$gene_id), " H8 tested among family:", sum(fam$ens %in% h8$gene_id), "\n")
fwrite(data.table(ens = fam$ens, symbol = fam$symbol, bulk_logFC = fam$logFC), paste0(S, "primary_family_genes.tsv"), sep = "\t")
fwrite(data.table(ens = pc$gene_id[keep[pc$gene_id]], symbol = pc$gene_name[keep[pc$gene_id]], h8_mu = pc$mu[keep[pc$gene_id]]), paste0(S, "pc_composite_genes.tsv"), sep = "\t")
