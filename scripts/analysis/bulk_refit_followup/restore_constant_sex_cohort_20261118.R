# Sex-adjusted GSE126848 sensitivity (n = 40: 23 Control M + 7 Disease M + 10 restored Disease F),
# using r8-recounted fragment counts for the 10 females and pre-voom keep_nz >= 3 filter.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(data.table)
})

base <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
out <- file.path(base, "figures/candidates/bulk-context-774-gse126848-f40-20261118")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

arm <- file.path(base, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
prior_path <- file.path(base, "figures/candidates/pi-stage-20261118/analysis/stage_extensions/cohort_disease_all_gene_results.tsv")
if (!file.exists(prior_path)) {
  prior_path <- file.path(base, "figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv")
}
f10_manifest <- fread(file.path(base, "scripts/analysis/bulk_refit_followup/library_manifest_gse126848_f10.tsv"))
f10_counts_root <- file.path(base, "figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_gse126848_f10_20261118/samples")

dge_all <- readRDS(file.path(arm, "merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(arm, "meta_matched.rds")))
m30 <- meta[dataset == "GSE126848"]
stopifnot(nrow(m30) == 30L, sum(m30$group_binary == "Control") == 23L, sum(m30$group_binary == "Disease") == 7L)

genes_24115 <- rownames(dge_all)
counts_30 <- as.matrix(dge_all$counts[, match(m30$sample_id, colnames(dge_all))])

counts_10_list <- lapply(f10_manifest$sample_id, function(sid) {
  fpath <- file.path(f10_counts_root, sid, "gene_counts.tsv")
  stopifnot(file.exists(fpath))
  dt <- fread(fpath)
  idx <- match(genes_24115, dt$gene_id)
  stopifnot(!anyNA(idx))
  as.numeric(dt[[sid]][idx])
})
counts_10 <- do.call(cbind, counts_10_list)
rownames(counts_10) <- genes_24115
colnames(counts_10) <- f10_manifest$sample_id

# Verify female sex markers XIST / RPS4Y1 on the 10 recounted female libraries
gmeta <- fread(file.path(base, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"))
xist_id <- gmeta[gene_name == "XIST"]$gene_id[1]
rps4y1_id <- gmeta[gene_name == "RPS4Y1"]$gene_id[1]
cat("XIST counts in 10 females:", paste(counts_10[xist_id, ], collapse = ","), "\n")
cat("RPS4Y1 counts in 10 females:", paste(counts_10[rps4y1_id, ], collapse = ","), "\n")
stopifnot(all(counts_10[xist_id, ] > 100), all(counts_10[rps4y1_id, ] < 20))

f10_meta <- data.table(
  sample_id = f10_manifest$sample_id,
  dataset = "GSE126848",
  condition = f10_manifest$condition,
  group_binary = f10_manifest$group_binary,
  inferred_sex = f10_manifest$inferred_sex,
  analysis_unit_id = paste0("ResourceGeneticCluster::", f10_manifest$sample_id),
  unit_evidence = "r8_recounted_female_sensitivity_restoration"
)
sel40 <- rbind(
  m30[, .(sample_id, dataset, condition, group_binary, inferred_sex, analysis_unit_id, unit_evidence)],
  f10_meta
)
setorder(sel40, sample_id)
stopifnot(nrow(sel40) == 40L,
          sum(sel40$group_binary == "Control") == 23L,
          sum(sel40$group_binary == "Disease") == 17L,
          sum(sel40$inferred_sex == "M") == 30L,
          sum(sel40$inferred_sex == "F") == 10L)

counts_40 <- cbind(counts_30, counts_10)[, sel40$sample_id, drop = FALSE]
dge_40 <- DGEList(counts = counts_40)
dge_40 <- calcNormFactors(dge_40)

info40 <- data.frame(
  group = factor(sel40$group_binary, levels = c("Control", "Disease")),
  inferred_sex = factor(sel40$inferred_sex, levels = c("F", "M")),
  row.names = sel40$sample_id
)
design40 <- model.matrix(~ group + inferred_sex, info40)
stopifnot(qr(design40)$rank == 3L, ncol(design40) == 3L)

keep_nz <- rowSums(dge_40$counts > 0) >= 3L
dge_sub <- dge_40[keep_nz, , keep.lib.sizes = FALSE]
v <- voomWithQualityWeights(dge_sub, design40, plot = FALSE)
stopifnot(all(is.finite(v$targets$sample.weights)), all(v$targets$sample.weights > 0))
fit <- eBayes(lmFit(v, design40))
tt <- topTable(fit, coef = "groupDisease", number = Inf, sort.by = "none")

prior <- fread(prior_path)
annotation <- unique(prior[, .(gene_id_versioned, gene_id_base, gene_name, gene_type, chromosome, is_canonical_chromosome, display_biotype)])
annotation <- annotation[match(genes_24115, gene_id_versioned)]
genes_sub <- rownames(dge_sub)
ok <- genes_24115 %in% genes_sub
n_bh <- length(genes_sub)
logfc <- rep(NA_real_, length(genes_24115)); logfc[ok] <- tt[genes_sub, "logFC"]
pval <- rep(NA_real_, length(genes_24115)); pval[ok] <- tt[genes_sub, "P.Value"]
fdr <- rep(NA_real_, length(genes_24115)); fdr[ok] <- p.adjust(tt[genes_sub, "P.Value"], "BH")
tstat <- rep(NA_real_, length(genes_24115)); tstat[ok] <- tt[genes_sub, "t"]

restored <- cbind(annotation, data.table(
  dataset = "GSE126848", contrast = "Disease_vs_Control",
  logFC = logfc, P.Value = pval, FDR = fdr, t = tstat,
  n_control = 23L, n_disease = 17L, n_participants = 40L,
  covariates = "group+inferred_sex", bh_family_size = n_bh
))
fwrite(restored, file.path(out, "gse126848_f40_all_gene_results.tsv"), sep = "\t")
fwrite(sel40, file.path(out, "GSE126848_selected_units.tsv"), sep = "\t")
fwrite(data.table(sample_id = sel40$sample_id, quality_weight = v$targets$sample.weights),
       file.path(out, "GSE126848_quality_weights.tsv"), sep = "\t")

# Direct comparison against n=30 male-only fit and 1,854 canonical DEGs
can <- fread(file.path(arm, "deg_results.csv"))[, .(gene_id_versioned = gene, full_logFC = logFC, canonical_deg = padj < .05 & abs(logFC) > .5)]
m30_prior <- prior[dataset == "GSE126848", .(gene_id_versioned, logFC_30 = logFC, FDR_30 = FDR)]
cmp <- merge(merge(restored, can, by = "gene_id_versioned"), m30_prior, by = "gene_id_versioned")

n_calls_40 <- sum(cmp$FDR < 0.05 & abs(cmp$logFC) > 0.5, na.rm = TRUE)
n_calls_30 <- sum(cmp$FDR_30 < 0.05 & abs(cmp$logFC_30) > 0.5, na.rm = TRUE)
can_sub <- cmp[canonical_deg == TRUE]
same_dir_40 <- sum(sign(can_sub$logFC) == sign(can_sub$full_logFC), na.rm = TRUE)
same_dir_frac_40 <- mean(sign(can_sub$logFC) == sign(can_sub$full_logFC), na.rm = TRUE)
cor_logfc_30_40 <- cor(cmp$logFC, cmp$logFC_30, use = "pairwise.complete.obs", method = "spearman")
cor_can_logfc <- cor(can_sub$logFC, can_sub$full_logFC, use = "pairwise.complete.obs", method = "spearman")

sum_dt <- data.table(
  model = "GSE126848_r8_n40_sex_adjusted",
  formula = "~ group + inferred_sex",
  n_control_M = 23L, n_disease_M = 7L, n_disease_F = 10L, n_total = 40L,
  n_tested = n_bh, n_calls_40 = n_calls_40, n_calls_30 = n_calls_30,
  canonical_genes = nrow(can_sub),
  canonical_same_dir_n_40 = same_dir_40,
  canonical_same_dir_frac_40 = same_dir_frac_40,
  spearman_logFC_vs_n30 = cor_logfc_30_40,
  spearman_canonical_logFC_vs_774 = cor_can_logfc
)
fwrite(sum_dt, file.path(out, "GSE126848_f40_summary.tsv"), sep = "\t")
print(sum_dt)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
