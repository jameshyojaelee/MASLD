#!/usr/bin/env Rscript
# F3ef repair: cohort-omission refits of the Figure 3E/3F fibrosis contrasts from the adopted
# 774-unit release. Model copied from scripts/figures/build_pi_stage_extensions.R:142-241
# (fit_contrast: ~ dataset + inferred_sex + group, voomWithQualityWeights, lmFit, eBayes,
# BH over all 24,115 genes, no fold-change filter). Single-cohort fits follow fit_cohort
# (~ inferred_sex + group). No library-complexity covariates (those belong to H82).
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(parallel)})
set.seed(20260923)
out <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F3ef/refit-20261015T171823Z"
rel <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
released <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv"
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

d <- readRDS(file.path(rel, "merged_dge.rds"))
m <- as.data.table(readRDS(file.path(rel, "meta_matched.rds")))
m <- m[sample_id %in% colnames(d) & !is.na(inferred_sex) & !is.na(dataset)]
stopifnot(nrow(m) == 774L, identical(dim(d), c(24115L, 774L)))

stage <- m[dataset %in% c("GSE130970", "GSE135251", "GSE162694") & fibrosis_stage %in% 0:4]
stage[, grp := paste0("F", fibrosis_stage)]
stopifnot(nrow(stage) == 387L)
fwrite(stage[, .N, by = .(dataset, grp, group_binary)][order(dataset, grp, group_binary)],
       file.path(out, "stage_roster_387_counts.tsv"), sep = "\t")

cohorts <- c("GSE130970", "GSE135251", "GSE162694")
contrasts <- rbind(
  data.table(contrast = c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4"),
             ref = c("F0", "F1", "F2", "F3"), cmp = c("F1", "F2", "F3", "F4")),
  data.table(contrast = c("F2_vs_F0", "F3_vs_F0", "F4_vs_F0"), ref = "F0", cmp = c("F2", "F3", "F4")))
specs <- list()
for (i in seq_len(nrow(contrasts))) {
  specs[[length(specs) + 1]] <- c(as.list(contrasts[i]), fit = "base", drop = NA, only = NA)
  for (cc in cohorts) specs[[length(specs) + 1]] <- c(as.list(contrasts[i]), fit = paste0("drop_", cc), drop = cc, only = NA)
}
for (cc in cohorts) specs[[length(specs) + 1]] <- list(contrast = "F0_to_F1", ref = "F0", cmp = "F1",
                                                      fit = paste0("alone_", cc), drop = NA, only = cc)

fit_one <- function(spec) {
  sel <- stage[grp %in% c(spec$ref, spec$cmp)]
  if (!is.na(spec$drop)) sel <- sel[dataset != spec$drop]
  if (!is.na(spec$only)) sel <- sel[dataset == spec$only]
  setorder(sel, sample_id)
  info <- data.frame(dataset = droplevels(factor(sel$dataset)), inferred_sex = droplevels(factor(sel$inferred_sex)),
                     group = factor(sel$grp, levels = c(spec$ref, spec$cmp)), row.names = sel$sample_id)
  terms <- c(if (nlevels(info$dataset) > 1) "dataset", if (nlevels(info$inferred_sex) > 1) "inferred_sex", "group")
  design <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), info)
  stopifnot(qr(design)$rank == ncol(design))
  v <- voomWithQualityWeights(d[, rownames(info)], design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = paste0("group", spec$cmp), number = Inf, sort.by = "none")
  res <- data.table(gene_id_versioned = rownames(d), logFC = tt[rownames(d), "logFC"], P.Value = tt[rownames(d), "P.Value"])
  res[, FDR := p.adjust(P.Value, "BH")]
  fwrite(res, file.path(out, "fits", sprintf("%s__%s.tsv.gz", spec$fit, spec$contrast)), sep = "\t")
  data.table(fit = spec$fit, contrast = spec$contrast, n_ref = sum(info$group == spec$ref),
             n_cmp = sum(info$group == spec$cmp), n_units = nrow(info), n_cohorts = nlevels(info$dataset),
             design = paste(terms, collapse = "+"), n_bh05 = sum(res$FDR < 0.05),
             n_up = sum(res$FDR < 0.05 & res$logFC > 0), n_down = sum(res$FDR < 0.05 & res$logFC < 0))
}
summ <- rbindlist(mclapply(specs, fit_one, mc.cores = max(1L, ncores), mc.preschedule = FALSE))
setorder(summ, contrast, fit)
fwrite(summ, file.path(out, "fit_counts.tsv"), sep = "\t")

# Reproduction of the released fits (base) against stage_extension_all_gene_results.tsv
relt <- fread(released, select = c("gene_id_versioned", "logFC", "FDR", "axis", "contrast"))
relt <- relt[axis %in% c("fibrosis_adjacent", "fibrosis")]
relt[contrast == "F2_vs_F0" | contrast == "F3_vs_F0" | contrast == "F4_vs_F0", contrast := contrast]
repro <- rbindlist(lapply(unique(contrasts$contrast), function(cn) {
  b <- fread(file.path(out, "fits", sprintf("base__%s.tsv.gz", cn)))
  r <- relt[contrast == cn]
  x <- merge(b, r, by = "gene_id_versioned")
  data.table(contrast = cn, n_join = nrow(x), max_abs_dlogFC = max(abs(x$logFC.x - x$logFC.y)),
             max_abs_dFDR = max(abs(x$FDR.x - x$FDR.y)), n_bh05_refit = sum(x$FDR.x < 0.05), n_bh05_release = sum(x$FDR.y < 0.05))
}))
fwrite(repro, file.path(out, "reproduction_vs_release.tsv"), sep = "\t")
print(repro)
print(summ)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
