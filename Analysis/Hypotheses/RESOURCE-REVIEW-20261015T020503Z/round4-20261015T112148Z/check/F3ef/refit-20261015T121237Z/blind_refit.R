#!/usr/bin/env Rscript
# Blind check F3ef: independent refit of the Fig 3E/3F/S3 stage and NAS models
# from the adopted 774-unit release (merged_dge.rds + meta_matched.rds).
# Model as in the released stage analysis: voomWithQualityWeights + lmFit + eBayes,
# ~ dataset + inferred_sex + group, BH over all 24,115 genes per contrast, no fold filter.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(parallel)})
set.seed(20260923)
out <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F3ef/refit-20261015T121237Z"
rel <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
ann_path <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

d <- readRDS(file.path(rel, "merged_dge.rds"))
m <- as.data.table(readRDS(file.path(rel, "meta_matched.rds")))
m <- m[sample_id %in% colnames(d) & !is.na(inferred_sex) & !is.na(dataset)]
stopifnot(nrow(m) == 774L, identical(dim(d), c(24115L, 774L)))
ann <- fread(ann_path)[, .(gene_id, gene_name, ensembl_base)]
ann <- ann[match(rownames(d), gene_id)]
stopifnot(!anyNA(ann$gene_id))

# Library complexity per unit, from the released counts.
gini <- apply(d$counts, 2, function(x) {
  x <- sort(as.numeric(x)); n <- length(x); sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
})
m[, `:=`(lib_size = d$samples[sample_id, "lib.size"], norm_factor = d$samples[sample_id, "norm.factors"],
         gini = gini[sample_id], control_status = as.integer(group_binary == "Control"))]

stage <- m[dataset %in% c("GSE130970", "GSE135251", "GSE162694") & fibrosis_stage %in% 0:4]
stage[, grp := paste0("F", fibrosis_stage)]
zc <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v)) / sd(v))
stage[, `:=`(z_lib = zc(log2(lib_size), dataset), z_gini = zc(gini, dataset), z_nf = zc(norm_factor, dataset))]
fwrite(stage[, .(sample_id, dataset, group_binary, fibrosis_stage, inferred_sex, lib_size, gini, norm_factor,
                 z_lib, z_gini, z_nf)], file.path(out, "stage_roster_387.tsv"), sep = "\t")

bin_nas <- function(x) fcase(x == 0, "NAS0", x <= 2, "NAS1-2", x <= 4, "NAS3-4", x <= 8, "NAS5-8", default = NA_character_)
nas <- m[!is.na(nas_score)]
nas[, grp := bin_nas(as.integer(nas_score))]
nas <- nas[!is.na(grp)]
fwrite(nas[, .N, by = .(dataset, group_binary, grp)][order(grp, dataset, group_binary)],
       file.path(out, "nas_group_composition.tsv"), sep = "\t")

fit_one <- function(spec) {
  tab <- if (spec$axis == "nas") nas else stage
  sel <- tab[grp %in% c(spec$ref, spec$cmp)]
  if (!is.na(spec$drop)) sel <- sel[dataset != spec$drop]
  if (!is.na(spec$only)) sel <- sel[dataset == spec$only]
  if (isTRUE(spec$disease_f0)) sel <- sel[!(grp == spec$ref & group_binary == "Control")]
  setorder(sel, sample_id)
  info <- data.frame(dataset = droplevels(factor(sel$dataset)), inferred_sex = droplevels(factor(sel$inferred_sex)),
                     group = factor(sel$grp, levels = c(spec$ref, spec$cmp)),
                     z_lib = sel$z_lib, z_gini = sel$z_gini, z_nf = sel$z_nf,
                     control_status = sel$control_status, row.names = sel$sample_id)
  terms <- c(if (nlevels(info$dataset) > 1) "dataset", "inferred_sex", spec$cov)
  if (spec$interaction) {
    fml <- as.formula("~ inferred_sex + dataset + dataset:group")
  } else {
    fml <- as.formula(paste("~", paste(c(terms, "group"), collapse = " + ")))
  }
  design <- model.matrix(fml, info)
  stopifnot(qr(design)$rank == ncol(design))
  v <- voomWithQualityWeights(d[, rownames(info)], design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  if (!spec$interaction) {
    cf <- paste0("group", spec$cmp)
    tt <- topTable(fit, coef = cf, number = Inf, sort.by = "none")
    res <- data.table(gene_id = rownames(d), logFC = tt[rownames(d), "logFC"], P = tt[rownames(d), "P.Value"])
    res[, FDR := p.adjust(P, "BH")]
  } else {
    gcols <- grep(":group", colnames(design), value = TRUE)
    stopifnot(length(gcols) == 3L)
    cm <- matrix(0, ncol(design), 2, dimnames = list(colnames(design), c("d2", "d3")))
    cm[gcols[2], 1] <- 1; cm[gcols[1], 1] <- -1
    cm[gcols[3], 2] <- 1; cm[gcols[1], 2] <- -1
    f2 <- eBayes(contrasts.fit(fit, cm))
    ft <- topTable(f2, number = Inf, sort.by = "none")
    res <- data.table(gene_id = rownames(d), P_int = ft[rownames(d), "P.Value"])
    res[, FDR_int := p.adjust(P_int, "BH")]
    for (k in seq_along(gcols)) set(res, j = paste0("b_", sub("dataset", "", sub(":group.*", "", gcols[k]))),
                                    value = fit$coefficients[rownames(d), gcols[k]])
  }
  res[, `:=`(fit = spec$label, contrast = spec$contrast, n_units = nrow(sel),
             n_ref = sum(sel$grp == spec$ref), n_cmp = sum(sel$grp == spec$cmp))]
  fwrite(res, file.path(out, "fits", paste0(spec$label, "__", spec$contrast, ".tsv.gz")), sep = "\t")
  data.table(fit = spec$label, contrast = spec$contrast, n_units = nrow(sel), n_ref = sum(sel$grp == spec$ref),
             n_cmp = sum(sel$grp == spec$cmp), n_cohorts = uniqueN(sel$dataset),
             n_bh05 = if (spec$interaction) sum(res$FDR_int < 0.05) else sum(res$FDR < 0.05),
             formula = paste(deparse(fml), collapse = ""))
}

mk <- function(label, ref, cmp, contrast, cov = character(0), drop = NA, only = NA, disease_f0 = FALSE,
               interaction = FALSE, axis = "stage")
  list(label = label, ref = ref, cmp = cmp, contrast = contrast, cov = cov, drop = drop, only = only,
       disease_f0 = disease_f0, interaction = interaction, axis = axis)
specs <- list()
adj <- list(c("F0", "F1"), c("F1", "F2"), c("F2", "F3"), c("F3", "F4"))
vsf0 <- paste0("F", 2:4)  # F1 vs F0 is identical to F0_to_F1
covs <- list(base = character(0), lib = "z_lib", gini = "z_gini", nf = "z_nf", joint = c("z_lib", "z_gini", "z_nf"))
for (cv in names(covs)) {
  for (p in adj) specs[[length(specs) + 1]] <- mk(cv, p[1], p[2], paste0(p[1], "_to_", p[2]), covs[[cv]])
  for (k in vsf0) specs[[length(specs) + 1]] <- mk(cv, "F0", k, paste0(k, "_vs_F0"), covs[[cv]])
}
for (dc in c("GSE130970", "GSE135251", "GSE162694")) {
  for (k in paste0("F", 1:4)) specs[[length(specs) + 1]] <- mk(paste0("drop_", dc), "F0", k, paste0(k, "_vs_F0"), drop = dc)
  specs[[length(specs) + 1]] <- mk(paste0("alone_", dc), "F0", "F1", "F1_vs_F0", only = dc)
}
for (k in paste0("F", 1:4)) {
  specs[[length(specs) + 1]] <- mk("armB_diseaseF0", "F0", k, paste0(k, "_vs_F0"), disease_f0 = TRUE)
  specs[[length(specs) + 1]] <- mk("armC_ctrlterm", "F0", k, paste0(k, "_vs_F0"), cov = "control_status")
}
for (k in c("F1", "F4")) specs[[length(specs) + 1]] <- mk("interaction", "F0", k, paste0(k, "_vs_F0"), interaction = TRUE)
for (p in list(c("NAS0", "NAS1-2"), c("NAS0", "NAS3-4"), c("NAS0", "NAS5-8"), c("NAS1-2", "NAS3-4"), c("NAS3-4", "NAS5-8")))
  specs[[length(specs) + 1]] <- mk("nas", p[1], p[2], paste0(p[2], "_vs_", p[1]), axis = "nas")

cat("fits:", length(specs), "cores:", ncores, "\n")
res <- mclapply(specs, function(s) tryCatch(fit_one(s), error = function(e)
  data.table(fit = s$label, contrast = s$contrast, error = conditionMessage(e))), mc.cores = ncores, mc.preschedule = FALSE)
summ <- rbindlist(res, fill = TRUE)
fwrite(summ, file.path(out, "fit_counts.tsv"), sep = "\t")
print(summ)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
