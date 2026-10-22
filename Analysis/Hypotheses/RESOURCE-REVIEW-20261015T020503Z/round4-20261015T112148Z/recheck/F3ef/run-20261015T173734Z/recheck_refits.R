#!/usr/bin/env Rscript
# Blind re-check F3ef: independent refit of fibrosis contrasts from the adopted 774-unit release
# inputs (merged_dge.rds, meta_matched.rds), following build_pi_stage_extensions.R fit_contrast.
# Adds: leave-one-cohort-out, single-cohort F0->F1, and cohort-by-stage interaction for F1/F4 vs F0.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(parallel)})
set.seed(20260923)
out <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F3ef/run-20261015T173734Z"
rel <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
released <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv"
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

dge <- readRDS(file.path(rel, "merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(rel, "meta_matched.rds")))
meta <- meta[sample_id %in% colnames(dge) & !is.na(inferred_sex) & !is.na(dataset)]
cat("units", nrow(meta), "genes", nrow(dge), "\n")
coh <- c("GSE130970", "GSE135251", "GSE162694")
st <- meta[dataset %in% coh & fibrosis_stage %in% 0:4]
st[, g := paste0("F", fibrosis_stage)]
cat("stage units", nrow(st), "\n")
fwrite(dcast(st[, .N, by = .(dataset, g)], dataset ~ g, value.var = "N"), file.path(out, "stage_counts.tsv"), sep = "\t")
fwrite(meta[, .N, by = .(dataset, fibrosis_stage, group_binary)][order(dataset, fibrosis_stage)],
       file.path(out, "fibrosis_codes_774.tsv"), sep = "\t")

run_fit <- function(lab, contrast, ref, cmp, keep) {
  s <- st[g %in% c(ref, cmp) & dataset %in% keep]
  setorder(s, sample_id)
  df <- data.frame(dataset = droplevels(factor(s$dataset)), sex = droplevels(factor(s$inferred_sex)),
                   group = factor(s$g, levels = c(ref, cmp)), row.names = s$sample_id)
  f <- "~ group"
  if (nlevels(df$sex) > 1) f <- paste("~ sex +", sub("~ ", "", f))
  if (nlevels(df$dataset) > 1) f <- paste("~ dataset +", sub("~ ", "", f))
  X <- model.matrix(as.formula(f), df)
  v <- voomWithQualityWeights(dge[, rownames(df)], X, plot = FALSE)
  fit <- eBayes(lmFit(v, X))
  tt <- topTable(fit, coef = paste0("group", cmp), number = Inf, sort.by = "none")
  r <- data.table(gene = rownames(dge), logFC = tt[rownames(dge), "logFC"], P = tt[rownames(dge), "P.Value"])
  r[, q := p.adjust(P, "BH")]
  fwrite(r, file.path(out, "fits", paste0(lab, "__", contrast, ".tsv.gz")), sep = "\t")
  data.table(fit = lab, contrast = contrast, n_ref = sum(df$group == ref), n_cmp = sum(df$group == cmp),
             n_units = nrow(df), formula = f, n_q05 = sum(r$q < 0.05), up = sum(r$q < 0.05 & r$logFC > 0),
             down = sum(r$q < 0.05 & r$logFC < 0))
}
cons <- data.table(contrast = c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4", "F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0"),
                   ref = c("F0", "F1", "F2", "F3", "F0", "F0", "F0", "F0"),
                   cmp = c("F1", "F2", "F3", "F4", "F1", "F2", "F3", "F4"))
jobs <- list()
for (i in seq_len(nrow(cons))) {
  jobs[[length(jobs) + 1]] <- list("all3", cons$contrast[i], cons$ref[i], cons$cmp[i], coh)
  for (c0 in coh) jobs[[length(jobs) + 1]] <- list(paste0("without_", c0), cons$contrast[i], cons$ref[i], cons$cmp[i], setdiff(coh, c0))
}
for (c0 in coh) jobs[[length(jobs) + 1]] <- list(paste0("only_", c0), "F0_to_F1", "F0", "F1", c0)
res <- rbindlist(mclapply(jobs, function(j) do.call(run_fit, j), mc.cores = max(1L, ncores), mc.preschedule = FALSE))
setorder(res, contrast, fit)
fwrite(res, file.path(out, "my_counts.tsv"), sep = "\t")
print(res)

# Reproduction of release
relt <- fread(released, select = c("gene_id_versioned", "gene_name", "logFC", "FDR", "axis", "contrast"))
rep <- rbindlist(lapply(cons$contrast, function(cn) {
  a <- fread(file.path(out, "fits", paste0("all3__", cn, ".tsv.gz")))
  ax <- if (grepl("_to_", cn)) "fibrosis_adjacent" else "fibrosis"
  r <- relt[axis == ax & contrast == cn]
  x <- merge(a, r, by.x = "gene", by.y = "gene_id_versioned")
  data.table(contrast = cn, n = nrow(x), max_dlogFC = max(abs(x$logFC.x - x$logFC.y)), max_dq = max(abs(x$q - x$FDR)),
             mine = sum(x$q < 0.05), release = sum(x$FDR < 0.05))
}))
fwrite(rep, file.path(out, "my_reproduction.tsv"), sep = "\t"); print(rep)

# Cohort-by-stage interaction model: ~ inferred_sex + dataset + dataset:group
inter <- rbindlist(lapply(c("F1", "F4"), function(cmp) {
  s <- st[g %in% c("F0", cmp)]; setorder(s, sample_id)
  df <- data.frame(dataset = factor(s$dataset, levels = coh), sex = factor(s$inferred_sex),
                   group = factor(s$g, levels = c("F0", cmp)), row.names = s$sample_id)
  X <- model.matrix(~ sex + dataset + dataset:group, df)
  colnames(X) <- make.names(colnames(X))
  v <- voomWithQualityWeights(dge[, rownames(df)], X, plot = FALSE)
  fit <- lmFit(v, X)
  gcols <- grep("group", colnames(X), value = TRUE)
  stopifnot(length(gcols) == 3L)
  C <- matrix(0, ncol(X), 2, dimnames = list(colnames(X), c("d2", "d3")))
  C[gcols[2], 1] <- 1; C[gcols[1], 1] <- -1; C[gcols[3], 2] <- 1; C[gcols[1], 2] <- -1
  f2 <- eBayes(contrasts.fit(fit, C))
  tt <- topTable(f2, number = Inf, sort.by = "none")
  fe <- eBayes(fit)
  co <- fe$coefficients[, gcols]; colnames(co) <- coh
  r <- data.table(gene = rownames(dge), F_P = tt[rownames(dge), "P.Value"], co[rownames(dge), ])
  r[, F_q := p.adjust(F_P, "BH")]
  r[, contrast := paste0(cmp, "_vs_F0")]
  r
}))
fwrite(inter, file.path(out, "interaction_F1_F4.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("DONE\n")
