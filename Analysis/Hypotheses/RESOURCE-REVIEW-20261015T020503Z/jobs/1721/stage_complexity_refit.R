#!/usr/bin/env Rscript
# Review G2, job 1721. Assay: bulk liver RNA-seq fragment counts; limma voomWithQualityWeights
# (native to counts), exactly as scripts/figures/build_pi_stage_extensions.R fit_contrast():
# design ~ dataset + inferred_sex + group, eBayes default, BH over all 24,115 genes per contrast,
# no fold-change filter (Fig 3E and 3F rule).
# Added here: within-cohort-standardized complexity covariates (log2 library size, Gini, TMM norm
# factor) one at a time and jointly; and F1-vs-F0 leave-one-cohort-out arms.
# Unit: genetic-identity cluster (not verified donor). Writes only inside this job directory.
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table); library(parallel)})
set.seed(20261015)
RUN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z"
REL <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
STG <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/stage_sample_manifest.tsv"
CX <- file.path(RUN, "scratch/review-R03/complexity_774.tsv")
OUT <- getwd(); FITS <- file.path(OUT, "fits"); dir.create(FITS, showWarnings = FALSE)
t0 <- Sys.time(); el <- function() as.numeric(difftime(Sys.time(), t0, units = "secs"))
BUDGET <- as.numeric(Sys.getenv("BUDGET_S", "455"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

d <- readRDS(file.path(REL, "merged_dge.rds"))
stopifnot(nrow(d) == 24115L, ncol(d) == 774L)
s <- fread(STG); cx <- fread(CX)
stopifnot(all(s$sample_id %in% cx$sample_id), all(s$sample_id %in% colnames(d)), !anyDuplicated(s$sample_id))
s <- merge(s, cx[, .(sample_id, lib_size, gini, norm_factor)], by = "sample_id")
zc <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v)) / sd(v))
s[, `:=`(z_loglib = zc(log2(lib_size), dataset), z_gini = zc(gini, dataset), z_nf = zc(norm_factor, dataset))]
fwrite(s, file.path(OUT, "stage_roster_with_complexity.tsv"), sep = "\t")
cat("roster", nrow(s), "units; load", round(el()), "s\n")

covsets <- list(base = character(0), lib = "z_loglib", gini = "z_gini", nf = "z_nf",
                joint = c("z_loglib", "z_gini", "z_nf"))
est_s <- function(n) 45 * (n / 128)^1.9  # job 1032 timings: 128 units 45 s, 205 units 111 s

fit_task <- function(tk) {
  sub <- copy(s[fib_group %in% c(tk$ref, tk$cmp)])
  if (!is.na(tk$drop)) sub <- sub[dataset != tk$drop]
  n <- nrow(sub); tag <- paste(tk$label, tk$contrast, sep = "__")
  if (el() + 1.2 * est_s(n) > BUDGET) {
    fwrite(data.table(label = tk$label, contrast = tk$contrast, status = "skipped_time_budget", n_units = n),
           file.path(FITS, paste0(tag, ".counts.tsv")), sep = "\t")
    return(NULL)
  }
  ts <- el()
  setorder(sub, sample_id)
  info <- data.frame(dataset = droplevels(factor(sub$dataset)), inferred_sex = droplevels(factor(sub$inferred_sex)),
                     group = factor(sub$fib_group, levels = c(tk$ref, tk$cmp)),
                     z_loglib = sub$z_loglib, z_gini = sub$z_gini, z_nf = sub$z_nf, row.names = sub$sample_id)
  rhs <- c("dataset", "inferred_sex", covsets[[tk$cov]], "group")
  X <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), info)
  stopifnot(qr(X)$rank == ncol(X))
  coef <- paste0("group", tk$cmp)
  dge <- d[, match(sub$sample_id, colnames(d))]
  fit <- eBayes(lmFit(voomWithQualityWeights(dge, X, plot = FALSE), X))
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  fdr <- p.adjust(tt$P.Value, method = "BH")
  fwrite(data.table(gene_id_versioned = rownames(tt), logFC = tt$logFC, P.Value = tt$P.Value, FDR = fdr),
         file.path(FITS, paste0(tag, ".genes.tsv.gz")), sep = "\t")
  # covariate coefficient summary (how many genes the complexity terms themselves reach at BH 0.05)
  cov_sig <- if (length(covsets[[tk$cov]])) paste(sapply(covsets[[tk$cov]], function(cc)
    paste0(cc, ":", sum(p.adjust(fit$p.value[, cc], "BH") < 0.05))), collapse = ";") else ""
  res <- data.table(label = tk$label, contrast = tk$contrast, status = "fitted", n_units = n,
                    n_ref = sum(info$group == tk$ref), n_cmp = sum(info$group == tk$cmp),
                    cohorts = paste(levels(info$dataset), collapse = ","), covariates = paste(rhs[-length(rhs)], collapse = "+"),
                    n_genes = nrow(tt), bh_sig = sum(fdr < 0.05), up = sum(fdr < 0.05 & tt$logFC > 0),
                    down = sum(fdr < 0.05 & tt$logFC < 0), cov_bh_sig = cov_sig,
                    start_s = round(ts), fit_s = round(el() - ts))
  fwrite(res, file.path(FITS, paste0(tag, ".counts.tsv")), sep = "\t")
  cat(tag, "n", n, "bh_sig", res$bh_sig, "fit_s", res$fit_s, "end", round(el()), "\n")
  res
}

mk <- function(cov, ref, cmp, contrast, drop = NA_character_, label = cov)
  list(cov = cov, ref = ref, cmp = cmp, contrast = contrast, drop = drop, label = label)
fvf0 <- function(cov, k) mk(cov, "F0", k, paste0(k, "_vs_F0"))
adj <- function(cov, a, b) mk(cov, a, b, paste0(a, "_to_", b))
tasks <- c(
  lapply(c("F1", "F2", "F3", "F4"), function(k) fvf0("joint", k)),
  list(adj("joint", "F1", "F2"), adj("joint", "F2", "F3"), adj("joint", "F3", "F4")),
  list(mk("base", "F0", "F1", "F1_vs_F0", drop = "GSE130970", label = "loco_minus_GSE130970")),
  unlist(lapply(c("F1", "F4", "F2", "F3"), function(k) lapply(c("lib", "gini", "nf"), function(cv) fvf0(cv, k))), recursive = FALSE),
  list(mk("base", "F0", "F1", "F1_vs_F0", drop = "GSE135251", label = "loco_minus_GSE135251"),
       mk("base", "F0", "F1", "F1_vs_F0", drop = "GSE162694", label = "loco_minus_GSE162694")),
  unlist(lapply(list(c("F1", "F2"), c("F2", "F3"), c("F3", "F4")), function(p)
    lapply(c("lib", "gini", "nf"), function(cv) adj(cv, p[1], p[2]))), recursive = FALSE)
)
cat("tasks", length(tasks), "\n")
res <- mclapply(tasks, function(tk) tryCatch(fit_task(tk), error = function(e) {
  fwrite(data.table(label = tk$label, contrast = tk$contrast, status = paste("error", conditionMessage(e))),
         file.path(FITS, paste0(tk$label, "__", tk$contrast, ".counts.tsv")), sep = "\t"); NULL }),
  mc.cores = 4, mc.preschedule = FALSE)
fwrite(rbindlist(lapply(list.files(FITS, "\\.counts\\.tsv$", full.names = TRUE), fread), fill = TRUE),
       file.path(OUT, "complexity_stage_counts.tsv"), sep = "\t")
cat("done", round(el()), "s\n")
