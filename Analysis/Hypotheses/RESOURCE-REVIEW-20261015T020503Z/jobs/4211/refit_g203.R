#!/usr/bin/env Rscript
# verify2-G2-03, job 4211. Independent check of Fig 3F adjacent-stage marks.
# Assay: bulk liver RNA-seq fragment counts. Model mirrors
# scripts/figures/build_pi_stage_extensions.R fit_contrast(): voomWithQualityWeights,
# design ~ dataset + inferred_sex [+ covariates] + group, eBayes default,
# BH over all 24,115 genes per contrast. Unit: genetic-identity cluster.
# Covariates computed here from merged_dge.rds: log2 library size, Gini over all
# 24,115 gene counts, stored TMM norm factor; each z-scored within cohort over the
# 387-unit stage roster.
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table); library(parallel)})
set.seed(20261015)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
DGE <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/merged_dge.rds")
STG <- file.path(ROOT, "figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_sample_manifest.tsv")
OUT <- getwd(); FITS <- file.path(OUT, "fits"); dir.create(FITS, showWarnings = FALSE)
t0 <- Sys.time(); el <- function() round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

d <- readRDS(DGE)
stopifnot(nrow(d) == 24115L, ncol(d) == 774L)
s <- fread(STG)
stopifnot(nrow(s) == 387L, !anyDuplicated(s$sample_id), all(s$sample_id %in% colnames(d)))
cnt <- d$counts[, s$sample_id]
gini <- apply(cnt, 2, function(x) { x <- sort(as.numeric(x)); n <- length(x); sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x)) })
s[, `:=`(lib_size = d$samples[sample_id, "lib.size"], norm_factor = d$samples[sample_id, "norm.factors"], gini = gini[sample_id])]
zc <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v)) / sd(v))
s[, `:=`(c_lib = zc(log2(lib_size), dataset), c_gini = zc(gini, dataset), c_nf = zc(norm_factor, dataset))]
fwrite(s, file.path(OUT, "stage_roster_covariates.tsv"), sep = "\t")
cat("roster", nrow(s), "load", el(), "s\n")

tasks <- list(
  list(label = "base",   ref = "F0", cmp = "F1", cov = character(0), drop = NA),
  list(label = "no130970", ref = "F0", cmp = "F1", cov = character(0), drop = "GSE130970"),
  list(label = "lib",    ref = "F0", cmp = "F1", cov = "c_lib", drop = NA),
  list(label = "gini",   ref = "F0", cmp = "F1", cov = "c_gini", drop = NA),
  list(label = "nf",     ref = "F0", cmp = "F1", cov = "c_nf", drop = NA),
  list(label = "joint",  ref = "F0", cmp = "F1", cov = c("c_lib", "c_gini", "c_nf"), drop = NA),
  list(label = "joint",  ref = "F1", cmp = "F2", cov = c("c_lib", "c_gini", "c_nf"), drop = NA),
  list(label = "joint",  ref = "F2", cmp = "F3", cov = c("c_lib", "c_gini", "c_nf"), drop = NA),
  list(label = "joint",  ref = "F3", cmp = "F4", cov = c("c_lib", "c_gini", "c_nf"), drop = NA)
)

fit_task <- function(tk) {
  sub <- copy(s[fib_group %in% c(tk$ref, tk$cmp)])
  if (!is.na(tk$drop)) sub <- sub[dataset != tk$drop]
  setorder(sub, sample_id)
  contrast <- paste0(tk$ref, "_to_", tk$cmp)
  info <- data.frame(dataset = droplevels(factor(sub$dataset)), inferred_sex = droplevels(factor(sub$inferred_sex)),
                     group = factor(sub$fib_group, levels = c(tk$ref, tk$cmp)), row.names = sub$sample_id)
  for (cv in tk$cov) info[[cv]] <- sub[[cv]]
  fml <- as.formula(paste("~", paste(c("dataset", "inferred_sex", tk$cov, "group"), collapse = " + ")))
  design <- model.matrix(fml, info)
  stopifnot(qr(design)$rank == ncol(design))
  coef <- paste0("group", tk$cmp)
  dge <- d[, match(sub$sample_id, colnames(d))]
  v <- voomWithQualityWeights(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  res <- data.table(gene_id_versioned = rownames(tt), logFC = tt$logFC, P.Value = tt$P.Value,
                    FDR = p.adjust(tt$P.Value, "BH"), label = tk$label, contrast = contrast,
                    n_units = nrow(sub), n_ref = sum(info$group == tk$ref), n_cmp = sum(info$group == tk$cmp))
  stopifnot(nrow(res) == 24115L)
  fwrite(res, file.path(FITS, paste0(tk$label, "__", contrast, ".tsv.gz")), sep = "\t")
  data.table(label = tk$label, contrast = contrast, formula = deparse(fml), n_units = nrow(sub),
             n_ref = sum(info$group == tk$ref), n_cmp = sum(info$group == tk$cmp),
             n_bh05 = sum(res$FDR < 0.05), finished_s = el())
}
res <- mclapply(tasks, fit_task, mc.cores = 4L, mc.preschedule = FALSE)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) print(res[bad])
summ <- rbindlist(res[!bad])
fwrite(summ, file.path(OUT, "fit_summary.tsv"), sep = "\t")
print(summ)
cat("done", el(), "s\n")
