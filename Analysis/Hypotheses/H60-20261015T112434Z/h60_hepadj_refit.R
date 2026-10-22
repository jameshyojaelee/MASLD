# H60 composition qualifier: R0 reproduction refit and R1 hepatocyte-fraction-adjusted refit.
# Recipe: refit_corrected_controls.R:186-218 (voomWithQualityWeights, lmFit, eBayes defaults,
# topTable coef group_binaryDisease, number=Inf, sort.by="none"). No random step; seed set anyway.
# Usage: Rscript h60_hepadj_refit.R <OUT>
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table)})
set.seed(20260923)
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
stopifnot(dir.exists(OUT))
RUN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H60-20261015T112434Z"
PRESPEC_SHA <- "420dbb486398504d2c9bcbdb59f4d4790784fa3edafa79c453942e1ae17bc8fb"
sha <- strsplit(system2("sha256sum", file.path(RUN, "PRESPEC_H60.md"), stdout = TRUE), " ")[[1]][1]
if (sha != PRESPEC_SHA) stop("PRESPEC hash mismatch: ", sha)
B <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
REL <- paste0(B, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/")

dge <- readRDS(paste0(REL, "merged_dge.rds"))
meta <- as.data.table(readRDS(paste0(REL, "meta_matched.rds")))
comp <- fread(paste0(B, "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz"))
deg <- fread(paste0(REL, "deg_results.csv"))
md <- fread(paste0(REL, "model_design.tsv"))
cat("dge dim", dim(dge), "\n")

ids <- colnames(dge)
sex <- meta$inferred_sex[match(ids, meta$sample_id)]
hep <- comp$Hepatocytes[match(ids, comp$sample_id)]
if (anyNA(sex)) stop("missing inferred_sex")
if (anyNA(hep)) stop("missing Hepatocytes")
info <- data.frame(
  group_binary = factor(dge$samples$group_binary, levels = c("Control", "Disease")),
  dataset = factor(dge$samples$dataset),
  inferred_sex = factor(sex),
  hepatocyte_fraction = as.numeric(hep)
)
rownames(info) <- ids

fit_one <- function(design, label) {
  stopifnot(qr(design)$rank == ncol(design), "group_binaryDisease" %in% colnames(design))
  v <- voomWithQualityWeights(dge, design)
  fit <- eBayes(lmFit(v, design))
  res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  cat(label, "genes", nrow(res), "\n")
  as.data.table(res)
}

d0 <- model.matrix(~ dataset + inferred_sex + group_binary, data = info)
# design matrix equality with the stored model_design.tsv
md_m <- as.matrix(md[match(ids, md$sample_id), -1, with = FALSE])
design_equal <- identical(colnames(md_m), colnames(d0)) && max(abs(md_m - d0)) == 0
cat("R0 design equals model_design.tsv:", design_equal, "\n")
r0 <- fit_one(d0, "R0")
same_order <- identical(r0$gene, deg$gene)
maxdiff <- if (same_order) max(abs(r0$t - deg$t)) else NA_real_
r0_pass <- isTRUE(same_order) && is.finite(maxdiff) && maxdiff < 1e-6
cat(sprintf("R0 same gene order %s max|t diff| %.3e pass %s\n", same_order, maxdiff, r0_pass))

d1 <- model.matrix(~ dataset + inferred_sex + hepatocyte_fraction + group_binary, data = info)
r1 <- fit_one(d1, "R1")
# VIF of hepatocyte fraction against the other design columns
X <- d1[, setdiff(colnames(d1), c("(Intercept)", "hepatocyte_fraction")), drop = FALSE]
r2 <- summary(lm(info$hepatocyte_fraction ~ X))$r.squared
vif <- 1 / (1 - r2)
r1_canon <- sum(r1$adj.P.Val < 0.05 & abs(r1$logFC) > 0.5)
cat(sprintf("R1 VIF(hepatocyte_fraction) %.4f canonical %d\n", vif, r1_canon))

fwrite(r0[, .(gene, logFC, AveExpr, t, P.Value, adj.P.Val)], file.path(OUT, "refit_R0_toptable.tsv.gz"), sep = "\t")
fwrite(r1[, .(gene, logFC, AveExpr, t, P.Value, adj.P.Val)], file.path(OUT, "refit_R1_hepadj_toptable.tsv.gz"), sep = "\t")
fwrite(data.table(
  item = c("R0_design_equals_model_design", "R0_same_gene_order", "R0_max_abs_t_diff", "R0_pass",
           "R1_VIF_hepatocyte_fraction", "R1_canonical_calls_padj005_abslogFC05", "R1_n_genes",
           "hep_median_control", "hep_median_disease"),
  value = c(design_equal, same_order, sprintf("%.6e", maxdiff), r0_pass, sprintf("%.6f", vif), r1_canon, nrow(r1),
            median(hep[info$group_binary == "Control"]), median(hep[info$group_binary == "Disease"]))
), file.path(OUT, "refit_checks.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_R.txt"))
cat("refit done\n")
