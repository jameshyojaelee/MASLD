# SS S3e reproduction check: compare an s3e_composition.R run on a release mock with the stored release values.
# Fractions vs Analysis/Deconvolution/corrected_774_20261009T141331Z (per cell type Pearson, max abs diff), ILR vs the
# ILR of the stored fractions, continuum position vs CG0 axis_raw, BF4 per-cohort P1 and H38 per-cohort b.
# Usage: Rscript s3e_verify_compare.R <dataset> <s3e out dir>
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
a <- commandArgs(trailingOnly = TRUE); ds <- a[1]; o <- a[2]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
ref <- read.delim(file.path(ROOT, "Analysis/Deconvolution/corrected_774_20261009T141331Z/results", ds, paste0(ds, "_bayesprism_proportions.tsv")), check.names = FALSE)
new <- fread(file.path(o, "fractions.tsv")); ct <- colnames(ref)
stopifnot(setequal(new$sample_alias, rownames(ref)))
R <- as.matrix(ref[new$sample_alias, ct]); N <- as.matrix(new[, ..ct])
tab <- data.table(cell_type = ct, pearson = sapply(ct, function(k) cor(R[, k], N[, k])), max_abs_diff = apply(abs(R - N), 2, max),
                  ref_mean = colMeans(R), new_mean = colMeans(N))
print(tab); fwrite(tab, file.path(o, "verify_fractions_vs_release.tsv"), sep = "\t")
cat("overall max abs diff", max(abs(R - N)), "\n")
u <- fread(file.path(ROOT, "Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z/units.tsv"))[dataset == ds]
b <- fread(file.path(o, "bf4_unit_scores.tsv"))
cg <- merge(u[, .(sample_alias = sample_id, y0 = as.numeric(scale(axis_raw)))], b[, .(sample_alias, y)], by = "sample_alias")
cat("continuum position vs CG0 axis_raw: n", nrow(cg), "pearson", cor(cg$y0, cg$y), "max abs diff", max(abs(cg$y0 - cg$y)), "\n")
r <- fromJSON(file.path(o, "results.json"))
bf4 <- fromJSON(file.path(ROOT, "Analysis/Hypotheses/BF4-20261024T010648Z/out-20261024T011020Z/results.json"))
cat("BF4 P1 per cohort: new", r$BF4_P1$share, "reference", bf4$P1$per_cohort[[ds]], "\n")
h <- fread(file.path(ROOT, "Analysis/Hypotheses/H38-20261009T180100Z/results/full-20261009T180646Z/t2_per_cohort.tsv"))[cohort == ds]
if (!is.null(r$H38)) cat("H38 b: new", r$H38$T2$b, "se", r$H38$T2$se, "persons", r$H38$T2$persons, "| reference", h$b, h$se, h$persons, "\n")
