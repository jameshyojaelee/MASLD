#!/usr/bin/env Rscript
# Quick login-node probe (chr22 only) of the sdY implied by summary statistics, for
# one study passed as argument. Full estimate: estimate_sdy.R (sbatch).
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
s <- commandArgs(TRUE)[1]
reg <- fread(file.path(FM, "config/gwas_registry.tsv"))[study_name == s]
afs <- fread(file.path(FM, "config/gwas_af_sources.tsv"))[study_name == s]
use_af <- afs$n_with_af > 0
path <- if (use_af) afs$af_sumstats_path else file.path(FM, reg$sumstats_path)
x <- fread(cmd = sprintf("awk -F'\\t' 'NR==1 || $1==\"22\"' %s", path))
if (!use_af) {
  pa <- fread(file.path(FM, "data/ld_ref/panel_af", tolower(reg$ancestry), "chr22.af.tsv.gz"),
              select = c("position", "bim_a1", "bim_a2", "af_a1"))
  pa <- unique(pa, by = "position")
  x <- merge(x, pa, by = "position")
  x <- x[(allele1 == bim_a1 & allele2 == bim_a2) | (allele1 == bim_a2 & allele2 == bim_a1)]
  x[, af := af_a1]
}
x[, maf := pmin(af, 1 - af)]
x <- x[is.finite(se) & se > 0 & maf >= 0.01]
nvx <- 2 * reg$N_tot * x$maf * (1 - x$maf)
cf <- coef(lm(nvx ~ I(1 / x$se^2) - 1))[[1]]
cat(sprintf("%s\t%s\tN=%d\tn_snps=%d\tsdY_reg=%.3f\tsdY_median=%.3f\n", s,
            if (use_af) "study_af" else "panel_af", reg$N_tot, nrow(x), sqrt(cf), median(x$se * sqrt(nvx))))
