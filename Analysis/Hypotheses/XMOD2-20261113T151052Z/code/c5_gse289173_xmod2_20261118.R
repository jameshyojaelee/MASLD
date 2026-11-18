#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CONF1_EXT <- file.path(ROOT, "Analysis/Hypotheses/CONF1-20261109T221846Z/out_c5_extract_20261118")
CONF1_TRIO <- file.path(ROOT, "Analysis/Hypotheses/CONF1-20261109T221846Z/out_c5_trio_20261118")
CONF1_SCRIPT <- file.path(ROOT, "Analysis/Hypotheses/CONF1-20261109T221846Z/scripts/c5_trio_gse289173_20261118.R")
XOUT <- file.path(ROOT, "Analysis/Hypotheses/XMOD2-20261113T151052Z/out_20261118")

out_trio <- file.path(XOUT, "trio_c5")
out_ex <- file.path(XOUT, "examples_c5")
dir.create(out_trio, recursive = TRUE, showWarnings = FALSE)
dir.create(out_ex, recursive = TRUE, showWarnings = FALSE)

rs_bin <- R.home("bin/Rscript")
st1 <- system2(rs_bin, c(CONF1_SCRIPT, shQuote(CONF1_EXT), shQuote(out_trio), "TRIO"))
stopifnot(st1 == 0L)
st2 <- system2(rs_bin, c(CONF1_SCRIPT, shQuote(CONF1_EXT), shQuote(out_ex), "EXAMPLES"))
stopifnot(st2 == 0L)

# Explicit regression check: XMOD2 trio_c5 must match CONF1 out_c5_trio_20261118 identically
if (file.exists(file.path(CONF1_TRIO, "primary.tsv"))) {
  a <- fread(file.path(out_trio, "primary.tsv"))
  b <- fread(file.path(CONF1_TRIO, "primary.tsv"))
  stopifnot(nrow(a) == nrow(b), identical(a$gene, b$gene),
            max(abs(a$p_perm - b$p_perm)) < 1e-12,
            max(abs(a$t_obs - b$t_obs)) < 1e-12)
}

ex_p <- fread(file.path(out_ex, "primary.tsv"))
ex_s <- fread(file.path(out_ex, "descriptive.tsv"))
stopifnot(sum(ex_p$holm < 0.05 & ex_p$ci_low <= 0 & ex_p$ci_high >= 0, na.rm = TRUE) == 0L)
stopifnot(sum(ex_s$holm < 0.05 & ex_s$ci_low <= 0 & ex_s$ci_high >= 0, na.rm = TRUE) == 0L)
cat("DONE XMOD2 C5 (trio_c5 + examples_c5) & CONF1 regression check passed\n")
