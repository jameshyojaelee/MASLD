#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CONF1_C2 <- file.path(ROOT, "Analysis/Hypotheses/CONF1-20261109T221846Z/out_c2_20261118")
CONF1_SCRIPT <- file.path(ROOT, "Analysis/Hypotheses/CONF1-20261109T221846Z/scripts/c2_methionine_hepatocyte_20261118.R")
XOUT <- file.path(ROOT, "Analysis/Hypotheses/XMOD2-20261113T151052Z/out_20261118")

out_trio <- file.path(XOUT, "trio_c2")
out_ex <- file.path(XOUT, "examples_c2")
dir.create(out_trio, recursive = TRUE, showWarnings = FALSE)
dir.create(out_ex, recursive = TRUE, showWarnings = FALSE)

rs_bin <- R.home("bin/Rscript")
st1 <- system2(rs_bin, c(CONF1_SCRIPT, shQuote(out_trio), "TRIO", "XMOD2_TRIO_C2"))
stopifnot(st1 == 0L)
st2 <- system2(rs_bin, c(CONF1_SCRIPT, shQuote(out_ex), "EXAMPLES", "XMOD2_EXAMPLES_C2"))
stopifnot(st2 == 0L)

# Explicit regression check: XMOD2 trio_c2 must match CONF1 out_c2_20261118 identically
if (file.exists(file.path(CONF1_C2, "combined.tsv"))) {
  a <- fread(file.path(out_trio, "combined.tsv"))
  b <- fread(file.path(CONF1_C2, "combined.tsv"))
  stopifnot(nrow(a) == nrow(b),
            max(abs(a$p_perm - b$p_perm)) < 1e-12,
            max(abs(a$T - b$T)) < 1e-12)
}

ex_pc <- fread(file.path(out_ex, "per_cohort.tsv"))
ex_cb <- fread(file.path(out_ex, "combined.tsv"))
ex_sn <- fread(file.path(out_ex, "sensitivity.tsv"))
stopifnot(sum(ex_pc$p_perm_holm_within_cohort < 0.05 & ex_pc$perm_ci_low <= 0 & ex_pc$perm_ci_high >= 0, na.rm = TRUE) == 0L)
stopifnot(sum(ex_cb$holm < 0.05 & ex_cb$ci_low <= 0 & ex_cb$ci_high >= 0, na.rm = TRUE) == 0L)
stopifnot(sum(ex_sn$holm < 0.05 & ex_sn$ci_low <= 0 & ex_sn$ci_high >= 0, na.rm = TRUE) == 0L)
cat("DONE XMOD2 C2 (trio_c2 + examples_c2) & CONF1 regression check passed\n")
