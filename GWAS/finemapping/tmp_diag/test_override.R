suppressPackageStartupMessages(library(data.table))
Sys.setenv(UNIFORM35_RUN_ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/runs/uniform35_v4_2026-08-04")
source("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src/seqfunc_v2/finemap/common.R")
cat("AF_PAL_MARGIN =", AF_PAL_MARGIN, " AF_CERT_OVERRIDE_MARGIN =", AF_CERT_OVERRIDE_MARGIN, "\n\n")
# the three real cases that decided the calibration
af    <- c(0.2765, 0.9985, 0.9969)          # rs12982412 (WRONG flip), rs6496572, rs6709525
panel <- c(0.6280, 0.0171, 0.0252)          # panel freq under the same-strand reading
r <- resolve_palindrome_af(af, panel)
m <- attr(r, "margin")
cat("verdicts :", paste(as.character(r), collapse=", "), "\n")
cat("margins  :", paste(round(m,3), collapse=", "), "\n")
cat("overrides certificate (margin >= 0.50)?", paste(m >= AF_CERT_OVERRIDE_MARGIN, collapse=", "), "\n\n")
stopifnot(all(as.character(r) == "opposite"))
stopifnot(!(m[1] >= AF_CERT_OVERRIDE_MARGIN))   # rs12982412 must NO LONGER override
stopifnot(m[2] >= AF_CERT_OVERRIDE_MARGIN, m[3] >= AF_CERT_OVERRIDE_MARGIN)
# attribute must survive the data.table := assignment used in 02
d <- data.table(af = af, panel_af_same = panel, palindromic = TRUE, match_class = "direct")
d[, c("af_call","af_margin") := { rr <- resolve_palindrome_af(af, panel_af_same)
                                  list(as.character(rr), as.numeric(attr(rr,"margin"))) }]
print(d[, .(af, panel_af_same, af_call, af_margin = round(af_margin,3))])
stopifnot(all(is.finite(d$af_margin)))
cat("\nRESULT: rs12982412 (the wrongly-flipped PIP=1 anchor) is now REJECTED;\n")
cat("        rs6496572 and rs6709525 (gnomAD-confirmed genuine flips) still override.\n")
