# Gate L10-02 pilot: positive-control locus (PNPLA3, chr22:44324727, hg19) fine-mapped with
# SuSiE-RSS exactly as the pipeline's run_susie() does it (ridge 1e-6, coverage 0.95, default L,
# set.seed(42); binary branch z + N_eff; quantitative branch bhat/shat + N), but without writing
# into the repo. LD read read-only through the pipeline's get_ld_per_locus() (copied file).
suppressPackageStartupMessages({library(susieR)})
source("finemapping_functions.R")
Sys.setenv(MASLD_PROJECT_ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
CHR <- 22; BP <- 44324727; START <- BP - 5e5; END <- BP + 5e5; LOCUS <- "22.44324727"
out <- "pilot_results.tsv"; cs_out <- "pilot_credible_sets.tsv"
cat("panel\tstudy\tn_variants\tlambda_s\tn_cs\tpip_44324727\tin_cs_44324727\tcs_size_44324727\ttop_variant\ttop_pip\ttop_abs_z\tz_44324727\tseconds\n", file = out)
cat("panel\tstudy\tcs\tsize\tpurity_min_abs_corr\tlead\tlead_pip\tlead_abs_z\n", file = cs_out)
run_one <- function(panel, study, binary, N_tot, N_cases) {
  t0 <- Sys.time()
  Sys.setenv(LD_PANEL = panel)
  ss <- fread(sprintf("ss_%s_chr22_pnpla3.tsv", study))
  ss <- ss[!is.na(beta) & !is.na(se) & se > 0]
  dat <- get_ld_per_locus(ss, LOCUS, CHR, START, END, ancestry = "EUR")
  if (is.null(dat)) { cat(panel, study, "LD extraction failed\n"); return(invisible(NULL)) }
  s <- dat[[1]]; R <- as.matrix(dat[[2]]) + 1e-6 * diag(nrow(dat[[2]]))
  z <- s$beta / s$se
  N_eff <- if (binary) 4 / (1 / N_cases + 1 / (N_tot - N_cases)) else N_tot
  lam <- tryCatch(as.numeric(susieR::estimate_s_rss(z = z, R = R, n = N_eff)), error = function(e) NA_real_)
  set.seed(42)
  fit <- if (binary) susie_rss(R = R, n = N_eff, z = z, coverage = 0.95)
         else susie_rss(R = R, n = N_eff, bhat = s$beta, shat = s$se, coverage = 0.95)
  pip <- fit$pip; cs <- fit$sets$cs
  i <- which(s$position == BP)
  in_cs <- if (length(i)) any(vapply(cs, function(x) i[1] %in% x, logical(1))) else NA
  cs_size <- if (length(i) && isTRUE(in_cs)) length(cs[[which(vapply(cs, function(x) i[1] %in% x, logical(1)))[1]]]) else NA
  top <- which.max(pip)
  cat(sprintf("%s\t%s\t%d\t%.4f\t%d\t%s\t%s\t%s\t%s\t%.4f\t%.2f\t%s\t%.0f\n", panel, study, nrow(s), lam, length(cs),
              if (length(i)) sprintf("%.4f", pip[i[1]]) else "absent", in_cs, cs_size, s$SNP[top], pip[top], abs(z[top]),
              if (length(i)) sprintf("%.2f", z[i[1]]) else "NA", as.numeric(difftime(Sys.time(), t0, units = "secs"))),
      file = out, append = TRUE)
  pur <- fit$sets$purity
  for (k in seq_along(cs)) {
    m <- cs[[k]][which.max(pip[cs[[k]]])]
    cat(sprintf("%s\t%s\t%s\t%d\t%.3f\t%s\t%.4f\t%.2f\n", panel, study, names(cs)[k], length(cs[[k]]),
                if (!is.null(pur)) pur[k, "min.abs.corr"] else NA, s$SNP[m], pip[m], abs(z[m])),
        file = cs_out, append = TRUE)
  }
  cat("done", panel, study, "\n")
}
run_one("polyfun", "MVP_NAFLD_EUR", TRUE, 432709, 29342)
run_one("ukbb", "MVP_NAFLD_EUR", TRUE, 432709, 29342)
run_one("polyfun", "UKBB_ALT", FALSE, 343850, NA)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
