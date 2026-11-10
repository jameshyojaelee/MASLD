#!/usr/bin/env Rscript
# U4 (c): estimate the trait SD (sdY) implied by each quantitative GWAS's summary
# statistics, and by the Broadaway eQTL table, to test the sdY = 1 assumption in
# 06_susie_coloc.R (d1$sdY <- 1 for quant GWAS; d2$sdY = 1 for eQTL).
#
# Model: var(beta_hat) ~= sdY^2 / (2 N maf (1 - maf)).
#   reg   : coloc:::sdY.est (lm(2 N maf(1-maf) ~ 1/varbeta - 1), sqrt of slope)
#   median: median(se * sqrt(2 N maf (1 - maf)))
# MAF: the study's own allele frequency (data/sumstats_af/*_af.tsv, column af)
# when the AF sidecar has one; otherwise 1000G panel AF of the study ancestry
# (data/ld_ref/panel_af/<pop>/chr*.af.tsv.gz), joined on position + allele pair.
# N: registry N_tot (the value coloc receives). The GWAS files carry no per-SNP N.
suppressPackageStartupMessages({ library(data.table) })
set.seed(20261023)
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
OUT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/sdy"
N_SAMPLE <- 300000L

reg <- fread(file.path(FM, "config/gwas_registry.tsv"))
afs <- fread(file.path(FM, "config/gwas_af_sources.tsv"), select = c("study_name", "tier", "mode", "af_sumstats_path", "n_with_af"))
reg <- merge(reg, afs, by = "study_name", all.x = TRUE)
q <- reg[trait_type == "quantitative"]
cat("quantitative studies:", nrow(q), "\n")

est <- function(se, maf, n) {
  ok <- is.finite(se) & se > 0 & is.finite(maf) & maf >= 0.01 & maf <= 0.5 & is.finite(n)
  se <- se[ok]; maf <- maf[ok]; n <- n[ok]
  nvx <- 2 * n * maf * (1 - maf)
  oneover <- 1 / se^2
  cf <- coef(lm(nvx ~ oneover - 1))[["oneover"]]
  ok5 <- maf >= 0.05
  cf5 <- coef(lm(nvx[ok5] ~ oneover[ok5] - 1))[[1]]
  c(n_snps = length(se),
    sdY_reg = if (cf > 0) sqrt(cf) else NA_real_,
    sdY_reg_maf05 = if (cf5 > 0) sqrt(cf5) else NA_real_,
    sdY_median = median(se * sqrt(nvx)),
    sdY_median_q25 = unname(quantile(se * sqrt(nvx), 0.25)),
    sdY_median_q75 = unname(quantile(se * sqrt(nvx), 0.75)))
}

panel_cache <- list()
panel_maf <- function(pop, ss) {
  # ss: data.table chromosome, position, allele1, allele2 ; returns maf vector
  ss[, rid := .I]
  res <- rep(NA_real_, nrow(ss))
  for (ch in sort(unique(ss$chromosome))) {
    key <- paste(pop, ch)
    p <- file.path(FM, "data/ld_ref/panel_af", tolower(pop), sprintf("chr%s.af.tsv.gz", ch))
    if (!file.exists(p)) next
    pa <- fread(p, select = c("position", "bim_a1", "bim_a2", "af_a1"), showProgress = FALSE)
    pa <- unique(pa, by = "position")
    sub <- ss[chromosome == ch]
    m <- merge(sub, pa, by = "position")
    m <- m[(toupper(allele1) == toupper(bim_a1) & toupper(allele2) == toupper(bim_a2)) |
           (toupper(allele1) == toupper(bim_a2) & toupper(allele2) == toupper(bim_a1))]
    res[m$rid] <- pmin(m$af_a1, 1 - m$af_a1)
  }
  res
}

rows <- list()
for (k in seq_len(nrow(q))) {
  r <- q[k]
  use_af <- !is.na(r$n_with_af) && r$n_with_af > 0 && file.exists(r$af_sumstats_path)
  path <- if (use_af) r$af_sumstats_path else file.path(FM, r$sumstats_path)
  cols <- c("chromosome", "position", "allele1", "allele2", "beta", "se", if (use_af) "af")
  t0 <- Sys.time()
  ss <- fread(path, select = cols, showProgress = FALSE)
  n_all <- nrow(ss)
  ss <- ss[is.finite(se) & se > 0 & chromosome %in% 1:22]
  idx <- sort(sample.int(nrow(ss), min(N_SAMPLE * 2L, nrow(ss))))
  ss <- ss[idx]
  if (use_af) {
    ss[, maf := pmin(af, 1 - af)]
    maf_src <- "study_af"
  } else {
    ss[, maf := panel_maf(r$ancestry, ss)]
    maf_src <- paste0("1000G_", tolower(r$ancestry), "_panel")
  }
  ss <- ss[is.finite(maf) & maf >= 0.01]
  if (nrow(ss) > N_SAMPLE) ss <- ss[sort(sample.int(nrow(ss), N_SAMPLE))]
  e <- est(ss$se, ss$maf, rep(r$N_tot, nrow(ss)))
  rows[[k]] <- data.table(study = r$study_name, tier = r$tier, ancestry = r$ancestry,
                          N_tot = r$N_tot, maf_source = maf_src, n_variants_file = n_all,
                          t(e), median_abs_beta = median(abs(ss$beta)), median_se = median(ss$se))
  cat(sprintf("%-28s %s N=%d sdY_reg=%.3f sdY_med=%.3f (%s, %.0fs)\n", r$study_name, maf_src, r$N_tot,
              e[["sdY_reg"]], e[["sdY_median"]], format(Sys.time()), as.numeric(Sys.time() - t0, units = "secs")))
  rm(ss); gc(verbose = FALSE)
}

# eQTL side (d2: N = EQTL_N = 1183, sdY = 1). Per-SNP N column is available.
# Sampled per chromosome (21 GB of text in total), 15,000 gene-SNP rows each.
eq_n_all <- 0
eq <- rbindlist(lapply(1:22, function(ch) {
  x <- fread(file.path(dirname(dirname(FM)), "data/broadaway_eqtl", sprintf("chr%d_marginal_summary_results.tsv", ch)),
             select = c("CHR", "POS", "EAF", "Beta", "SE", "N", "ENSG"), showProgress = FALSE)
  eq_n_all <<- eq_n_all + nrow(x)
  x <- x[is.finite(SE) & SE > 0 & is.finite(EAF)]
  x <- x[sort(sample.int(nrow(x), min(15000L, nrow(x))))]
  gc(verbose = FALSE)
  x
}))
eq[, maf := pmin(EAF, 1 - EAF)]
e1 <- est(eq$SE, eq$maf, eq$N)
e2 <- est(eq$SE, eq$maf, rep(1183, nrow(eq)))
rows[[length(rows) + 1]] <- data.table(study = "Broadaway_eQTL(perSNP_N)", tier = NA, ancestry = "EUR",
  N_tot = median(eq$N), maf_source = "eqtl_EAF", n_variants_file = eq_n_all, t(e1),
  median_abs_beta = median(abs(eq$Beta)), median_se = median(eq$SE))
rows[[length(rows) + 1]] <- data.table(study = "Broadaway_eQTL(N=1183)", tier = NA, ancestry = "EUR",
  N_tot = 1183, maf_source = "eqtl_EAF", n_variants_file = eq_n_all, t(e2),
  median_abs_beta = median(abs(eq$Beta)), median_se = median(eq$SE))

out <- rbindlist(rows, fill = TRUE)
out[, outside_0.7_1.4 := sdY_reg < 0.7 | sdY_reg > 1.4]
fwrite(out, file.path(OUT, "sdy_by_study_raw.tsv"), sep = "\t")
print(out)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("seed 20261023; N_SAMPLE", N_SAMPLE, "\n")
