#!/usr/bin/env Rscript
# =============================================================================
# mvp_power_mde.R
#
# Per-ancestry MINIMUM DETECTABLE EFFECT (MDE) / power framework for the MVP
# multi-ancestry liver/NAFLD GWAS strata.
#
# PURPOSE
#   The MVP EUR/AFR/AMR/EAS strata differ by >50x in case count (NAFLD: EUR
#   29,342 vs EAS 505). For the cross-ancestry COLOC interpretation we must NOT
#   read "no signal in AFR/AMR/EAS" as biological absence when the stratum was
#   simply underpowered. This script computes, for each stratum, the smallest
#   effect detectable at 80% power and the (ancestry-specific Kanai 2016)
#   genome-wide significance threshold, as a function of risk-allele frequency
#   (RAF). An observed null is then interpretable as a TRUE negative only for
#   effects at or above that stratum's MDE; below it, the null is a power
#   ceiling, not absence.
#
# STATISTICAL MODEL  (standard single-variant additive score/Wald test)
#   Under an additive (per-allele) model the asymptotic non-centrality
#   parameter (NCP) of the 1-df chi-square association statistic is
#
#       NCP = beta^2 * Var(G) * N_eff
#
#   where beta is the per-allele effect on the relevant scale, Var(G) =
#   2*f*(1-f) is the genotype variance for a SNP at allele frequency f under
#   Hardy-Weinberg, and N_eff is the effective sample size.
#     - QUANTITATIVE (inverse-normal-transformed) trait, beta in trait-SD units:
#           N_eff = N_tot
#           NCP   = beta^2 * 2*f*(1-f) * N_tot
#       (Sham & Purcell 2014, Nat Rev Genet; the IVNT makes residual variance ~1
#        so beta is already in SD units and the 1-beta^2*Var(G) term ~= 1.)
#     - BINARY case-control trait, beta = per-allele log-odds-ratio:
#           N_eff = N_cases * N_controls / N_tot         (balanced-equivalent N)
#           NCP   = beta^2 * 2*f*(1-f) * (N_cases*N_controls/N_tot)
#       (logistic score-test NCP; the N_cases*N_controls/N_tot factor is the
#        standard case-control effective N, e.g. Pirinen 2013, Willer METAL
#        weighting; equivalent to N * phi * (1-phi) with phi = case fraction.)
#
#   To reach power = 0.80 at two-sided significance alpha the required NCP is
#       NCP_req = ( qnorm(1 - alpha/2) + qnorm(power) )^2
#   Solving NCP = NCP_req for beta gives the minimum detectable effect:
#       beta_min = sqrt( NCP_req / (2*f*(1-f) * N_eff) )
#       MDE_OR   = exp(beta_min)     (binary)
#       MDE_beta = beta_min          (quantitative, in trait-SD units)
#
# ALPHA (Kanai et al. 2016, J Hum Genet — ancestry-specific genome-wide
#        significance from LD-pruned effective number of tests):
#       EUR 9.26e-8 | AFR 3.24e-8 | EAS 1.61e-7 | AMR 1.83e-7
#   (AFR is the most stringent because its lower LD => more independent tests;
#    EAS/AMR are more lenient.)
#
# OUTPUT
#   GWAS/MR_Data/MVP/mvp_power_mde.csv  — tidy (stratum x RAF) MDE table
#   + per-ancestry MDE summary at RAF=0.2 printed to stdout.
#
# Author: cross-ancestry power module
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
proj_root <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
mvp_dir  <- file.path(proj_root, "GWAS", "MR_Data", "MVP")
in_tsv   <- file.path(mvp_dir, "mvp_sample_sizes.tsv")
out_csv  <- file.path(mvp_dir, "mvp_power_mde.csv")

stopifnot(file.exists(in_tsv))

# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
power_target <- 0.80

# Ancestry-specific genome-wide significance thresholds (Kanai 2016).
alpha_kanai <- c(
  EUR = 9.26e-8,
  AFR = 3.24e-8,
  EAS = 1.61e-7,
  AMR = 1.83e-7
)

# Risk-allele-frequency grid spanning low-frequency to common variants.
raf_grid <- c(0.01, 0.05, 0.10, 0.20, 0.30, 0.50)

# ----------------------------------------------------------------------------
# Core MDE functions
# ----------------------------------------------------------------------------

# Required non-centrality parameter for a 1-df chi-square test to reach
# `power` at two-sided significance `alpha`.
ncp_required <- function(alpha, power) {
  (qnorm(1 - alpha / 2) + qnorm(power))^2
}

# Minimum detectable per-allele effect (beta) on whatever scale N_eff/Var(G)
# imply. Returns beta in the natural scale (log-OR for binary, SD for quant).
mde_beta <- function(N_eff, raf, alpha, power) {
  varG    <- 2 * raf * (1 - raf)
  ncp_req <- ncp_required(alpha, power)
  sqrt(ncp_req / (varG * N_eff))
}

# ----------------------------------------------------------------------------
# Load strata
# ----------------------------------------------------------------------------
dt <- fread(in_tsv, na.strings = c("", "NA"))

# Effective sample size:
#   quantitative -> N_tot
#   binary       -> N_cases * N_controls / N_tot
dt[, N_eff := ifelse(
  trait_type == "binary",
  as.numeric(N_cases) * as.numeric(N_controls) / as.numeric(N_tot),
  as.numeric(N_tot)
)]

# Map each stratum to its ancestry-specific Kanai alpha.
dt[, alpha_kanai := alpha_kanai[ancestry]]
if (anyNA(dt$alpha_kanai)) {
  stop("Unmapped ancestry in input (no Kanai alpha): ",
       paste(unique(dt$ancestry[is.na(dt$alpha_kanai)]), collapse = ", "))
}

# ----------------------------------------------------------------------------
# Build the (stratum x RAF) MDE table
# ----------------------------------------------------------------------------
grid <- CJ(study = dt$study, RAF = raf_grid, sorted = FALSE, unique = TRUE)
res  <- merge(grid, dt, by = "study", all.x = TRUE)

res[, beta_min := mapply(mde_beta, N_eff, RAF, alpha_kanai, power_target)]

# MDE on the reportable scale: OR for binary, SD-beta for quantitative.
res[, MDE_OR_or_beta := ifelse(trait_type == "binary",
                               exp(beta_min),  # odds ratio
                               beta_min)]      # per-allele SD

res[, power_target := power_target]

out <- res[, .(
  stratum = study,
  ancestry,
  trait_type,
  N_eff = round(N_eff, 1),
  RAF,
  alpha_kanai,
  MDE_OR_or_beta = round(MDE_OR_or_beta, 4),
  power_target
)]
setorder(out, stratum, RAF)

fwrite(out, out_csv)
cat("Wrote MDE table:", out_csv, "  (", nrow(out), "rows )\n\n")

# ----------------------------------------------------------------------------
# Sanity check: EUR NAFLD at RAF=0.3 should give a modest OR (~1.04-1.07)
# ----------------------------------------------------------------------------
chk <- out[stratum == "MVP_NAFLD_EUR" & RAF == 0.30]
if (nrow(chk) == 1) {
  cat(sprintf(
    "SANITY CHECK | EUR NAFLD (N_eff=%.0f) at RAF=0.30, alpha=%.2e:\n  MDE OR = %.4f  (expected ~1.04-1.07 for a study of this size) -> %s\n\n",
    chk$N_eff, chk$alpha_kanai, chk$MDE_OR_or_beta,
    ifelse(chk$MDE_OR_or_beta >= 1.04 & chk$MDE_OR_or_beta <= 1.07,
           "PASS", "CHECK")
  ))
}

# ----------------------------------------------------------------------------
# Per-ancestry summary at a "typical" common-variant RAF = 0.20
# ----------------------------------------------------------------------------
cat("================================================================\n")
cat("PER-STRATUM MDE @ 80% power, RAF=0.20, ancestry-specific Kanai alpha\n")
cat("  binary  -> smallest detectable per-allele ODDS RATIO\n")
cat("  quant   -> smallest detectable per-allele effect (trait SD units)\n")
cat("================================================================\n")

summ <- out[RAF == 0.20]
setorder(summ, trait_type, ancestry, stratum)
summ[, label := ifelse(trait_type == "binary",
                       sprintf("OR>=%.3f", MDE_OR_or_beta),
                       sprintf("beta>=%.3f SD", MDE_OR_or_beta))]
print(summ[, .(stratum, ancestry, trait_type,
               N_eff, alpha_kanai, MDE = label)],
      nrows = nrow(summ))

# Compact NAFLD-only readout (the headline cross-ancestry trait).
cat("\n---- NAFLD cross-ancestry MDE (RAF=0.20) ----\n")
naf <- summ[grepl("^MVP_NAFLD_", stratum)]
setorder(naf, MDE_OR_or_beta)
for (i in seq_len(nrow(naf))) {
  cat(sprintf("  %-4s NAFLD: N_cases-equiv N_eff=%7.0f  ->  detects OR>=%.3f at RAF 0.20\n",
              naf$ancestry[i], naf$N_eff[i], naf$MDE_OR_or_beta[i]))
}
cat("\nDone.\n")
