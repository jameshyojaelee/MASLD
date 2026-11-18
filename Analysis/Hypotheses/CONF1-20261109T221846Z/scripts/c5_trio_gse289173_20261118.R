#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
EXT_DIR <- if (length(args) >= 1L) args[1] else file.path(ROOT, "Analysis/Hypotheses/CONF1-20261109T221846Z/out_c5_extract_20261118")
OUT_DIR <- if (length(args) >= 2L) args[2] else file.path(ROOT, "Analysis/Hypotheses/CONF1-20261109T221846Z/out_c5_trio_20261118")
GENES_ARG <- if (length(args) >= 3L) args[3] else "TRIO"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (GENES_ARG == "TRIO") {
  GENES <- c("GNMT", "MAT1A", "CYP2C19")
  DIG1_EFFECT <- c(GNMT = 1.89, MAT1A = 0.82, CYP2C19 = 1.72)
  IS_TRIO <- TRUE
} else {
  ex_tab <- fread(file.path(ROOT, "scripts/figures/example_genes.tsv"))
  GENES <- as.character(ex_tab[order(order), symbol])
  DIG1_EFFECT <- setNames(rep(NA_real_, length(GENES)), GENES)
  IS_TRIO <- FALSE
}
NONHEP <- c("COL1A1", "COL3A1", "DCN", "PECAM1", "STAB2", "PTPRC", "CD163", "LYZ")
SEED <- 20261109L
B <- 100000L

to_bool <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")

donors_all <- fread(file.path(EXT_DIR, "donors.tsv"))
donors_all[, sex_concordant := to_bool(sex_concordant)]
donors_all[, exclude_sex_discordant := to_bool(exclude_sex_discordant)]
donors_all[, disease := ifelse(group == "ctrl", 0L, ifelse(group %in% c("eMASH", "aMASH"), 1L, NA_integer_))]
donors_all[, series_fl := as.integer(series == "FL_series")]

pb_dt <- fread(cmd = sprintf("zcat %s", shQuote(file.path(EXT_DIR, "pseudobulk_hep.tsv.gz"))), data.table = FALSE)
syms <- as.character(pb_dt[[1]])
stopifnot(!anyDuplicated(syms))
pb_mat <- as.matrix(pb_dt[, -1, drop = FALSE])
storage.mode(pb_mat) <- "numeric"
rownames(pb_mat) <- syms

# Eligible donors >= 30 hep nuclei (primary excludes sex-discordant GSM8785704)
elig_all47 <- donors_all[n_hep >= 30L]
elig_conc46 <- donors_all[n_hep >= 30L & !exclude_sex_discordant]

compute_logcpm <- function(don_dt) {
  sub_mat <- pb_mat[, don_dt$donor, drop = FALSE]
  dge <- DGEList(counts = sub_mat)
  dge <- calcNormFactors(dge, method = "TMM")
  cpm(dge, log = TRUE, prior.count = 1)
}

lcpm_conc46 <- compute_logcpm(elig_conc46)
lcpm_all47 <- compute_logcpm(elig_all47)

# Gate check on primary concordant donors (excl GSM8785704)
gate_rows <- lapply(GENES, function(g) {
  pres <- g %in% rownames(pb_mat)
  df <- if (pres) mean(pb_mat[g, elig_conc46$donor] > 0) else 0
  data.table(gene = g, present = pres, detect_frac = df, testable = pres && df >= 0.80)
})
gate_dt <- rbindlist(gate_rows)
fwrite(gate_dt, file.path(OUT_DIR, "gene_gate.tsv"), sep = "\t")
fwrite(elig_conc46, file.path(OUT_DIR, "donors_eligible.tsv"), sep = "\t")

# Save long donor_logcpm for primary cohort & testable genes
dl_rows <- rbindlist(lapply(gate_dt[testable == TRUE, gene], function(g) {
  data.table(gene = g, donor = elig_conc46$donor, group = elig_conc46$group,
             series = elig_conc46$series, logcpm = as.numeric(lcpm_conc46[g, elig_conc46$donor]))
}))
fwrite(dl_rows, file.path(OUT_DIR, "donor_logcpm.tsv"), sep = "\t")

# HC3 + Freedman-Lane studentized t-permutation core (verified against V5 05_e6.R)
fit_hc3_obs <- function(y, X, j = 2L) {
  if (nrow(X) <= ncol(X) + 1L || qr(X)$rank < ncol(X)) return(NULL)
  XtXi <- solve(crossprod(X))
  Cm <- XtXi %*% t(X)
  b <- as.vector(Cm %*% y)
  e <- as.vector(y - X %*% b)
  h <- rowSums((X %*% XtXi) * X)
  if (any(1.0 - h <= 1e-8)) return(NULL)
  W <- (e^2) / ((1.0 - h)^2)
  se <- sqrt(sum((Cm[j, ]^2) * W))
  if (!is.finite(se) || se <= 0) return(NULL)
  df <- nrow(X) - ncol(X)
  tc <- qt(0.975, df = df)
  list(beta = b[j], se = se, t = b[j] / se, df = df,
       ci_low = b[j] - tc * se, ci_high = b[j] + tc * se,
       Cm = Cm, h = h, max_h = max(h))
}

fl_null_t <- function(y, X, obs, P_idx, j = 2L, chunk = 5000L) {
  Z <- X[, -j, drop = FALSE]
  gamma0 <- solve(crossprod(Z), crossprod(Z, y))
  f0 <- as.vector(Z %*% gamma0)
  r0 <- as.vector(y - f0)
  n <- length(y)
  B_n <- ncol(P_idx)
  cm_j2 <- obs$Cm[j, ]^2
  inv_1mh2 <- 1.0 / ((1.0 - obs$h)^2)
  b_null <- numeric(B_n)
  se_null <- numeric(B_n)
  t_null <- numeric(B_n)
  for (s in seq.int(1L, B_n, by = chunk)) {
    e_i <- min(B_n, s + chunk - 1L)
    idx <- s:e_i
    Y_star <- matrix(r0[P_idx[, idx, drop = FALSE]], nrow = n, ncol = length(idx)) + f0
    B_star <- obs$Cm %*% Y_star
    E_star <- Y_star - X %*% B_star
    W_star <- (E_star^2) * inv_1mh2
    se_v <- sqrt(as.vector(crossprod(cm_j2, W_star)))
    b_v <- as.vector(B_star[j, ])
    b_null[idx] <- b_v
    se_null[idx] <- se_v
    t_null[idx] <- b_v / se_v
  }
  list(beta = b_null, se = se_null, t = t_null)
}

run_arm <- function(don_sub, lcpm_mat, x_var, formula_str, gene_vec, arm_name, seed = SEED, B_draws = B) {
  # Build model matrix with x_var forced into column 2 (after Intercept)
  rhs_vars <- trimws(strsplit(formula_str, "\\+")[[1]])
  stopifnot(rhs_vars[1] == x_var)
  cov_mat <- if (length(rhs_vars) > 1L) as.matrix(don_sub[, rhs_vars[-1], with = FALSE]) else NULL
  X <- cbind(Intercept = 1.0, target = as.numeric(don_sub[[x_var]]))
  if (!is.null(cov_mat)) X <- cbind(X, cov_mat)
  stopifnot(qr(X)$rank == ncol(X))

  set.seed(seed)
  n <- nrow(don_sub)
  P_idx <- replicate(B_draws, sample.int(n))

  n0 <- if (x_var == "disease") sum(don_sub[[x_var]] == 0L) else sum(don_sub[[x_var]] == 0L)
  n1 <- if (x_var == "disease") sum(don_sub[[x_var]] == 1L) else sum(don_sub[[x_var]] > 0L)

  out_rows <- lapply(gene_vec, function(g) {
    if (!(g %in% rownames(lcpm_mat))) return(NULL)
    y <- as.numeric(lcpm_mat[g, don_sub$donor])
    obs <- fit_hc3_obs(y, X, j = 2L)
    if (is.null(obs)) return(NULL)
    nl <- fl_null_t(y, X, obs, P_idx, j = 2L)
    ok <- is.finite(nl$t)
    p_perm <- (1.0 + sum(abs(nl$t[ok]) >= abs(obs$t) * (1.0 - 1e-12))) / (sum(ok) + 1.0)
    nt_sd <- sd(nl$t[ok])
    nb_sd <- sd(nl$beta[ok])
    data.table(
      analysis = arm_name,
      gene = g,
      n = n,
      n0 = n0,
      n1 = n1,
      beta = obs$beta,
      se = obs$se,
      t_obs = obs$t,
      df = obs$df,
      max_leverage = obs$max_h,
      ci_low = obs$ci_low,
      ci_high = obs$ci_high,
      null_t_sd = nt_sd,
      null_beta_sd = nb_sd,
      n_draws = sum(ok),
      p_perm = p_perm,
      mde = 2.8 * obs$se * nt_sd,
      mde_t = 2.8 * nt_sd
    )
  })
  dt <- rbindlist(Filter(Negate(is.null), out_rows))
  dt[, holm := p.adjust(p_perm, method = "holm")]
  # Assert 0 genes have Holm < 0.05 while 95% HC3 CI includes 0
  bad <- sum(dt$holm < 0.05 & dt$ci_low <= 0 & dt$ci_high >= 0, na.rm = TRUE)
  stopifnot(bad == 0L)
  dt
}

# ---- Dilution guard across all genome-wide testable genes on primary arm ----
prim_don <- elig_conc46[!is.na(disease)]
test_all_g <- rownames(pb_mat)[rowMeans(pb_mat[, elig_conc46$donor, drop = FALSE] > 0) >= 0.80]
X_prim <- cbind(Intercept = 1.0, disease = as.numeric(prim_don$disease),
                female = as.numeric(prim_don$female), logn = log10(as.numeric(prim_don$n_hep)))
XtXi_p <- solve(crossprod(X_prim))
Cm_p <- XtXi_p %*% t(X_prim)
h_p <- rowSums((X_prim %*% XtXi_p) * X_prim)
Y_all <- t(lcpm_conc46[test_all_g, prim_don$donor, drop = FALSE]) # n x G
B_all <- Cm_p %*% Y_all
E_all <- Y_all - X_prim %*% B_all
W_all <- (E_all^2) / ((1.0 - h_p)^2)
se_all <- sqrt(as.vector(crossprod(Cm_p[2L, ]^2, W_all)))
ctrl_mean <- rowMeans(lcpm_conc46[test_all_g, prim_don[disease == 0L, donor], drop = FALSE])
all_beta_dt <- data.table(
  gene = test_all_g,
  beta = as.numeric(B_all[2L, ]),
  se = se_all,
  t_obs = as.numeric(B_all[2L, ]) / se_all,
  ctrl_mean_logcpm = as.numeric(ctrl_mean)
)
all_beta_dt[, decile := as.integer(cut(ctrl_mean_logcpm,
                                       breaks = quantile(ctrl_mean_logcpm, probs = seq(0, 1, 0.1), na.rm = TRUE),
                                       include.lowest = TRUE, labels = FALSE))]
fwrite(all_beta_dt, file.path(OUT_DIR, "all_gene_betas.tsv.gz"), sep = "\t")
med_all <- median(all_beta_dt$beta, na.rm = TRUE)
dec_q20_b <- all_beta_dt[, .(decile_p20 = quantile(beta, 0.20, na.rm = TRUE),
                             decile_p20_t = quantile(t_obs, 0.20, na.rm = TRUE)), by = decile]

# ---- Primary arm (excl sex-discordant GSM8785704, n=40) ----
prep_covs <- function(dt, use_corrected_sex = FALSE) {
  d <- copy(dt)
  d[, log10_n_hep := log10(as.numeric(n_hep))]
  if (use_corrected_sex) d[, female := female_corrected]
  d
}
elig_conc46 <- prep_covs(elig_conc46, FALSE)
elig_all47 <- prep_covs(elig_all47, TRUE)

prim <- run_arm(elig_conc46[!is.na(disease)], lcpm_conc46, "disease", "disease + female + log10_n_hep",
                gate_dt[testable == TRUE, gene], "primary_excl_GSM8785704")
# Attach dilution guard & decision
prim[, decile := all_beta_dt[match(prim$gene, all_beta_dt$gene), decile]]
prim[, decile_p20 := dec_q20_b[match(prim$decile, dec_q20_b$decile), decile_p20]]
prim[, decile_p20_t := dec_q20_b[match(prim$decile, dec_q20_b$decile), decile_p20_t]]
prim[, frac_decile_below := vapply(seq_len(.N), function(i) {
  sub <- all_beta_dt[decile == prim$decile[i], beta]
  mean(sub <= prim$beta[i], na.rm = TRUE)
}, 0.0)]
prim[, frac_decile_below_t := vapply(seq_len(.N), function(i) {
  sub <- all_beta_dt[decile == prim$decile[i], t_obs]
  mean(sub <= prim$t_obs[i], na.rm = TRUE)
}, 0.0)]
prim[, passes_guard := beta <= decile_p20]
prim[, dig1_effect := DIG1_EFFECT[gene]]
prim[, median_beta_all := med_all]
prim[, decision := if (IS_TRIO) ifelse(beta < 0 & holm < 0.05 & passes_guard, "replicated",
                                       ifelse(mde <= dig1_effect, "not_replicated", "uninformative")) else "not_applicable"]
stopifnot(sum(prim$holm < 0.05 & prim$ci_low <= 0 & prim$ci_high >= 0, na.rm = TRUE) == 0L)
fwrite(prim, file.path(OUT_DIR, "primary.tsv"), sep = "\t")

# ---- Sensitivity arms (descriptive + series + sex correction) ----
tgenes <- gate_dt[testable == TRUE, gene]
sens_list <- list(
  run_arm(elig_all47[!is.na(disease)], lcpm_all47, "disease", "disease + female + log10_n_hep",
          tgenes, "GSM8785704_corrected_sex"),
  run_arm(elig_conc46[series == "FL_series" & !is.na(disease)], lcpm_conc46, "disease", "disease + female + log10_n_hep",
          tgenes, "FL_series_only_excl_GSM8785704"),
  run_arm(elig_all47[series == "FL_series" & !is.na(disease)], lcpm_all47, "disease", "disease + female + log10_n_hep",
          tgenes, "FL_series_only_corrected_fl723"),
  run_arm(elig_conc46[!is.na(disease)], lcpm_conc46, "disease", "disease + series_fl + female + log10_n_hep",
          tgenes, "series_adjusted_excl_GSM8785704"),
  run_arm(elig_all47[!is.na(disease)], lcpm_all47, "disease", "disease + series_fl + female + log10_n_hep",
          tgenes, "series_adjusted_corrected_fl723"),
  run_arm(elig_conc46, lcpm_conc46, "stage_ord", "stage_ord + female + log10_n_hep",
          tgenes, "ordinal_stage_all_excl_GSM8785704"),
  run_arm(elig_conc46[series == "FL_series"], lcpm_conc46, "stage_ord", "stage_ord + female + log10_n_hep",
          tgenes, "ordinal_stage_FL_series_only"),
  run_arm(elig_conc46[!is.na(disease) & n_hep >= 100L], lcpm_conc46, "disease", "disease + female + log10_n_hep",
          tgenes, "min_100_hep_excl_GSM8785704"),
  run_arm(elig_conc46[!is.na(disease)], lcpm_conc46, "disease", "disease + female",
          tgenes, "no_nuclei_cov_excl_GSM8785704"),
  run_arm(elig_conc46[!is.na(disease)], lcpm_conc46, "disease", "disease + female + log10_n_hep",
          NONHEP[NONHEP %in% rownames(lcpm_conc46)], "nonhep_markers_in_hep_pseudobulk")
)
sens_dt <- rbindlist(sens_list, fill = TRUE)
# Verify 0 discordant Holm < 0.05 vs 95% HC3 CI across every arm
stopifnot(sum(sens_dt$holm < 0.05 & sens_dt$ci_low <= 0 & sens_dt$ci_high >= 0, na.rm = TRUE) == 0L)
fwrite(sens_dt, file.path(OUT_DIR, "descriptive.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("DONE c5_trio_gse289173_20261118\n")
