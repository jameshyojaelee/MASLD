#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT_DIR <- if (length(args) >= 1L) args[1] else file.path(ROOT, "Analysis/Hypotheses/CONF1-20261109T221846Z/out_c2_20261118")
GENES_ARG <- if (length(args) >= 2L) args[2] else "TRIO"
MODE_LABEL <- if (length(args) >= 3L) args[3] else "CONF1_C2"
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

SEED <- 20261109L
B <- 100000L
to_bool <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")

# ---- Load exact cohort inputs ----
p_v5_cnt <- file.path(ROOT, "Analysis/Hypotheses/V5-20261021T122723Z/outputs/job3/pseudobulk_hep.tsv.gz")
p_v5_don <- file.path(ROOT, "Analysis/Hypotheses/V5-20261021T122723Z/outputs/donor_qc.tsv")
p_v5_cel <- file.path(ROOT, "Analysis/Hypotheses/V5-20261021T122723Z/outputs/job3/cells_per_donor_set.tsv")
p_v5_feat <- file.path(ROOT, "Analysis/Hypotheses/V5-20261021T122723Z/downloads_untrusted/matrices/GSM9859944_FRPBiopsy_01_Male_Healthy_features.tsv.gz")
p_v3_cnt <- file.path(ROOT, "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk/pseudobulk_hepatocyte.tsv.gz")
p_v3_don <- file.path(ROOT, "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk/donor_metadata.tsv")
p_v3_cel <- file.path(ROOT, "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk/nuclei_per_donor_lineage.tsv")
p_uc_cnt <- file.path(ROOT, "Analysis/Hypotheses/NUCLEUS-FIX-20261107/unique_counts/pseudobulk/Hepatocytes_pseudobulk.csv")
p_uc_don <- file.path(ROOT, "Analysis/Hypotheses/NUCLEUS-FIX-20261107/unique_counts/metadata/donor_metadata.tsv")
p_uc_lib <- file.path(ROOT, "Analysis/Hypotheses/NUCLEUS-FIX-20261107/unique_counts/metadata/library_metadata.tsv")
p_rv1_sex <- file.path(ROOT, "Analysis/Hypotheses/RV1-20261031T155620Z/out-20261102T130905Z/donors.tsv")
p_docF <- file.path(ROOT, "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/metadata/documented_F_source.tsv")
p_ann <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")

# V5 GSE298719
v5_dq <- fread(p_v5_don)[to_bool(included)]
v5_cl <- fread(p_v5_cel)[set == "hep" & n_cells >= 30L]
v5_m <- merge(v5_dq[, .(donor = gsm, female = as.integer(sex_f), disease = as.integer(disease), fibrosis = as.numeric(fibrosis))],
              v5_cl[, .(donor = gsm, n_hep = as.integer(n_cells))], by = "donor")
v5_m[, `:=`(cohort = "GSE298719", stage_ord = ifelse(disease == 0L, 0L, ifelse(fibrosis >= 4, 3L, 2L)), eligible = TRUE)]
v5_raw <- fread(cmd = sprintf("zcat %s", shQuote(p_v5_cnt)), data.table = FALSE)
v5_feat <- fread(cmd = sprintf("zcat %s", shQuote(p_v5_feat)), header = FALSE, col.names = c("ens", "sym", "tp"))
v5_syms <- v5_feat$sym[match(as.character(v5_raw[[1]]), v5_feat$ens)]
keep_v5 <- !is.na(v5_syms) & nzchar(v5_syms)
v5_mat0 <- as.matrix(v5_raw[keep_v5, -1, drop = FALSE]); storage.mode(v5_mat0) <- "numeric"
v5_mat <- rowsum(v5_mat0, group = v5_syms[keep_v5], reorder = FALSE)[, v5_m$donor, drop = FALSE]

# V3 GSE256398
v3_don <- fread(p_v3_don)
v3_cel <- fread(p_v3_cel)
v3_m <- merge(v3_don[, .(donor, stage_ord = as.integer(stage_ordinal), female = as.integer(tolower(sex) == "female"))],
              v3_cel[, .(donor, n_hep = as.integer(hepatocyte))], by = "donor")
v3_m <- v3_m[n_hep >= 30L]
v3_m[, `:=`(cohort = "GSE256398", disease = ifelse(stage_ord >= 2L, 1L, ifelse(stage_ord == 0L, 0L, NA_integer_)), eligible = TRUE)]
v3_raw <- fread(cmd = sprintf("zcat %s", shQuote(p_v3_cnt)), data.table = FALSE)
v3_syms <- as.character(v3_raw[[2]])
keep_v3 <- !is.na(v3_syms) & nzchar(v3_syms)
v3_mat0 <- as.matrix(v3_raw[keep_v3, -(1:2), drop = FALSE]); storage.mode(v3_mat0) <- "numeric"
v3_mat <- rowsum(v3_mat0, group = v3_syms[keep_v3], reorder = FALSE)[, v3_m$donor, drop = FALSE]

# Atlas GSE202379 & GSE244832 (unique_counts collapsed cleanly per donor)
uc_don <- fread(p_uc_don)
uc_lib <- fread(p_uc_lib)
uc_don[, exclude_flag := to_bool(exclude)]
uc_don[, binary_flag := to_bool(binary_eligible)]
uc_lib[, prep_flag := to_bool(prep_eligible)]
rv1_sex <- fread(p_rv1_sex)[, .(donor, female = as.integer(tolower(sex) == "female"))]

uc_raw <- fread(p_uc_cnt, data.table = FALSE)
uc_genes <- as.character(uc_raw[[1]])
uc_mat0 <- as.matrix(uc_raw[, -1, drop = FALSE]); storage.mode(uc_mat0) <- "numeric"
rownames(uc_mat0) <- uc_genes

load_atlas_cohort <- function(ds_name) {
  d_sub <- uc_don[dataset == ds_name & !exclude_flag & binary_flag & condition_harmonized %in% c("Healthy", "MASLD")]
  l_sub <- uc_lib[dataset == ds_name & prep_flag & sample %in% colnames(uc_mat0) & donor %in% d_sub$donor]
  dons <- unique(l_sub$donor)
  cmat <- sapply(dons, function(d) {
    cs <- l_sub[donor == d, sample]
    if (length(cs) == 1L) uc_mat0[, cs] else rowSums(uc_mat0[, cs, drop = FALSE])
  })
  if (is.null(dim(cmat))) cmat <- matrix(cmat, nrow = nrow(uc_mat0), ncol = length(dons))
  rownames(cmat) <- uc_genes; colnames(cmat) <- dons
  nhep_dt <- l_sub[, .(n_hep = as.integer(sum(n_Hepatocytes))), by = donor]
  m <- merge(nhep_dt, d_sub[, .(donor, stage_ord = as.integer(round(as.numeric(stage_ordinal))))], by = "donor")
  m <- merge(m, rv1_sex, by = "donor")
  m <- m[n_hep >= 30L]
  m[, `:=`(cohort = ds_name, disease = ifelse(stage_ord >= 2L, 1L, ifelse(stage_ord == 0L, 0L, NA_integer_)), eligible = TRUE)]
  list(meta = m, counts = cmat[, m$donor, drop = FALSE])
}
a202 <- load_atlas_cohort("GSE202379")
a244 <- load_atlas_cohort("GSE244832")

COHORT_ORDER <- c("GSE298719", "GSE256398", "GSE202379", "GSE244832")
cohorts <- list(
  GSE298719 = list(meta = v5_m, counts = v5_mat),
  GSE256398 = list(meta = v3_m, counts = v3_mat),
  GSE202379 = a202,
  GSE244832 = a244
)

# Normalize within each cohort across all >=30-hep eligible donors
for (cname in COHORT_ORDER) {
  cm <- cohorts[[cname]]$counts
  dge <- calcNormFactors(DGEList(counts = cm), method = "TMM")
  cohorts[[cname]]$logcpm <- cpm(dge, log = TRUE, prior.count = 1)
  cohorts[[cname]]$meta[, log10_n_hep := log10(as.numeric(n_hep))]
}

donors_out <- rbindlist(lapply(COHORT_ORDER, function(c) cohorts[[c]]$meta[, .(donor, female, stage_ord, n_hep, cohort, disease, eligible)]))
fwrite(donors_out, file.path(OUT_DIR, "donors.tsv"), sep = "\t")

# Gate check per cohort
gate_dt <- rbindlist(lapply(COHORT_ORDER, function(c) {
  cm <- cohorts[[c]]$counts
  rbindlist(lapply(GENES, function(g) {
    pres <- g %in% rownames(cm)
    df <- if (pres) mean(cm[g, ] > 0) else NA_real_
    data.table(cohort = c, gene = g, present = pres, detect_frac = df, testable = pres && df >= 0.80)
  }))
}))
fwrite(gate_dt, file.path(OUT_DIR, "gene_gate.tsv"), sep = "\t")

# Long donor_logcpm
d_lcpm <- rbindlist(lapply(COHORT_ORDER, function(c) {
  lc <- cohorts[[c]]$logcpm
  tg <- gate_dt[cohort == c & testable == TRUE, gene]
  rbindlist(lapply(tg, function(g) {
    data.table(cohort = c, gene = g, donor = cohorts[[c]]$meta$donor, logcpm = as.numeric(lc[g, cohorts[[c]]$meta$donor]))
  }))
}))
fwrite(d_lcpm, file.path(OUT_DIR, "donor_logcpm.tsv"), sep = "\t")

# ---- HC3 + Freedman-Lane studentized permutation helpers ----
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

fl_null_draws <- function(y, X, obs, P_idx, j = 2L, chunk = 5000L) {
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

run_meta_suite <- function(cohort_specs, gene_vec, seed = SEED, B_draws = B) {
  # cohort_specs: named list of data.table(donor, lcpm, X, n0, n1)
  set.seed(seed)
  cnames <- names(cohort_specs)
  P_list <- lapply(cnames, function(c) replicate(B_draws, sample.int(nrow(cohort_specs[[c]]$X))))
  names(P_list) <- cnames

  per_rows <- list()
  null_store <- list()

  for (c in cnames) {
    sp <- cohort_specs[[c]]
    n <- nrow(sp$X)
    tg <- gate_dt[cohort == c & testable == TRUE & gene %in% gene_vec, gene]
    null_store[[c]] <- list()
    for (g in tg) {
      y <- as.numeric(sp$lcpm[g, sp$donors])
      obs <- fit_hc3_obs(y, sp$X, j = 2L)
      if (is.null(obs)) next
      nl <- fl_null_draws(y, sp$X, obs, P_list[[c]], j = 2L)
      null_store[[c]][[g]] <- list(obs = obs, nl = nl, n0 = sp$n0, n1 = sp$n1, n = n)
      ok <- is.finite(nl$t)
      p_p <- (1.0 + sum(abs(nl$t[ok]) >= abs(obs$t) * (1.0 - 1e-12))) / (sum(ok) + 1.0)
      per_rows[[paste0(c, "_", g)]] <- data.table(
        cohort = c, gene = g, n = n, n0 = sp$n0, n1 = sp$n1,
        beta = obs$beta, se = obs$se, t_obs = obs$t, df = obs$df, max_leverage = obs$max_h,
        ci_low = obs$ci_low, ci_high = obs$ci_high,
        null_t_sd = sd(nl$t[ok]), null_beta_sd = sd(nl$beta[ok]),
      n_draws = sum(ok), n_draws_dropped = sum(!ok), p_perm = p_p,
      perm_q95_abs_t = as.numeric(quantile(abs(nl$t[ok]), 0.95)),
      perm_ci_low = obs$beta - as.numeric(quantile(abs(nl$t[ok]), 0.95)) * obs$se,
      perm_ci_high = obs$beta + as.numeric(quantile(abs(nl$t[ok]), 0.95)) * obs$se
    )
  }
  }
  per_dt <- rbindlist(per_rows)
  per_dt[, p_perm_holm_within_cohort := p.adjust(p_perm, method = "holm"), by = cohort]

  combine_set <- function(cset, gene, label) {
    avail <- cset[vapply(cset, function(c) !is.null(null_store[[c]][[gene]]), TRUE)]
    if (length(avail) == 0L) return(NULL)
    bs <- vapply(avail, function(c) null_store[[c]][[gene]]$obs$beta, 0.0)
    ses <- vapply(avail, function(c) null_store[[c]][[gene]]$obs$se, 0.0)
    ns <- vapply(avail, function(c) null_store[[c]][[gene]]$n, 0L)
    dfs <- vapply(avail, function(c) null_store[[c]][[gene]]$obs$df, 0L)
    ws <- 1.0 / (ses^2)
    b_fe <- sum(ws * bs) / sum(ws)
    se_fe <- 1.0 / sqrt(sum(ws))
    T_obs <- b_fe / se_fe

    # Freedman-Lane IVW studentized permutation null across cohorts in avail
    num_b <- numeric(B_draws)
    den_b <- numeric(B_draws)
    for (c in avail) {
      nl <- null_store[[c]][[gene]]$nl
      w_b <- 1.0 / (nl$se^2)
      num_b <- num_b + w_b * nl$beta
      den_b <- den_b + w_b
    }
    b_fe_null <- num_b / den_b
    se_fe_null <- 1.0 / sqrt(den_b)
    T_null <- b_fe_null / se_fe_null
    ok <- is.finite(T_null)
    p_p <- (1.0 + sum(abs(T_null[ok]) >= abs(T_obs) * (1.0 - 1e-12))) / (sum(ok) + 1.0)
    nt_sd <- sd(T_null[ok])
    zc <- qnorm(0.975)
    tc_min <- qt(0.975, df = min(dfs))
    data.table(
      gene = gene,
      cohorts = paste(avail, collapse = ";"),
      n_cohorts = length(avail),
      n_donors = sum(ns),
      df_min = min(dfs),
      beta_fe = b_fe,
      se_fe = se_fe,
      T = T_obs,
      ci_low = b_fe - zc * se_fe,
      ci_high = b_fe + zc * se_fe,
      ci_low_tdf = b_fe - tc_min * se_fe,
      ci_high_tdf = b_fe + tc_min * se_fe,
      n_negative = sum(bs < 0),
      null_sd = nt_sd,
      null_beta_sd = sd(b_fe_null[ok]),
      mde = 2.8 * se_fe * nt_sd,
      mde_t = 2.8 * nt_sd,
      p_perm = p_p,
      n_draws = sum(ok),
      analysis = label
    )
  }

  list(per_cohort = per_dt, null_store = null_store, combine_set = combine_set)
}

build_spec <- function(cnames, x_var = "disease", covs = c("female", "log10_n_hep"), filter_expr = NULL) {
  res <- list()
  for (c in cnames) {
    m <- copy(cohorts[[c]]$meta)
    if (!is.null(filter_expr)) m <- m[eval(filter_expr)]
    if (x_var == "disease") m <- m[!is.na(disease)]
    rhs <- c(x_var, covs)
    X <- cbind(Intercept = 1.0, as.matrix(m[, rhs, with = FALSE]))
    n0 <- sum(m[[x_var]] == 0L)
    n1 <- sum(m[[x_var]] > 0L)
    res[[c]] <- list(donors = m$donor, lcpm = cohorts[[c]]$logcpm, X = X, n0 = n0, n1 = n1)
  }
  res
}

# ---- Run primary 4-cohort HC3 Freedman-Lane suite ----
specs_prim <- build_spec(COHORT_ORDER, "disease", c("female", "log10_n_hep"))
suite_prim <- run_meta_suite(specs_prim, GENES)
fwrite(suite_prim$per_cohort, file.path(OUT_DIR, "per_cohort.tsv"), sep = "\t")

# Combined arms
comb_rep <- rbindlist(lapply(GENES, function(g) suite_prim$combine_set(c("GSE202379", "GSE244832"), g, "primary_replication")))
comb_4 <- rbindlist(lapply(GENES, function(g) suite_prim$combine_set(COHORT_ORDER, g, "four_cohort_descriptive")))
comb_dig1 <- rbindlist(lapply(GENES, function(g) suite_prim$combine_set(c("GSE298719", "GSE256398"), g, "dig1_cohorts_descriptive")))

if (IS_TRIO) {
  comb_rep[, holm := p.adjust(p_perm, method = "holm", n = length(GENES))]
  comb_rep[, dig1_effect := DIG1_EFFECT[gene]]
  comb_rep[, decision := ifelse(beta_fe < 0 & holm < 0.05 & n_negative == n_cohorts, "replicated",
                                ifelse(mde <= dig1_effect, "not_replicated", "uninformative"))]
  comb_4[, holm := NA_real_]; comb_4[, dig1_effect := DIG1_EFFECT[gene]]; comb_4[, decision := NA_character_]
  comb_dig1[, holm := NA_real_]; comb_dig1[, dig1_effect := DIG1_EFFECT[gene]]; comb_dig1[, decision := NA_character_]
} else {
  # XMOD2: four-cohort Holm over family of 8, plus two-cohort replication Holm over 8
  comb_4[, holm := p.adjust(p_perm, method = "holm", n = length(GENES))]
  comb_rep[, holm := p.adjust(p_perm, method = "holm", n = length(GENES))]
  comb_dig1[, holm := NA_real_]
  comb_4[, `:=`(dig1_effect = NA_real_, decision = "not_applicable")]
  comb_rep[, `:=`(dig1_effect = NA_real_, decision = "not_applicable")]
  comb_dig1[, `:=`(dig1_effect = NA_real_, decision = NA_character_)]
}

comb_all <- rbindlist(list(comb_rep, comb_4, comb_dig1), fill = TRUE)
pc_disc <- suite_prim$per_cohort[p_perm < 0.05 & ci_low <= 0 & ci_high >= 0]
if (nrow(pc_disc) > 0L) {
  cat("NOTE per_cohort p_perm < 0.05 with t_df CI touching 0 (df small-sample residual perm):",
      paste(sprintf("%s:%s(p=%.5f,ci=[%.4f,%.4f])", pc_disc$cohort, pc_disc$gene, pc_disc$p_perm, pc_disc$ci_low, pc_disc$ci_high), collapse = "; "), "\n")
}
stopifnot(sum(suite_prim$per_cohort$p_perm_holm_within_cohort < 0.05 &
              suite_prim$per_cohort$perm_ci_low <= 0 & suite_prim$per_cohort$perm_ci_high >= 0, na.rm = TRUE) == 0L)
stopifnot(sum(comb_all$holm < 0.05 & comb_all$ci_low <= 0 & comb_all$ci_high >= 0, na.rm = TRUE) == 0L)
fwrite(comb_all, file.path(OUT_DIR, "combined.tsv"), sep = "\t")

# Heterogeneity (fixed + DerSimonian-Laird random effects on observed HC3 SE)
het_rows <- lapply(GENES, function(g) {
  sub <- suite_prim$per_cohort[gene == g]
  k <- nrow(sub)
  if (k < 2L) return(NULL)
  w <- 1.0 / (sub$se^2)
  fe <- sum(w * sub$beta) / sum(w)
  se_fe <- sqrt(1.0 / sum(w))
  Q <- sum(w * (sub$beta - fe)^2)
  Qp <- pchisq(Q, df = k - 1L, lower.tail = FALSE)
  I2 <- if (Q > k - 1L) (Q - (k - 1L)) / Q else 0.0
  C_dl <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0.0, (Q - (k - 1L)) / C_dl)
  wr <- 1.0 / (sub$se^2 + tau2)
  re <- sum(wr * sub$beta) / sum(wr)
  se_re <- sqrt(1.0 / sum(wr))
  data.table(gene = g, k = k, fixed = fe, fixed_se = se_fe, Q = Q, Q_p = Qp, I2 = I2, tau2 = tau2, random = re, random_se = se_re)
})
fwrite(rbindlist(Filter(Negate(is.null), het_rows)), file.path(OUT_DIR, "heterogeneity.tsv"), sep = "\t")

# Leave-one-out
loo_rows <- rbindlist(lapply(COHORT_ORDER, function(drop_c) {
  keep_c <- setdiff(COHORT_ORDER, drop_c)
  dt <- rbindlist(lapply(GENES, function(g) suite_prim$combine_set(keep_c, g, "leave_one_out")))
  dt[, left_out := drop_c]
  dt
}))
fwrite(loo_rows, file.path(OUT_DIR, "leave_one_out.tsv"), sep = "\t")

# ---- Sensitivity suite (all HC3 Freedman-Lane) ----
sens_list <- list()
# 1) Ordinal stage (all eligible donors)
sp_ord <- build_spec(COHORT_ORDER, "stage_ord", c("female", "log10_n_hep"))
st_ord <- run_meta_suite(sp_ord, GENES)
s_ord_rep <- rbindlist(lapply(GENES, function(g) st_ord$combine_set(c("GSE202379", "GSE244832"), g, "ordinal_stage_rep")))
s_ord_4 <- rbindlist(lapply(GENES, function(g) st_ord$combine_set(COHORT_ORDER, g, "ordinal_stage_all")))
s_ord_rep[, set := "replication"]; s_ord_4[, set := "four_cohort"]
sens_list[["ord_rep"]] <- s_ord_rep; sens_list[["ord_4"]] <- s_ord_4

# 2) No nuclei covariate on replication cohorts
sp_nonuc <- build_spec(c("GSE202379", "GSE244832"), "disease", c("female"))
st_nonuc <- run_meta_suite(sp_nonuc, GENES)
s_nn <- rbindlist(lapply(GENES, function(g) st_nonuc$combine_set(c("GSE202379", "GSE244832"), g, "no_nuclei_cov")))
s_nn[, set := "replication"]
sens_list[["no_nuc"]] <- s_nn

# 3) min_50_hep on replication cohorts
sp_50 <- build_spec(c("GSE202379", "GSE244832"), "disease", c("female", "log10_n_hep"), quote(n_hep >= 50L))
st_50 <- run_meta_suite(sp_50, GENES)
s_50 <- rbindlist(lapply(GENES, function(g) st_50$combine_set(c("GSE202379", "GSE244832"), g, "min_50_hep")))
s_50[, set := "replication"]
sens_list[["min50"]] <- s_50

# 4) Steatosis counted as disease
for (c in c("GSE202379", "GSE244832")) cohorts[[c]]$meta[, disease_steat := as.integer(stage_ord >= 1L)]
sp_steat <- lapply(c("GSE202379", "GSE244832"), function(c) {
  m <- cohorts[[c]]$meta
  X <- cbind(Intercept = 1.0, disease_steat = as.numeric(m$disease_steat),
             female = as.numeric(m$female), log10_n_hep = m$log10_n_hep)
  list(donors = m$donor, lcpm = cohorts[[c]]$logcpm, X = X, n0 = sum(m$disease_steat == 0L), n1 = sum(m$disease_steat == 1L))
})
names(sp_steat) <- c("GSE202379", "GSE244832")
st_steat <- run_meta_suite(sp_steat, GENES)
s_st <- rbindlist(lapply(GENES, function(g) st_steat$combine_set(c("GSE202379", "GSE244832"), g, "steatosis_as_disease")))
s_st[, set := "replication"]
sens_list[["steat"]] <- s_st

sens_dt <- rbindlist(sens_list, fill = TRUE)
sens_dt[, holm := p.adjust(p_perm, method = "holm"), by = .(analysis, set)]
stopifnot(sum(sens_dt$holm < 0.05 & sens_dt$ci_low <= 0 & sens_dt$ci_high >= 0, na.rm = TRUE) == 0L)
fwrite(sens_dt, file.path(OUT_DIR, "sensitivity.tsv"), sep = "\t")

# ---- Documented F stage in GSE202379 ----
docF <- fread(p_docF)[, .(donor, F_stage = as.numeric(F_stage_documented))]
m202 <- merge(cohorts$GSE202379$meta, docF, by = "donor")
m202 <- m202[!is.na(F_stage)]
X_F <- cbind(Intercept = 1.0, F_stage = m202$F_stage, female = as.numeric(m202$female), log10_n_hep = m202$log10_n_hep)
set.seed(SEED)
P_F <- replicate(B, sample.int(nrow(m202)))
f_rows <- lapply(GENES, function(g) {
  y <- as.numeric(cohorts$GSE202379$logcpm[g, m202$donor])
  obs <- fit_hc3_obs(y, X_F, j = 2L)
  nl <- fl_null_draws(y, X_F, obs, P_F, j = 2L)
  ok <- is.finite(nl$t)
  p_p <- (1.0 + sum(abs(nl$t[ok]) >= abs(obs$t) * (1.0 - 1e-12))) / (sum(ok) + 1.0)
  data.table(cohort = "GSE202379", gene = g, n = nrow(m202), beta = obs$beta, se = obs$se, t_obs = obs$t,
             df = obs$df, ci_low = obs$ci_low, ci_high = obs$ci_high,
             n0 = sum(m202$F_stage == 0), n1 = sum(m202$F_stage > 0), p_perm = p_p)
})
f_dt <- rbindlist(f_rows)
f_dt[, holm := p.adjust(p_perm, method = "holm")]
stopifnot(sum(f_dt$holm < 0.05 & f_dt$ci_low <= 0 & f_dt$ci_high >= 0, na.rm = TRUE) == 0L)
fwrite(f_dt, file.path(OUT_DIR, "f_stage_GSE202379.tsv"), sep = "\t")

# ---- Scale check across common genes passing detection gate in GSE202379 & GSE244832 ----
g202 <- rownames(cohorts$GSE202379$counts)[rowMeans(cohorts$GSE202379$counts > 0) >= 0.80]
g244 <- rownames(cohorts$GSE244832$counts)[rowMeans(cohorts$GSE244832$counts > 0) >= 0.80]
common_g <- intersect(g202, g244)
vec_hc3 <- function(cname) {
  sp <- specs_prim[[cname]]
  X <- sp$X
  XtXi <- solve(crossprod(X))
  Cm <- XtXi %*% t(X)
  h <- rowSums((X %*% XtXi) * X)
  Y <- t(sp$lcpm[common_g, sp$donors, drop = FALSE])
  Bm <- Cm %*% Y
  Em <- Y - X %*% Bm
  Wm <- (Em^2) / ((1.0 - h)^2)
  se <- sqrt(as.vector(crossprod(Cm[2L, ]^2, Wm)))
  list(beta = as.numeric(Bm[2L, ]), se = se)
}
r202 <- vec_hc3("GSE202379"); r244 <- vec_hc3("GSE244832")
w1 <- 1.0 / (r202$se^2); w2 <- 1.0 / (r244$se^2)
b_fe_all <- (w1 * r202$beta + w2 * r244$beta) / (w1 + w2)
se_fe_all <- 1.0 / sqrt(w1 + w2)
T_fe_all <- b_fe_all / se_fe_all
h_mean <- 0.5 * (rowMeans(cohorts$GSE202379$logcpm[common_g, cohorts$GSE202379$meta[disease == 0L, donor], drop = FALSE]) +
                 rowMeans(cohorts$GSE244832$logcpm[common_g, cohorts$GSE244832$meta[disease == 0L, donor], drop = FALSE]))
sc_dt <- data.table(gene = common_g, beta_fe = b_fe_all, se_fe = se_fe_all, T = T_fe_all, healthy_logcpm = h_mean)
sc_dt[, decile := as.integer(cut(healthy_logcpm, breaks = quantile(healthy_logcpm, probs = seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE))]
med_T <- median(sc_dt$T, na.rm = TRUE)
sc_genes <- rbindlist(lapply(GENES, function(g) {
  row <- sc_dt[gene == g]
  if (nrow(row) == 0L) return(NULL)
  d_sub <- sc_dt[decile == row$decile[1]]
  data.table(gene = g, set = "GSE202379;GSE244832", n_genes_common = length(common_g),
             beta_fe = row$beta_fe[1], T = row$T[1], healthy_logcpm = row$healthy_logcpm[1],
             decile = row$decile[1], n_in_decile = nrow(d_sub),
             frac_decile_below_T = mean(d_sub$T <= row$T[1]),
             frac_all_below_T = mean(sc_dt$T <= row$T[1]),
             median_T_all = med_T)
}))
fwrite(sc_dt, file.path(OUT_DIR, "scale_check_all_genes.tsv.gz"), sep = "\t")
fwrite(sc_genes, file.path(OUT_DIR, "scale_check.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("DONE ", MODE_LABEL, " -> ", OUT_DIR, "\n")
