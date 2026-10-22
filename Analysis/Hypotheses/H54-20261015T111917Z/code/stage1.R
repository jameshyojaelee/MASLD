#!/usr/bin/env Rscript
# H54 stage 1 (job 2). Controls, prespecified diagnostics, power gate G1, stored
# permutations and bootstrap indices, stop_rule.tsv.
#
# THIS STAGE COMPUTES NO NEW-COHORT S1 OR E ON REAL DATA. On GSE281797 and GSE225740
# it builds design blocks, testability, leave-shared-out eligibility (no correlation),
# stage leakage, residual-component diagnostics, and statistics on GENE-PERMUTED data
# only (NC1, NC3). t4_one_pass() is called on real data only for the GSE193066 holdout.
#
# Usage: stage1.R [out_dir] [composition_dir]
#   env H54_DEBUG_DIVISOR=k divides the draw counts by k (debug runs only; never the
#   real run, which must write <run>/stage1 with full counts).

source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z/code/h54_lib.R")
args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) >= 1L) args[[1]] else file.path(H54_RUN, "stage1")
COMP <- if (length(args) >= 2L) args[[2]] else file.path(H54_RUN, "composition")
DIV <- as.integer(Sys.getenv("H54_DEBUG_DIVISOR", "1"))
if (DIV != 1L && normalizePath(dirname(OUT)) == normalizePath(H54_RUN) && basename(OUT) == "stage1")
  fail("debug divisor not allowed for the real stage1 directory")
n_perm <- N_PERM %/% DIV; n_boot <- N_BOOT %/% DIV; n_nc1 <- N_NC1 %/% DIV; n_g1 <- N_G1 %/% DIV
n_perm_e <- min(N_PERM_E, n_perm)
t0 <- Sys.time()

verify_prespec_and_inputs()
if (!file.exists(file.path(COMP, "code_sha256.txt"))) fail("composition output incomplete: ", COMP)
refuse_if_exists(OUT)
w <- function(x, f) fwrite(x, file.path(OUT, f), sep = "\t", na = "NA")
msg("stage1 -> %s; draws: perm %d, boot %d, NC1 %d, G1 %d (divisor %d)", OUT, n_perm, n_boot, n_nc1, n_g1, DIV)

frozen <- load_frozen()
pr <- frozen$pairs
pidx <- pair_index_matrix(frozen)
features <- frozen$features
msg("frozen: %d programs, %d pairs, %d coordinated, %d C3, %d informative_null",
    length(features), nrow(pr), sum(pr$is_coord), sum(pr$is_c3), sum(pr$is_infnull))

# ------------------------------------------------------------------ sets
sm <- build_set_meta()
countsA <- read_counts(file.path(H54_ROOT, "data/GSE281797/GSE281797_counts_ensembl_stable.tsv.gz"))
countsB <- read_counts(file.path(H54_ROOT, "data/GSE225740/GSE225740_counts_ensembl_stable.tsv.gz"))
new_sets <- c("A91", "B87", "A66", "B53", "A94", "A37")
S <- list()
for (s in new_sets) {
  meta <- sm$sets[[s]]
  ex <- build_set_expression(if (SET_COHORT[[s]] == "A") countsA else countsB, meta, frozen)
  S[[s]] <- c(list(meta = meta), ex)
}
rm(countsA, countsB); invisible(gc())
H <- build_holdout(frozen, pooled = TRUE)
S[["H105"]] <- H
all_sets <- c(new_sets, "H105")

# ------------------------------------------------------------ composition
read_props <- function(ds) {
  p <- file.path(COMP, "results", ds, paste0(ds, "_bayesprism_proportions.tsv"))
  x <- read.delim(p, row.names = 1, check.names = FALSE)
  dt <- data.table(sample_id = rownames(x), as.data.table(x))
  dt
}
props <- list(A = read_props("GSE281797"), B = read_props("GSE225740"))
if (!setequal(props$A$sample_id, sm$A_all$sample_id)) fail("GSE281797 composition sample ids differ")
if (!setequal(props$B$sample_id, sm$B_all$sample_id)) fail("GSE225740 composition sample ids differ")
lineages_comp <- setdiff(names(props$A), "sample_id")
if (!setequal(lineages_comp, setdiff(names(props$B), "sample_id"))) fail("lineage sets differ between cohorts")
comp_qc <- rbindlist(lapply(names(props), function(cn) {
  M <- as.matrix(props[[cn]][, ..lineages_comp])
  data.table(cohort = cn, n = nrow(M), n_lineages = ncol(M),
             hepatocyte_median = median(M[, "Hepatocytes"]),
             max_lineage_zero_fraction = max(colMeans(M == 0)),
             lineage_max_zero = colnames(M)[which.max(colMeans(M == 0))],
             row_sum_range = paste(signif(range(rowSums(M)), 6), collapse = "-"))
}))
comp_qc[, passes_hepatocyte_range := hepatocyte_median >= composition_hepatocyte_median_range[[1]] &
          hepatocyte_median <= composition_hepatocyte_median_range[[2]]]
comp_qc[, passes_zero_fraction := max_lineage_zero_fraction <= composition_max_zero_fraction]
comp_qc[, composition_arm_readable := passes_hepatocyte_range & passes_zero_fraction]
comp_lineage <- rbindlist(lapply(names(props), function(cn) {
  M <- as.matrix(props[[cn]][, ..lineages_comp])
  data.table(cohort = cn, lineage = colnames(M), median = apply(M, 2, median),
             zero_fraction = colMeans(M == 0), p05 = apply(M, 2, quantile, 0.05))
}))

# CLR principal components within each analysis set (build_composition_scores, frozen v9 code).
comp_cov <- list(); comp_var <- list()
for (s in new_sets) {
  meta <- S[[s]]$meta
  pc <- props[[SET_COHORT[[s]]]][match(meta$sample_id, sample_id)]
  pc[, dataset := s]
  cs <- build_composition_scores(pc, data.table(sample_id = meta$sample_id, dataset = s),
                                 lineages_comp, composition_pseudocount,
                                 n_components = composition_n_clr_components)
  qc <- S[[s]]$qc[match(meta$sample_id, sample_id)]
  comp_cov[[s]] <- data.table(set = s, sample_id = meta$sample_id,
                              composition_pc1 = cs$pcs$composition_pc1,
                              composition_pc2 = cs$pcs$composition_pc2,
                              composition_pc3 = cs$pcs$composition_pc3,
                              log_library_size = log(qc$library_size_universe),
                              top100_share = qc$top100_share)
  comp_var[[s]] <- cbind(set = s, cs$pc_variance)
}
comp_cov_dt <- rbindlist(comp_cov)
w(comp_qc, "diag11_composition_qc.tsv"); w(comp_lineage, "diag11_composition_by_lineage.tsv")
w(rbindlist(comp_var), "diag11_composition_clr_variance.tsv"); w(comp_cov_dt, "composition_covariates_by_set.tsv")

COMP_FORMULA <- ~ sex_final + fibrosis_factor + nas_factor + composition_pc1 + composition_pc2 +
  composition_pc3 + log_library_size + top100_share

# ------------------------------------------- diagnostics 1-6, G2, G3 per set
blocks <- list(); opsl <- list(); okpairs <- list()
d1 <- list(); d1tab <- list(); d2 <- list(); d2drop <- list(); d3 <- list(); d4 <- list()
d4full <- list(); d5 <- list(); d6 <- list(); d6cor <- list()
disc_cov <- fread(file.path(T4_DIR, "global_component_covariates.tsv"))
for (s in all_sets) {
  meta <- S[[s]]$meta
  dz <- design_ok(meta)
  dzc <- if (s != "H105") {
    mm <- merge(meta, comp_cov[[s]][, !"set"], by = "sample_id", sort = FALSE)
    design_ok(mm, COMP_FORMULA)
  } else c(ok = NA, cols = NA, rank = NA, df = NA)
  d1[[s]] <- data.table(set = s, role = SET_ROLE[[s]], n = nrow(meta),
                        design_columns = dz[["cols"]], design_rank = dz[["rank"]], residual_df = dz[["df"]],
                        composition_arm_columns = dzc[["cols"]], composition_arm_rank = dzc[["rank"]],
                        composition_arm_df = dzc[["df"]])
  tb <- meta[, .N, by = .(fibrosis_stage, nas_bin = as.character(t4_nas_bin(nas_score)), sex_final)]
  d1tab[[s]] <- cbind(set = s, tb[order(fibrosis_stage, nas_bin, sex_final)])

  block <- make_block(S[[s]]$PM, meta, seq_len(nrow(meta)), frozen)
  if (is.null(block)) fail("design not estimable for full set ", s)
  blocks[[s]] <- block
  ops <- ops_for(block, frozen)
  opsl[[s]] <- ops
  okp <- eligible_pairs(block, ops, frozen)
  okpairs[[s]] <- okp
  cv <- block$coverage
  tst <- cv[cohort_testable == TRUE, feature_id]
  d2[[s]] <- data.table(set = s, universe_rows_present = S[[s]]$n_universe_present,
                        panel_symbols = length(S[[s]]$panel), panel_symbols_frozen = uniqueN(frozen$weights$gene_symbol),
                        variable_genes = length(block$genes), programs_testable = length(tst),
                        min_retained_weight = min(cv$cohort_weight_coverage),
                        min_retained_weight_testable = min(cv[cohort_testable == TRUE, cohort_weight_coverage]))
  if (any(!cv$cohort_testable)) {
    dd <- cv[cohort_testable == FALSE, .(set = s, feature_id, cohort_weight_coverage)]
    dd[, lineage := frozen$lineage[feature_id]]
    dd[, community := frozen$communities$community[match(feature_id, frozen$communities$feature_id)]]
    d2drop[[s]] <- dd
  }
  tst_f <- features %in% tst
  both <- tst_f[pr$i] & tst_f[pr$j]
  d3[[s]] <- data.table(set = s,
                        coord_both_testable = sum(pr$is_coord & both),
                        coord_eligible_k1 = sum(pr$is_coord & okp),
                        c3_both_testable = sum(pr$is_c3 & both),
                        c3_eligible_k3 = sum(pr$is_c3 & okp),
                        noncoord_eligible = sum(pr$is_noncoord & okp),
                        infnull_eligible = sum(pr$is_infnull & okp),
                        same_lineage_coord_eligible = sum(pr$is_coord & pr$same_lineage & okp),
                        cross_lineage_coord_eligible = sum(pr$is_coord & !pr$same_lineage & okp))
  lk <- t4_stage_leakage(list(block))
  d4full[[s]] <- cbind(set = s, lk)
  d4[[s]] <- data.table(set = s, max_raw_stage_r2 = max(lk$raw_stage_r2),
                        median_stage_r2 = median(lk$stage_r2), max_stage_r2 = max(lk$stage_r2),
                        chance_stage_r2 = lk$chance_stage_r2[[1]],
                        frac_stage_r2_above_chance = mean(lk$stage_r2 > lk$chance_stage_r2))
  qc <- S[[s]]$qc[match(meta$sample_id, sample_id)]
  d5[[s]] <- data.table(set = s, lib_min = min(qc$library_size_universe), lib_median = median(qc$library_size_universe),
                        lib_max = max(qc$library_size_universe), top100_min = min(qc$top100_share),
                        top100_median = median(qc$top100_share), top100_max = max(qc$top100_share))
  rc <- residual_components(block, 3L)
  d6[[s]] <- data.table(set = s, leading_component_share = rc$share, discovery_mean_share = 0.402511574178318)
  covs <- data.table(log_library_size = log(qc$library_size_universe), top100_share = qc$top100_share)
  if (s != "H105") {
    pp <- props[[SET_COHORT[[s]]]][match(meta$sample_id, sample_id), ..lineages_comp]
    covs <- cbind(covs, pp, comp_cov[[s]][, .(composition_pc1, composition_pc2, composition_pc3)])
  }
  d6cor[[s]] <- rbindlist(lapply(names(covs), function(cn) {
    data.table(set = s, covariate = cn, component = 1:3,
               spearman = suppressWarnings(apply(rc$u, 2, function(u) cor(u, covs[[cn]], method = "spearman"))))
  }))
}
D1 <- rbindlist(d1); D2 <- rbindlist(d2); D3 <- rbindlist(d3); D4 <- rbindlist(d4)
w(D1, "diag01_design.tsv"); w(rbindlist(d1tab), "diag01_fibrosis_nas_sex.tsv")
w(D2, "diag02_coverage.tsv"); w(rbindlist(d2drop, fill = TRUE), "diag02_dropped_programs.tsv")
w(rbindlist(lapply(names(blocks), function(s) cbind(set = s, blocks[[s]]$coverage))), "diag02_program_coverage.tsv")
w(D3, "diag03_edge_eligibility.tsv"); w(D4, "diag04_stage_leakage.tsv")
w(rbindlist(d4full), "diag04_stage_leakage_by_program.tsv.gz")
w(rbindlist(d5), "diag05_library.tsv"); w(rbindlist(d6), "diag06_leading_component_share.tsv")
w(rbindlist(d6cor), "diag06_component_covariate_spearman.tsv")
w(disc_cov, "diag06_discovery_reference_global_component_covariates.tsv")
print(D1); print(D2); print(D3); print(D4)

# Prespec counts check (n, design columns, df, testable programs, coordinated edges).
expect <- data.table(set = c("A91", "A66", "A94", "B87", "B53"), cols = c(7L, 5L, 7L, 8L, 5L),
                     df = c(84L, 61L, 87L, 79L, 48L), prog = c(113L, 113L, 113L, 110L, 110L),
                     edges = c(2426L, 2426L, 2426L, 2324L, 2324L))
chk <- merge(expect, merge(D1[, .(set, design_columns, residual_df)],
                           merge(D2[, .(set, programs_testable)], D3[, .(set, coord_both_testable)], by = "set"),
                           by = "set"), by = "set")
chk[, matches := cols == design_columns & df == residual_df & prog == programs_testable & edges == coord_both_testable]
w(chk, "prespec_count_check.tsv"); print(chk)
if (!all(chk$matches)) fail("set counts differ from the prespec table")
if (!all(D1[set %in% c("A91", "B87")]$composition_arm_columns == c(12, 13)) ||
    !all(D1[set %in% c("A91", "B87")]$composition_arm_df == c(79, 74)))
  message("NOTE: composition-arm design columns/df differ from prespec (A91 12/79, B87 13/74)")

# ----------------------------------------------- label permutations (stored)
testable_of <- function(s) blocks[[s]]$coverage[cohort_testable == TRUE, feature_id]
cohort_primary <- c(A = "A91", B = "B87", H = "H105")
for (cn in names(cohort_primary)) {
  for (s in names(SET_COHORT)[SET_COHORT == cn])
    if (!setequal(testable_of(s), testable_of(cohort_primary[[cn]])))
      message("NOTE: testable programs in ", s, " differ from ", cohort_primary[[cn]])
}
perms <- list(); perms_size <- list(); size_strata <- list()
for (cn in names(cohort_primary)) {
  tst <- testable_of(cohort_primary[[cn]])
  if (length(tst) != EXPECTED_TESTABLE[[cn]]) fail("cohort ", cn, " testable ", length(tst))
  set.seed(H54_SEED)
  perms[[cn]] <- draw_lineage_perms(n_perm, tst, frozen$lineage, features)
}
for (cn in names(cohort_primary)) {
  tst <- testable_of(cohort_primary[[cn]])
  # Size-matched strata: lineage x tertile of program gene count (nonzero frozen weight),
  # tertiles taken within lineage over the cohort's testable programs.
  st <- rep(NA_character_, length(features)); names(st) <- features
  for (L in unique(frozen$lineage[tst])) {
    f <- tst[frozen$lineage[tst] == L]
    ng <- frozen$n_genes[f]
    ter <- floor(3 * (rank(ng, ties.method = "first") - 1) / length(f)) + 1
    st[f] <- paste(L, ter, sep = "_t")
  }
  size_strata[[cn]] <- data.table(cohort = cn, feature_id = features, n_genes = frozen$n_genes, stratum = st)
  set.seed(H54_SEED)
  perms_size[[cn]] <- draw_lineage_perms(n_perm, tst, st, features)
}
for (cn in names(perms)) {
  stopifnot(all(apply(perms[[cn]], 1, function(p) all(frozen$lineage[p] == frozen$lineage))),
            all(apply(perms[[cn]], 1, function(p) !anyDuplicated(p))))
}
saveRDS(list(perms = perms, perms_size_matched = perms_size, features = features, seed = H54_SEED,
             n_perm = n_perm, n_perm_e = n_perm_e), file.path(OUT, "label_permutations.rds"))
w(rbindlist(size_strata), "size_matched_strata.tsv")

d10 <- rbindlist(lapply(all_sets, function(s) {
  cn <- SET_COHORT[[s]]
  tst <- testable_of(s)
  tab <- table(frozen$lineage[tst])
  data.table(set = s, cohort = cn, n_lineages = length(tab),
             programs_per_lineage = paste(names(tab), as.integer(tab), sep = ":", collapse = ";"),
             mean_nonfinite_mapped_c1 = perm_nonfinite_fraction(okpairs[[s]], pr$is_coord, perms[[cn]], pidx, frozen),
             mean_nonfinite_mapped_c3 = perm_nonfinite_fraction(okpairs[[s]], pr$is_c3, perms[[cn]], pidx, frozen),
             identity_fraction_mean = mean(apply(perms[[cn]], 1, function(p) mean(p == seq_along(p)))))
}))
w(d10, "diag10_label_permutations.tsv"); print(d10)

# ----------------------------------------------- bootstrap indices (stored)
boot <- list(); d9 <- list()
for (s in new_sets) {
  meta <- S[[s]]$meta
  set.seed(H54_SEED)
  idx <- draw_stratified_boot(n_boot, meta$fibrosis_stage)
  ok <- apply(idx, 1, function(r) design_ok(meta[r])[["ok"]])
  boot[[s]] <- idx
  d9[[s]] <- data.table(set = s, n_draws = nrow(idx), n_design_failed = sum(ok == 0),
                        fraction_failed = mean(ok == 0))
  stopifnot(all(apply(idx, 1, function(r) all(table(meta$fibrosis_stage[r]) == table(meta$fibrosis_stage)))))
}
D9 <- rbindlist(d9)
saveRDS(list(indices = boot, sample_order = lapply(S[new_sets], function(x) x$meta$sample_id), seed = H54_SEED),
        file.path(OUT, "bootstrap_indices.rds"))
w(D9, "diag09_bootstrap_design_failures.tsv"); print(D9)

# ----------------------------------------- PC1, PC1b, PC2, PC3 on H105
hb <- blocks[["H105"]]
hz <- t4_one_pass(list(hb), pr, list(opsl[["H105"]]))
h_eval <- evaluate_z(hz$z_partial_k1, hz$z_partial_k3, perms$H, pidx, frozen, n_perm_e, all_perms = TRUE)
pc1_s1_diff <- abs(h_eval$s1_k1_signmatch - PC1_S1)
pc1_base_diff <- abs(h_eval$base_noncoord - PC1_BASE)
pc1_pass <- pc1_s1_diff < 1e-9 && pc1_base_diff < 1e-9
msg("PC1: S1 %.15f (diff %.2e), comparator %.15f (diff %.2e) -> %s",
    h_eval$s1_k1_signmatch, pc1_s1_diff, h_eval$base_noncoord, pc1_base_diff, pc1_pass)
p2 <- p_perm_upper(h_eval$s1_k1, h_eval$null_k1)
p3 <- p_perm_upper(h_eval$s1_k3, h_eval$null_k3)
Hb <- build_holdout(frozen, pooled = FALSE)
hbb <- make_block(Hb$PM, Hb$meta, seq_len(nrow(Hb$meta)), frozen)
hbz <- t4_one_pass(list(hbb), pr, list(ops_for(hbb, frozen)))
hb_eval <- evaluate_z(hbz$z_partial_k1, hbz$z_partial_k3, perms$H, pidx, frozen, n_perm_e)
f01 <- which(H$meta$fibrosis_stage <= 1L)
stopifnot(length(f01) == 43L)
hf <- make_block(H$PM, H$meta, f01, frozen)
hfz <- t4_one_pass(list(hf), pr, list(ops_for(hf, frozen)))
hf_eval <- evaluate_z(hfz$z_partial_k1, hfz$z_partial_k3, perms$H, pidx, frozen, n_perm_e)
# NC2 on the holdout only (descriptive here; the new-cohort NC2 needs real z and is in stage 2).
nc2_h <- p_perm_upper(s1_frac(hz$z_partial_k1, pr$is_infnull)[["s1"]],
                      perm_s1(hz$z_partial_k1, pr$is_infnull, perms$H, pidx, frozen))
pcs <- data.table(
  control = c("PC1", "PC1", "PC1b", "PC2", "PC2", "PC2", "PC2", "PC3", "PC3", "PC3", "PC3",
              "H105_E", "H105_E", "H105_F0F1", "H105_F0F1", "H105_F0F1", "H105_F0F1", "NC2_H105"),
  quantity = c("S1_k1_signmatch", "noncoord_comparator", "S1_k1_withincohort_TMM",
               "S1_k1", "label_null_mean_1000", "label_null_sd_1000", "p_perm_k1",
               "S1_k3", "label_null_mean_1000_k3", "label_null_sd_1000_k3", "p_perm_k3",
               "E_k1", "E_k3", "S1_k1", "S1_k3", "E_k1", "E_k3", "p_perm_k1_informative_null"),
  value = c(h_eval$s1_k1_signmatch, h_eval$base_noncoord, hb_eval$s1_k1,
            h_eval$s1_k1, mean(h_eval$null_k1), sd(h_eval$null_k1), p2,
            h_eval$s1_k3, mean(h_eval$null_k3), sd(h_eval$null_k3), p3,
            h_eval$e_k1, h_eval$e_k3, hf_eval$s1_k1, hf_eval$s1_k3, hf_eval$e_k1, hf_eval$e_k3, nc2_h),
  reference = c(PC1_S1, PC1_BASE, PC1_S1, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
pcs[, abs_diff := abs(value - reference)]
w(pcs, "controls_PC_H105.tsv"); print(pcs)
pc1b_shift <- hb_eval$s1_k1 - h_eval$s1_k1
pc2_validated <- p2 < 0.05
w(data.table(pair_index = pr$pair_index, z_k1 = hz$z_partial_k1, z_k3 = hz$z_partial_k3),
  "H105_edge_z.tsv.gz")

# ------------------------------------------------- NC1 and NC3
nc1 <- list(); nc3 <- list(); nc1_z <- list()
for (s in c("A91", "B87", "H105")) {
  cn <- SET_COHORT[[s]]
  block <- blocks[[s]]; ops <- opsl[[s]]
  set.seed(H54_SEED)
  Z1 <- matrix(NA_real_, nrow(pr), n_nc1); Z3 <- Z1
  rows <- vector("list", n_nc1)
  for (d in seq_len(n_nc1)) {
    zz <- t4_one_pass(list(block), pr, list(ops), Zlist = list(t4_permute_genes(block$Z)))
    Z1[, d] <- zz$z_partial_k1; Z3[, d] <- zz$z_partial_k3
    a1 <- s1_frac(Z1[, d], pr$is_coord); a3 <- s1_frac(Z3[, d], pr$is_c3)
    okn <- is.finite(Z1[, d]) & pr$is_noncoord
    rows[[d]] <- data.table(set = s, draw = d, s1_k1 = a1[["s1"]], s1_k3 = a3[["s1"]], n_c1 = a1[["n"]],
                            n_c3 = a3[["n"]], base_noncoord = mean(sign(Z1[okn, d]) == sign(pr$z_disc_all[okn])))
  }
  nc1[[s]] <- rbindlist(rows)
  nc1_z[[s]] <- list(z1 = Z1, z3 = Z3)
  msg("NC1 %s: mean S1_k1 %.4f sd %.4f", s, mean(nc1[[s]]$s1_k1), sd(nc1[[s]]$s1_k1))
  # NC3: label test on every gene-permuted dataset with the stored permutations.
  nc3[[s]] <- rbindlist(lapply(seq_len(n_nc1), function(d) {
    z1 <- Z1[, d]; z3 <- Z3[, d]
    o1 <- s1_frac(z1, pr$is_coord)[["s1"]]; o3 <- s1_frac(z3, pr$is_c3)[["s1"]]
    q1 <- p_perm_upper(o1, perm_s1(z1, pr$is_coord, perms[[cn]], pidx, frozen))
    q3 <- p_perm_upper(o3, perm_s1(z3, pr$is_c3, perms[[cn]], pidx, frozen))
    data.table(set = s, draw = d, p_perm_k1 = q1, p_perm_k3 = q3, p_max = max(q1, q3),
               reject = max(q1, q3) < 0.05)
  }))
  msg("NC3 %s: %d of %d rejections (max p < 0.05)", s, sum(nc3[[s]]$reject), n_nc1)
}
NC1 <- rbindlist(nc1); NC3 <- rbindlist(nc3)
w(NC1, "NC1_gene_permutation_draws.tsv"); w(NC3, "NC3_label_test_under_gene_permutation.tsv")
saveRDS(nc1_z, file.path(OUT, "NC1_gene_permuted_z.rds"))
NC1s <- NC1[, .(n = .N, s1_k1_mean = mean(s1_k1), s1_k1_sd = sd(s1_k1), s1_k1_max = max(s1_k1),
                s1_k3_mean = mean(s1_k3), base_mean = mean(base_noncoord)), by = set]
NC3s <- NC3[, .(n = .N, n_reject = sum(reject), rate = mean(reject),
                n_reject_k1 = sum(p_perm_k1 < 0.05), n_reject_k3 = sum(p_perm_k3 < 0.05)), by = set]
w(NC1s, "NC1_summary.tsv"); w(NC3s, "NC3_summary.tsv"); print(NC1s); print(NC3s)

# ------------------------------------------------- G1 power on the holdout
scen <- data.table(scenario = c("all_m91", "all_m87", "all_m66", "all_m53", "F0F1_m66", "F0F1_m53"),
                   pool = c("all", "all", "all", "all", "F0F1", "F0F1"),
                   m = c(91L, 87L, 66L, 53L, 66L, 53L))
set.seed(H54_SEED)
g1_idx <- lapply(seq_len(nrow(scen)), function(k) {
  pool <- if (scen$pool[[k]] == "all") seq_len(nrow(H$meta)) else f01
  t(vapply(seq_len(n_g1), function(b) sample(pool, scen$m[[k]], replace = TRUE), integer(scen$m[[k]])))
})
cache <- new.env()
g1 <- list()
for (k in seq_len(nrow(scen))) {
  idx <- g1_idx[[k]]
  g1[[k]] <- rbindlist(lapply(seq_len(nrow(idx)), function(b) {
    blk <- make_block(H$PM, H$meta, idx[b, ], frozen)
    if (is.null(blk)) return(data.table(scenario = scen$scenario[[k]], m = scen$m[[k]], b = b, design_ok = FALSE))
    ops <- ops_for(blk, frozen, H$panel, cache)
    zz <- t4_one_pass(list(blk), pr, list(ops))
    ev <- evaluate_z(zz$z_partial_k1, zz$z_partial_k3, perms$H, pidx, frozen, n_perm_e)
    data.table(scenario = scen$scenario[[k]], m = scen$m[[k]], b = b, design_ok = TRUE,
               s1_k1 = ev$s1_k1, s1_k3 = ev$s1_k3, null_mean_k1 = ev$null_mean_k1,
               null_mean_k3 = ev$null_mean_k3, e_k1 = ev$e_k1, e_k3 = ev$e_k3,
               base_noncoord = ev$base_noncoord)
  }), fill = TRUE)
  msg("G1 %s done", scen$scenario[[k]])
}
G1 <- rbindlist(g1, fill = TRUE)
w(G1, "G1_holdout_draws.tsv.gz")
G1s <- G1[, .(n_draws = .N, n_design_failed = sum(!design_ok),
              e_k1_mean = mean(e_k1, na.rm = TRUE), e_k1_sd = sd(e_k1, na.rm = TRUE),
              e_k3_mean = mean(e_k3, na.rm = TRUE), e_k3_sd = sd(e_k3, na.rm = TRUE),
              s1_k1_sd = sd(s1_k1, na.rm = TRUE)), by = .(scenario, m)]
G1s[, `:=`(mde_k1 = MDE_MULT * e_k1_sd, mde_k3 = MDE_MULT * e_k3_sd)]
w(G1s, "G1_summary.tsv"); print(G1s)

E_all <- c(k1 = h_eval$e_k1, k3 = h_eval$e_k3)
E_f01 <- c(k1 = hf_eval$e_k1, k3 = hf_eval$e_k3)
g1_sd <- function(scen_name, col) {
  v <- G1s[G1s$scenario == scen_name][[col]]
  if (length(v) != 1L) fail("G1 summary lookup failed for ", scen_name)
  v
}
g1_dec <- rbindlist(lapply(c("A91", "B87", "A66", "B53"), function(s) {
  m_s <- EXPECTED_N[[s]]
  rbindlist(lapply(c("k1", "k3"), function(k) {
    sdcol <- paste0("e_", k, "_sd")
    if (SET_ROLE[[s]] == "primary") {
      sdv <- g1_sd(paste0("all_m", m_s), sdcol)
      comp <- E_all[[k]]; rule <- "MDE(all-stage pool) < holdout all-stage E"
    } else {
      sdv <- max(g1_sd(paste0("all_m", m_s), sdcol), g1_sd(paste0("F0F1_m", m_s), sdcol))
      comp <- min(E_all[[k]], E_f01[[k]]); rule <- "MDE(larger SD of two pools) < min(holdout all-stage E, holdout F0-F1 E)"
    }
    data.table(set = s, k = k, m = m_s, sd_e = sdv, mde = MDE_MULT * sdv, comparator_e = comp,
               passes = MDE_MULT * sdv < comp, rule = rule)
  }))
}))
stopifnot(nrow(g1_dec) == 8L)
w(g1_dec, "G1_decision.tsv"); print(g1_dec)

# ------------------------------------------------- stop rules
sr <- list()
add <- function(rule, scope, quantity, value, threshold, fired, note = "")
  sr[[length(sr) + 1L]] <<- data.table(rule = rule, scope = scope, quantity = quantity,
                                       value = as.character(value), threshold = threshold,
                                       fired = fired, note = note)
add(1, "all", "abs(S1_k1 - 0.897361912613355)", signif(pc1_s1_diff, 6), "< 1e-9", pc1_s1_diff >= 1e-9)
add(1, "all", "abs(comparator - 0.709472008490316)", signif(pc1_base_diff, 6), "< 1e-9", pc1_base_diff >= 1e-9)
for (s in all_sets)
  add(2, s, "max raw_stage_r2", signif(D4[set == s, max_raw_stage_r2], 6), "<= 1e-8",
      D4[set == s, max_raw_stage_r2] > 1e-8)
for (s in c("A91", "B87")) {
  v <- NC1s[set == s, s1_k1_mean]
  add(3, paste0("cohort_", SET_COHORT[[s]]), paste0("NC1 mean S1_k1 (", s, ")"), round(v, 6), "0.44-0.54",
      v < 0.44 || v > 0.54)
}
g2min <- c(primary = 60L, secondary = 40L)
for (s in c("A91", "B87", "A66", "B53")) {
  n <- D1[set == s, n]; df <- D1[set == s, residual_df]
  g2 <- n >= g2min[[SET_ROLE[[s]]]] && df >= 30
  g3 <- D2[set == s, programs_testable] >= 100 && D3[set == s, coord_eligible_k1] >= 2000
  add(4, s, "G2 n / residual df", paste0(n, " / ", df), paste0("n >= ", g2min[[SET_ROLE[[s]]]], ", df >= 30"), !g2)
  add(4, s, "G3 programs testable / coordinated edges eligible at k1",
      paste0(D2[set == s, programs_testable], " / ", D3[set == s, coord_eligible_k1]), ">= 100 / >= 2000", !g3,
      "edges counted by testability and leave-shared-out eligibility (no z); stage 2 re-checks finite count")
  add(5, s, "G1 k1: MDE vs comparator E", paste0(round(g1_dec[g1_dec$set == s & g1_dec$k == "k1"]$mde, 4), " vs ",
                                                  round(g1_dec[g1_dec$set == s & g1_dec$k == "k1"]$comparator_e, 4)),
      "MDE < E", !g1_dec[g1_dec$set == s & g1_dec$k == "k1"]$passes, "fired at k1 = set untestable")
  add(5, s, "G1 k3: MDE vs comparator E", paste0(round(g1_dec[g1_dec$set == s & g1_dec$k == "k3"]$mde, 4), " vs ",
                                                  round(g1_dec[g1_dec$set == s & g1_dec$k == "k3"]$comparator_e, 4)),
      "MDE < E", !g1_dec[g1_dec$set == s & g1_dec$k == "k3"]$passes, "fired at k3 only = set tested at k1 only")
}
for (s in new_sets)
  add(6, s, "fraction of stored bootstrap draws with rank-deficient design",
      signif(D9[set == s, fraction_failed], 6), "<= 0.05", D9[set == s, fraction_failed] > 0.05)
add(7, "A66", "primary A91 passed", "evaluated in stage 2", "A91 passes", NA, "secondary opened only if own primary passed")
add(7, "B53", "primary B87 passed", "evaluated in stage 2", "B87 passes", NA, "secondary opened only if own primary passed")
for (s in c("A91", "B87", "H105")) {
  v <- NC3s[set == s, n_reject]
  add(8, if (s == "H105") "H105" else paste0("cohort_", SET_COHORT[[s]]),
      paste0("NC3 rejections of ", n_nc1, " (", s, ")"), v, paste0("<= ", floor(0.08 * n_nc1)), v > floor(0.08 * n_nc1),
      if (s == "H105") "holdout calibration, reported only" else "")
}
add("PC2", "primary_verdict", "H105 p_perm_k1 (label null)", signif(p2, 6), "< 0.05", !pc2_validated,
    "fired = primary verdict 'null not validated'")
add("PC2", "H105", "label null mean (1000 perms)", round(mean(h_eval$null_k1), 6), "0.45-0.70 predicted", NA, "prediction only")
add("PC3", "H105", "p_perm_k3", signif(p3, 6), "< 0.05 predicted", NA, "reported, not a stop rule")
add("PC1b", "H105", "S1_k1 within-cohort TMM minus pooled", signif(pc1b_shift, 6), "|shift| <= 0.02 else printed next to every S1", NA,
    "descriptive")
add("diag11", "cohort_A", "composition arm readable (hepatocyte median, zero fraction)",
    paste0(round(comp_qc[cohort == "A", hepatocyte_median], 4), " / ", round(comp_qc[cohort == "A", max_lineage_zero_fraction], 4)),
    "0.60-0.95 / <= 0.05", !comp_qc[cohort == "A", composition_arm_readable], "fired = composition arm not read")
add("diag11", "cohort_B", "composition arm readable (hepatocyte median, zero fraction)",
    paste0(round(comp_qc[cohort == "B", hepatocyte_median], 4), " / ", round(comp_qc[cohort == "B", max_lineage_zero_fraction], 4)),
    "0.60-0.95 / <= 0.05", !comp_qc[cohort == "B", composition_arm_readable], "fired = composition arm not read")
SR <- rbindlist(sr)
w(SR, "stop_rule.tsv")
print(SR)

# Per-set clearance for stage 2 (stage 2 recomputes this from stop_rule.tsv).
clear <- rbindlist(lapply(new_sets, function(s) {
  cs <- paste0("cohort_", SET_COHORT[[s]])
  f <- function(r, sc) any(SR$fired[SR$rule == as.character(r) & SR$scope == sc] %in% TRUE)
  r1 <- f(1, "all"); r2 <- f(2, s); r3 <- f(3, cs); r4 <- f(4, s); r6 <- f(6, s); r8 <- f(8, cs)
  g1k1 <- if (s %in% g1_dec$set) g1_dec[g1_dec$set == s & g1_dec$k == "k1"]$passes else NA
  g1k3 <- if (s %in% g1_dec$set) g1_dec[g1_dec$set == s & g1_dec$k == "k3"]$passes else NA
  stopped <- r1 || r2 || r3 || (SET_ROLE[[s]] != "descriptive" && (r4 || isFALSE(g1k1)))
  ks <- if (stopped) "none" else if (SET_ROLE[[s]] == "descriptive") "k1,k3 (descriptive)" else
    if (isTRUE(g1k3)) "k1,k3" else "k1 only"
  data.table(set = s, role = SET_ROLE[[s]], rule1 = r1, rule2 = r2, rule3 = r3, rule4 = r4,
             g1_k1 = g1k1, g1_k3 = g1k3, rule6_cannot_pass = r6, rule8_null_not_validated = r8,
             pc2_null_not_validated = !pc2_validated, cleared_for_stage2 = !stopped, tested_k = ks,
             secondary_requires_primary_pass = SET_ROLE[[s]] == "secondary")
}))
w(clear, "stage2_clearance.tsv"); print(clear)

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
fwrite(data.table(file = list.files(file.path(H54_RUN, "code"), full.names = TRUE))[
  , sha256 := vapply(file, sha256_of, "")], file.path(OUT, "code_sha256.tsv"), sep = "\t")
fwrite(data.table(file = list.files(OUT, full.names = TRUE))[, sha256 := vapply(file, sha256_of, "")],
       file.path(OUT, "..", paste0(basename(OUT), "_outputs_sha256.tsv")), sep = "\t")
msg("STAGE1 COMPLETE in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
