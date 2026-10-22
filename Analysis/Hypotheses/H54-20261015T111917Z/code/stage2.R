#!/usr/bin/env Rscript
# H54 stage 2 (job 3, the outcome stage). Submitted by hand only after
# stage1/stop_rule.tsv has been read. Implements PRESPEC_H54.md "Primary statistic and
# test", "Testing family and correction", the descriptive outputs and "Reading".
#
# Usage: stage2.R [out_dir] [stage1_dir]
# Debug (never the real run): H54_DEBUG_GENE_PERMUTED=1 permutes participants
# independently per gene in every set before anything is computed, so no real
# outcome is seen; H54_DEBUG_DIVISOR=k uses the first 1/k of the stored bootstrap
# draws; H54_DEBUG_FORCE_OPEN_SECONDARY=1 opens the secondary sets regardless of
# the primary result so their code path is exercised.

source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z/code/h54_lib.R")
args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) >= 1L) args[[1]] else file.path(H54_RUN, "stage2")
S1DIR <- if (length(args) >= 2L) args[[2]] else file.path(H54_RUN, "stage1")
DEBUG_PERM <- Sys.getenv("H54_DEBUG_GENE_PERMUTED", "0") == "1"
DIV <- as.integer(Sys.getenv("H54_DEBUG_DIVISOR", "1"))
FORCE_SEC <- Sys.getenv("H54_DEBUG_FORCE_OPEN_SECONDARY", "0") == "1"
REAL <- normalizePath(dirname(OUT), mustWork = TRUE) == normalizePath(H54_RUN) && basename(OUT) == "stage2"
if (REAL && (DEBUG_PERM || DIV != 1L || FORCE_SEC)) fail("debug switches are not allowed for the real stage2 directory")
if (!REAL && !DEBUG_PERM) fail("a non-real stage 2 output directory is allowed only on gene-permuted data")
RESULT <- if (REAL) file.path(H54_RUN, "RESULT_H54.md") else file.path(OUT, "RESULT_H54_DEBUG.md")
t0 <- Sys.time()

verify_prespec_and_inputs()
# Stage 1 outputs must be the ones stage 1 wrote.
s1sha <- fread(file.path(dirname(S1DIR), paste0(basename(S1DIR), "_outputs_sha256.tsv")))
for (k in seq_len(nrow(s1sha))) if (sha256_of(s1sha$file[[k]]) != s1sha$sha256[[k]])
  fail("stage 1 output changed: ", s1sha$file[[k]])
msg("stage 1 outputs verified (%d files)", nrow(s1sha))
if (file.exists(RESULT)) fail("Refusing to overwrite ", RESULT)
refuse_if_exists(OUT)
w <- function(x, f) fwrite(x, file.path(OUT, f), sep = "\t", na = "NA")

SR <- fread(file.path(S1DIR, "stop_rule.tsv"), colClasses = list(character = c("rule", "value")))
CL <- fread(file.path(S1DIR, "stage2_clearance.tsv"))
G1D <- fread(file.path(S1DIR, "G1_decision.tsv"))
PCS <- fread(file.path(S1DIR, "controls_PC_H105.tsv"))
CQC <- fread(file.path(S1DIR, "diag11_composition_qc.tsv"))
LP <- readRDS(file.path(S1DIR, "label_permutations.rds"))
BI <- readRDS(file.path(S1DIR, "bootstrap_indices.rds"))
COV <- fread(file.path(S1DIR, "composition_covariates_by_set.tsv"))
if (any(SR[rule == "1", fired] %in% TRUE)) {
  writeLines(c("# H54 result", "", "Verdict: UNTESTABLE. Stop rule 1 fired in stage 1 (holdout code path not reproduced).",
               "No new-cohort statistic was computed."), RESULT)
  fail("stop rule 1 fired; nothing computed")
}
n_boot_use <- nrow(BI$indices[[1]]) %/% DIV
n_perm_e <- LP$n_perm_e
if (REAL && (LP$n_perm != N_PERM || nrow(BI$indices[[1]]) != N_BOOT))
  fail("stage 1 was not run with the full draw counts")

frozen <- load_frozen()
pr <- frozen$pairs
pidx <- pair_index_matrix(frozen)
features <- frozen$features
sm <- build_set_meta()
countsA <- read_counts(file.path(H54_ROOT, "data/GSE281797/GSE281797_counts_ensembl_stable.tsv.gz"))
countsB <- read_counts(file.path(H54_ROOT, "data/GSE225740/GSE225740_counts_ensembl_stable.tsv.gz"))
new_sets <- c("A91", "B87", "A66", "B53", "A94", "A37")
S <- list()
for (s in new_sets) {
  meta <- sm$sets[[s]]
  if (!identical(meta$sample_id, BI$sample_order[[s]])) fail("sample order differs from stage 1 for ", s)
  S[[s]] <- c(list(meta = meta), build_set_expression(if (SET_COHORT[[s]] == "A") countsA else countsB, meta, frozen))
}
rm(countsA, countsB); invisible(gc())
if (DEBUG_PERM) {
  set.seed(H54_SEED + 54L)
  for (s in new_sets) {
    PM <- S[[s]]$PM
    S[[s]]$PM <- t(apply(PM, 1L, function(v) v[sample.int(length(v))]))
    dimnames(S[[s]]$PM) <- dimnames(PM)
  }
  msg("DEBUG: every set gene-permuted (participants shuffled independently per gene)")
}

COMP_FORMULA <- ~ sex_final + fibrosis_factor + nas_factor + composition_pc1 + composition_pc2 +
  composition_pc3 + log_library_size + top100_share
with_formula <- function(f, expr) {
  old <- T4_RESIDUAL_FORMULA
  assign("T4_RESIDUAL_FORMULA", f, envir = globalenv())
  on.exit(assign("T4_RESIDUAL_FORMULA", old, envir = globalenv()))
  force(expr)
}

# --------------------------------------------------------- core per-set work
observe <- function(s, idx = NULL, formula = NULL, meta_override = NULL) {
  meta <- if (is.null(meta_override)) S[[s]]$meta else meta_override
  if (is.null(idx)) idx <- seq_len(nrow(meta))
  blk <- if (is.null(formula)) make_block(S[[s]]$PM, meta, idx, frozen) else
    with_formula(formula, make_block(S[[s]]$PM, meta, idx, frozen))
  if (is.null(blk)) return(NULL)
  zz <- t4_one_pass(list(blk), pr, list(ops_for(blk, frozen)))
  list(block = blk, z1 = zz$z_partial_k1, z3 = zz$z_partial_k3)
}

boot_set <- function(s) {
  cn <- SET_COHORT[[s]]
  idx <- BI$indices[[s]][seq_len(n_boot_use), , drop = FALSE]
  cache <- new.env()
  rbindlist(lapply(seq_len(nrow(idx)), function(b) {
    blk <- make_block(S[[s]]$PM, S[[s]]$meta, idx[b, ], frozen)
    if (is.null(blk)) return(data.table(set = s, b = b, design_ok = FALSE))
    zz <- t4_one_pass(list(blk), pr, list(ops_for(blk, frozen, S[[s]]$panel, cache)))
    ev <- evaluate_z(zz$z_partial_k1, zz$z_partial_k3, LP$perms[[cn]], pidx, frozen, n_perm_e)
    data.table(set = s, b = b, design_ok = TRUE, s1_k1 = ev$s1_k1, s1_k3 = ev$s1_k3,
               null_mean_k1 = ev$null_mean_k1, null_mean_k3 = ev$null_mean_k3,
               e_k1 = ev$e_k1, e_k3 = ev$e_k3, n_c1 = ev$n_c1, n_c3 = ev$n_c3)
  }), fill = TRUE)
}

run_set <- function(s) {
  cn <- SET_COHORT[[s]]
  ob <- observe(s)
  if (is.null(ob)) fail("full-set design failed for ", s)
  ev <- evaluate_z(ob$z1, ob$z3, LP$perms[[cn]], pidx, frozen, n_perm_e, all_perms = TRUE)
  bd <- boot_set(s)
  ok <- bd[design_ok == TRUE]
  st <- data.table(
    set = s, role = SET_ROLE[[s]], n = nrow(S[[s]]$meta),
    n_c1 = ev$n_c1, n_c3 = ev$n_c3,
    s1_k1 = ev$s1_k1, s1_k3 = ev$s1_k3,
    null_mean_k1_1000 = mean(ev$null_k1), null_sd_k1 = sd(ev$null_k1), null_q95_k1 = quantile(ev$null_k1, 0.95),
    null_mean_k3_1000 = mean(ev$null_k3), null_sd_k3 = sd(ev$null_k3), null_q95_k3 = quantile(ev$null_k3, 0.95),
    null_mean_k1_200 = ev$null_mean_k1, null_mean_k3_200 = ev$null_mean_k3,
    e_k1 = ev$e_k1, e_k3 = ev$e_k3,
    p_perm_k1 = p_perm_upper(ev$s1_k1, ev$null_k1), p_perm_k3 = p_perm_upper(ev$s1_k3, ev$null_k3),
    b_total = nrow(bd), b_ok = nrow(ok),
    p_boot_k1 = (1 + sum(ok$e_k1 <= 0)) / (nrow(ok) + 1), p_boot_k3 = (1 + sum(ok$e_k3 <= 0)) / (nrow(ok) + 1),
    e_k1_ci_lo = quantile(ok$e_k1, 0.025), e_k1_ci_hi = quantile(ok$e_k1, 0.975),
    e_k3_ci_lo = quantile(ok$e_k3, 0.025), e_k3_ci_hi = quantile(ok$e_k3, 0.975),
    s1_k1_ci_lo = quantile(ok$s1_k1, 0.025), s1_k1_ci_hi = quantile(ok$s1_k1, 0.975),
    s1_k3_ci_lo = quantile(ok$s1_k3, 0.025), s1_k3_ci_hi = quantile(ok$s1_k3, 0.975),
    base_noncoord = ev$base_noncoord, n_noncoord = ev$n_noncoord)
  list(stat = st, boot = bd, z1 = ob$z1, z3 = ob$z3, null_k1 = ev$null_k1, null_k3 = ev$null_k3,
       block = ob$block)
}

clear_row <- function(s) CL[set == s]

# --------------------------------------------------------- primary family
R <- list()
fam_row <- function(s, opened = TRUE) {
  cr <- clear_row(s)
  base <- data.table(set = s, role = SET_ROLE[[s]], opened = opened, cleared = cr$cleared_for_stage2,
                     tested_k = cr$tested_k, cannot_pass = cr$rule6_cannot_pass | cr$rule8_null_not_validated)
  if (!opened || !cr$cleared_for_stage2) {
    return(cbind(base, data.table(p_perm_k1 = NA_real_, p_perm_k3 = NA_real_, p_boot_k1 = NA_real_,
                                  p_boot_k3 = NA_real_, p_star_raw = NA_real_, p_star = 1, g3_finite_c1 = NA)))
  }
  R[[s]] <<- run_set(s)
  st <- R[[s]]$stat
  g3 <- st$n_c1 >= 2000L
  pr_raw <- if (cr$tested_k == "k1,k3") max(st$p_perm_k1, st$p_perm_k3, st$p_boot_k1, st$p_boot_k3) else
    max(st$p_perm_k1, st$p_boot_k1)
  pstar <- if (base$cannot_pass || !g3) 1 else pr_raw
  cbind(base, data.table(p_perm_k1 = st$p_perm_k1, p_perm_k3 = st$p_perm_k3, p_boot_k1 = st$p_boot_k1,
                         p_boot_k3 = st$p_boot_k3, p_star_raw = pr_raw, p_star = pstar, g3_finite_c1 = g3))
}

holm2 <- function(fam) {
  fam <- copy(fam)
  o <- order(fam$p_star, fam$set)
  fam[, `:=`(holm_rank = NA_integer_, alpha_level = NA_real_, rejected = FALSE)]
  fam[o[1], `:=`(holm_rank = 1L, alpha_level = 0.025, rejected = p_star <= 0.025)]
  first <- fam$rejected[o[1]]
  fam[o[2], `:=`(holm_rank = 2L, alpha_level = if (first) 0.05 else 0.025,
                 rejected = first && p_star <= 0.05)]
  fam
}

bound_at <- function(fam) {
  fam[, `:=`(e_k1_lower_bound = NA_real_, e_k3_lower_bound = NA_real_)]
  for (r in seq_len(nrow(fam))) {
    s <- fam$set[[r]]
    if (!is.null(R[[s]])) {
      ok <- R[[s]]$boot[design_ok == TRUE]
      fam[r, `:=`(e_k1_lower_bound = quantile(ok$e_k1, alpha_level),
                  e_k3_lower_bound = quantile(ok$e_k3, alpha_level))]
    }
  }
  fam
}

prim <- rbindlist(lapply(c("A91", "B87"), fam_row))
prim <- bound_at(holm2(prim))
w(prim, "primary_family.tsv"); print(prim)

# --------------------------------------------------------- secondary family
own_primary <- c(A66 = "A91", B53 = "B87")
sec <- rbindlist(lapply(c("A66", "B53"), function(s) {
  opened <- FORCE_SEC || isTRUE(prim[set == own_primary[[s]], rejected])
  fam_row(s, opened)
}))
sec <- bound_at(holm2(sec))
w(sec, "secondary_family.tsv"); print(sec)

# --------------------------------------------------------- descriptive sets
# A94 is the sensitivity version of A91. A37 is read only beside A66, so it is
# computed only when A66 was opened.
desc_sets <- c("A94", if (!is.null(R[["A66"]])) "A37")
for (s in desc_sets) if (clear_row(s)$cleared_for_stage2) R[[s]] <- run_set(s)
STAT <- rbindlist(lapply(R, `[[`, "stat"), fill = TRUE)
w(STAT, "set_statistics.tsv"); w(rbindlist(lapply(R, `[[`, "boot"), fill = TRUE), "bootstrap_draws.tsv.gz")
for (s in names(R)) w(data.table(pair_index = pr$pair_index, program_a = pr$program_a, program_b = pr$program_b,
                                 z_k1 = R[[s]]$z1, z_k3 = R[[s]]$z3), paste0("edge_z_", s, ".tsv.gz"))
print(STAT)

# --------------------------------------------------------- per-set descriptives
desc <- list(); lin_split <- list(); lin_pair <- list(); size <- list(); comp <- list(); nc2 <- list()
s23 <- list()
part <- setNames(frozen$communities$community, frozen$communities$feature_id)
mod_of <- function(edges, assign) {
  e <- edges[program_a %in% names(assign) & program_b %in% names(assign)]
  if (nrow(e) < 10L) return(NA_real_)
  ca <- assign[e$program_a]; cb <- assign[e$program_b]
  wv <- e$z; tot <- sum(wv)
  ki <- tapply(c(wv, wv), c(ca, cb), sum)[as.character(sort(unique(assign)))]
  ki[is.na(ki)] <- 0
  sum(wv[ca == cb]) / tot - sum((ki / (2 * tot))^2)
}
rewire_once <- function(e, n_swaps) {
  a <- e$program_a; b <- e$program_b; m <- length(a)
  for (k in seq_len(n_swaps)) {
    ij <- sample.int(m, 2L); i <- ij[[1L]]; j <- ij[[2L]]
    if (a[i] == a[j] || b[i] == b[j] || a[i] == b[j] || b[i] == a[j]) next
    tmp <- b[i]; b[i] <- b[j]; b[j] <- tmp
  }
  data.table(program_a = a, program_b = b, z = e$z)
}
for (s in names(R)) {
  cn <- SET_COHORT[[s]]
  z1 <- R[[s]]$z1; z3 <- R[[s]]$z3
  c1 <- pr$is_coord & is.finite(z1)
  # Comparators, two numbers each, never the difference alone.
  nonc <- which(pr$is_noncoord & is.finite(z1))
  conc <- sign(z1[nonc]) == sign(pr$z_disc_all[nonc])
  set.seed(H54_SEED)
  edge_null <- replicate(1000L, mean(conc[sample.int(length(nonc), min(sum(c1), length(nonc)))]))
  desc[[s]] <- data.table(set = s, s1_k1 = mean(z1[c1] > 0), noncoord_concordance = mean(conc),
                          n_noncoord = length(nonc), edge_draw_null_mean = mean(edge_null),
                          edge_draw_null_sd = sd(edge_null), edge_draw_size = min(sum(c1), length(nonc)),
                          edge_draw_with_replacement = FALSE)
  lin_split[[s]] <- data.table(set = s, edges = c("same_lineage", "cross_lineage"),
                               n = c(sum(c1 & pr$same_lineage), sum(c1 & !pr$same_lineage)),
                               s1_k1 = c(mean(z1[c1 & pr$same_lineage] > 0), mean(z1[c1 & !pr$same_lineage] > 0)))
  lp <- pr[c1, .(lineage_pair = paste(pmin(lineage_a, lineage_b), pmax(lineage_a, lineage_b), sep = "|"), z = z1[c1])]
  lin_pair[[s]] <- cbind(set = s, lp[, .(n = .N, s1_k1 = mean(z > 0)), by = lineage_pair][order(lineage_pair)])
  # Size-matched label null.
  ns1 <- perm_s1(z1, pr$is_coord, LP$perms_size_matched[[cn]], pidx, frozen)
  ns3 <- perm_s1(z3, pr$is_c3, LP$perms_size_matched[[cn]], pidx, frozen)
  size[[s]] <- data.table(set = s, s1_k1 = R[[s]]$stat$s1_k1, size_null_mean_k1 = mean(ns1),
                          p_size_k1 = p_perm_upper(R[[s]]$stat$s1_k1, ns1),
                          s1_k3 = R[[s]]$stat$s1_k3, size_null_mean_k3 = mean(ns3),
                          p_size_k3 = p_perm_upper(R[[s]]$stat$s1_k3, ns3))
  # NC2: informative-null pairs through the same statistic and label null (k = 1).
  inf_obs <- s1_frac(z1, pr$is_infnull)
  inf_null <- perm_s1(z1, pr$is_infnull, LP$perms[[cn]], pidx, frozen)
  nc2[[s]] <- data.table(set = s, n_informative_null = inf_obs[["n"]], s1_k1 = inf_obs[["s1"]],
                         null_mean = mean(inf_null), p_perm_k1 = p_perm_upper(inf_obs[["s1"]], inf_null),
                         nc2_passes_test = p_perm_upper(inf_obs[["s1"]], inf_null) < 0.05)
  # Composition and library-complexity arm.
  cv <- COV[set == s, !"set"]
  mm <- merge(S[[s]]$meta, cv, by = "sample_id", sort = FALSE)
  mm <- mm[match(S[[s]]$meta$sample_id, sample_id)]
  ob <- observe(s, formula = COMP_FORMULA, meta_override = mm)
  if (is.null(ob)) {
    comp[[s]] <- data.table(set = s, design_ok = FALSE)
  } else {
    evc <- evaluate_z(ob$z1, ob$z3, LP$perms[[cn]], pidx, frozen, n_perm_e, all_perms = TRUE)
    comp[[s]] <- data.table(set = s, design_ok = TRUE, design_columns = ob$block$n_par, residual_df = ob$block$df,
                            s1_k1 = evc$s1_k1, s1_k3 = evc$s1_k3, e_k1 = evc$e_k1, e_k3 = evc$e_k3,
                            p_perm_k1 = p_perm_upper(evc$s1_k1, evc$null_k1),
                            p_perm_k3 = p_perm_upper(evc$s1_k3, evc$null_k3),
                            composition_qc_readable = CQC[cohort == cn, composition_arm_readable])
  }
  # S2 modularity (frozen partition) and S3 community-3 density.
  tst <- R[[s]]$block$coverage[cohort_testable == TRUE, feature_id]
  g_edges <- pr[is.finite(z1) & z1 > 0, .(program_a, program_b)][, z := z1[is.finite(z1) & z1 > 0]]
  ptst <- part[names(part) %in% tst]
  s2 <- mod_of(g_edges, ptst)
  set.seed(H54_SEED)
  null_rank <- replicate(200L, mod_of(g_edges, setNames(sample(ptst), names(ptst))))
  set.seed(H54_SEED)
  null_rew <- replicate(200L, mod_of(rewire_once(g_edges, 10L * nrow(g_edges)), ptst))
  c3 <- intersect(frozen$communities[community == 3, feature_id], tst)
  dens <- function(set) {
    e <- g_edges[program_a %in% set & program_b %in% set]
    nrow(e) / (length(set) * (length(set) - 1) / 2)
  }
  s3 <- dens(c3)
  lin <- frozen$lineage[tst]
  c3_lin <- table(lin[c3])
  set.seed(H54_SEED)
  null_c3 <- replicate(1000L, dens(unlist(lapply(names(c3_lin), function(L) {
    pool <- names(lin)[lin == L]; sample(pool, min(c3_lin[[L]], length(pool)))
  }))))
  s23[[s]] <- data.table(set = s,
    statistic = c("S2a_modularity_vs_label_permutation", "S2b_modularity_vs_rewiring", "S3_community3_density"),
    observed = c(s2, s2, s3),
    null_mean = c(mean(null_rank, na.rm = TRUE), mean(null_rew, na.rm = TRUE), mean(null_c3)),
    null_sd = c(sd(null_rank, na.rm = TRUE), sd(null_rew, na.rm = TRUE), sd(null_c3)),
    p_value = c((1 + sum(null_rank >= s2, na.rm = TRUE)) / (1 + sum(is.finite(null_rank))),
                (1 + sum(null_rew >= s2, na.rm = TRUE)) / (1 + sum(is.finite(null_rew))),
                (1 + sum(null_c3 >= s3)) / 1001),
    n_programs_community3 = length(c3),
    note = if (length(c3) != 15L) "reduced community 3; not compared with the holdout 0.876" else "")
}
w(rbindlist(desc), "comparators.tsv"); w(rbindlist(lin_split), "lineage_split.tsv")
w(rbindlist(lin_pair), "lineage_pair.tsv"); w(rbindlist(size), "size_matched_null.tsv")
NC2 <- rbindlist(nc2); w(NC2, "NC2_informative_null.tsv")
COMPA <- rbindlist(comp, fill = TRUE); w(COMPA, "composition_arm.tsv")
w(rbindlist(s23), "S2_S3_structure.tsv")

# F2-or-higher strata, stage pattern only (columns of the primary set; no bootstrap).
f2 <- rbindlist(lapply(intersect(c("A91", "B87"), names(R)), function(s) {
  idx <- which(S[[s]]$meta$fibrosis_stage >= 2L)
  ob <- observe(s, idx = idx)
  if (is.null(ob)) return(data.table(set = s, stratum = "F2plus", n = length(idx), design_ok = FALSE))
  ev <- evaluate_z(ob$z1, ob$z3, LP$perms[[SET_COHORT[[s]]]], pidx, frozen, n_perm_e)
  data.table(set = s, stratum = "F2plus", n = length(idx), design_ok = TRUE, s1_k1 = ev$s1_k1, s1_k3 = ev$s1_k3,
             null_mean_k1 = ev$null_mean_k1, null_mean_k3 = ev$null_mean_k3, e_k1 = ev$e_k1, e_k3 = ev$e_k3,
             base_noncoord = ev$base_noncoord)
}), fill = TRUE)
w(f2, "F2plus_strata.tsv")

# --------------------------------------------------------- verdict (computed)
pc2_ok <- !any(SR[rule == "PC2" & scope == "primary_verdict", fired] %in% TRUE)
pc1b_shift <- as.numeric(SR[rule == "PC1b", value])
nc2_fail <- function(s) isTRUE(NC2[set == s, nc2_passes_test])
set_reading <- function(fam, s) {
  r <- fam[set == s]
  cr <- clear_row(s)
  if (!r$opened) return("not opened (stop rule 7: own primary did not pass)")
  if (!cr$cleared_for_stage2) {
    why <- c(if (cr$rule1) "1", if (cr$rule2) "2", if (cr$rule3) "3", if (cr$rule4) "4",
             if (isFALSE(cr$g1_k1)) "5 (G1 at k = 1)")
    return(paste0("untestable (stop rule ", paste(why, collapse = ", "), ")"))
  }
  if (!isTRUE(r$g3_finite_c1)) return("untestable (G3: fewer than 2,000 finite coordinated edges at k = 1)")
  if (isTRUE(cr$rule8_null_not_validated)) return("null not validated (stop rule 8); cannot pass")
  if (isTRUE(cr$rule6_cannot_pass)) return("participant condition untestable (stop rule 6); cannot pass")
  if (!pc2_ok) return(if (r$rejected) "null not validated (PC2); the p would pass but neither reading applies" else
    "null not validated (PC2)")
  if (r$rejected && nc2_fail(s)) return("pass, but NC2 passed in this cohort: not interpretable")
  if (r$rejected && cr$tested_k == "k1,k3") return("pass at k = 1 and k = 3")
  if (r$rejected) return("pass at k = 1 only (set reduced by stop rule 5)")
  st <- R[[s]]$stat
  k1_ok <- max(st$p_perm_k1, st$p_boot_k1) <= r$alpha_level
  if (cr$tested_k == "k1,k3" && k1_ok) return("fail: k = 1 parts pass at the set's Holm level, k = 3 parts do not")
  "fail"
}
prim[, reading := vapply(set, function(s) set_reading(prim, s), "")]
sec[, reading := vapply(set, function(s) set_reading(sec, s), "")]
pass_both <- prim[reading == "pass at k = 1 and k = 3", set]
pass_k1only <- prim[grepl("^pass at k = 1 only|k = 1 parts pass", reading), set]
untestable <- prim[grepl("^untestable|cannot pass|not validated", reading), set]
overall <- if (!pc2_ok) {
  "NULL NOT VALIDATED (PC2 failed): the primary was computed and is reported, but neither the pass nor the null reading applies."
} else if (length(pass_both) == 2L) {
  "REPLICATED IN BOTH PRIMARY SETS at k = 1 and k = 3 (A91 and B87)."
} else if (length(pass_both) == 1L) {
  paste0("REPLICATED IN ONE PRIMARY SET (", pass_both, ") at k = 1 and k = 3; the other set did not pass.")
} else if (length(pass_k1only)) {
  paste0("RESIDUAL COVARIANCE REPLICATES AT k = 1 ONLY (", paste(pass_k1only, collapse = ", "),
         "); this is composition-linked covariance, not coordination.")
} else if (length(untestable) == 2L) {
  "UNTESTABLE: stop rules ended both primary sets."
} else if (length(untestable) == 1L) {
  paste0("NULL in the testable primary set; ", untestable, " untestable.")
} else {
  "NULL: neither primary set passed (G1 passed and controls passed); bounds only an excess E_1 above the MDE."
}
fmt <- function(x, d = 4) ifelse(is.na(x), "NA", trimws(formatC(x, digits = d, format = "f")))
fmtp <- function(x) ifelse(is.na(x), "NA", trimws(formatC(x, digits = 6, format = "g")))
line_set <- function(s) {
  st <- R[[s]]$stat
  if (is.null(st)) return(paste0("| ", s, " | not computed |  |  |  |  |  |"))
  paste0("| ", s, " | ", st$n, " | ", fmt(st$s1_k1), " / ", fmt(st$s1_k3), " | ",
         fmt(st$null_mean_k1_200), " / ", fmt(st$null_mean_k3_200), " | ",
         fmt(st$e_k1), " [", fmt(st$e_k1_ci_lo), ", ", fmt(st$e_k1_ci_hi), "] | ",
         fmt(st$e_k3), " [", fmt(st$e_k3_ci_lo), ", ", fmt(st$e_k3_ci_hi), "] | ",
         fmtp(st$p_perm_k1), " / ", fmtp(st$p_perm_k3), " / ", fmtp(st$p_boot_k1), " / ", fmtp(st$p_boot_k3), " |")
}
fam_lines <- function(fam) c("| set | p_perm k1 / k3 | p_boot k1 / k3 | p* | Holm rank, level | rejected | E_k1 / E_k3 lower bound at level | reading |",
  "|---|---|---|---|---|---|---|---|",
  vapply(seq_len(nrow(fam)), function(r) with(fam[r], paste0(
    "| ", set, " | ", fmtp(p_perm_k1), " / ", fmtp(p_perm_k3), " | ", fmtp(p_boot_k1), " / ", fmtp(p_boot_k3),
    " | ", fmtp(p_star), " | ", holm_rank, ", ", alpha_level, " | ", rejected, " | ",
    fmt(e_k1_lower_bound), " / ", fmt(e_k3_lower_bound), " | ", reading, " |")), ""), "")
comp_reading <- vapply(intersect(pass_both, COMPA$set), function(s) {
  r <- COMPA[set == s]
  if (!isTRUE(r$composition_qc_readable)) return(paste0(s, ": composition arm fails the frozen acceptance rule; not read."))
  if (isTRUE(r$design_ok) && r$p_perm_k1 < 0.05 && r$p_perm_k3 < 0.05)
    paste0(s, ": the result held after adding three estimated composition components, log library size and top-100 count share to the design (descriptive).")
  else paste0(s, ": part of this covariance follows estimated cell-type abundance or library complexity.")
}, "")
size_reading <- vapply(intersect(pass_both, rbindlist(size)$set), function(s) {
  r <- rbindlist(size)[set == s]
  if (r$p_size_k1 >= 0.05 || r$p_size_k3 >= 0.05) paste0(s, ": size-matched null p >= 0.05; the excess may follow program size within lineage.") else
    paste0(s, ": size-matched null also rejected (p_k1 ", fmtp(r$p_size_k1), ", p_k3 ", fmtp(r$p_size_k3), ").")
}, "")
res <- c(
  if (!REAL) c("# DEBUG RUN ON GENE-PERMUTED DATA. NOT A RESULT.", ""),
  "# H54 result",
  "",
  paste0("Computed ", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), " by code/stage2.R. Prespec sha256 ", H54_PRESPEC_SHA, "."),
  "",
  paste0("**Verdict (primary): ", overall, "**"),
  "",
  "Unit of replication: participant (one biopsy each). Design: ~ sex_final + fibrosis_factor + nas_factor (NAS bins le2/3to4/ge5).",
  "Primary family: A91, B87, Holm at family-wise 0.05 on p* = max(p_perm_k1, p_perm_k3, p_boot_k1, p_boot_k3) (k1 parts only for a set reduced by stop rule 5).",
  "Secondary family: A66, B53, same rule; an unopened test has p* = 1. Everything else is descriptive. A set that cannot pass (stop rules 6, 8, or G3 at stage 2) enters Holm with p* = 1.",
  if (abs(pc1b_shift) > 0.02) paste0("PC1b: within-cohort TMM moved the holdout S1_k1 by ", fmt(pc1b_shift),
                                    " (more than 0.02); every new-cohort S1 below uses only that path.") else
    paste0("PC1b: within-cohort TMM moved the holdout S1_k1 by ", fmt(pc1b_shift), "."),
  "",
  "## Primary family", "", fam_lines(prim),
  "## Secondary family (cross-sectional F0-F1 strata)", "", fam_lines(sec),
  "## Per-set statistics (S1_k1 / S1_k3; label-null mean over permutations 1-200; E with participant-bootstrap 95% percentile interval; p_perm k1 / k3 / p_boot k1 / k3)", "",
  "| set | n | S1 k1 / k3 | null mean k1 / k3 | E_k1 [95%] | E_k3 [95%] | p |",
  "|---|---|---|---|---|---|---|",
  vapply(c("A91", "B87", "A66", "B53", "A94", "A37"), line_set, ""), "",
  "## Comparators (two numbers each)", "",
  vapply(seq_len(nrow(rbindlist(desc))), function(r) with(rbindlist(desc)[r], paste0(
    "- ", set, ": non-coordinated sign concordance ", fmt(noncoord_concordance), " (", n_noncoord,
    " pairs); edge-draw null mean ", fmt(edge_draw_null_mean), " (sd ", fmt(edge_draw_null_sd), ")")), ""), "",
  "## Controls computed in stage 2", "",
  vapply(seq_len(nrow(NC2)), function(r) with(NC2[r], paste0(
    "- NC2 ", set, ": informative-null pairs ", n_informative_null, ", S1_k1 ", fmt(s1_k1), ", label-null mean ",
    fmt(null_mean), ", p ", fmtp(p_perm_k1), if (nc2_passes_test) " (PASSES: null too lenient in this cohort)" else "")), ""), "",
  "## Composition arm and size-matched null (descriptive)", "",
  if (length(comp_reading)) paste0("- ", comp_reading) else "- No primary set passed; composition arm reported in stage2/composition_arm.tsv only.",
  if (length(size_reading)) paste0("- ", size_reading) else NULL, "",
  "## Limits that apply to every reading", "",
  "- Participant overlap with discovery cohorts was not checked by genotype in either cohort. GSE225740 identities cannot be checked (no reads, no deposited sex, age or batch).",
  "- Bulk program covariance mixes cell-type abundance, biopsy content and within-cell programs.",
  "- F0-F1 strata are cross-sectional; 29 of 66 A66 participants have diagnosis 'No Pathology'.",
  "",
  "Tables: stage2/*.tsv. Stop rules: stage1/stop_rule.tsv.")
writeLines(res, RESULT)
file.copy(RESULT, file.path(OUT, basename(RESULT)))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
fwrite(data.table(file = list.files(file.path(H54_RUN, "code"), full.names = TRUE))[
  , sha256 := vapply(file, sha256_of, "")], file.path(OUT, "code_sha256.tsv"), sep = "\t")
msg("STAGE2 COMPLETE in %.1f min: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), overall)
