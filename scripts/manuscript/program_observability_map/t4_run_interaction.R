#!/usr/bin/env Rscript

# T4 (SYS-I): residual program interaction structure on the sealed discovery set.
#
# THE QUESTION
# Which of the 117 frozen programs move together in a donor AFTER their common
# histologic drive is removed? Pairs that covary because both track fibrosis are
# removed by construction. What is left is coordination among donors carrying
# the SAME recorded histology.
#
# WHAT THIS IS NOT
# Not a co-expression network. No module is discovered, no gene graph is built,
# no external outcome is read, and the 117 programs are used exactly as frozen.
# Residual coordination between fixed programs is a different object from a
# co-expression module graph and needs no new data.
#
# FOUR ARMS, EACH REMOVING ONE MORE UNINTERESTING REASON FOR AN EDGE
#   full     : residual correlation of the programs as frozen. Reported for
#              completeness and never on its own, because shared genes alone
#              produce it.
#   disjoint : PRIMARY. The genes the two programs share are deleted from BOTH
#              before correlating. A pair that cannot give up its shared genes
#              and keep 80 percent of each program's weight is refused as
#              membership-entangled rather than called.
#   partial  : the leading component of the program residual space is projected
#              out of both members as well, at k = 1 and again at k = 3. Scores
#              from one bulk library share a large amount of variance for reasons
#              specific to no pair, so this arm asks which edges are more than
#              those global axes. Only POSITIVE partial edges are interpretable;
#              see the note above the support assignment.
#
# Design rationale and the failure modes each piece exists to block are in
# t4_interaction_lib.R; the synthetic checks are in t4_run_tests.R and the
# submission script runs them first.

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
  library(Matrix)
  library(parallel)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(script_dir, "calibration_lib.R"))
source(file.path(script_dir, "systems_contract.R"))
source(file.path(script_dir, "t4_interaction_lib.R"))

set.seed(seed)
n_ref <- as.integer(Sys.getenv("T4_REFERENCE_PERMUTATIONS", "4000"))
n_cal <- as.integer(Sys.getenv("T4_CALIBRATION_PERMUTATIONS", "2000"))
n_emp <- as.integer(Sys.getenv("T4_EMPIRICAL_PERMUTATIONS", "150000"))
n_grp <- as.integer(Sys.getenv("T4_GROUP_PERMUTATIONS", "1000"))
n_cores <- as.integer(Sys.getenv("T4_CORES", "16"))
alpha <- 0.05
# Prespecified, before any edge was looked at. A pair must keep 80 percent of
# each program's weight after its shared genes are deleted, which is the same
# coverage floor score_programs() uses to call a program testable at all.
min_retained_weight <- program_weight_coverage
# Prespecified scale at which a non-call becomes informative rather than merely
# underpowered: an excess correlation of 0.20 beyond what shared membership
# produces.
interpretable_excess <- 0.20

# Defaults to the system's own sealed directory. The override exists only so a
# smoke run can be pointed at scratch; a real run never sets it.
output_root <- Sys.getenv("T4_OUTPUT_ROOT",
                          file.path(workstream_root, "systems", "T4_interaction"))
if (dir.exists(output_root)) fail("Refusing to overwrite: ", output_root)
if (!file.exists(nonholdout_dge)) fail("Sealed inputs absent; run partition_inputs.R first")

# ----------------------------------------------------------------- [1] inputs

message("[1/10] Loading sealed non-holdout input")
assert_frozen_programs(program_registry, program_membership)
dge <- readRDS(nonholdout_dge)
meta <- as.data.table(readRDS(nonholdout_meta))
assert_holdout_sealed(meta, holdout_cohort, gse193066_crosswalk)
library_size <- data.table(sample_id = colnames(dge),
                           lib_size = if (!is.null(dge$samples$lib.size))
                             dge$samples$lib.size else colSums(dge$counts))
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))
symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
rm(logcpm, dge); gc()

registry <- fread(program_registry)
membership <- fread(program_membership)
scored <- score_programs(symbol_matrix, meta, membership, registry,
                         program_weight_coverage)
features <- scored$coverage[testable == TRUE, feature_id]
assert_true(length(features) == expected_testable_programs,
            sprintf("Expected %d testable programs; found %d",
                    expected_testable_programs, length(features)))

discovery_meta <- meta[dataset %in% discovery_cohorts &
                         !is.na(fibrosis_stage) & !is.na(nas_score)]
assert_true(nrow(discovery_meta) == 469L, "Discovery census drift")

# Biological unit, stated non-vacuously. The sealed metadata carries no
# participant column, so any "one sample per participant" assertion built from it
# compares sample_id to itself and is worthless. The substantive evidence is
# external and is asserted instead: no discovery cohort has a donor-pairing
# table, whereas every dataset in this project that does have repeat sampling
# does have one, and the single bulk cohort with repeat biopsies (GSE193066) is
# the sealed holdout.
assert_biological_unit(discovery_meta, "sample_id",
                       cohorts = discovery_cohorts,
                       pairing_root = file.path(project_root, "data"))
holdout_participants <- fread(gse193066_crosswalk)[, unique(participant_token)]
assert_true(!any(discovery_meta$sample_id %in% holdout_participants),
            "HOLDOUT LEAK: a discovery sample matches a GSE193066 participant token")

# ------------------------------------------------------------------ [2] blocks

message("[2/10] Building per-cohort saturated-residual blocks")
weights <- t4_membership_weights(membership, features)
panel <- intersect(sort(unique(weights$gene_symbol)), rownames(symbol_matrix))
assert_true(length(panel) > 1000L, "Program gene panel collapsed")
panel_matrix <- symbol_matrix[panel, discovery_meta$sample_id, drop = FALSE]
rm(symbol_matrix); gc()

blocks <- lapply(discovery_cohorts, function(cc)
  t4_build_block(panel_matrix, discovery_meta[dataset == cc], weights, features,
                 program_weight_coverage))
names(blocks) <- discovery_cohorts
block_summary <- rbindlist(lapply(blocks, function(b) data.table(
  cohort = b$cohort, n_donors = b$n, n_model_columns = b$n_par, residual_df = b$df,
  n_panel_genes = length(b$genes),
  n_programs_scored = sum(b$coverage$cohort_testable))))
print(block_summary)
assert_true(all(block_summary$n_programs_scored == length(features)),
            "A program is not scorable in every discovery cohort; the edge family would differ by cohort")
assert_true(all(block_summary$residual_df > 40L), "A cohort has too little residual df")

pairs <- t4_pair_table(features, weights)
ops <- lapply(blocks, t4_disjoint_operators, pairs = pairs, features = features,
              min_retained_weight = min_retained_weight)
eligible_all <- Reduce(`&`, lapply(discovery_cohorts, function(cc) {
  e <- rep(TRUE, nrow(pairs))
  e[ops[[cc]]$pair_index] <- ops[[cc]]$eligible
  e
}))
pairs[, membership_entangled := !eligible_all]
message("      ", nrow(pairs), " pairs; ", pairs[n_shared > 0L, .N],
        " share a gene; ", sum(!eligible_all), " are membership-entangled")

# ------------------------------------------------------------- [3] observed

message("[3/10] Observed pass")
observed <- t4_one_pass(blocks, pairs, ops, want_matrix = TRUE)
print(round(observed$globals, 4))

message("[4/10] Stage-leakage gate")
leak <- t4_stage_leakage(blocks)
leak_summary <- leak[, .(max_abs_spearman = max(abs(spearman_stage)),
                         max_raw_stage_r2 = max(raw_stage_r2),
                         mean_stage_r2 = mean(stage_r2),
                         max_stage_r2 = max(stage_r2),
                         chance_stage_r2 = chance_stage_r2[[1L]]), by = cohort]
print(leak_summary)
# Saturation makes the unscaled residual exactly orthogonal to the stage factor.
# Anything above machine precision here means the model is not saturated, which
# is the exact artifact that forced a retraction elsewhere in this workstream.
assert_true(max(leak$raw_stage_r2) < 1e-16,
            sprintf("Stage structure survives the saturated residual model (max raw R2 = %.3g)",
                    max(leak$raw_stage_r2)))
# The leverage-scaled residual the correlations use cannot be exactly orthogonal,
# so it is held to the chance level for a 4-df regressor instead.
assert_true(all(leak_summary$mean_stage_r2 < leak_summary$chance_stage_r2),
            "Leverage-scaled residuals retain more stage structure than chance")

# What the global residual component is made of. It is projected out in the
# partial arm regardless, but an axis that turns out to be library depth or
# estimated composition should be named as such rather than left mysterious.
composition <- fread(nonholdout_composition)
comp_cols <- setdiff(names(composition), c("sample_id", "dataset"))
component_diag <- rbindlist(lapply(seq_along(blocks), function(ci) {
  b <- blocks[[ci]]
  U <- observed$global_components[[ci]]
  covariates <- c(list(log_library_size = log(library_size[match(b$samples, sample_id),
                                                           lib_size])),
                  lapply(setNames(comp_cols, comp_cols), function(cn)
                    suppressWarnings(as.numeric(composition[match(b$samples, sample_id)][[cn]]))))
  rbindlist(lapply(names(covariates), function(cn) {
    v <- covariates[[cn]]
    if (sum(is.finite(v)) < 20L || !isTRUE(stats::sd(v, na.rm = TRUE) > 0)) return(NULL)
    data.table(cohort = b$cohort, covariate = cn,
               component = seq_len(ncol(U)),
               spearman = suppressWarnings(as.numeric(
                 stats::cor(U, v, method = "spearman", use = "pairwise.complete.obs"))))
  }))
}))
print(component_diag[component == 1L,
                     .(max_abs_spearman = max(abs(spearman), na.rm = TRUE),
                       covariate = covariate[which.max(abs(spearman))]), by = cohort])

# ------------------------------------------------------------ [5] permutations

# Reproducible regardless of how many cores the scheduler gives us.
ref_seeds <- sample.int(.Machine$integer.max, n_ref)
cal_seeds <- sample.int(.Machine$integer.max, n_cal)
emp_seeds <- sample.int(.Machine$integer.max, n_emp)
grp_seeds <- sample.int(.Machine$integer.max, n_grp)
permute_genes <- function(b) t4_permute_genes(b$Z)
run_replicate <- function(s, permute) {
  set.seed(s)
  t4_one_pass(blocks, pairs, ops, Zlist = lapply(blocks, permute))
}
# Workers are given at most 250 replicates each and the phase runs in waves, so
# a phase reports progress rather than going silent for hours and a stalled
# worker is visible in the wave timings instead of being invisible.
#
# A worker whose result is large SPILLS it to a file and returns the path.
# Returning it through the fork pipe instead is what made an early version of
# this script unusable: the 4,000-replicate block, whose workers return only
# accumulator vectors, finished in 144 seconds, while the 2,000-replicate block,
# whose workers returned 200 MB of per-replicate statistics, had not finished
# forty minutes later on the same node.
t4_spill_dir <- file.path(tempdir(), "t4_spill")
dir.create(t4_spill_dir, showWarnings = FALSE, recursive = TRUE)
spilled <- function(value) {
  f <- tempfile(tmpdir = t4_spill_dir, fileext = ".rds")
  saveRDS(value, f, compress = FALSE)
  f
}
chunked <- function(seeds, fn, label = "") {
  per_worker <- max(1L, min(250L, ceiling(length(seeds) / n_cores)))
  nchunk <- ceiling(length(seeds) / per_worker)
  chunks <- split(seq_along(seeds), cut(seq_along(seeds), nchunk, labels = FALSE))
  waves <- split(seq_len(nchunk), ceiling(seq_len(nchunk) / n_cores))
  out <- vector("list", nchunk)
  started <- Sys.time()
  for (w in seq_along(waves)) {
    got <- mclapply(chunks[waves[[w]]], function(idx) fn(seeds[idx]),
                    mc.cores = n_cores, mc.preschedule = FALSE)
    bad <- vapply(got, inherits, logical(1), what = "try-error")
    if (any(bad)) fail("A permutation worker failed: ",
                       as.character(got[[which(bad)[[1L]]]]))
    out[waves[[w]]] <- got
    done <- sum(lengths(chunks[unlist(waves[seq_len(w)])]))
    elapsed <- as.numeric(Sys.time() - started, units = "secs")
    message(sprintf("      %s wave %d/%d: %d/%d replicates, %.0f s elapsed, %.1f rep/s",
                    label, w, length(waves), done, length(seeds), elapsed,
                    done / max(elapsed, 1e-6)))
  }
  out
}
ARMS <- c("z_full", "z_disjoint", "z_partial_k1", "z_partial_k3")

message("[5/10] Null A reference block: ", n_ref, " per-gene permutations on ",
        n_cores, " cores")
ref_parts <- chunked(ref_seeds, label = "reference", fn = function(ss) {
  acc <- lapply(ARMS, function(a) list(s = numeric(nrow(pairs)),
                                       q = numeric(nrow(pairs)),
                                       c = numeric(nrow(pairs))))
  names(acc) <- ARMS
  n <- 0L
  for (s in ss) {
    p <- run_replicate(s, permute_genes)
    for (a in ARMS) {
      v <- p[[a]]; ok <- is.finite(v)
      acc[[a]]$s <- acc[[a]]$s + ifelse(ok, v, 0)
      acc[[a]]$q <- acc[[a]]$q + ifelse(ok, v^2, 0)
      acc[[a]]$c <- acc[[a]]$c + ok
    }
    n <- n + 1L
  }
  c(acc, list(n = n))
})
assert_true(sum(vapply(ref_parts, `[[`, integer(1), "n")) == n_ref,
            "Reference permutation count drift")
null_ref <- lapply(ARMS, function(a) {
  s <- Reduce(`+`, lapply(ref_parts, function(p) p[[a]]$s))
  q <- Reduce(`+`, lapply(ref_parts, function(p) p[[a]]$q))
  cnt <- Reduce(`+`, lapply(ref_parts, function(p) p[[a]]$c))
  mu <- s / cnt
  sd <- sqrt(pmax(q / cnt - mu^2, 0) * cnt / pmax(cnt - 1, 1))
  list(mean = ifelse(cnt > 2, mu, NA_real_),
       sd = ifelse(cnt > 2 & sd > 0, sd, NA_real_))
})
names(null_ref) <- ARMS
rm(ref_parts); gc()

standardized <- lapply(ARMS, function(a)
  (observed[[a]] - null_ref[[a]]$mean) / null_ref[[a]]$sd)
names(standardized) <- ARMS
parametric_p <- lapply(standardized, function(z) 2 * stats::pnorm(-abs(z)))

# ------------------------------------------------------- [6] calibration block

message("[6/10] Null A calibration block: ", n_cal, " held-out permutations")
cal_parts <- chunked(cal_seeds, label = "calibration", fn = function(ss) {
  out <- lapply(ARMS, function(a) matrix(NA_real_, length(ss), nrow(pairs)))
  names(out) <- ARMS
  gl <- matrix(NA_real_, length(ss), length(observed$globals))
  for (k in seq_along(ss)) {
    p <- run_replicate(ss[[k]], permute_genes)
    for (a in ARMS) out[[a]][k, ] <- p[[a]]
    gl[k, ] <- p$globals
  }
  spilled(c(out, list(gl = gl)))
})
cal_parts <- lapply(cal_parts, readRDS)
cal_z <- lapply(ARMS, function(a) do.call(rbind, lapply(cal_parts, `[[`, a)))
names(cal_z) <- ARMS
cal_gl <- do.call(rbind, lapply(cal_parts, `[[`, "gl"))
colnames(cal_gl) <- names(observed$globals)
rm(cal_parts); gc()
assert_true(nrow(cal_z$z_disjoint) == n_cal, "Calibration permutation count drift")

families <- lapply(ARMS, function(a) {
  ok <- is.finite(parametric_p[[a]])
  if (a != "z_full") ok <- ok & !pairs$membership_entangled
  which(ok)
})
names(families) <- ARMS
view_label <- c(z_full = "edge_full_membership",
                z_disjoint = "edge_disjoint_membership",
                z_partial_k1 = "edge_partial_global_component_k1",
                z_partial_k3 = "edge_partial_global_component_k3")

message("[7/10] Calibration gate")
calibrations <- lapply(ARMS, function(a) {
  fam <- families[[a]]
  cp <- 2 * stats::pnorm(-abs(sweep(sweep(cal_z[[a]][, fam, drop = FALSE], 2L,
                                          null_ref[[a]]$mean[fam], "-"), 2L,
                                    null_ref[[a]]$sd[fam], "/")))
  calibrate_view(observed_p = parametric_p[[a]][fam],
                 permuted_p_fn = function(i) cp[i, ],
                 n_reps = n_cal, n_tests = length(fam), alpha = alpha,
                 label = view_label[[a]], observed_p_is_empirical = FALSE)
})
names(calibrations) <- ARMS
calibration <- rbindlist(lapply(calibrations, calibration_row))
print(calibration)
gates <- lapply(calibrations, gate_view, strict = FALSE)
primary_passed <- gates$z_disjoint$passed && gates$z_partial_k1$passed
if (!primary_passed) {
  message("      GATE FAILED: ",
          paste(c(gates$z_disjoint$reasons, gates$z_partial_k1$reasons), collapse = "; "))
}

# ------------------------------------------------------- [8] empirical block

# The parametric reference above has no resolution floor and its false-call rate
# is measured, which is what licenses the counts. This block additionally gives
# each edge an empirical p from a permutation family large enough to resolve BH
# on ~6,200 pairs, so the primary calls do not rest on the normal reference
# being exactly right in the far tail. Only exceedance counts are accumulated;
# the replicate statistics are never stored.
message("[8/10] Empirical block: ", n_emp, " permutations (resolution floor ",
        signif(1 / (n_emp + 1), 3), " against a BH threshold of ",
        signif(alpha / length(families$z_disjoint), 3), ")")
emp_arms <- c("z_disjoint", "z_partial_k1", "z_partial_k3")
# An unevaluable edge gets an unreachable target so it never accumulates a hit.
emp_target <- lapply(emp_arms, function(a) {
  d <- abs(observed[[a]] - null_ref[[a]]$mean)
  ifelse(is.finite(d), d, Inf)
})
names(emp_target) <- emp_arms
emp_parts <- chunked(emp_seeds, label = "empirical", fn = function(ss) {
  acc <- lapply(emp_arms, function(a) list(hit = numeric(nrow(pairs)),
                                           n = numeric(nrow(pairs))))
  names(acc) <- emp_arms
  for (s in ss) {
    p <- run_replicate(s, permute_genes)
    for (a in emp_arms) {
      d <- abs(p[[a]] - null_ref[[a]]$mean)
      ok <- is.finite(d)
      acc[[a]]$hit <- acc[[a]]$hit + (ok & d >= emp_target[[a]])
      acc[[a]]$n <- acc[[a]]$n + ok
    }
    rm(p)
  }
  acc
})
empirical_p <- lapply(emp_arms, function(a) {
  hit <- Reduce(`+`, lapply(emp_parts, function(p) p[[a]]$hit))
  cnt <- Reduce(`+`, lapply(emp_parts, function(p) p[[a]]$n))
  # An edge with no observed statistic was given an unreachable target, so it
  # would otherwise come back with a spuriously tiny p in the full edge table.
  ifelse(cnt > 0 & is.finite(observed[[a]]), (1 + hit) / (cnt + 1), NA_real_)
})
names(empirical_p) <- emp_arms
rm(emp_parts); gc()
empirical_resolvable <- (1 / (n_emp + 1)) < (alpha / length(families$z_disjoint))

# ------------------------------------------------------------- [9] robustness

message("[9/10] Null B: ", n_grp, " program-block permutations")
groups <- t4_gene_groups(blocks, weights, features)
permute_groups <- function(b) t4_permute_groups(b$Z, groups[[b$cohort]])
grp_parts <- chunked(grp_seeds, label = "group-null", fn = function(ss) {
  out <- lapply(emp_arms, function(a) matrix(NA_real_, length(ss), nrow(pairs)))
  names(out) <- emp_arms
  gl <- matrix(NA_real_, length(ss), length(observed$globals))
  for (k in seq_along(ss)) {
    p <- run_replicate(ss[[k]], permute_groups)
    for (a in emp_arms) out[[a]][k, ] <- p[[a]]
    gl[k, ] <- p$globals
  }
  spilled(c(out, list(gl = gl)))
})
grp_parts <- lapply(grp_parts, readRDS)
grp_z <- lapply(emp_arms, function(a) do.call(rbind, lapply(grp_parts, `[[`, a)))
names(grp_z) <- emp_arms
grp_gl <- do.call(rbind, lapply(grp_parts, `[[`, "gl"))
colnames(grp_gl) <- names(observed$globals)
rm(grp_parts); gc()
grp_q <- lapply(emp_arms, function(a) {
  mu <- colMeans(grp_z[[a]], na.rm = TRUE)
  sdv <- apply(grp_z[[a]], 2L, stats::sd, na.rm = TRUE)
  z <- (observed[[a]] - mu) / sdv
  q <- rep(NA_real_, nrow(pairs))
  fam <- families[[a]]
  q[fam] <- p.adjust(2 * stats::pnorm(-abs(z[fam])), method = "BH")
  q
})
names(grp_q) <- emp_arms
grp_fam <- families$z_disjoint
grp_mu <- colMeans(grp_z$z_disjoint, na.rm = TRUE)[grp_fam]
grp_sd <- apply(grp_z$z_disjoint[, grp_fam, drop = FALSE], 2L, stats::sd, na.rm = TRUE)
grp_false_calls <- vapply(seq_len(n_grp), function(b) {
  z <- (grp_z$z_disjoint[b, grp_fam] - grp_mu) / grp_sd
  sum(p.adjust(2 * stats::pnorm(-abs(z)), method = "BH") < alpha, na.rm = TRUE)
}, numeric(1))

# -------------------------------------------------------------- [10] assemble

message("[10/10] Assembling, gating and sealing")
bh <- function(p, fam) { q <- rep(NA_real_, length(p)); q[fam] <- p.adjust(p[fam], "BH"); q }
q_par <- lapply(ARMS, function(a) bh(parametric_p[[a]], families[[a]]))
names(q_par) <- ARMS
q_emp <- lapply(emp_arms, function(a) bh(empirical_p[[a]], families[[a]]))
names(q_emp) <- emp_arms

zc <- observed$zc_disjoint
w_c <- vapply(blocks, function(b) b$df - 3L, numeric(1))
wmat <- matrix(w_c, nrow(pairs), length(blocks), byrow = TRUE) * is.finite(zc)
zbar <- rowSums(replace(zc, !is.finite(zc), 0) * wmat) / rowSums(wmat)
Qhet <- rowSums(wmat * (replace(zc, !is.finite(zc), 0) - zbar)^2)
dfhet <- pmax(rowSums(is.finite(zc)) - 1L, 1L)

reg <- registry[, .(feature_id = program_uid, cell_type, module_name)]
edges <- copy(pairs)
edges <- merge(edges, setnames(copy(reg), paste0(names(reg), "_a")),
               by.x = "program_a", by.y = "feature_id_a", sort = FALSE)
edges <- merge(edges, setnames(copy(reg), paste0(names(reg), "_b")),
               by.x = "program_b", by.y = "feature_id_b", sort = FALSE)
setorder(edges, pair_index)
edges[, `:=`(
  same_cell_type = cell_type_a == cell_type_b,
  r_full = tanh(observed$z_full),
  r_full_null_mean = tanh(null_ref$z_full$mean),
  q_full = q_par$z_full,
  r_disjoint = tanh(observed$z_disjoint),
  r_disjoint_null_mean = tanh(null_ref$z_disjoint$mean),
  excess_correlation = tanh(observed$z_disjoint) - tanh(null_ref$z_disjoint$mean),
  z_disjoint_standardized = standardized$z_disjoint,
  p_disjoint = parametric_p$z_disjoint, q_disjoint = q_par$z_disjoint,
  p_disjoint_empirical = empirical_p$z_disjoint,
  q_disjoint_empirical = q_emp$z_disjoint,
  q_disjoint_group_null = grp_q$z_disjoint,
  r_partial = tanh(observed$z_partial_k1),
  r_partial_null_mean = tanh(null_ref$z_partial_k1$mean),
  excess_correlation_partial = tanh(observed$z_partial_k1) - tanh(null_ref$z_partial_k1$mean),
  z_partial_standardized = standardized$z_partial_k1,
  q_partial = q_par$z_partial_k1,
  q_partial_empirical = q_emp$z_partial_k1,
  q_partial_group_null = grp_q$z_partial_k1,
  r_partial_k3 = tanh(observed$z_partial_k3),
  excess_correlation_partial_k3 = tanh(observed$z_partial_k3) -
    tanh(null_ref$z_partial_k3$mean),
  q_partial_k3 = q_par$z_partial_k3,
  q_partial_k3_empirical = q_emp$z_partial_k3,
  n_cohorts = observed$n_cohorts_disjoint,
  n_same_sign = observed$n_same_sign_disjoint,
  n_same_sign_partial = observed$n_same_sign_k1,
  heterogeneity_Q = Qhet, heterogeneity_df = dfhet,
  heterogeneity_I2 = pmax(0, (Qhet - dfhet) / pmax(Qhet, 1e-12)),
  mde_excess_correlation = t4_minimum_detectable_excess(
    null_ref$z_disjoint$mean, null_ref$z_disjoint$sd, alpha = alpha))]
for (k in seq_along(discovery_cohorts)) {
  edges[[paste0("r_disjoint_", discovery_cohorts[[k]])]] <- tanh(zc[, k])
}
# Only POSITIVE partial edges may be called. Projecting k dimensions out of p
# programs forces the leftover covariances to sum down, so a pair that was
# nothing but the global axis reappears with a negative partial correlation.
# That is an artifact of the projection, not evidence of anticoordination, and
# t4_run_tests.R demonstrates it on planted data.
partial_positive <- edges$excess_correlation_partial > 0
edges[, support := fifelse(
  membership_entangled, "membership_entangled",
  fifelse(!is.finite(q_disjoint), "unevaluable",
          fifelse(q_disjoint < alpha & q_partial < alpha & partial_positive &
                    n_same_sign_partial == n_cohorts,
                  "coordinated_beyond_global_component",
                  fifelse(q_disjoint < alpha & q_partial < alpha & partial_positive,
                          "coordinated_beyond_global_component_heterogeneous",
                          fifelse(q_disjoint < alpha, "coordinated_global_component_only",
                                  fifelse(is.finite(mde_excess_correlation) &
                                            mde_excess_correlation <= interpretable_excess,
                                          "informative_null", "indeterminate"))))))]
assert_negatives_have_mde(edges[support %in% c("informative_null", "indeterminate")],
                          "support", "mde_excess_correlation")
coordinated_classes <- c("coordinated_beyond_global_component",
                         "coordinated_beyond_global_component_heterogeneous")

# Global structure statistics, each against the held-out Null A block.
global_p <- function(obs, null) (1 + sum(null[is.finite(null)] >= obs)) /
  (sum(is.finite(null)) + 1)
globals <- data.table(
  statistic = names(observed$globals),
  observed = as.numeric(observed$globals),
  null_mean = colMeans(cal_gl, na.rm = TRUE),
  null_p95 = apply(cal_gl, 2L, stats::quantile, 0.95, na.rm = TRUE),
  null_max = apply(cal_gl, 2L, max, na.rm = TRUE),
  empirical_p = vapply(seq_along(observed$globals), function(k)
    global_p(observed$globals[[k]], cal_gl[, k]), numeric(1)),
  group_null_mean = colMeans(grp_gl, na.rm = TRUE),
  group_null_p = vapply(seq_along(observed$globals), function(k)
    global_p(observed$globals[[k]], grp_gl[, k]), numeric(1)),
  n_permutations = n_cal,
  empirical_p_floor = 1 / (n_cal + 1))
print(globals)

# Community structure on the surviving edges, with the same statistic recomputed
# on null replicates so that "it forms blocks" is a comparison, not an assertion.
called <- edges[support %in% coordinated_classes]
edge_weight <- ifelse(edges$support %in% coordinated_classes &
                        edges$excess_correlation_partial > 0,
                      edges$excess_correlation_partial, NA_real_)
comm <- t4_modularity(pairs, edge_weight, features)
# The null graph is RANK-MATCHED, not threshold-matched. Under the permutation
# null almost nothing survives BH, so a threshold-matched null graph is empty and
# its modularity is undefined rather than low, which would compare an observed
# number against nothing. Taking the same NUMBER of strongest positive null edges
# holds graph size fixed and isolates the question actually being asked: given a
# graph this dense, is the observed one organised into blocks?
n_called <- nrow(called)
null_modularity <- vapply(seq_len(min(n_cal, 200L)), function(b) {
  ex <- tanh(cal_z$z_partial_k1[b, ]) - tanh(null_ref$z_partial_k1$mean)
  ex[pairs$membership_entangled | !is.finite(ex) | ex <= 0] <- NA_real_
  if (!n_called || sum(!is.na(ex)) < n_called) return(NA_real_)
  cutoff <- sort(ex, decreasing = TRUE)[n_called]
  t4_modularity(pairs, ifelse(!is.na(ex) & ex >= cutoff, ex, NA_real_),
                features)$modularity
}, numeric(1))
# And a second null that keeps the observed graph and only rewires it, so the
# degree sequence cannot be what produces the blocks.
rewired_modularity <- t4_rewired_modularity(comm$graph, 200L)

summary_rows <- data.table(
  metric = c("n_donors", "n_cohorts", "n_programs", "n_pairs",
             "n_pairs_sharing_genes", "n_pairs_membership_entangled",
             "n_pairs_in_family", "min_retained_weight", "interpretable_excess",
             "n_reference_permutations", "n_calibration_permutations",
             "n_empirical_permutations", "n_group_permutations",
             "gate_passed_disjoint", "gate_passed_partial",
             "null_false_call_mean_disjoint", "null_false_call_max_disjoint",
             "null_false_call_mean_partial", "group_null_false_call_mean",
             "empirical_per_edge_resolvable", "empirical_p_floor",
             "bh_threshold_needed",
             "n_edges_disjoint_significant", "n_edges_partial_k1_positive",
             "n_edges_partial_k3_positive",
             "n_coordinated_beyond_global_component",
             "n_coordinated_all_cohorts", "n_coordinated_same_cell_type",
             "n_coordinated_cross_cell_type",
             "n_coordinated_empirical_confirmed",
             "n_coordinated_group_null_confirmed",
             "n_informative_null", "n_indeterminate",
             "median_mde_excess_correlation",
             "median_abs_r_full", "median_abs_r_disjoint", "median_abs_r_partial",
             "global_component_share", "max_raw_stage_r2",
             "modularity_observed", "modularity_rank_matched_null_mean",
             "modularity_rank_matched_null_p95", "modularity_rewired_mean",
             "modularity_rewired_p95", "n_communities"),
  value = c(nrow(discovery_meta), length(discovery_cohorts), length(features),
            nrow(pairs), pairs[n_shared > 0L, .N], sum(pairs$membership_entangled),
            length(families$z_disjoint), min_retained_weight, interpretable_excess,
            n_ref, n_cal, n_emp, n_grp,
            as.numeric(gates$z_disjoint$passed), as.numeric(gates$z_partial_k1$passed),
            calibrations$z_disjoint$null_call_mean,
            calibrations$z_disjoint$null_call_max,
            calibrations$z_partial_k1$null_call_mean, mean(grp_false_calls),
            as.numeric(empirical_resolvable), 1 / (n_emp + 1),
            alpha / length(families$z_disjoint),
            edges[q_disjoint < alpha, .N],
            edges[q_partial < alpha & excess_correlation_partial > 0, .N],
            edges[q_partial_k3 < alpha & excess_correlation_partial_k3 > 0, .N],
            nrow(called),
            called[support == "coordinated_beyond_global_component", .N],
            called[same_cell_type == TRUE, .N], called[same_cell_type == FALSE, .N],
            called[q_partial_empirical < alpha, .N],
            called[q_partial_group_null < alpha, .N],
            edges[support == "informative_null", .N],
            edges[support == "indeterminate", .N],
            stats::median(edges$mde_excess_correlation, na.rm = TRUE),
            stats::median(abs(edges$r_full), na.rm = TRUE),
            stats::median(abs(edges$r_disjoint), na.rm = TRUE),
            stats::median(abs(edges$r_partial), na.rm = TRUE),
            observed$globals[["global_component_share"]], max(leak$raw_stage_r2),
            comm$modularity, mean(null_modularity, na.rm = TRUE),
            as.numeric(stats::quantile(null_modularity, 0.95, na.rm = TRUE)),
            mean(rewired_modularity, na.rm = TRUE),
            as.numeric(stats::quantile(rewired_modularity, 0.95, na.rm = TRUE)),
            comm$n_communities))

# No count leaves this script without its calibration row, and the primary views
# have to have passed. If they did not, refusing is the intended behaviour: the
# counts are blanked and the diagnostics are still sealed so the failure is on
# the record.
if (!primary_passed) {
  withheld <- grep("^n_(coordinated|edges|informative|indeterminate)|^modularity_observed",
                   summary_rows$metric, value = TRUE)
  summary_rows[metric %in% withheld, value := NA_real_]
  edges[, `:=`(q_disjoint = NA_real_, q_partial = NA_real_,
               support = "withheld_gate_failed")]
  called <- called[0L]
  message("      COUNTS WITHHELD")
} else {
  assert_counts_calibrated(
    summary_rows, calibration[view %in% view_label[c("z_disjoint", "z_partial_k1")]])
}
# The k = 3 sensitivity carries its own gate. If it fails, its count is withheld
# on its own terms; it does not withhold the primary, and the primary is not
# allowed to borrow its number.
if (!gates$z_partial_k3$passed) {
  summary_rows[metric == "n_edges_partial_k3_positive", value := NA_real_]
  edges[, `:=`(q_partial_k3 = NA_real_, q_partial_k3_empirical = NA_real_)]
  message("      k = 3 sensitivity WITHHELD: ",
          paste(gates$z_partial_k3$reasons, collapse = "; "))
}
assert_inclusion_criterion(c("evidence_interpretation", "experiment_routing"))
assert_language(c(readLines(file.path(script_dir, "t4_interaction_lib.R")),
                  readLines(script_arg)))

tmp <- atomic_dir(output_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)
fwrite(edges, file.path(tmp, "program_interaction_edges.tsv.gz"), sep = "\t",
       quote = FALSE, na = "NA", compress = "gzip")
write_tsv(called[order(q_partial)][, .(
  program_a, program_b, cell_type_a, cell_type_b, module_name_a, module_name_b,
  same_cell_type, n_shared, jaccard, weight_cosine,
  r_full, r_disjoint, r_partial, excess_correlation, excess_correlation_partial,
  z_partial_standardized, q_disjoint, q_partial, q_partial_empirical,
  q_partial_group_null, q_partial_k3, r_partial_k3,
  n_same_sign_partial, n_cohorts, heterogeneity_I2, support)],
  file.path(tmp, "coordinated_edges.tsv"))
write_tsv(globals, file.path(tmp, "global_structure.tsv"))
write_tsv(calibration, file.path(tmp, "calibration.tsv"))
write_tsv(summary_rows, file.path(tmp, "interaction_summary.tsv"))
write_tsv(block_summary, file.path(tmp, "cohort_blocks.tsv"))
write_tsv(leak, file.path(tmp, "stage_leakage_diagnostic.tsv"))
write_tsv(component_diag, file.path(tmp, "global_component_covariates.tsv"))
write_tsv(data.table(index = seq_along(observed$eigenvalues),
                     eigenvalue = observed$eigenvalues,
                     share = observed$eigenvalues / sum(observed$eigenvalues)),
          file.path(tmp, "residual_spectrum.tsv"))
if (!is.null(comm$membership)) write_tsv(comm$membership,
                                         file.path(tmp, "program_communities.tsv"))
write_tsv(rbind(
  data.table(null = "rank_matched_permutation", replicate = seq_along(null_modularity),
             modularity = null_modularity),
  data.table(null = "degree_preserving_rewiring", replicate = seq_along(rewired_modularity),
             modularity = rewired_modularity)),
  file.path(tmp, "modularity_null.tsv"))
saveRDS(list(null_reference = null_ref, calibration_globals = cal_gl,
             group_globals = grp_gl, correlation_matrix = observed$correlation_matrix,
             global_components = observed$global_components),
        file.path(tmp, "null_reference.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))

inputs <- c(nonholdout_dge, nonholdout_meta, nonholdout_composition, gene_annotation,
            program_registry, program_membership)
write_tsv(data.table(path = inputs, sha256 = vapply(inputs, sha256_file, character(1))),
          file.path(tmp, "input_manifest.tsv"))
code <- file.path(script_dir, c("t4_interaction_lib.R", "t4_run_interaction.R",
                                "t4_run_tests.R", "config.R", "analysis_lib.R",
                                "calibration_lib.R", "systems_contract.R"))
write_tsv(data.table(script = basename(code),
                     sha256 = vapply(code, sha256_file, character(1))),
          file.path(tmp, "code_manifest.tsv"))
artifacts <- setdiff(list.files(tmp), "artifact_manifest.tsv")
write_tsv(data.table(artifact = artifacts,
                     sha256 = vapply(file.path(tmp, artifacts), sha256_file, character(1))),
          file.path(tmp, "artifact_manifest.tsv"))
jsonlite::write_json(list(
  state = if (primary_passed) "T4_INTERACTION_COMPLETE_HOLDOUT_UNOPENED"
          else "T4_INTERACTION_WITHHELD_GATE_FAILED",
  release_id = release_id, workstream_id = workstream_id, system = "T4_interaction",
  n_donors = nrow(discovery_meta), n_programs = length(features),
  n_pairs = nrow(pairs), n_pairs_in_family = length(families$z_disjoint),
  n_permutations_reference = n_ref, n_permutations_calibration = n_cal,
  n_permutations_empirical = n_emp, n_permutations_group_null = n_grp,
  primary_view = "edge_partial_global_component_k1",
  primary_gate_passed = primary_passed,
  empirical_per_edge_resolvable = empirical_resolvable,
  null_false_call_rate_mean = calibrations$z_disjoint$null_call_mean,
  n_coordinated = nrow(called),
  holdout_accessed = FALSE, external_outcomes_read = FALSE
), file.path(tmp, "T4_INTERACTION_READY.json"), pretty = TRUE, auto_unbox = TRUE)

publish_dir(tmp, output_root)
Sys.chmod(list.files(output_root, full.names = TRUE), mode = "0440")

cat("\n=== T4 INTERACTION ===\n")
print(summary_rows)
cat("T4_COMPLETE\t", output_root, "\n", sep = "")
