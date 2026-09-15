#!/usr/bin/env Rscript
# 05: the first outcome read.
#
# Everything before this point ran without seeing a fibrosis stage, a NAS score
# or a diagnosis. This script reads them for the first time, and it may not
# change a single module. Its jobs are to give each frozen module a direction
# and an endpoint family, to label which modules the discovery cohorts support,
# and to place that support against the nulls that decide what it is worth.
#
# WHAT CHANGED IN v2 (see 00_contract.json, amendments_v2). Endpoints are
# z-scored within cohort so every effect is a standardized slope; the two axes
# are corrected jointly (320 tests) rather than OR-ed after separate BH; the
# expression-matched null bins mean log-CPM rather than the centred z it binned
# in v1, and matches each module's own profile; the competitive draws are made
# once here at 10,000 per size, saved for every later assay, and carry their
# within-set coherence so a coherence-matched sensitivity is possible; and the
# competitive machinery is run under 20 label permutations so its rejection
# rate under no association is on record.
#
# WHY TWO NULLS. Bulk liver expression is broadly correlated with fibrosis on
# this substrate, so a gene set of the right size is not expected to score
# zero. The label permutation asks whether the association exceeds chance for
# THIS gene set. The competitive draw asks whether a connected, co-expressed
# set of the same size drawn from the same graph would have done as well.
# The two answer different questions and are reported as different states.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))
suppressPackageStartupMessages(library(igraph))

contract <- cam_contract()
out <- cam_dir("discovery")
nulls_dir <- cam_dir("nulls")
mod_dir <- file.path(cam_out_root(), "modules")
cam_assert(file.exists(file.path(mod_dir, "MODULES_FROZEN.ok")),
           "Modules are not frozen; run 04 first")

# The freeze must predate this script's outputs. The validator re-checks it.
freeze_time <- file.info(file.path(mod_dir, "MODULES_FROZEN.ok"))$mtime
cam_say("modules frozen at ", format(freeze_time), "; reading outcomes now")

set.seed(contract$seed)
z <- readRDS(file.path(cam_out_root(), "graph", "cohort_z_expression.rds"))
g <- readRDS(file.path(cam_out_root(), "graph", "coexpression_graph.rds"))
rho_unclipped <- readRDS(file.path(cam_out_root(), "graph", "spearman_rho_unclipped.rds"))
mean_logcpm <- fread(file.path(cam_out_root(), "graph", "gene_mean_logcpm.tsv"))
membership <- fread(file.path(mod_dir, "module_membership.tsv"))
registry <- fread(file.path(mod_dir, "module_registry.tsv"))
manifest <- fread(cam_input("bulk_manifest", contract))
cam_assert(identical(colnames(z), manifest$sample_id),
           "Expression columns and manifest rows are not aligned")

modules <- split(membership$gene_symbol, membership$module_id)
module_ids <- names(modules)
cohorts <- sort(unique(manifest$dataset))
n_family <- length(module_ids)
cam_say(n_family, " modules, ", length(cohorts), " discovery cohorts")

# Equal-weight score, then standardised within cohort so that a cohort with a
# wider severity range cannot dominate the pooled estimate through its variance.
# No direction is applied here or in any assay script; direction enters once,
# in step 15.
score_sets <- function(sets) {
  S <- t(vapply(sets, function(M) {
    idx <- intersect(M, rownames(z))
    if (!length(idx)) return(rep(NA_real_, ncol(z)))
    colMeans(z[idx, , drop = FALSE])
  }, numeric(ncol(z))))
  for (ds in cohorts) {
    j <- which(manifest$dataset == ds)
    S[, j] <- t(apply(S[, j, drop = FALSE], 1L, standardize_vector))
  }
  S
}
S_obs <- score_sets(modules)
rownames(S_obs) <- module_ids
cam_write_rds(S_obs, file.path(out, "module_scores_discovery.rds"))

axes <- list(fibrosis = "fibrosis_stage", nas = "nas_score")

# Per cohort, the endpoint is z-scored over the participants entering the fit,
# so the coefficient is SD of score per SD of endpoint. Participants missing
# the grade are dropped from that axis only.
cohort_design <- function(axis_col) {
  lapply(cohorts, function(ds) {
    j <- which(manifest$dataset == ds)
    y <- manifest[[axis_col]][j]; sex <- manifest$inferred_sex[j]
    yz <- cam_endpoint_z(y, data.frame(sex = sex))
    keep <- is.finite(yz) & !is.na(sex)
    if (sum(keep) < 20L || length(unique(y[keep])) < 2L) return(NULL)
    X <- if (length(unique(sex[keep])) > 1L)
      stats::model.matrix(~ sex[keep])[, -1, drop = FALSE] else NULL
    list(ds = ds, j = j[keep], y = yz[keep], strata = sex[keep], X = X)
  })
}
designs <- lapply(axes, cohort_design)

fit_axis <- function(S, des, y_override = NULL) {
  beta <- se <- matrix(NA_real_, nrow(S), length(cohorts),
                       dimnames = list(rownames(S), cohorts))
  n_used <- setNames(integer(length(cohorts)), cohorts)
  for (ci in seq_along(des)) {
    d <- des[[ci]]; if (is.null(d)) next
    y <- if (is.null(y_override)) d$y else y_override[[ci]]
    fitres <- cam_fast_assoc(S[, d$j, drop = FALSE], y, d$X)
    beta[, d$ds] <- fitres$beta; se[, d$ds] <- fitres$se
    n_used[d$ds] <- length(d$j)
  }
  list(beta = beta, se = se, n_used = n_used)
}

obs_rows <- list(); meta_store <- list()
for (ax in names(axes)) {
  f <- fit_axis(S_obs, designs[[ax]])
  m <- cam_meta_rows(f$beta, f$se)
  meta_store[[ax]] <- list(fit = f, meta = m)
  concordant <- apply(sign(f$beta), 1L, function(v) {
    v <- v[is.finite(v)]; length(v) > 0L && length(unique(v)) == 1L
  })
  obs_rows[[ax]] <- data.table(
    module_id = module_ids, axis = ax,
    meta_beta = m$beta, meta_se = m$se, meta_z = m$z, meta_p = m$p_two_sided,
    n_cohorts = m$n_cohorts, direction_concordant = concordant,
    n_participants = sum(f$n_used)
  )
}
observed <- rbindlist(obs_rows)
# One family for the two axes: BH over 320 tests.
observed[, meta_q_joint := ml_complete_bh(meta_p, 2L * n_family)]
for (ax in names(axes)) cam_say(
  "axis ", ax, ": ", observed[axis == ax, sum(meta_q_joint < 0.05, na.rm = TRUE)],
  " of ", n_family, " at joint BH q<0.05, ",
  observed[axis == ax, sum(meta_q_joint < 0.05 & direction_concordant, na.rm = TRUE)],
  " of those direction-concordant")

# Cross-check the vectorised path against lm() on real fits.
check_rows <- list()
ds_big <- cohorts[[which.max(table(manifest$dataset))]]
dd <- designs$fibrosis[[which(cohorts == ds_big)]]
for (i in head(seq_along(module_ids), 5L)) {
  d <- data.table(score = S_obs[i, dd$j], y = dd$y, sex = dd$strata)
  ref <- summary(stats::lm(score ~ y + factor(sex), data = d))$coefficients["y", ]
  check_rows[[length(check_rows) + 1L]] <- data.table(
    module_id = module_ids[i], cohort = ds_big,
    fast_beta = meta_store$fibrosis$fit$beta[i, ds_big], lm_beta = ref[[1]],
    abs_difference = abs(meta_store$fibrosis$fit$beta[i, ds_big] - ref[[1]]))
}
vec_check <- rbindlist(check_rows)
cam_assert(max(vec_check$abs_difference) < 1e-8, "Vectorised association disagrees with lm()")
cam_write_tsv(vec_check, file.path(out, "vectorised_fit_check.tsv"))

# --- direction, endpoint family and the discovery label ----------------------
# A module's family is the axis that supports it (the smaller q if both do).
# Modules with no support are assigned fibrosis and form an unselected
# comparison class; they are not negative controls, because nothing selected
# them to be null.
wide <- dcast(observed, module_id ~ axis,
              value.var = c("meta_beta", "meta_q_joint", "direction_concordant"))
wide[, fibrosis_supported := is.finite(meta_q_joint_fibrosis) & meta_q_joint_fibrosis < 0.05 &
       direction_concordant_fibrosis]
wide[, nas_supported := is.finite(meta_q_joint_nas) & meta_q_joint_nas < 0.05 &
       direction_concordant_nas]
wide[, discovery_supported := fibrosis_supported | nas_supported]
wide[, endpoint_family := fifelse(fibrosis_supported & nas_supported,
        fifelse(meta_q_joint_fibrosis <= meta_q_joint_nas, "fibrosis", "activity"),
        fifelse(nas_supported, "activity", "fibrosis"))]
wide[, family_beta := fifelse(endpoint_family == "fibrosis", meta_beta_fibrosis, meta_beta_nas)]
wide[, direction := sign(family_beta)]
wide[!is.finite(direction) | direction == 0, direction := 1]
wide[, discovery_effect := abs(family_beta)]
wide[, discovery_label := fifelse(discovery_supported, "discovery_supported",
                                  "discovery_unsupported_comparison_class")]
cam_say("discovery-supported modules: ", sum(wide$discovery_supported), " of ", n_family,
        " (fibrosis family ", sum(wide$endpoint_family == "fibrosis"),
        ", activity family ", sum(wide$endpoint_family == "activity"), ")")

# --- null 1: label permutation ----------------------------------------------
n_perm <- contract$nulls$label_permutation$draws
cam_say("label permutation, ", n_perm, " draws")
perm_exceed <- matrix(0L, n_family, length(axes), dimnames = list(module_ids, names(axes)))
permute_design <- function(des, seed_offset) {
  lapply(seq_along(des), function(ci) {
    d <- des[[ci]]; if (is.null(d)) return(NULL)
    cam_permute_within(d$y, d$strata, seed_offset + ci)
  })
}
for (ax in names(axes)) {
  obs_z <- abs(meta_store[[ax]]$meta$z)
  for (b in seq_len(n_perm)) {
    yp <- permute_design(designs[[ax]], contract$seed + b * 1000L)
    fr <- fit_axis(S_obs, designs[[ax]], yp)
    mz <- abs(cam_meta_rows(fr$beta, fr$se)$z)
    perm_exceed[, ax] <- perm_exceed[, ax] + as.integer(is.finite(mz) & is.finite(obs_z) & mz >= obs_z)
  }
  cam_say("  ", ax, " permutation complete")
}
perm_p <- (1 + perm_exceed) / (1 + n_perm)

# --- null 2: competitive connected draws, drawn ONCE for every assay ---------
n_comp <- contract$nulls$competitive_connected$draws
sizes <- sort(unique(registry$n_genes))
size_of <- setNames(registry$n_genes, registry$module_id)[module_ids]
cam_say("competitive null: ", n_comp, " connected draws for each of ", length(sizes), " sizes")
sets_by_size <- list(); coh_by_size <- list(); comp_diag <- list()
module_coherence <- vapply(modules, function(M) cam_set_coherence(rho_unclipped, M), numeric(1))
for (s in sizes) {
  draws <- cam_connected_draws(g, s, n_comp, contract$seed + s)
  comp_diag[[length(comp_diag) + 1L]] <- data.table(
    module_size = s, n_requested = n_comp, n_drawn = length(draws$sets),
    n_failed = draws$n_failed, failure_rate = draws$failure_rate,
    n_eligible_nodes = draws$n_eligible_nodes)
  sets_by_size[[as.character(s)]] <- draws$sets
  coh_by_size[[as.character(s)]] <- vapply(draws$sets, function(M)
    cam_set_coherence(rho_unclipped, M), numeric(1))
}
comp_diagnostics <- rbindlist(comp_diag)
cam_write_tsv(comp_diagnostics, file.path(out, "competitive_draw_diagnostics.tsv"))
cam_write_rds(sets_by_size, file.path(nulls_dir, "competitive_sets_by_size.rds"))
cam_write_rds(coh_by_size, file.path(nulls_dir, "competitive_set_coherence_by_size.rds"))
coh_tab <- rbindlist(lapply(sizes, function(s) {
  mc <- module_coherence[size_of == s]; dc <- coh_by_size[[as.character(s)]]
  data.table(module_size = s, n_modules = length(mc),
             module_median_coherence = stats::median(mc),
             draw_median_coherence = stats::median(dc, na.rm = TRUE),
             draw_q05 = stats::quantile(dc, 0.05, na.rm = TRUE),
             draw_q95 = stats::quantile(dc, 0.95, na.rm = TRUE),
             fraction_modules_above_draw_q95 = mean(mc > stats::quantile(dc, 0.95, na.rm = TRUE)))
}))
cam_write_tsv(coh_tab, file.path(nulls_dir, "coherence_by_size.tsv"))
cam_write_tsv(data.table(module_id = module_ids, within_set_mean_spearman = unname(module_coherence)),
              file.path(nulls_dir, "module_coherence.tsv"))
cam_say("coherence: modules above the draws' 95th percentile in ",
        round(100 * mean(coh_tab$fraction_modules_above_draw_q95 * coh_tab$n_modules) /
                mean(coh_tab$n_modules) , 1), " percent of cases (size-weighted)")

# Closest draws by coherence, per module: the coherence-matched sensitivity.
n_matched <- 2000L
matched_idx <- lapply(module_ids, function(mid) {
  s <- as.character(size_of[[mid]]); dc <- coh_by_size[[s]]
  o <- order(abs(dc - module_coherence[[mid]]))[seq_len(min(n_matched, sum(is.finite(dc))))]
  list(idx = o, gap = stats::median(dc[o], na.rm = TRUE) - module_coherence[[mid]])
})
names(matched_idx) <- module_ids
cam_write_rds(matched_idx, file.path(nulls_dir, "competitive_coherence_matched_index.rds"))

# Observed competitive p, unmatched and coherence-matched, and the calibration
# run: the same machinery under 20 within-stratum label permutations.
n_cal <- contract$nulls$calibration$label_permutations
cal_designs <- lapply(seq_len(n_cal), function(b)
  permute_design(designs$fibrosis, contract$seed + 777000L + b * 1000L))
comp_p <- comp_p_matched <- matrix(NA_real_, n_family, length(axes),
                                   dimnames = list(module_ids, names(axes)))
comp_gap <- setNames(vapply(matched_idx, function(x) x$gap, numeric(1)), module_ids)
cal_hits <- matrix(NA_real_, n_family, n_cal, dimnames = list(module_ids, NULL))
for (s in sizes) {
  sets <- sets_by_size[[as.character(s)]]
  if (!length(sets)) next
  Sd <- score_sets(sets)
  idx <- which(size_of == s)
  nullz <- list()
  for (ax in names(axes)) {
    fd <- fit_axis(Sd, designs[[ax]])
    nullz[[ax]] <- abs(cam_meta_rows(fd$beta, fd$se)$z)
  }
  for (i in idx) for (ax in names(axes)) {
    obs <- abs(meta_store[[ax]]$meta$z[i])
    comp_p[i, ax] <- cam_empirical_p_ge(obs, nullz[[ax]])
    comp_p_matched[i, ax] <- cam_empirical_p_ge(obs, nullz[[ax]][matched_idx[[module_ids[i]]]$idx])
  }
  for (b in seq_len(n_cal)) {
    fo <- fit_axis(S_obs[idx, , drop = FALSE], designs$fibrosis, cal_designs[[b]])
    oz <- abs(cam_meta_rows(fo$beta, fo$se)$z)
    fd <- fit_axis(Sd, designs$fibrosis, cal_designs[[b]])
    nz <- abs(cam_meta_rows(fd$beta, fd$se)$z)
    cal_hits[idx, b] <- vapply(oz, function(o) cam_empirical_p_ge(o, nz), numeric(1)) < 0.05
  }
  cam_say("  size ", s, ": ", length(sets), " sets scored; calibration done")
}
cal_rate <- mean(cal_hits, na.rm = TRUE)
cal_by_perm <- colMeans(cal_hits, na.rm = TRUE)
cal_ok <- cal_rate >= contract$nulls$calibration$acceptable_rejection_rate_at_0.05[1] &&
  cal_rate <= contract$nulls$calibration$acceptable_rejection_rate_at_0.05[2]
cam_write_tsv(data.table(permutation = seq_len(n_cal), rejection_rate_at_0p05 = cal_by_perm),
              file.path(nulls_dir, "competitive_calibration_by_permutation.tsv"))
cam_write_json(list(
  label_permutations = n_cal, axis = "fibrosis",
  rejection_rate_at_0p05 = cal_rate,
  acceptable_range = contract$nulls$calibration$acceptable_rejection_rate_at_0.05,
  calibrated = cal_ok,
  consequence_if_not = "specificity states are reported as descriptive only"
), file.path(nulls_dir, "competitive_calibration.json"))
cam_say("competitive calibration: rejection rate ", round(cal_rate, 4),
        if (cal_ok) " (within range)" else " (OUTSIDE range; specificity descriptive only)")

# --- null 3: expression-decile-matched sets, per module ---------------------
n_dec <- contract$nulls$expression_decile_matched$draws
cam_say("expression-decile-matched null, ", n_dec, " draws per module")
mu <- setNames(mean_logcpm$mean_logcpm, mean_logcpm$gene_symbol)[rownames(z)]
bins <- t3_expression_deciles(matrix(mu, ncol = 1, dimnames = list(rownames(z), NULL)))
pool_by_bin <- split(rownames(z), bins[rownames(z)])
dec_p <- matrix(NA_real_, n_family, length(axes), dimnames = list(module_ids, names(axes)))
for (i in seq_len(n_family)) {
  M <- modules[[i]]
  set.seed(contract$seed + 7771L + i)
  counts <- tabulate(bins[M], nbins = length(pool_by_bin))
  sets <- lapply(seq_len(n_dec), function(b) {
    unlist(lapply(seq_along(counts), function(bb) {
      if (!counts[[bb]]) return(character(0))
      pool <- setdiff(pool_by_bin[[bb]], M)
      if (length(pool) < counts[[bb]]) return(character(0))
      sample(pool, counts[[bb]])
    }), use.names = FALSE)
  })
  Sd <- score_sets(sets)
  for (ax in names(axes)) {
    fd <- fit_axis(Sd, designs[[ax]])
    nullz <- abs(cam_meta_rows(fd$beta, fd$se)$z)
    dec_p[i, ax] <- cam_empirical_p_ge(abs(meta_store[[ax]]$meta$z[i]), nullz)
  }
}

# --- composition sensitivity ------------------------------------------------
# The 15 orthonormal log-ratio coordinates of the accepted 16-lineage
# composition are added to the fit. Attenuation is reported per axis as the
# ratio of adjusted to unadjusted slope; a count alone hides that the two axes
# behave differently.
cam_say("composition-adjusted sensitivity")
comp_tab <- fread(cam_input("bulk_composition", contract))
comp_tab <- comp_tab[match(manifest$sample_id, sample_id)]
cam_assert(!any(is.na(comp_tab$sample_id)), "Composition table does not cover the manifest")
prop_cols <- setdiff(names(comp_tab), c("sample_id", "dataset"))
ilr <- cam_ilr(as.matrix(comp_tab[, ..prop_cols]))
comp_adj_rows <- list()
for (ax in names(axes)) {
  beta <- se <- matrix(NA_real_, n_family, length(cohorts))
  for (ci in seq_along(cohorts)) {
    d <- designs[[ax]][[ci]]; if (is.null(d) || length(d$j) < 30L) next
    X <- cbind(d$X, ilr[d$j, , drop = FALSE])
    fr <- cam_fast_assoc(S_obs[, d$j, drop = FALSE], d$y, X)
    beta[, ci] <- fr$beta; se[, ci] <- fr$se
  }
  m <- cam_meta_rows(beta, se)
  comp_adj_rows[[ax]] <- data.table(
    module_id = module_ids, axis = ax,
    composition_adjusted_beta = m$beta,
    composition_adjusted_q = ml_complete_bh(m$p_two_sided, n_family),
    composition_attenuation_ratio = m$beta / meta_store[[ax]]$meta$beta)
}
composition_adjusted <- rbindlist(comp_adj_rows)
cam_write_tsv(composition_adjusted, file.path(out, "composition_adjusted.tsv"))

# --- assemble ---------------------------------------------------------------
null_long <- rbindlist(lapply(names(axes), function(ax) data.table(
  module_id = module_ids, axis = ax,
  label_permutation_p = perm_p[, ax],
  competitive_connected_p = comp_p[, ax],
  competitive_coherence_matched_p = comp_p_matched[, ax],
  competitive_coherence_gap = unname(comp_gap[module_ids]),
  expression_matched_p = dec_p[, ax]
)))
observed <- merge(observed, null_long, by = c("module_id", "axis"), sort = FALSE)
observed <- merge(observed, composition_adjusted, by = c("module_id", "axis"), sort = FALSE)
observed[, competitive_connected_q := ml_complete_bh(competitive_connected_p, n_family), by = axis]
observed[, competitive_coherence_matched_q := ml_complete_bh(competitive_coherence_matched_p, n_family), by = axis]
cam_assert_no_prohibited_columns(observed, contract)
cam_write_tsv(observed, file.path(out, "discovery_axis_results.tsv"))

labels <- wide[, .(module_id, direction, endpoint_family, discovery_label, discovery_supported,
                   discovery_effect, fibrosis_supported, nas_supported,
                   meta_beta_fibrosis, meta_q_joint_fibrosis, meta_beta_nas, meta_q_joint_nas)]
labels <- merge(labels, registry[, .(module_id, n_genes, min_loco_heldout_score_spearman,
                                     score_stability_flag)], by = "module_id")
cam_assert_no_prohibited_columns(labels, contract)
cam_write_tsv(labels, file.path(out, "module_discovery_labels.tsv"))

summ_axis <- function(ax) {
  o <- observed[axis == ax]; sig <- o$meta_q_joint < 0.05
  list(
    n_joint_bh_q05 = sum(sig, na.rm = TRUE),
    n_label_permutation_p05 = sum(o$label_permutation_p < 0.05, na.rm = TRUE),
    n_expression_matched_p05 = sum(o$expression_matched_p < 0.05, na.rm = TRUE),
    n_competitive_p05_uncorrected = sum(o$competitive_connected_p < 0.05, na.rm = TRUE),
    n_competitive_bh_q05 = sum(o$competitive_connected_q < 0.05, na.rm = TRUE),
    n_competitive_coherence_matched_bh_q05 = sum(o$competitive_coherence_matched_q < 0.05, na.rm = TRUE),
    n_composition_adjusted_q05_among_bh = sum(sig & o$composition_adjusted_q < 0.05, na.rm = TRUE),
    median_composition_attenuation_ratio_among_bh =
      stats::median(o$composition_attenuation_ratio[sig], na.rm = TRUE),
    n_participants = o$n_participants[1], n_cohorts = max(o$n_cohorts, na.rm = TRUE))
}
cam_write_json(list(
  n_modules = n_family,
  n_discovery_supported = sum(wide$discovery_supported),
  n_fibrosis_supported = sum(wide$fibrosis_supported),
  n_nas_supported = sum(wide$nas_supported),
  n_endpoint_family_fibrosis = sum(wide$endpoint_family == "fibrosis"),
  n_endpoint_family_activity = sum(wide$endpoint_family == "activity"),
  n_direction_negative = sum(wide$direction < 0),
  fibrosis = summ_axis("fibrosis"), nas = summ_axis("nas"),
  label_permutation_draws = n_perm, competitive_draws = n_comp,
  coherence_matched_draws = n_matched, expression_matched_draws = n_dec,
  competitive_calibration_rejection_rate = cal_rate, competitive_calibrated = cal_ok,
  max_competitive_failure_rate = max(comp_diagnostics$failure_rate),
  modules_frozen_at = format(freeze_time, tz = "UTC", "%Y-%m-%dT%H:%M:%SZ"),
  outcome_read_at = format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ")
), file.path(out, "discovery_summary.json"))

writeLines("discovery labelled", file.path(out, "READY"))
cam_say("05 complete")
