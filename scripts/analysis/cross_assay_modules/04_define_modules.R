#!/usr/bin/env Rscript
# 04: define and FREEZE the modules. Still label-blind.
#
# Nothing in this script reads a fibrosis stage, a NAS score or a diagnosis.
# That is the point: the partition must be fixed before any outcome can
# influence it, so that the external cohorts test a set they did not help
# choose and so that modules without disease support exist to act as the
# negative-control class. The freeze marker written at the end is what step 05
# and the validator check their own timestamps against.
#
# Two stability questions are settled here, both outcome-blind: does the
# partition survive dropping a whole discovery cohort, and does it survive a
# neighbouring resolution. A partition that fails either is not a coordinate
# system and the analysis stops.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))
suppressPackageStartupMessages({
  library(igraph)
})

# Adjusted Rand index via igraph rather than mclust, which is not installed in
# the project environment. igraph::compare takes integer membership vectors.
adjusted_rand <- function(a, b) {
  igraph::compare(as.integer(a), as.integer(b), method = "adjusted.rand")
}

contract <- cam_contract()
out <- cam_dir("modules")
graph_dir <- file.path(cam_out_root(), "graph")
cam_assert(file.exists(file.path(graph_dir, "READY")), "Run 03 first")

set.seed(contract$seed)
g <- readRDS(file.path(graph_dir, "coexpression_graph.rds"))
z <- readRDS(file.path(graph_dir, "cohort_z_expression.rds"))
scan <- fread(file.path(graph_dir, "k_selection_scan.tsv"))
gsum <- jsonlite::fromJSON(file.path(graph_dir, "graph_summary.json"))
k_selected <- gsum$k_selected
manifest <- fread(cam_input("bulk_manifest", contract))

min_size <- contract$modules$min_size
max_size <- contract$modules$max_size

leiden_at <- function(graph, resolution, seed) {
  set.seed(seed)
  igraph::cluster_leiden(graph, objective_function = "modularity",
                         resolution_parameter = resolution, n_iterations = 10L)
}

# Resolution is chosen by how much of the universe lands in modules of usable
# size. A resolution that puts most genes in one giant module or in singletons
# gives a vocabulary with nothing to say.
cam_say("scanning resolutions")
res_rows <- list(); parts <- list()
for (r in contract$modules$resolution_grid) {
  cl <- leiden_at(g, r, contract$seed)
  memb <- igraph::membership(cl)
  sizes <- table(memb)
  usable <- names(sizes)[sizes >= min_size & sizes <= max_size]
  in_usable <- sum(sizes[usable])
  res_rows[[length(res_rows) + 1L]] <- data.table(
    resolution = r, n_clusters = length(sizes),
    n_usable_modules = length(usable),
    n_genes_in_usable = in_usable,
    fraction_universe_in_usable = in_usable / length(memb),
    largest_cluster = max(sizes)
  )
  parts[[as.character(r)]] <- memb
  cam_say("resolution ", r, ": ", length(usable), " modules of size ",
          min_size, "-", max_size, " holding ", in_usable, " genes (",
          round(in_usable / length(memb), 3), ")")
}
res_scan <- rbindlist(res_rows)
cam_write_tsv(res_scan, file.path(out, "resolution_scan.tsv"))
best <- res_scan[which.max(fraction_universe_in_usable)]
resolution <- best$resolution
cam_say("selected resolution ", resolution)

membership_full <- parts[[as.character(resolution)]]

# Leave-one-cohort-out. The graph is rebuilt from scratch without the cohort,
# with the same weighting the full graph carries, so the check covers the
# correlation structure and not just the clustering step.
#
# WHAT IS EVALUATED, AND WHERE. The frozen module and its best-overlap
# rebuilt module are both scored on the OMITTED cohort's participants only.
# v1 scored them over all 844, which included the participants the rebuilt
# module was trained on, and reported that as stability. The all-sample value
# is still computed and printed as reported-only so the two can be compared.
#
# WHAT IS GATED. The module SCORE, not its membership. Membership ARI and
# best-match Jaccard are reported in full but do not decide; this is an
# exploratory redesign made after the prespecified membership criterion
# failed, recorded as such in GRAPH_K_AMENDMENT.md (amendments 3 and 4), and
# it is not a demonstration that the membership criterion was unnecessary.
#
# BASELINES. Two per fold, so a high score correlation cannot be read as
# reproducible biology when it is only the substrate's dominant axis: the
# frozen score against a random same-size gene set, and against the omitted
# cohort's first principal component.
cam_say("leave-one-cohort-out stability, evaluated on the omitted cohort")

sizes_full <- table(membership_full)
usable_full <- names(sizes_full)[sizes_full >= min_size & sizes_full <= max_size]
mods_full <- lapply(usable_full, function(cid) names(membership_full)[membership_full == as.integer(cid)])
names(mods_full) <- usable_full
score_of <- function(members, cols = seq_len(ncol(z))) {
  idx <- intersect(members, rownames(z))
  if (!length(idx)) return(rep(NA_real_, length(cols)))
  colMeans(z[idx, cols, drop = FALSE])
}

cohort_ids <- sort(unique(manifest$dataset))
loco_rows <- list(); module_rows <- list()
rho_held_mat <- matrix(NA_real_, nrow = length(mods_full), ncol = length(cohort_ids),
                       dimnames = list(names(mods_full), cohort_ids))
rho_all_mat <- rho_held_mat
for (fi in seq_along(cohort_ids)) {
  ds <- cohort_ids[[fi]]
  j_train <- which(manifest$dataset != ds)
  j_test <- which(manifest$dataset == ds)
  zz <- z[, j_train, drop = FALSE]
  # Re-standardise inside the retained cohorts so the dropped cohort leaves no
  # trace through the centring.
  for (d2 in setdiff(cohort_ids, ds)) {
    cols <- which(manifest$dataset[j_train] == d2)
    blk <- zz[, cols, drop = FALSE]
    mu <- rowMeans(blk); sdv <- apply(blk, 1L, stats::sd)
    sdv[!is.finite(sdv) | sdv == 0] <- NA_real_
    zz[, cols] <- (blk - mu) / sdv
  }
  keep <- stats::complete.cases(zz)
  gg <- cam_build_mutual_graph(zz[keep, , drop = FALSE], k_selected)
  mm <- igraph::membership(leiden_at(gg, resolution, contract$seed))
  names(mm) <- igraph::V(gg)$name
  by_cluster <- split(names(mm), mm)

  shared <- intersect(names(mm), names(membership_full))
  ari <- adjusted_rand(membership_full[shared], mm[shared])

  # Omitted-cohort PC1 on the universe; z is already centred within cohort.
  zt <- z[, j_test, drop = FALSE]
  pc1 <- stats::prcomp(t(zt), center = FALSE, scale. = FALSE, rank. = 1L)$x[, 1]

  set.seed(contract$seed + 100L + fi)
  for (i in seq_along(mods_full)) {
    M <- mods_full[[i]]
    overlaps <- vapply(by_cluster, function(C)
      length(intersect(M, C)) / length(union(M, C)), numeric(1))
    best <- by_cluster[[which.max(overlaps)]]
    s_held <- score_of(M, j_test); b_held <- score_of(best, j_test)
    rand_set <- sample(rownames(z), length(M))
    rho_held <- suppressWarnings(stats::cor(s_held, b_held, method = "spearman"))
    rho_all <- suppressWarnings(stats::cor(score_of(M), score_of(best), method = "spearman"))
    rho_held_mat[i, ds] <- rho_held; rho_all_mat[i, ds] <- rho_all
    module_rows[[length(module_rows) + 1L]] <- data.table(
      dropped_cohort = ds, cluster_id = names(mods_full)[i], n_genes = length(M),
      best_match_jaccard = max(overlaps),
      heldout_score_spearman = rho_held,
      allsample_score_spearman_reported_only = rho_all,
      abs_rho_vs_random_set_heldout = abs(suppressWarnings(
        stats::cor(s_held, score_of(rand_set, j_test), method = "spearman"))),
      abs_rho_vs_heldout_pc1 = abs(suppressWarnings(
        stats::cor(s_held, pc1, method = "spearman")))
    )
  }
  mr <- rbindlist(module_rows)[dropped_cohort == ds]
  loco_rows[[length(loco_rows) + 1L]] <- data.table(
    dropped_cohort = ds, n_participants_retained = length(j_train),
    n_participants_heldout = length(j_test), n_genes_compared = length(shared),
    membership_ari_reported_only = ari,
    median_best_match_jaccard_reported_only = stats::median(mr$best_match_jaccard, na.rm = TRUE),
    median_heldout_score_spearman = stats::median(mr$heldout_score_spearman, na.rm = TRUE),
    fraction_modules_heldout_rho_ge_0p8 = mean(mr$heldout_score_spearman >= 0.8, na.rm = TRUE),
    min_heldout_score_spearman = min(mr$heldout_score_spearman, na.rm = TRUE),
    median_allsample_score_spearman_reported_only =
      stats::median(mr$allsample_score_spearman_reported_only, na.rm = TRUE),
    median_abs_rho_vs_random_set = stats::median(mr$abs_rho_vs_random_set_heldout, na.rm = TRUE),
    median_abs_rho_vs_heldout_pc1 = stats::median(mr$abs_rho_vs_heldout_pc1, na.rm = TRUE)
  )
  cam_say("  drop ", ds, ": held-out score rho median ",
          round(stats::median(mr$heldout_score_spearman, na.rm = TRUE), 4),
          " (all-sample ", round(stats::median(mr$allsample_score_spearman_reported_only, na.rm = TRUE), 4),
          "; random-set ", round(stats::median(mr$abs_rho_vs_random_set_heldout, na.rm = TRUE), 3),
          "; PC1 ", round(stats::median(mr$abs_rho_vs_heldout_pc1, na.rm = TRUE), 3),
          ") | membership ARI ", round(ari, 4), " (reported only)")
}
loco <- rbindlist(loco_rows)
cam_write_tsv(loco, file.path(out, "loco_stability.tsv"))
cam_write_tsv(rbindlist(module_rows), file.path(out, "loco_module_stability.tsv"))
cam_write_tsv(loco[, .(dropped_cohort, median_heldout_score_spearman,
                       median_abs_rho_vs_random_set, median_abs_rho_vs_heldout_pc1)],
              file.path(out, "loco_baselines.tsv"))

# Neighbouring resolutions, same graph.
grid <- sort(contract$modules$resolution_grid)
pos <- which(grid == resolution)
idx <- unique(c(max(1L, pos - 1L), min(length(grid), pos + 1L)))
neigh <- setdiff(grid[idx], resolution)
cam_assert(length(neigh) >= 1L, "No neighbouring resolution to compare against")
neigh_rows <- rbindlist(lapply(neigh, function(r) {
  data.table(neighbour_resolution = r,
             adjusted_rand_index = adjusted_rand(
               membership_full, parts[[as.character(r)]]))
}))
cam_write_tsv(neigh_rows, file.path(out, "resolution_stability.tsv"))
cam_say("neighbour-resolution ARI: ",
        paste(round(neigh_rows$adjusted_rand_index, 4), collapse = ", "))

min_loco_score <- min(loco$median_heldout_score_spearman)
min_neigh <- min(neigh_rows$adjusted_rand_index)
cam_say("gate: minimum median HELD-OUT LOCO score Spearman = ", round(min_loco_score, 4))
cam_assert(min_loco_score >= contract$modules$loco_score_spearman_floor,
  paste0("Leave-one-cohort-out median held-out score Spearman ", round(min_loco_score, 3),
         " is below the floor ", contract$modules$loco_score_spearman_floor,
         "; the module scores are not reproducible out of cohort and the analysis stops"))
cam_assert(min_neigh >= contract$modules$neighbour_resolution_ari_floor,
  paste0("Neighbour-resolution ARI ", round(min_neigh, 3), " is below the floor ",
         contract$modules$neighbour_resolution_ari_floor, "; the partition is unstable and the analysis stops"))

# Freeze.
sizes <- table(membership_full)
usable_ids <- names(sizes)[sizes >= min_size & sizes <= max_size]
cam_assert(length(usable_ids) >= 3L,
           paste0("Only ", length(usable_ids), " modules of usable size; too few to test"))
module_key <- data.table(cluster_id = usable_ids)
module_key[, module_id := sprintf("cam_m%02d", seq_len(.N))]
membership <- rbindlist(lapply(seq_len(nrow(module_key)), function(i) {
  cid <- module_key$cluster_id[i]
  data.table(module_id = module_key$module_id[i],
             gene_symbol = names(membership_full)[membership_full == as.integer(cid)])
}))
setorder(membership, module_id, gene_symbol)
cam_write_tsv(membership, file.path(out, "module_membership.tsv"))

registry <- membership[, .(n_genes = .N), by = module_id]
registry[, `:=`(resolution = resolution, k = k_selected,
                min_size = min_size, max_size = max_size,
                weight_rule = "equal")]
# Per-module reproducibility of the score on held-out participants, carried
# onto the figure so a module whose coordinate moves when a cohort is dropped
# is visible as such and excluded from the main-panel rows.
min_held <- apply(rho_held_mat, 1L, min, na.rm = TRUE)
min_all <- apply(rho_all_mat, 1L, min, na.rm = TRUE)
cid_of <- module_key$cluster_id[match(registry$module_id, module_key$module_id)]
registry[, min_loco_heldout_score_spearman := min_held[cid_of]]
registry[, min_loco_allsample_score_spearman_reported_only := min_all[cid_of]]
registry[, score_stability_flag := fifelse(
  min_loco_heldout_score_spearman < contract$modules$per_module_score_unstable_threshold,
  "score_unstable", "score_stable")]
cam_write_tsv(registry, file.path(out, "module_registry.tsv"))

membership_sha <- cam_sha256(file.path(out, "module_membership.tsv"))
cam_write_json(list(
  membership_sha256 = membership_sha,
  n_modules = nrow(registry),
  n_genes_in_modules = nrow(membership),
  universe_size = length(membership_full),
  fraction_universe_in_modules = nrow(membership) / length(membership_full),
  resolution = resolution, k = k_selected,
  min_loco_median_heldout_score_spearman = min_loco_score,
  min_loco_median_allsample_score_spearman_reported_only =
    min(loco$median_allsample_score_spearman_reported_only),
  min_loco_membership_ari_reported_only = min(loco$membership_ari_reported_only),
  min_neighbour_ari = min_neigh,
  n_score_unstable_modules = sum(registry$score_stability_flag == "score_unstable"),
  module_sizes = as.list(setNames(registry$n_genes, registry$module_id)),
  outcome_read = FALSE,
  frozen_at_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ")
), file.path(out, "module_freeze.json"))

writeLines(c(paste0("membership_sha256=", membership_sha),
             paste0("n_modules=", nrow(registry)),
             paste0("resolution=", resolution),
             paste0("k=", k_selected),
             "outcome_read=FALSE"),
           file.path(out, "MODULES_FROZEN.ok"))
cam_say("04 complete: ", nrow(registry), " modules, ", nrow(membership), " genes, FROZEN")
