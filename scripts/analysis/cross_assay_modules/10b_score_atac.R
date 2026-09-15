#!/usr/bin/env Rscript
# 10b: GSE296875 promoter accessibility, and the matched RNA of the same nuclei.
#
# THE WELL IS A BATCH. These donors were processed in eight wells, RNA and ATAC
# from the same nuclei, and the well alone predicts fibrosis. v1 computed a
# within-well permutation but reported an unadjusted pooled coefficient as the
# primary estimate. Here the primary fit carries the well as a fixed effect,
# the association p used by step 15 is the within-well permutation p on that
# adjusted coefficient, and the leave-one-well-out range is a sensitivity, not
# a confidence interval. Deleting one well at a time does not remove
# confounding among the retained wells; the endpoint-by-well table is saved so
# a reader can see how much within-well contrast exists at all.
#
# EXPECTATION, RECORDED IN ADVANCE. With 39 donors, 8 wells and a binary
# fibrosis label, this arm is expected to return `indeterminate` for almost
# every module, with its upper bound printed.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))

contract <- cam_contract()
out <- cam_dir("assays")
ap_dir <- file.path(cam_out_root(), "chromatin", "atac_promoter")
cam_assert(file.exists(file.path(ap_dir, "gene_promoter_counts.npy")), "Run 10a first")

membership <- fread(file.path(cam_out_root(), "modules", "module_membership.tsv"))
modules <- split(membership$gene_symbol, membership$module_id)
module_ids <- names(modules)

py <- "/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
load_npy <- function(path) {
  tsv <- paste0(path, ".tsv")
  if (!file.exists(tsv)) system2(py, c("-c", shQuote(sprintf(
    "import numpy as np; a=np.load('%s'); np.savetxt('%s', a, delimiter='\\t', fmt='%%.6f')",
    path, tsv))))
  as.matrix(fread(tsv))
}

donor_dir <- cam_input("gse296875_donor_dir", contract)
axis <- readLines(file.path(ap_dir, "gene_axis.txt"))
A <- load_npy(file.path(ap_dir, "gene_promoter_counts.npy"))
lib <- as.vector(load_npy(file.path(ap_dir, "library_size_all_peaks.npy")))
donors <- fread(file.path(donor_dir, "donor_axis.tsv"))
endpoints <- fread(cam_input("gse296875_endpoints", contract))
cam_assert(nrow(A) == nrow(donors), "ATAC matrix rows and donor axis disagree")

atac <- log2(t(A) / rep(lib, each = ncol(A)) * 1e6 + 1)
rownames(atac) <- axis
colnames(atac) <- as.character(donors$donor_id)

rna_genes <- fread(file.path(donor_dir, "rna_genes.tsv"))
R <- load_npy(file.path(donor_dir, "donor_rna_all.npy"))
rna <- log2(t(R) / rep(rowSums(R), each = ncol(R)) * 1e6 + 1)
rownames(rna) <- rna_genes$symbol
rna <- rna[!is.na(rownames(rna)) & rownames(rna) != "", , drop = FALSE]
rna <- rowsum(rna, group = rownames(rna)) /
  as.vector(table(rownames(rna))[sort(unique(rownames(rna)))])
colnames(rna) <- as.character(donors$donor_id)

em <- endpoints[match(donors$donor_id, donor_id)]
fib <- suppressWarnings(as.numeric(em$fibrosis_any == "true" | em$fibrosis_any == TRUE))
fib[em$fibrosis_observed %in% c("false", FALSE)] <- NA_real_
well <- donors$well_id
cam_say("ATAC donors ", ncol(atac), " in ", length(unique(well)), " wells; fibrosis observed ",
        sum(is.finite(fib)))

# Endpoint by well: how much within-well contrast the design offers.
by_well <- data.table(well = well, fibrosis = fib)[, .(
  n_donors = .N, n_fibrosis_observed = sum(is.finite(fibrosis)),
  n_positive = sum(fibrosis == 1, na.rm = TRUE), n_negative = sum(fibrosis == 0, na.rm = TRUE)),
  by = well][order(well)]
by_well[, both_classes := n_positive > 0 & n_negative > 0]
cam_write_tsv(by_well, file.path(out, "atac_endpoint_by_well.tsv"))
n_wells_both <- sum(by_well$both_classes)
cam_say("wells carrying both fibrosis classes: ", n_wells_both, " of ", nrow(by_well))

score_all <- function(mat) {
  z <- cam_zscore_rows(mat)
  S <- t(vapply(seq_along(module_ids), function(i) {
    idx <- intersect(modules[[i]], rownames(z))
    if (!length(idx)) return(rep(NA_real_, ncol(z)))
    standardize_vector(colMeans(z[idx, , drop = FALSE], na.rm = TRUE))
  }, numeric(ncol(mat))))
  rownames(S) <- module_ids
  list(scores = S, measured = rownames(z))
}
atac_sc <- score_all(atac)
rna_sc <- score_all(rna)
cam_write_rds(atac_sc$scores, file.path(ap_dir, "atac_module_scores.rds"))
cam_write_rds(rna_sc$scores, file.path(ap_dir, "gse296875_rna_module_scores.rds"))
writeLines(atac_sc$measured, file.path(ap_dir, "atac_measured_genes.txt"))
writeLines(rna_sc$measured, file.path(ap_dir, "gse296875_rna_measured_genes.txt"))
mod_genes <- unique(membership$gene_symbol)
za <- cam_zscore_rows(atac); za <- za[rownames(za) %in% mod_genes, , drop = FALSE]
zr <- cam_zscore_rows(rna); zr <- zr[rownames(zr) %in% mod_genes, , drop = FALSE]
cam_write_rds(za, file.path(ap_dir, "atac_module_gene_z.rds"))
cam_write_rds(zr, file.path(ap_dir, "gse296875_rna_module_gene_z.rds"))

# Well-adjusted primary fit on the donors with a recorded endpoint. All 160
# modules share one design, so every fit is one matrix operation.
n_perm <- 10000L
yz <- cam_endpoint_z(fib)
ok <- is.finite(yz)
S <- atac_sc$scores[, ok, drop = FALSE]
y <- yz[ok]; well_ok <- well[ok]
S[!is.finite(S)] <- NA_real_
usable_row <- stats::complete.cases(S)
Sfit <- S; Sfit[!usable_row, ] <- 0
W <- stats::model.matrix(~ factor(well_ok))[, -1, drop = FALSE]

obs_fit <- cam_fast_assoc(Sfit, y, W)
unadj_fit <- cam_fast_assoc(Sfit, y)
cam_say("observed well-adjusted fit on ", sum(ok), " donors with recorded fibrosis")

lowo_beta <- matrix(NA_real_, nrow(S), length(unique(well_ok)),
                    dimnames = list(module_ids, unique(well_ok)))
for (w in unique(well_ok)) {
  keep <- well_ok != w
  if (sum(keep) < 12L) next
  Wk <- stats::model.matrix(~ factor(well_ok[keep]))[, -1, drop = FALSE]
  lowo_beta[, w] <- cam_fast_assoc(Sfit[, keep, drop = FALSE], y[keep], Wk)$beta
}

null_abs <- matrix(NA_real_, nrow(S), n_perm)
for (b in seq_len(n_perm)) {
  yp <- cam_permute_within(y, well_ok, contract$seed + b)
  null_abs[, b] <- abs(cam_fast_assoc(Sfit, yp, W)$beta)
}
cam_say("within-well permutation complete (", n_perm, " draws)")

rows <- list()
for (i in seq_along(module_ids)) {
  tb <- cam_testability(modules[[i]], atac_sc$measured, contract)
  estimable <- tb$testable && usable_row[i]
  obs <- if (estimable) abs(obs_fit$beta[i]) else NA_real_
  nullb <- null_abs[i, ]
  lw <- lowo_beta[i, ]
  rows[[length(rows) + 1L]] <- data.table(
    assay = "atac_gse296875", module_id = module_ids[i],
    endpoint = "fibrosis_binary", endpoint_family = "fibrosis",
    endpoint_note = "any fibrosis vs none, z-scored; well as fixed effect",
    unit = "biological_donor",
    n_units = if (estimable) sum(ok) else 0L,
    df = if (estimable) obs_fit$df else NA_integer_,
    testable = tb$testable, n_members = tb$n_members, n_measured = tb$n_measured,
    fraction_measured = tb$fraction_measured,
    beta = if (estimable) obs_fit$beta[i] else NA_real_,
    se = if (estimable) obs_fit$se[i] else NA_real_,
    p_two_sided = if (estimable) obs_fit$p_two_sided[i] else NA_real_,
    unadjusted_beta_reported_only = if (estimable) unadj_fit$beta[i] else NA_real_,
    lowo_min_beta = if (all(is.na(lw))) NA_real_ else min(lw, na.rm = TRUE),
    lowo_max_beta = if (all(is.na(lw))) NA_real_ else max(lw, na.rm = TRUE),
    lowo_all_same_sign = if (all(is.na(lw))) NA else length(unique(sign(lw[!is.na(lw)]))) == 1L,
    within_well_null_mean_abs = mean(nullb, na.rm = TRUE),
    within_well_exceedance_count = sum(nullb >= obs, na.rm = TRUE),
    within_well_draws = sum(is.finite(nullb)),
    within_well_p = cam_empirical_p_ge(obs, nullb)
  )
}
res <- rbindlist(rows)
cam_assert_no_prohibited_columns(res, contract)
cam_write_tsv(res, file.path(out, "atac_results.tsv"))
cam_write_json(list(
  donors = ncol(atac), wells = length(unique(well)),
  fibrosis_observed = sum(is.finite(fib)),
  fibrosis_positive = sum(fib == 1, na.rm = TRUE),
  wells_with_both_classes = n_wells_both,
  min_wells_with_both_classes_required = contract$assays$atac_gse296875$min_wells_with_both_classes,
  genes_with_promoter_peak = nrow(atac),
  modules_testable = sum(res$testable),
  within_well_draws = n_perm,
  within_well_null_mean_of_mean_abs_beta = mean(res$within_well_null_mean_abs, na.rm = TRUE),
  direction_applied_at_scoring = FALSE,
  note = "primary beta is well-adjusted; the association p consumed by 15 is within_well_p"
), file.path(out, "atac_summary.json"))
writeLines("atac scored", file.path(out, "READY_atac"))
cam_say("10b complete")
