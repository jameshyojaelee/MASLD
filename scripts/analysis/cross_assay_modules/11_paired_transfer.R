#!/usr/bin/env Rscript
# 11: within-participant transfer.
#
# Every other column compares an association measured in one cohort with an
# association measured in another. Here the same people are measured twice, so
# the question becomes whether a participant whose RNA module score is high
# also has a high chromatin score, holding the histology constant.
#
# CONDITIONING IS THE POINT, AND IT IS NOT MECHANISM. Without conditioning both
# scores track severity and correlate for that reason alone. A partial
# correlation given histology removes that route; it does not establish a
# shared regulatory mechanism, because cell mixture, residual disease variation,
# genotype and technical variation remain possible common causes.
#
# v2 corrections. A module is tested only if it meets the coverage rule in BOTH
# assays of the pair. The bootstrap resamples participant triples (RNA score,
# chromatin score, covariates) and refits the partial statistic; v1 bootstrapped
# the raw rank correlation and reported that interval beside a partial estimate.
# The GSE296875 leave-one-well-out range is labelled as a sensitivity, not an
# interval. The sign of the coupling is a reported column and a figure aesthetic.
# Three sensitivities are added: both scores recomputed on the members measured
# in both assays; for GSE296875 the correlation given well AND fibrosis; and for
# GSE296875 the correlation given well AND the donor's hepatocyte nuclear
# fraction. The last one matters most: RNA and ATAC come from the SAME NUCLEI,
# so a donor with more hepatocyte nuclei carries both a higher hepatocyte-gene
# RNA score and a higher hepatocyte-promoter ATAC score, and that common cause
# would produce coupling with no regulatory relationship at all.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))

contract <- cam_contract()
out <- cam_dir("assays")
chrom <- file.path(cam_out_root(), "chromatin")
cam_assert(file.exists(file.path(out, "READY_h3k27ac")), "Run 09b first")
cam_assert(file.exists(file.path(out, "READY_atac")), "Run 10b first")

membership <- fread(file.path(cam_out_root(), "modules", "module_membership.tsv"))
modules <- split(membership$gene_symbol, membership$module_id)
labels <- fread(file.path(cam_out_root(), "discovery", "module_discovery_labels.tsv"))
module_ids <- labels$module_id
n_family <- length(module_ids)
n_boot <- contract$paired_transfer$gse267145$bootstrap_replicates
n_null <- 2000L

partial_cor <- function(x, y, engine) {
  suppressWarnings(stats::cor(engine$resid(rank(x)), engine$resid(rank(y))))
}
shared_member_scores <- function(zx, zy, members) {
  g <- intersect(intersect(members, rownames(zx)), rownames(zy))
  if (length(g) < 2L) return(list(x = NULL, y = NULL, n = length(g)))
  list(x = standardize_vector(colMeans(zx[g, , drop = FALSE], na.rm = TRUE)),
       y = standardize_vector(colMeans(zy[g, , drop = FALSE], na.rm = TRUE)), n = length(g))
}
empty_row <- function(pairing, mid, tb_r, tb_c, n) data.table(
  pairing = pairing, module_id = mid, n_participants = n,
  testable_rna = tb_r$testable, testable_chromatin = tb_c$testable,
  testable_both = tb_r$testable && tb_c$testable,
  fraction_measured_rna = tb_r$fraction_measured, fraction_measured_chromatin = tb_c$fraction_measured,
  raw_spearman = NA_real_, partial_spearman = NA_real_, coupling_sign = NA_character_,
  boot_ci_low = NA_real_, boot_ci_high = NA_real_, lowo_min = NA_real_, lowo_max = NA_real_,
  n_shared_members = NA_integer_, shared_member_partial_spearman = NA_real_,
  partial_given_well_and_fibrosis = NA_real_, partial_given_well_and_hepatocyte_fraction = NA_real_,
  p_given_well_and_fibrosis = NA_real_, p_given_well_and_hepatocyte_fraction = NA_real_,
  null_mean = NA_real_, exceedance_count = NA_integer_, n_null_draws = 0L, empirical_p = NA_real_)

rows <- list()

# --- GSE267145: RNA vs H3K27ac, 99 participants ------------------------------
h3 <- readRDS(file.path(chrom, "h3k27ac_module_scores.rds"))
rn <- readRDS(file.path(chrom, "gse267145_rna_module_scores.rds"))
zh <- readRDS(file.path(chrom, "h3k27ac_module_gene_z.rds"))
zr1 <- readRDS(file.path(chrom, "gse267145_rna_module_gene_z.rds"))
meas_h <- readLines(file.path(chrom, "h3k27ac_measured_genes.txt"))
meas_r1 <- readLines(file.path(chrom, "gse267145_rna_measured_genes.txt"))
shared <- intersect(colnames(h3), colnames(rn))
cam_say("GSE267145 paired participants: ", length(shared))
h3 <- h3[, shared, drop = FALSE]; rn <- rn[, shared, drop = FALSE]
zh <- zh[, shared, drop = FALSE]; zr1 <- zr1[, shared, drop = FALSE]
join <- fread(cam_input("gse267145_join", contract))
jm <- join[match(shared, participant_id)]
Z <- as.matrix(data.table(
  steatosis = suppressWarnings(as.numeric(jm$steatosis)),
  ballooning = suppressWarnings(as.numeric(jm$ballooning)),
  lobular_inflammation = suppressWarnings(as.numeric(jm$lobular_inflammation)),
  fibrosis = suppressWarnings(as.numeric(jm$fibrosis)),
  sex = as.numeric(factor(jm$sex))))
ok_all <- stats::complete.cases(Z)
engine <- cam_partial_engine(Z[ok_all, , drop = FALSE])
# Bootstrap designs are module-independent: one engine per resample, reused.
set.seed(contract$seed + 500000L)
n_ok <- sum(ok_all)
boot_idx <- matrix(sample.int(n_ok, n_ok * n_boot, replace = TRUE), nrow = n_ok)
boot_engines <- lapply(seq_len(n_boot), function(b)
  cam_partial_engine(Z[ok_all, , drop = FALSE][boot_idx[, b], , drop = FALSE]))
for (i in seq_len(n_family)) {
  mid <- module_ids[i]
  tb_r <- cam_testability(modules[[mid]], meas_r1, contract)
  tb_c <- cam_testability(modules[[mid]], meas_h, contract)
  x <- rn[mid, ]; y <- h3[mid, ]
  ok <- is.finite(x) & is.finite(y) & ok_all
  if (!(tb_r$testable && tb_c$testable) || sum(ok) < 20L) {
    rows[[length(rows) + 1L]] <- empty_row("gse267145_rna_vs_h3k27ac", mid, tb_r, tb_c, sum(ok)); next
  }
  same <- identical(ok, ok_all)
  eng <- if (same) engine else cam_partial_engine(Z[ok, , drop = FALSE])
  xo <- x[ok]; yo <- y[ok]
  raw <- suppressWarnings(stats::cor(xo, yo, method = "spearman"))
  ex <- eng$resid(rank(xo)); ey <- eng$resid(rank(yo))
  partial <- suppressWarnings(stats::cor(ex, ey))
  # Freedman-Lane: only the covariate-free residual of x is permuted.
  set.seed(contract$seed + i)
  EX <- vapply(seq_len(n_null), function(b) sample(ex), numeric(length(ex)))
  EXr <- EX - eng$design %*% qr.coef(qr(eng$design), EX)
  nullv <- as.vector(stats::cor(EXr, ey))
  # Bootstrap of participant triples with the partial statistic refitted.
  boot <- if (same) {
    vapply(seq_len(n_boot), function(b) {
      j <- boot_idx[, b]; partial_cor(xo[j], yo[j], boot_engines[[b]]) }, numeric(1))
  } else {
    set.seed(contract$seed + 600000L + i)
    vapply(seq_len(n_boot), function(b) {
      j <- sample.int(length(xo), replace = TRUE)
      partial_cor(xo[j], yo[j], cam_partial_engine(Z[ok, , drop = FALSE][j, , drop = FALSE])) }, numeric(1))
  }
  sm <- shared_member_scores(zr1, zh, modules[[mid]])
  sm_partial <- if (is.null(sm$x)) NA_real_ else partial_cor(sm$x[ok], sm$y[ok], eng)
  rows[[length(rows) + 1L]] <- data.table(
    pairing = "gse267145_rna_vs_h3k27ac", module_id = mid, n_participants = sum(ok),
    testable_rna = TRUE, testable_chromatin = TRUE, testable_both = TRUE,
    fraction_measured_rna = tb_r$fraction_measured, fraction_measured_chromatin = tb_c$fraction_measured,
    raw_spearman = raw, partial_spearman = partial,
    coupling_sign = if (partial >= 0) "positive" else "inverse",
    boot_ci_low = unname(stats::quantile(boot, 0.025, na.rm = TRUE)),
    boot_ci_high = unname(stats::quantile(boot, 0.975, na.rm = TRUE)),
    lowo_min = NA_real_, lowo_max = NA_real_,
    n_shared_members = sm$n, shared_member_partial_spearman = sm_partial,
    partial_given_well_and_fibrosis = NA_real_, partial_given_well_and_hepatocyte_fraction = NA_real_,
    p_given_well_and_fibrosis = NA_real_, p_given_well_and_hepatocyte_fraction = NA_real_,
    null_mean = mean(nullv, na.rm = TRUE),
    exceedance_count = sum(abs(nullv) >= abs(partial), na.rm = TRUE),
    n_null_draws = sum(is.finite(nullv)),
    empirical_p = cam_empirical_p_ge(abs(partial), abs(nullv)))
}

# --- GSE296875: RNA vs ATAC, same nuclei, 39 donors --------------------------
# Same nuclei means the two scores share every technical source, including the
# well. The correlation is computed given well, and the null permutes within
# well; the sensitivity adds fibrosis so the reader can see what severity
# explains.
ap <- file.path(chrom, "atac_promoter")
at <- readRDS(file.path(ap, "atac_module_scores.rds"))
ar <- readRDS(file.path(ap, "gse296875_rna_module_scores.rds"))
za <- readRDS(file.path(ap, "atac_module_gene_z.rds"))
zr2 <- readRDS(file.path(ap, "gse296875_rna_module_gene_z.rds"))
meas_a <- readLines(file.path(ap, "atac_measured_genes.txt"))
meas_r2 <- readLines(file.path(ap, "gse296875_rna_measured_genes.txt"))
donors <- fread(file.path(cam_input("gse296875_donor_dir", contract), "donor_axis.tsv"))
ep <- fread(cam_input("gse296875_endpoints", contract))[match(donors$donor_id, donor_id)]
fib <- suppressWarnings(as.numeric(ep$fibrosis_any == "true" | ep$fibrosis_any == TRUE))
fib[ep$fibrosis_observed %in% c("false", FALSE)] <- NA_real_
well <- donors$well_id
cam_say("GSE296875 paired donors: ", ncol(at), " in ", length(unique(well)), " wells")
Wd <- stats::model.matrix(~ factor(well))[, -1, drop = FALSE]
eng_w <- cam_partial_engine(Wd)
okf <- is.finite(fib)
eng_wf <- cam_partial_engine(cbind(Wd[okf, , drop = FALSE], fibrosis = fib[okf]))
# Hepatocyte nuclear fraction: the composition common cause of the two scores.
hep <- donors$hepatocyte_nuclei / donors$nuclei
cam_say("hepatocyte nuclear fraction range ",
        paste(round(range(hep, na.rm = TRUE), 3), collapse = "-"))
okh <- is.finite(hep)
eng_wh <- cam_partial_engine(cbind(Wd[okh, , drop = FALSE], hepatocyte_fraction = hep[okh]))
for (i in seq_len(n_family)) {
  mid <- module_ids[i]
  tb_r <- cam_testability(modules[[mid]], meas_r2, contract)
  tb_c <- cam_testability(modules[[mid]], meas_a, contract)
  x <- ar[mid, ]; y <- at[mid, ]
  ok <- is.finite(x) & is.finite(y)
  if (!(tb_r$testable && tb_c$testable) || sum(ok) < 20L) {
    rows[[length(rows) + 1L]] <- empty_row("gse296875_rna_vs_atac", mid, tb_r, tb_c, sum(ok)); next
  }
  eng <- if (all(ok)) eng_w else cam_partial_engine(Wd[ok, , drop = FALSE])
  raw <- suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
  ey <- eng$resid(rank(y[ok]))
  partial <- suppressWarnings(stats::cor(eng$resid(rank(x[ok])), ey))
  nullv <- vapply(seq_len(n_null), function(b) {
    xp <- cam_permute_within(x[ok], well[ok], contract$seed + 900000L + b * 173L + i)
    suppressWarnings(stats::cor(eng$resid(rank(xp)), ey))
  }, numeric(1))
  lowo <- vapply(unique(well), function(w) {
    keep <- ok & well != w
    if (sum(keep) < 12L) return(NA_real_)
    e2 <- cam_partial_engine(Wd[keep, , drop = FALSE])
    partial_cor(x[keep], y[keep], e2)
  }, numeric(1))
  sm <- shared_member_scores(zr2, za, modules[[mid]])
  sm_partial <- if (is.null(sm$x)) NA_real_ else partial_cor(sm$x[ok], sm$y[ok], eng)
  # Each sensitivity is a test in its own right: the same within-well
  # permutation of the RNA score, evaluated with the extended adjustment.
  sens <- function(okx, extra, seed_off) {
    if (sum(okx) < 20L) return(list(r = NA_real_, p = NA_real_))
    e <- cam_partial_engine(cbind(Wd[okx, , drop = FALSE], extra[okx]))
    eyx <- e$resid(rank(y[okx]))
    r <- suppressWarnings(stats::cor(e$resid(rank(x[okx])), eyx))
    nv <- vapply(seq_len(n_null), function(b) {
      xp <- cam_permute_within(x[okx], well[okx], contract$seed + seed_off + b * 173L + i)
      suppressWarnings(stats::cor(e$resid(rank(xp)), eyx))
    }, numeric(1))
    list(r = r, p = cam_empirical_p_ge(abs(r), abs(nv)))
  }
  sf <- sens(ok & okf, fib, 910000L)
  sh <- sens(ok & okh, hep, 920000L)
  pwf <- sf$r; pwh <- sh$r
  rows[[length(rows) + 1L]] <- data.table(
    pairing = "gse296875_rna_vs_atac", module_id = mid, n_participants = sum(ok),
    testable_rna = TRUE, testable_chromatin = TRUE, testable_both = TRUE,
    fraction_measured_rna = tb_r$fraction_measured, fraction_measured_chromatin = tb_c$fraction_measured,
    raw_spearman = raw, partial_spearman = partial,
    coupling_sign = if (partial >= 0) "positive" else "inverse",
    boot_ci_low = NA_real_, boot_ci_high = NA_real_,
    lowo_min = if (all(is.na(lowo))) NA_real_ else min(lowo, na.rm = TRUE),
    lowo_max = if (all(is.na(lowo))) NA_real_ else max(lowo, na.rm = TRUE),
    n_shared_members = sm$n, shared_member_partial_spearman = sm_partial,
    partial_given_well_and_fibrosis = pwf, partial_given_well_and_hepatocyte_fraction = pwh,
    p_given_well_and_fibrosis = sf$p, p_given_well_and_hepatocyte_fraction = sh$p,
    null_mean = mean(nullv, na.rm = TRUE),
    exceedance_count = sum(abs(nullv) >= abs(partial), na.rm = TRUE),
    n_null_draws = sum(is.finite(nullv)),
    empirical_p = cam_empirical_p_ge(abs(partial), abs(nullv)))
}

res <- rbindlist(rows, fill = TRUE)
# BH within each pairing over the 160-module family; only modules testable on
# both sides carry a p.
res[, q_value := ml_complete_bh(fifelse(testable_both, empirical_p, NA_real_), n_family), by = pairing]
res[, q_given_well_and_fibrosis := ml_complete_bh(fifelse(testable_both, p_given_well_and_fibrosis, NA_real_), n_family), by = pairing]
res[, q_given_well_and_hepatocyte_fraction := ml_complete_bh(fifelse(testable_both, p_given_well_and_hepatocyte_fraction, NA_real_), n_family), by = pairing]
cam_assert_no_prohibited_columns(res, contract)
cam_write_tsv(res, file.path(out, "paired_transfer_results.tsv"))
summ <- res[, .(
  n_testable_both = sum(testable_both), n_untestable_rna = sum(!testable_rna),
  n_untestable_chromatin = sum(!testable_chromatin),
  n_q05 = sum(q_value < 0.05, na.rm = TRUE),
  n_q05_positive = sum(q_value < 0.05 & coupling_sign == "positive", na.rm = TRUE),
  n_q05_inverse = sum(q_value < 0.05 & coupling_sign == "inverse", na.rm = TRUE),
  median_null_mean = stats::median(null_mean, na.rm = TRUE),
  max_abs_partial = max(abs(partial_spearman), na.rm = TRUE),
  median_attenuation_given_hepatocyte_fraction = stats::median(
    partial_given_well_and_hepatocyte_fraction[q_value < 0.05] /
      partial_spearman[q_value < 0.05], na.rm = TRUE),
  n_q05_retaining_half_under_hepatocyte_fraction = sum(
    q_value < 0.05 & abs(partial_given_well_and_hepatocyte_fraction) >= 0.5 * abs(partial_spearman),
    na.rm = TRUE),
  n_q05_and_q05_given_fibrosis = sum(q_value < 0.05 & q_given_well_and_fibrosis < 0.05, na.rm = TRUE),
  n_q05_and_q05_given_hepatocyte_fraction = sum(q_value < 0.05 & q_given_well_and_hepatocyte_fraction < 0.05, na.rm = TRUE)), by = pairing]
cam_write_tsv(summ, file.path(out, "paired_transfer_summary_by_pairing.tsv"))
cam_write_json(list(
  gse267145_participants = res[pairing == "gse267145_rna_vs_h3k27ac", max(n_participants)],
  gse296875_donors = res[pairing == "gse296875_rna_vs_atac", max(n_participants)],
  by_pairing = summ, bh_family = n_family, bootstrap_replicates = n_boot, null_draws = n_null,
  adjustment_gse267145 = "steatosis, ballooning, lobular inflammation, fibrosis, sex",
  adjustment_gse296875 = "well (fixed effect); sensitivities add fibrosis, and hepatocyte nuclear fraction",
  hepatocyte_fraction_note = "RNA and ATAC come from the same nuclei; donor hepatocyte fraction is a common cause of both scores"
), file.path(out, "paired_transfer_summary.json"))
writeLines("paired transfer done", file.path(out, "READY_paired"))
cam_say("11 complete")
