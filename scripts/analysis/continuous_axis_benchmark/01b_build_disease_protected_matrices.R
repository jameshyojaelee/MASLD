#!/usr/bin/env Rscript
# 01b: build the disease-protected ComBat matrices for arm A2b.
#
# WHY THIS EXISTS, AND WHAT IT IS NOT
#
# The ComBat damage audit in 01 refit F4-vs-F0 with `dataset` already in the
# model. That is a WITHIN-cohort contrast, so by construction it cannot detect
# ComBat removing BETWEEN-cohort disease signal. Cohort is partly disease here
# (GSE126848 = 53 entirely unstaged; GSE213621 = 361 entirely staged), and the
# published ComBat drove PC1's cohort eta-squared to zero, i.e. it removed ALL
# between-cohort variance including whatever part was genuinely disease. The Q1
# null could therefore be an over-correction artifact. This script builds the
# matrix that tests that.
#
# THIS ARM IS POST-HOC AND CANNOT PASS THE Q1 GATE. The gate in
# 00_prespecification.json was sealed over arms A1-A6 and evaluated on
# 2026-08-14. A2b was added on 2026-08-17 in response to a named hole in the
# audit above. It is a diagnostic on the Q1 null, not a competitor for it, and
# `q1/GATE.json` is not reopened.
#
# THE CIRCULARITY, STATED UP FRONT. Putting fibrosis stage in ComBat's `mod`
# preserves fibrosis-associated variance BY CONSTRUCTION. An axis that then
# correlates with fibrosis is partly measuring what was protected (Nygaard,
# Rodland & Hovig 2016, Biostatistics 17:29-39, on ComBat with the group term in
# mod inflating downstream differential signal). That is why 02i runs a
# permuted-protection null: protect a within-cohort shuffle of the same labels,
# which leaves the same AMOUNT of variance unremoved but attaches it to the wrong
# donors. Only the excess of the real arm over that null is interpretable.
#
# Two matrices are built:
#   E_cbd       844 donors, stage as a factor with an explicit `unstaged` level
#   E_cbd_popB  760 staged donors only, no `unstaged` level
#
# E_cbd matches A1/A2's population so the axis is directly comparable, but its
# `unstaged` coefficient is identified only through the 31 unstaged GSE162694
# donors. E_cbd_popB removes that fragility at the cost of a different
# population, and is carried at PC1-diagnostic level only.

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(sva)
})

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
ARMS <- file.path(OUT, "arms")
assert_true(dir.exists(ARMS), paste0("No arms/ directory in ", OUT, "; run 01 first"))

pre <- read_prespec()
set.seed(pre$seeds$master)

meta <- load_manifest()
E_qn <- readRDS(file.path(ARMS, "E_qn.rds"))
assert_true(identical(colnames(E_qn), meta$sample_id),
            "E_qn columns are not in manifest order")
pca_genes <- fread(file.path(ARMS, "pca_input_genes.tsv"))$gene_id

save_rds_once <- function(x, path) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite an existing matrix: ", path))
  saveRDS(x, path)
  invisible(path)
}

# ------------------------------------------------------- protected covariate --
# `unstaged` is an explicit level, never NA and never dropped. Dropping the 84
# unstaged donors silently would change the population out from under the
# comparison with A1 and A2.
meta[, stage_protected := factor(
  ifelse(is.na(fibrosis_stage), "unstaged", paste0("F", fibrosis_stage)),
  levels = c("F0", "F1", "F2", "F3", "F4", "unstaged"))]
assert_true(!anyNA(meta$stage_protected), "stage_protected contains NA")

xtab <- as.data.table(table(dataset = meta$dataset, stage = meta$stage_protected))
write_tsv_once(xtab, file.path(ARMS, "A2b_stage_by_cohort.tsv"))

# --------------------------------------------------------------- rank audit --
# ComBat errors if a mod column is confounded with batch. Check it here so the
# failure is a legible statement about identifiability rather than an sva
# traceback, and record HOW MUCH data identifies each protected coefficient.
rank_audit <- function(info, label) {
  bmod <- model.matrix(~ -1 + droplevels(info$dataset))
  mmod <- model.matrix(~ inferred_sex + stage_protected, data = droplevels(info))
  design <- cbind(bmod, mmod[, -1, drop = FALSE])
  n_unique_cohorts_per_level <- info[, .(n = .N, n_cohorts = uniqueN(dataset)),
                                     by = stage_protected][order(stage_protected)]
  list(
    audit = data.table(
      population = label,
      n_donors = nrow(info),
      n_design_columns = ncol(design),
      design_rank = qr(design)$rank,
      full_rank = qr(design)$rank == ncol(design),
      n_batches = uniqueN(info$dataset)),
    levels = cbind(population = label, n_unique_cohorts_per_level))
}

ra_all <- rank_audit(meta, "POP-A_844")
popB <- meta[!is.na(fibrosis_stage)]
ra_b <- rank_audit(popB, "POP-B_760")

write_tsv_once(rbind(ra_all$audit, ra_b$audit),
               file.path(ARMS, "A2b_design_rank_audit.tsv"))
write_tsv_once(rbind(ra_all$levels, ra_b$levels),
               file.path(ARMS, "A2b_protected_level_support.tsv"))
print(ra_all$audit); print(ra_b$audit); print(ra_all$levels)
assert_true(ra_all$audit$full_rank,
            "The disease-protected design is rank deficient on all 844; A2b is not estimable as specified")

# --------------------------------------------------------------- the matrices --
log_step("ComBat: batch = cohort, mod = ~ inferred_sex + stage_protected (844)")
E_cbd <- ComBat(dat = E_qn,
                batch = as.character(meta$dataset),
                mod = model.matrix(~ inferred_sex + stage_protected, data = meta),
                par.prior = TRUE)
dimnames(E_cbd) <- dimnames(E_qn)
save_rds_once(E_cbd, file.path(ARMS, "E_cbd.rds"))
log_step("E_cbd persisted")

log_step("ComBat on POP-B only (760 staged, no unstaged level)")
mb <- droplevels(popB)
E_cbd_popB <- ComBat(dat = E_qn[, mb$sample_id],
                     batch = as.character(mb$dataset),
                     mod = model.matrix(~ inferred_sex + stage_protected, data = mb),
                     par.prior = TRUE)
dimnames(E_cbd_popB) <- list(rownames(E_qn), mb$sample_id)
save_rds_once(E_cbd_popB, file.path(ARMS, "E_cbd_popB.rds"))
log_step("E_cbd_popB persisted")

# ------------------------------------------------------------- PC1 diagnostics --
# Same three quantities 01 reported for E_raw / E_qn / E_cb, so the new rows drop
# straight into that table. eta-squared formula is copied from 01 unchanged.
eta2_kw <- function(x, g) {
  kw <- kruskal.test(x ~ g)
  (unname(kw$statistic) - length(unique(g)) + 1) / (length(x) - length(unique(g)))
}

pc_row <- function(mat, info, label) {
  p <- prcomp(t(mat[pca_genes, info$sample_id]), center = TRUE, scale. = FALSE)
  pc1 <- p$x[, 1]
  ok <- !is.na(info$fibrosis_stage)
  data.table(matrix_id = label,
             n_donors = nrow(info),
             pc1_pct_variance = 100 * p$sdev[1]^2 / sum(p$sdev^2),
             kruskal_eta2_pc1_on_cohort = eta2_kw(pc1, droplevels(info$dataset)),
             spearman_pc1_fibrosis = cor(pc1[ok], info$fibrosis_stage[ok], method = "spearman"))
}

log_step("PC1 diagnostics")
# E_cb restricted to POP-B is the like-for-like comparator for E_cbd_popB: same
# donors, so the only difference is what was protected.
E_cb <- readRDS(file.path(ARMS, "E_cb.rds"))
pc_diag <- rbind(
  pc_row(E_cbd, meta, "E_cbd"),
  pc_row(E_cbd_popB, popB, "E_cbd_popB"),
  pc_row(E_cb, popB, "E_cb_restricted_to_popB"),
  pc_row(E_qn, popB, "E_qn_restricted_to_popB"))
write_tsv_once(pc_diag, file.path(ARMS, "A2b_pc1_diagnostics.tsv"))
print(pc_diag)

# ----------------------------------------------------------- circularity meter --
# Refit F4-vs-F0 on the protected matrix exactly as 01 did on E_cb. Protection
# should push the attenuation slope toward 1. That is NOT a result: it is the
# quantitative statement of the circularity, and it is written next to the number
# so the two are never read apart.
log_step("F4 vs F0 attenuation on the protected matrix")
E_raw <- readRDS(file.path(ARMS, "E_raw.rds"))
fib <- population(meta, "POP-C")[fibrosis_stage %in% c(0L, 4L)]
fib[, grp := factor(ifelse(fibrosis_stage == 4L, "F4", "F0"), levels = c("F0", "F4"))]
fit_on <- function(mat, info) {
  d <- model.matrix(~ dataset + inferred_sex + grp, data = droplevels(info))
  d <- d[, qr(d)$pivot[seq_len(qr(d)$rank)], drop = FALSE]
  fit <- eBayes(lmFit(mat[, info$sample_id], d))
  topTable(fit, coef = "grpF4", number = Inf, sort.by = "none")$logFC
}
lfc_raw <- fit_on(E_raw, fib)
lfc_cbd <- fit_on(E_cbd, fib)
lfc_cb <- fit_on(E_cb, fib)
write_tsv_once(data.table(
  quantity = c("n_F0", "n_F4",
               "attenuation_slope_cb_on_raw", "attenuation_slope_cbd_on_raw",
               "pearson_r_cb_vs_raw", "pearson_r_cbd_vs_raw",
               "median_abs_lfc_raw", "median_abs_lfc_cb", "median_abs_lfc_cbd"),
  value = c(sum(fib$grp == "F0"), sum(fib$grp == "F4"),
            coef(lm(lfc_cb ~ 0 + lfc_raw))[[1]], coef(lm(lfc_cbd ~ 0 + lfc_raw))[[1]],
            cor(lfc_cb, lfc_raw), cor(lfc_cbd, lfc_raw),
            median(abs(lfc_raw)), median(abs(lfc_cb)), median(abs(lfc_cbd))),
  note = c("", "", "published ComBat", "disease-protected ComBat",
           "published ComBat", "disease-protected ComBat",
           "", "", "a slope near 1 is BY CONSTRUCTION, not evidence")),
  file.path(ARMS, "A2b_combat_damage_audit.tsv"))

writeLines(c(
  paste0("workstream\t", WORKSTREAM_ID),
  paste0("exemption_id\t", EXEMPTION_ID),
  "arm\tA2b",
  "status\tpost-hoc diagnostic; sealed 2026-08-14 gate is NOT reopened",
  paste0("seed_master\t", pre$seeds$master),
  paste0("R_version\t", R.version.string),
  paste0("run_utc\t", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))),
  file.path(OUT, "A2b_run_parameters.txt"))

log_step("A2B_SUBSTRATE_COMPLETE")
