#!/usr/bin/env Rscript
# F5ef: protein-histology association for all 117 frozen programs (v2 registry and membership).
#
# Same assay, normalization, scoring rule and test as the superseded 22-program Fig 5E column:
#   - PXD051911 liver DIA-MS, single-gene protein groups, log2, <=50% missing, quantile normalized,
#     duplicate symbols collapsed by median, gene-wise z across 58 samples
#     (copied from Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R).
#   - Program score = L1-weighted mean of measured member-protein z; testable if >=8 measured members
#     and retained L1 weight >= 0.20 (same file).
#   - Histologic burden = PC1 of steatosis, ballooning, inflammation, fibrosis in MASLD participants;
#     rank-residual partial Spearman adjusted for acquisition batch, age, BMI, sex;
#     batch-stratified Freedman-Lane permutation p (9,999); patient bootstrap CI within batch (4,999);
#     BH over all testable programs; robust = q < 0.05 and equal-weight, leave-top-gene and
#     both batch signs agree (copied from 02b_histology_burden.R).
# Unit of inference: participant (46 MASLD participants with complete histology).
# Change from the superseded column: program universe is the complete 117-program v2 registry,
# not the 22 programs selected on the superseded July disease-stage fit; BH family = testable of 117.
#
# Gates (stop on failure):
#   G-A0 derived metadata equals results/proteomics/panel4c_metadata.tsv
#   G-A  v1 22-program scores reproduce results/proteomics/module_protein_scores.tsv (max |d| < 1e-10)
#   G-B  the copied 02b test on the stored 22-program scores reproduces
#        results/proteomics/module_protein_histology_burden_masld.tsv (rho, p, CI, padj; |d| < 1e-10)
# Seeds: 117-program run uses 20260923 + program index (permutation) and 30260923 + index (bootstrap).

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(parallel)
})

P <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
WS <- file.path(P, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef")
OUT <- file.path(WS, "out")
PR <- file.path(P, "Analysis/Multimodal_Program_Projection/results/proteomics")
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
N_PERM <- 9999L
N_BOOT <- 4999L

RAW_FILE <- file.path(P, "data/PXD051911/liver_protein_quant.txt")
META_FILE <- file.path(P, "data/PXD051911/meta_data.txt")
V1_REG <- file.path(P, "Analysis/Multimodal_Program_Projection/results/frozen_programs.tsv")
V1_MEM <- file.path(P, "Analysis/Multimodal_Program_Projection/results/frozen_program_membership.tsv")
V2_REG <- file.path(P, "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv")
V2_MEM <- file.path(P, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv")

log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)

# ---------------------------------------------------------------- data (copied from 02)
raw <- fread(RAW_FILE)
meta <- fread(META_FILE)
sample_cols <- setdiff(names(raw), c("ProteinAccessions", "Genes", "ProteinDescriptions"))
meta <- meta[!is.na(liver_proteomics_filename) & liver_proteomics_filename != "" &
               liver_proteomics_filename %in% sample_cols]
meta <- unique(meta, by = "liver_proteomics_filename")
stopifnot(nrow(meta) == 58L)
meta[, sample_id := liver_proteomics_filename]
meta[, group := factor(fifelse(saf_diagnosis == "No_MASLD", "Control", "MASLD"), levels = c("Control", "MASLD"))]
meta[, acquisition_batch := factor(fifelse(grepl("^2019", sample_id), "2019", "2020"), levels = c("2019", "2020"))]
meta[, `:=`(
  age = as.numeric(alder), bmi_num = as.numeric(bmi), sex = factor(gender),
  Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)),
  Steatosis = as.numeric(steatosis_score),
  Ballooning = as.numeric(hepatocellular_ballooning_score),
  Inflammation = as.numeric(lobular_inflammation_score),
  NAS = as.numeric(nafld_activity_score)
)]
meta[, `:=`(age_z = as.numeric(scale(age)), bmi_z = as.numeric(scale(bmi_num)))]
meta[, sample_order := match(sample_id, sample_cols)]
setorder(meta, sample_order)
samples <- meta$sample_id

raw <- raw[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE)]
expr <- as.matrix(raw[, ..samples])
storage.mode(expr) <- "double"
expr[!is.finite(expr) | expr <= 0] <- NA_real_
expr <- log2(expr)
rownames(expr) <- raw$Genes
keep <- rowMeans(is.na(expr)) <= 0.50
expr <- expr[keep, , drop = FALSE]
expr <- normalizeBetweenArrays(expr, method = "quantile")
collapse_gene <- function(idx) {
  z <- expr[idx, , drop = FALSE]
  if (nrow(z) == 1L) return(as.numeric(z[1, ]))
  apply(z, 2, median, na.rm = TRUE)
}
idx_by_gene <- split(seq_len(nrow(expr)), rownames(expr))
gene_expr <- t(vapply(idx_by_gene, collapse_gene, numeric(ncol(expr))))
colnames(gene_expr) <- colnames(expr)
gene_expr[!is.finite(gene_expr)] <- NA_real_
z_all <- t(scale(t(gene_expr[, samples, drop = FALSE])))
z_all[!is.finite(z_all)] <- NA_real_
log_msg("proteins after filter: ", nrow(gene_expr))

# G-A0
stored_meta <- fread(file.path(PR, "panel4c_metadata.tsv"))
chk <- meta[, .(sample_id, group = as.character(group), acquisition_batch = as.character(acquisition_batch),
                sex = as.character(sex), age, bmi = bmi_num, Steatosis, Ballooning, Inflammation, Fibrosis, NAS)]
stored_meta[, `:=`(group = as.character(group), acquisition_batch = as.character(acquisition_batch), sex = as.character(sex))]
ga0 <- isTRUE(all.equal(as.data.frame(chk), as.data.frame(stored_meta[, names(chk), with = FALSE]), check.attributes = FALSE))
log_msg("G-A0 metadata identical: ", ga0)
stopifnot(ga0)

# ---------------------------------------------------------------- scoring (copied from 02)
score_program <- function(m) {
  m <- m[, .(original_l1_weight = sum(original_l1_weight)), by = gene_symbol]
  measured <- m[gene_symbol %in% rownames(z_all)]
  n_measured <- nrow(measured)
  retained_weight <- sum(measured$original_l1_weight)
  testable <- n_measured >= 8L && retained_weight >= 0.20
  out <- list(n_measured = n_measured, retained_l1_weight = retained_weight, testable = testable, scores = NULL)
  if (!testable) return(out)
  measured[, weight := original_l1_weight / sum(original_l1_weight)]
  zz <- z_all[measured$gene_symbol, samples, drop = FALSE]
  observed_weight <- colSums(is.finite(zz) * measured$weight)
  weighted_score <- colSums(zz * measured$weight, na.rm = TRUE) / observed_weight
  weighted_score[!is.finite(weighted_score)] <- NA_real_
  equal_score <- colMeans(zz, na.rm = TRUE)
  equal_score[!is.finite(equal_score)] <- NA_real_
  top_gene <- measured$gene_symbol[which.max(measured$weight)]
  leave <- measured[gene_symbol != top_gene]
  leave[, weight := original_l1_weight / sum(original_l1_weight)]
  leave_zz <- z_all[leave$gene_symbol, samples, drop = FALSE]
  leave_observed_weight <- colSums(is.finite(leave_zz) * leave$weight)
  leave_score <- colSums(leave_zz * leave$weight, na.rm = TRUE) / leave_observed_weight
  leave_score[!is.finite(leave_score)] <- NA_real_
  out$scores <- data.table(sample_id = samples, group = as.character(meta$group),
                           acquisition_batch = as.character(meta$acquisition_batch),
                           score = weighted_score, equal_score = equal_score, leave_top_score = leave_score)
  out
}

# G-A: v1 22-program rescoring
v1_reg <- fread(V1_REG)
v1_mem <- fread(V1_MEM)
v1_scores <- rbindlist(lapply(v1_reg$program_id, function(pid) {
  m <- v1_mem[program_id == pid & mapped_symbol == TRUE & !is.na(gene_symbol)]
  s <- score_program(m)
  if (is.null(s$scores)) return(NULL)
  cbind(program_id = pid, s$scores)
}))
stored_scores <- fread(file.path(PR, "module_protein_scores.tsv"))
cmp <- merge(stored_scores, v1_scores, by = c("program_id", "sample_id"), suffixes = c(".s", ".n"))
ga <- nrow(cmp) == nrow(stored_scores) && nrow(v1_scores) == nrow(stored_scores)
dmax <- max(abs(c(cmp$score.s - cmp$score.n, cmp$equal_score.s - cmp$equal_score.n,
                  cmp$leave_top_score.s - cmp$leave_top_score.n)), na.rm = TRUE)
na_agree <- all(is.na(cmp$score.s) == is.na(cmp$score.n))
log_msg("G-A rows ", nrow(cmp), "/", nrow(stored_scores), " max|d| ", signif(dmax, 3), " NA pattern agree ", na_agree)
stopifnot(ga, dmax < 1e-10, na_agree)

# v2 117-program scoring
v2_reg <- fread(V2_REG, select = c("cell_type", "module", "program_uid", "module_name"))
stopifnot(nrow(v2_reg) == 117L)
v2_mem <- fread(V2_MEM)
v2_mem <- v2_mem[mapped_symbol_status == "gencode_v49_unique_symbol_confirmed"]
v2_mem[, gene_symbol := trimws(mapped_symbol)]
setorder(v2_reg, cell_type, module)
v2_cov <- list()
v2_scores <- rbindlist(lapply(seq_len(nrow(v2_reg)), function(i) {
  pid <- v2_reg$program_uid[i]
  s <- score_program(v2_mem[program_uid == pid, .(gene_symbol, original_l1_weight)])
  v2_cov[[pid]] <<- data.table(program_uid = pid, n_measured = s$n_measured,
                               retained_l1_weight = s$retained_l1_weight, testable_scoring = s$testable)
  if (is.null(s$scores)) return(NULL)
  cbind(program_id = pid, s$scores)
}))
v2_cov <- rbindlist(v2_cov)
log_msg("v2 programs testable for protein scoring: ", sum(v2_cov$testable_scoring), "/117")
fwrite(v2_scores, file.path(OUT, "protein_program_scores_117.tsv"), sep = "\t")

# ---------------------------------------------------------------- histology test (copied from 02b)
components <- c("Steatosis", "Ballooning", "Inflammation", "Fibrosis")
pmeta <- fread(file.path(PR, "panel4c_metadata.tsv"))  # fresh read: batch stays integer as in 02b
meta_masld <- pmeta[group == "MASLD" & complete.cases(pmeta[, ..components])]
derive_burden <- function(d) {
  z <- scale(as.matrix(d[, ..components]))
  if (any(!is.finite(z))) return(NULL)
  fit <- prcomp(z, center = FALSE, scale. = FALSE)
  loading <- fit$rotation[, 1L]
  score <- fit$x[, 1L]
  if (sum(loading) < 0) { loading <- -loading; score <- -score }
  list(score = as.numeric(score), loading = loading, variance_explained = summary(fit)$importance[2L, 1L])
}
axis_fit <- derive_burden(meta_masld)
stopifnot(!is.null(axis_fit), all(axis_fit$loading > 0))
meta_masld[, histology_burden := axis_fit$score]
stored_axis <- fread(file.path(PR, "histology_burden_axis.tsv"))
log_msg("axis loadings max|d| vs stored: ", signif(max(abs(stored_axis$loading - axis_fit$loading[stored_axis$feature])), 3),
        "; n MASLD ", nrow(meta_masld))

rank_residualize <- function(y, d, include_batch = TRUE) {
  rhs <- c("age", "bmi")
  if (include_batch && uniqueN(d$acquisition_batch[!is.na(d$acquisition_batch)]) > 1L) rhs <- c("acquisition_batch", rhs)
  if (uniqueN(d$sex[!is.na(d$sex)]) > 1L) rhs <- c(rhs, "sex")
  nuisance <- model.matrix(reformulate(rhs), data = d)
  ok <- is.finite(y) & apply(nuisance, 1L, function(z) all(is.finite(z)))
  out <- rep(NA_real_, length(y))
  if (sum(ok) > ncol(nuisance) + 2L) {
    out[ok] <- residuals(lm.fit(nuisance[ok, , drop = FALSE], rank(y[ok], ties.method = "average")))
  }
  out
}
partial_spearman <- function(d, score_col, include_batch = TRUE) {
  score_resid <- rank_residualize(d[[score_col]], d, include_batch)
  burden_resid <- rank_residualize(d$histology_burden, d, include_batch)
  ok <- is.finite(score_resid) & is.finite(burden_resid)
  if (sum(ok) < 8L) return(c(n = sum(ok), rho = NA_real_, pvalue_asymptotic = NA_real_))
  tst <- suppressWarnings(cor.test(score_resid[ok], burden_resid[ok], method = "pearson"))
  c(n = sum(ok), rho = unname(tst$estimate), pvalue_asymptotic = tst$p.value)
}
freedman_lane_pvalue <- function(d, score_col, n_perm, seed) {
  rhs <- c("acquisition_batch", "age", "bmi")
  if (uniqueN(d$sex[!is.na(d$sex)]) > 1L) rhs <- c(rhs, "sex")
  nuisance <- model.matrix(reformulate(rhs), data = d)
  ok <- is.finite(d[[score_col]]) & is.finite(d$histology_burden) & apply(nuisance, 1L, function(z) all(is.finite(z)))
  if (sum(ok) < 8L) return(NA_real_)
  nuisance <- nuisance[ok, , drop = FALSE]
  x_rank <- rank(d[[score_col]][ok], ties.method = "average")
  y_rank <- rank(d$histology_burden[ok], ties.method = "average")
  batch <- droplevels(factor(d$acquisition_batch[ok]))
  nuisance_qr <- qr(nuisance)
  x_resid <- qr.resid(nuisance_qr, x_rank)
  y_resid <- qr.resid(nuisance_qr, y_rank)
  observed <- cor(x_resid, y_resid)
  fitted_reduced <- y_rank - y_resid
  by_batch <- split(seq_along(y_resid), batch)
  set.seed(seed)
  exceed <- 0L
  for (i in seq_len(n_perm)) {
    perm_index <- seq_along(y_resid)
    for (ii in by_batch) perm_index[ii] <- sample(ii, length(ii), replace = FALSE)
    y_star <- fitted_reduced + y_resid[perm_index]
    perm_stat <- cor(x_resid, qr.resid(nuisance_qr, y_star))
    if (is.finite(perm_stat) && abs(perm_stat) >= abs(observed)) exceed <- exceed + 1L
  }
  (exceed + 1) / (n_perm + 1)
}
bootstrap_rho <- function(d, score_col) {
  index_by_batch <- split(seq_len(nrow(d)), d$acquisition_batch)
  idx <- unlist(lapply(index_by_batch, function(ii) sample(ii, length(ii), replace = TRUE)), use.names = FALSE)
  b <- copy(d[idx])
  burden <- derive_burden(b)
  if (is.null(burden) || any(burden$loading <= 0)) return(NA_real_)
  b[, histology_burden := burden$score]
  unname(partial_spearman(b, score_col)["rho"])
}
run_family <- function(scores, perm_base, boot_base) {
  scores <- copy(scores)
  scores[, acquisition_batch := as.integer(acquisition_batch)]
  analysis_data <- merge(
    scores[group == "MASLD"],
    meta_masld[, .(sample_id, acquisition_batch, sex, age, bmi,
                   Steatosis, Ballooning, Inflammation, Fibrosis, histology_burden)],
    by = c("sample_id", "acquisition_batch"), all = FALSE, sort = FALSE)
  program_ids <- unique(analysis_data$program_id)
  rows <- rbindlist(mclapply(seq_along(program_ids), function(program_index) {
    pid <- program_ids[[program_index]]
    d <- analysis_data[program_id == pid]
    primary <- partial_spearman(d, "score")
    equal <- partial_spearman(d, "equal_score")
    leave_top <- partial_spearman(d, "leave_top_score")
    batch_rho <- vapply(c("2019", "2020"), function(batch_id) {
      unname(partial_spearman(d[acquisition_batch == batch_id], "score", include_batch = FALSE)["rho"])
    }, numeric(1))
    permutation_p <- freedman_lane_pvalue(d, "score", N_PERM, seed = perm_base + program_index)
    set.seed(boot_base + program_index)
    boot <- replicate(N_BOOT, bootstrap_rho(d, "score"))
    ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    data.table(program_id = pid, n = as.integer(primary["n"]), rho = unname(primary["rho"]),
               pvalue = permutation_p, pvalue_asymptotic = unname(primary["pvalue_asymptotic"]),
               ci_low = ci[1L], ci_high = ci[2L], n_boot = N_BOOT, n_boot_finite = sum(is.finite(boot)),
               n_permutations = N_PERM, equal_rho = unname(equal["rho"]), leave_top_rho = unname(leave_top["rho"]),
               rho_2019 = batch_rho[1L], rho_2020 = batch_rho[2L], perm_seed = perm_base + program_index,
               boot_seed = boot_base + program_index)
  }, mc.cores = N_CORES, mc.preschedule = FALSE))
  rows[, padj := p.adjust(pvalue, method = "BH")]
  rows[, sensitivity_sign_agree := is.finite(rho) & is.finite(equal_rho) & is.finite(leave_top_rho) &
         sign(rho) == sign(equal_rho) & sign(rho) == sign(leave_top_rho)]
  rows[, batch_sign_agree := is.finite(rho) & is.finite(rho_2019) & is.finite(rho_2020) &
         sign(rho) == sign(rho_2019) & sign(rho) == sign(rho_2020)]
  rows[, robust := is.finite(padj) & padj < 0.05 & sensitivity_sign_agree & batch_sign_agree]
  rows
}

# G-B: reproduce the stored 22-program column with the original seeds and stored scores
log_msg("G-B start")
gb <- run_family(stored_scores, 20260806L, 30260806L)
stored_b <- fread(file.path(PR, "module_protein_histology_burden_masld.tsv"))
cb <- merge(stored_b[testable == TRUE], gb, by = "program_id", suffixes = c(".s", ".n"))
gb_d <- max(abs(c(cb$rho.s - cb$rho.n, cb$pvalue.s - cb$pvalue.n, cb$ci_low.s - cb$ci_low.n,
                  cb$ci_high.s - cb$ci_high.n, cb$padj.s - cb$padj.n)))
gb_ok <- nrow(cb) == sum(stored_b$testable) && nrow(gb) == sum(stored_b$testable) && gb_d < 1e-10 &&
  all(cb$robust.s == cb$robust.n)
log_msg("G-B rows ", nrow(cb), " max|d| ", signif(gb_d, 3), " robust agree ", all(cb$robust.s == cb$robust.n))
fwrite(cb, file.path(OUT, "gateB_22_reproduction.tsv"), sep = "\t")
stopifnot(gb_ok)

# 117-program family
log_msg("117 family start")
fam <- run_family(v2_scores, 20260923L, 30260923L)
res <- merge(v2_reg, v2_cov, by = "program_uid", sort = FALSE)
res <- merge(res, fam, by.x = "program_uid", by.y = "program_id", all.x = TRUE, sort = FALSE)
res[, testable := is.finite(rho)]
res[is.na(robust), robust := FALSE]
res[, `:=`(correction_family = paste0("BH over ", sum(res$testable), " testable of 117 frozen programs"),
           inference_unit = "participant (PXD051911 MASLD, complete histology)",
           inference_method = "rank-residual partial Spearman vs histologic-burden PC1; batch-stratified Freedman-Lane permutation")]
setorder(res, cell_type, module)
fwrite(res, file.path(OUT, "protein_histology_117.tsv"), sep = "\t")
stopifnot(all(res[testable == TRUE, n_boot_finite] >= 0.95 * N_BOOT))
log_msg("117: testable ", sum(res$testable), "; q<0.05 ", sum(res$padj < 0.05, na.rm = TRUE),
        "; robust ", sum(res$robust))
print(res[order(padj)][1:15, .(program_uid, module_name, n, rho, ci_low, ci_high, pvalue, padj,
                               sensitivity_sign_agree, batch_sign_agree, robust)])
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_protein.txt"))
log_msg("done")
