# Corrected source labels; fixed scores, weights and 117-program family.
suppressPackageStartupMessages({library(data.table); library(nlme); library(parallel)})
setDTthreads(1)
ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
TASK <- file.path(ROOT, "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026")
source(file.path(ROOT, "Analysis/Hypotheses/H59-20261015T112432Z/scripts/h59_lib.R"))
MODE <- Sys.getenv("SINGLECELL_PROGRAM_MODE", "primary")
stopifnot(MODE %in% c("primary", "qc", "tcell20", "ambient"))
OUT <- Sys.getenv("SINGLECELL_PROGRAM_OUT", file.path(TASK, if (MODE == "primary") "programs" else paste0("programs_", MODE)))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(OUT, "results_primary_family.tsv"))) stop("Refusing to overwrite completed program results")
CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
DRAWS <- 19999L
set.seed(20260923L)
meta <- fread(file.path(TASK, "metadata/donor_metadata.tsv"))
stopifnot(!anyDuplicated(meta$donor))
scores_path <- file.path(ROOT, H, "donor_program_scores_primary.tsv")
scores <- fread(scores_path)
scores[, c("dataset", "disease_stage_coarse", "stage_ordinal", "exclude_stage_analysis") := NULL]
all_scores <- merge(scores, meta[, .(donor, dataset, disease_stage_coarse, exclude_stage_analysis)], by = "donor")
all_scores[, stage := match(disease_stage_coarse, STAGES) - 1L]
eligible <- all_scores[!exclude_stage_analysis & !is.na(stage)]
fit <- eligible[dataset %in% FIT_DATASETS]
eq <- fread(file.path(ROOT, H, "donor_program_scores_equal_run.tsv"))
eq[, c("dataset", "disease_stage_coarse", "stage_ordinal", "exclude_stage_analysis") := NULL]
eq <- merge(eq, meta[, .(donor, dataset, disease_stage_coarse, exclude_stage_analysis)], by = "donor")
eq[, stage := match(disease_stage_coarse, STAGES) - 1L]
eq <- eq[!exclude_stage_analysis & !is.na(stage) & dataset %in% FIT_DATASETS]
reg <- fread(inp("registry"))
programs <- unique(reg[, .(cell_type, module = as.integer(module), module_name, stability_median)])
setorder(programs, cell_type, module)
stopifnot(nrow(programs) == 117L, uniqueN(eligible$donor) == 61L,
          !"GSE202379_P98" %in% fit$donor,
          all(fit[donor == "GSE244832_D11", stage] == 1L),
          all(fit[donor == "GSE244832_D14", stage] == 2L))
fwrite(unique(fit[, .(donor, dataset, disease_stage_coarse, stage)]), file.path(OUT, "fitted_donor_labels.tsv"), sep = "\t")
fwrite(eligible, file.path(OUT, "donor_program_scores_primary.tsv"), sep = "\t")
AN7 <- file.path(ROOT, "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN7/20261021T124200Z/out")
if (MODE == "qc") {
  cached <- fread(file.path(AN7, "donor_scores_v2qc_cells.tsv"))
  fit <- merge(cached[, .(cell_type, module=as.integer(module), donor, score)],
               unique(fit[, .(donor, dataset, disease_stage_coarse, stage)]), by="donor")
}
if (MODE == "ambient") {
  cached <- fread(inp("ambient_scores"))
  fit <- merge(cached[is.finite(corrected), .(cell_type, module=as.integer(module), donor, score=corrected)],
               unique(fit[, .(donor, dataset, disease_stage_coarse, stage)]), by="donor")
  programs <- programs[cell_type != "tcells"]
  stopifnot(nrow(programs)==106L)
}
if (MODE == "tcell20") {
  cached <- fread(file.path(AN7, "donor_scores_all_cells.tsv"))
  fit <- merge(fit, cached[, .(cell_type, module=as.integer(module), donor, n_cells)], by=c("cell_type","module","donor"))
  fit <- fit[!(cell_type == "tcells" & n_cells < 20L)]
  programs <- programs[cell_type == "tcells"]
}

# Reproduce the adopted named effects before using the same model functions.
old_meta <- fread(inp("metadata"))
setnames(old_meta, "sample", "donor")
old_meta[, stage := match(disease_stage_coarse, STAGES) - 1L]
old_scores <- merge(scores, old_meta[!exclude_stage_analysis & !is.na(stage), .(donor, dataset, stage)], by = "donor")
old_results <- fread(file.path(ROOT, "Analysis/Hypotheses/H59-20261015T112432Z/out/outcome/results_primary_family.tsv"))
for (mod in c(8L, 20L)) {
  d <- prog_data(old_scores, "hepatocytes", mod, FIT_DATASETS)
  previous <- fit_gls(d$score, d$stage, d$dataset)
  expected <- old_results[cell_type == "hepatocytes" & module == mod]
  stopifnot(nrow(d) == expected$n_fit, abs(previous[["beta"]] - expected$gls_beta) < 1e-8)
}
donor_order <- sort(unique(eligible$donor), method = "radix")
keys <- matrix(runif(DRAWS * length(donor_order)), nrow = DRAWS)
one <- function(i) {
  ct <- programs$cell_type[i]; mod <- programs$module[i]
  cache <- file.path(OUT,paste0("checkpoint_",ct,"_",mod,".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  d <- prog_data(fit, ct, mod, FIT_DATASETS)
  observed <- fit_gls(d$score, d$stage, d$dataset)
  equal <- prog_data(eq, ct, mod, FIT_DATASETS)
  equal_fit <- fit_gls(equal$score, equal$stage, equal$dataset)
  permuted <- perm_stage(d, keys, donor_order)
  identity <- perm_stage(d, matrix(seq_along(donor_order), nrow=1L), donor_order)
  stopifnot(identical(as.numeric(identity[1,]), as.numeric(d$stage)))
  for (ds in unique(d$dataset)) {
    idx <- which(d$dataset == ds)
    stopifnot(all(rowSums(permuted[,idx,drop=FALSE]) == sum(d$stage[idx])))
  }
  null <- t(vapply(seq_len(DRAWS), function(b) fit_gls(d$score, permuted[b, ], d$dataset)[c("t", "t3")], c(t = 0, t3 = 0)))
  finite <- is.finite(null[, "t"])
  if (!all(finite)) stop("Incomplete permutation null: ", ct, " ", mod)
  stage_ok <- all(vapply(FIT_DATASETS, function(ds) {
    count <- table(d[dataset == ds, stage])
    sum(count) >= 8L && sum(count >= 2L) >= 2L
  }, TRUE))
  p <- (1 + sum(abs(null[, "t"]) >= abs(observed[["t"]]) * (1 - 1e-12))) / (DRAWS + 1)
  p3 <- if (all(is.finite(null[, "t3"]))) (1 + sum(abs(null[, "t3"]) >= abs(observed[["t3"]]) * (1 - 1e-12))) / (DRAWS + 1) else NA_real_
  rdf <- nrow(d) - 3L
  strict <- d[donor != "GSE202379_P30"]
  fs <- fit_gls(strict$score, strict$stage, strict$dataset)
  s1 <- within_slope(d, FIT_DATASETS[1]); s2 <- within_slope(d, FIT_DATASETS[2])
  row <- data.table(program = paste(ct, mod, sep = "_"), cell_type = ct, module = mod,
                    module_name = programs$module_name[i], stability_median = programs$stability_median[i],
                    n_fit = nrow(d), n_Healthy = sum(d$stage == 0), n_Steatosis = sum(d$stage == 1), n_Steatohepatitis = sum(d$stage == 2),
                    gls_beta = observed[["beta"]], gls_se = observed[["se"]], gls_t = observed[["t"]],
                    equal_run_gls_beta = equal_fit[["beta"]],
                    gls_ci_low = observed[["beta"]] - qt(.975, rdf) * observed[["se"]],
                    gls_ci_high = observed[["beta"]] + qt(.975, rdf) * observed[["se"]],
                    perm_p = if (stage_ok) p else 1, perm_p_t3 = if (stage_ok) p3 else 1,
                    gls_param_p = 2 * pt(abs(observed[["t"]]), rdf, lower.tail = FALSE),
                    null_critical_abs_t = unname(quantile(abs(null[, "t"]), 1 - .05 / 117)),
                    stage_count_eligible = stage_ok, permutation_draws = DRAWS,
                    sigma_GSE202379 = observed[["s1"]], sigma_GSE244832 = observed[["s2"]],
                    slope_GSE202379 = s1[["slope"]], se_GSE202379 = s1[["se"]],
                    slope_GSE244832 = s2[["slope"]], se_GSE244832 = s2[["se"]],
                    strict_n = nrow(strict), strict_beta = fs[["beta"]], strict_se = fs[["se"]],
                    strict_param_p = 2 * pt(abs(fs[["t"]]), nrow(strict) - 3, lower.tail = FALSE))
  fwrite(row, file.path(OUT, paste0("checkpoint_", ct, "_", mod, ".tsv")), sep = "\t")
  result <- list(row = row, null_t = null[, "t"], null_t3 = null[, "t3"])
  saveRDS(result,cache)
  result
}
results <- mclapply(seq_len(nrow(programs)), one, mc.cores = CORES, mc.preschedule = FALSE)
if (any(vapply(results, inherits, TRUE, "try-error"))) stop("A program fit failed; no family result written")
family <- rbindlist(lapply(results, `[[`, "row"))
family[, `:=`(perm_q = bh117(perm_p), perm_q_t3 = bh117(perm_p_t3), gls_param_q = bh117(gls_param_p), strict_param_q = bh117(strict_param_p))]
if (MODE == "tcell20") {
  combined <- fread(file.path(TASK,"programs/results_primary_family.tsv"))
  combined[match(family$program, program), perm_p := family$perm_p]
  combined[match(family$program, program), perm_p_t3 := family$perm_p_t3]
  combined[, `:=`(perm_q = bh117(perm_p), perm_q_t3 = bh117(perm_p_t3))]
  family[, `:=`(perm_q = combined$perm_q[match(program, combined$program)],
               perm_q_t3 = combined$perm_q_t3[match(program, combined$program)])]
  fwrite(combined[,.(program,perm_p,perm_p_t3,perm_q,perm_q_t3)],file.path(OUT,"full_117_primary_p_with_tcell20.tsv"),sep="\t")
}
family[, selected_H59 := perm_q < .05 & stability_median >= .5 & stage_count_eligible]
family[, robust_display_H59 := selected_H59 & is.finite(perm_q_t3) & perm_q_t3 < .05 &
         is.finite(equal_run_gls_beta) & sign(gls_beta) == sign(equal_run_gls_beta)]
family[, analysis_release := "source-diagnosis-correction-20261026"]
fwrite(family, file.path(OUT, "results_primary_family.tsv"), sep = "\t")
fwrite(family[cell_type == "hepatocytes" & module %in% c(8, 20)], file.path(OUT, "named_programs.tsv"), sep = "\t")
null_t <- do.call(cbind, lapply(results, `[[`, "null_t"))
null_t3 <- do.call(cbind, lapply(results, `[[`, "null_t3"))
colnames(null_t) <- colnames(null_t3) <- family$program
saveRDS(list(t = null_t, hc3_t = null_t3, donor_order = donor_order, seed = 20260923L), file.path(OUT, "permutation_nulls.rds"))
comparison <- merge(old_results[, .(program, previous_beta = gls_beta, previous_q = perm_q, previous_n = n_fit)],
                    family[, .(program, corrected_beta = gls_beta, corrected_q = perm_q, corrected_n = n_fit)], by = "program")
fwrite(comparison, file.path(OUT, "old_vs_corrected.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
writeLines(sprintf("%d/117 programs selected; %d permutations; BH over 117; seed 20260923.", sum(family$selected_H59), DRAWS), file.path(OUT, "COMPLETED"))
