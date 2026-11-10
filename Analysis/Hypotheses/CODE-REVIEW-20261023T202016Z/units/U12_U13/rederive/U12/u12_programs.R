#!/usr/bin/env Rscript
# U12 program-level ATAC check. Starts from the saved donor program scores
# (scores themselves are not re-derived here). Refits MASH vs normal per cohort
# with robust eBayes, trend = FALSE (corrected) and trend = TRUE (to reproduce the
# pre-correction values), BH within cohort over contrast-testable programs,
# and counts evidence states under the producer rule and under a rule where
# source_dependent also requires the three score variants to agree in sign.

suppressPackageStartupMessages({ library(data.table); library(limma) })
set.seed(1)  # no random draws below
args <- commandArgs(trailingOnly = TRUE); OUT <- args[1]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CAND <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates"
R1 <- file.path(CAND, "atac-context-v3-candidate-2026-08-11-r1")
CORR <- file.path(CAND, "atac-context-v3-statistics-correction-20261022")

scores <- fread(file.path(R1, "programs/program_donor_scores.tsv.gz"))
res_r1 <- fread(file.path(R1, "programs/program_atac_results.tsv"))
res_co <- fread(file.path(CORR, "programs/program_atac_results.tsv"))
elig <- rbindlist(lapply(c("GSE244832", "GSE281367"), function(co) {
  rbindlist(lapply(c("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"), function(lin) {
    d <- fread(file.path(R1, "counts", co, paste0(lin, ".donors.tsv")))
    d[condition %in% c("NORMAL", "MASH") & n_cells >= 20, .(cohort, lineage, donor_id, condition)]
  }))
}))

fit_one <- function(sc, trend) {
  out <- list()
  for (co in unique(sc$cohort)) for (lin in unique(sc[cohort == co, lineage])) {
    d <- elig[cohort == co & lineage == lin]
    if (sum(d$condition == "NORMAL") < 4 || sum(d$condition == "MASH") < 4) next
    s <- sc[cohort == co & lineage == lin & donor_id %in% d$donor_id]
    w <- dcast(s, program_uid ~ donor_id, value.var = "score")
    M <- as.matrix(w[, -1]); rownames(M) <- w$program_uid
    grp <- factor(d$condition[match(colnames(M), d$donor_id)], levels = c("NORMAL", "MASH"))
    X <- model.matrix(~ grp)
    f <- eBayes(lmFit(M, X), robust = TRUE, trend = trend)
    eff_eq <- dcast(s, program_uid ~ donor_id, value.var = "equal_weight_score")
    eff_lv <- dcast(s, program_uid ~ donor_id, value.var = "leave_top_gene_score")
    mdiff <- function(w2) {
      A <- as.matrix(w2[, -1]); rownames(A) <- w2$program_uid
      A <- A[rownames(M), colnames(M), drop = FALSE]
      rowMeans(A[, grp == "MASH", drop = FALSE]) - rowMeans(A[, grp == "NORMAL", drop = FALSE])
    }
    out[[paste(co, lin)]] <- data.table(
      cohort = co, lineage = lin, program_uid = rownames(M),
      n_normal = sum(grp == "NORMAL"), n_mash = sum(grp == "MASH"),
      effect = f$coefficients[, 2], se = f$stdev.unscaled[, 2] * sqrt(f$s2.post),
      p = f$p.value[, 2], eq_effect = mdiff(eff_eq), lv_effect = mdiff(eff_lv))
  }
  r <- rbindlist(out)
  r[, q := p.adjust(p, "BH"), by = cohort]
  r[, agree := sign(effect) == sign(eq_effect) & sign(effect) == sign(lv_effect)]
  r
}

states <- function(r, require_agree_for_source_dependent) {
  w <- dcast(r, program_uid ~ cohort, value.var = c("q", "effect", "agree"))
  for (nm in c("q_GSE244832", "q_GSE281367")) if (!nm %in% names(w)) w[, (nm) := NA_real_]
  s244 <- !is.na(w$q_GSE244832) & w$q_GSE244832 < 0.05
  s281 <- !is.na(w$q_GSE281367) & w$q_GSE281367 < 0.05
  a244 <- !is.na(w$agree_GSE244832) & w$agree_GSE244832
  a281 <- !is.na(w$agree_GSE281367) & w$agree_GSE281367
  both_testable <- !is.na(w$q_GSE244832) & !is.na(w$q_GSE281367)
  sd <- if (require_agree_for_source_dependent) (s244 & !s281 & a244) | (s281 & !s244 & a281) else xor(s244, s281)
  w[, state := fcase(
    both_testable & s244 & s281 & w$effect_GSE244832 * w$effect_GSE281367 > 0 & a244 & a281, "supported",
    both_testable & s244 & s281 & w$effect_GSE244832 * w$effect_GSE281367 < 0, "discordant",
    both_testable & sd, "source_dependent",
    !both_testable, "untestable",
    default = "indeterminate")]
  w
}

fit_corr <- fit_one(scores, trend = FALSE)
fit_old <- fit_one(scores, trend = TRUE)
cmp <- function(mine, prod, tag) {
  m <- merge(mine, prod[is.finite(pvalue), .(cohort, program_uid, pvalue, qvalue, standard_error, effect_prod = effect)],
             by = c("cohort", "program_uid"), all = TRUE)
  data.table(comparison = tag, n_mine = nrow(mine), n_prod = prod[is.finite(pvalue), .N], n_joined = nrow(m),
             family_244 = mine[cohort == "GSE244832", .N], family_281 = mine[cohort == "GSE281367", .N],
             max_abs_effect = m[, max(abs(effect - effect_prod), na.rm = TRUE)],
             max_abs_p = m[, max(abs(p - pvalue), na.rm = TRUE)],
             max_abs_q = m[, max(abs(q - qvalue), na.rm = TRUE)],
             max_abs_se = m[, max(abs(se - standard_error), na.rm = TRUE)],
             n_q_lt_0.05_mine = mine[q < 0.05, .N], n_q_lt_0.05_prod = prod[qvalue < 0.05, .N],
             min_q_mine_244 = mine[cohort == "GSE244832", min(q)], min_q_mine_281 = mine[cohort == "GSE281367", min(q)])
}
summary <- rbind(cmp(fit_corr, res_co, "trend=FALSE vs corrected 2026-10-22"),
                 cmp(fit_old, res_r1, "trend=TRUE vs r1 2026-08-11"),
                 cmp(fit_corr, res_r1, "trend=FALSE vs r1 2026-08-11"))
fwrite(summary, file.path(OUT, "program_compare.tsv"), sep = "\t")

st <- rbindlist(list(
  states(fit_corr, FALSE)[, .N, by = state][, rule := "corrected, producer rule"],
  states(fit_corr, TRUE)[, .N, by = state][, rule := "corrected, source_dependent requires sign agreement"],
  states(fit_old, FALSE)[, .N, by = state][, rule := "pre-correction, producer rule"],
  states(fit_old, TRUE)[, .N, by = state][, rule := "pre-correction, source_dependent requires sign agreement"]))
fwrite(st, file.path(OUT, "program_state_counts.tsv"), sep = "\t")
prod_states <- res_co[, .N, by = .(cohort, cross_cohort_state, within_source_state)]
fwrite(prod_states, file.path(OUT, "program_states_producer_corrected.tsv"), sep = "\t")
# Programs within 0.05 of the q = 0.05 line, either fit.
near <- merge(fit_corr[, .(cohort, program_uid, q_corrected = q, p_corrected = p)],
              fit_old[, .(cohort, program_uid, q_old = q, p_old = p)], by = c("cohort", "program_uid"))
fwrite(near[order(pmin(q_corrected, q_old))][1:20], file.path(OUT, "program_lowest_q.tsv"), sep = "\t")
print(summary); print(st); print(near[order(pmin(q_corrected, q_old))][1:10])
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_programs.txt"))
