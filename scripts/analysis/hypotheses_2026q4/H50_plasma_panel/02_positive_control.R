#!/usr/bin/env Rscript
# H50 positive control (data adequacy): does the published 5-protein panel (C7, ALDOB, ICAM1, CTSD, IGF2;
# Boel 2025, PMID 40032974) discriminate Kleiner F2 or higher in PXD051911 initial-visit plasma?
# Modality: DIA-MS log2 intensities from the deposited filtered, batch-corrected matrix; logistic regression.
# Unit: the patient (one initial-visit plasma sample each). Not independent: the panel was derived in this cohort.
# This script uses no Ductular-program protein.
# Usage: Rscript 02_positive_control.R <new_output_dir>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
out <- args[[1]]
if (dir.exists(out)) stop("refusing to overwrite ", out)
dir.create(out, recursive = TRUE)
SEED <- 20260923L
N_REPEAT <- 20L; N_FOLD <- 5L; N_BOOT <- 2000L
PX <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/PXD051911"
PANEL <- c(C7 = "P10643", ALDOB = "P05062", ICAM1 = "P05362", CTSD = "P07339", IGF2 = "P01344")

meta <- fread(file.path(PX, "meta_data.txt"), na.strings = c("NA", ""))
filt <- fread(file.path(PX, "plasma_protein_quant_filtered_and_batch_corrected.txt"), na.strings = c("NA", "", "NaN"))
raw  <- fread(file.path(PX, "plasma_protein_quant.txt"), select = c("ProteinAccessions", "Genes"))
# the accession for each panel gene must be the only protein group carrying that gene symbol
for (g in names(PANEL)) {
  hit <- raw[vapply(strsplit(Genes, ";", fixed = TRUE), function(x) g %in% x, logical(1)), ProteinAccessions]
  stopifnot(identical(hit, unname(PANEL[g])))
}
stopifnot(all(PANEL %in% filt$ProteinAccessions))

s1 <- meta[sample_group == "initial_sample" & !is.na(plasma_proteomics_filename) &
           !is.na(kleiner_fibrosis_grade) & unique_identifier %in% names(filt)]
stopifnot(!anyDuplicated(s1$patient_name))
s1[, y := as.integer(kleiner_fibrosis_grade %in% c("F2", "F3", "F4"))]
s1[, `:=`(age = as.numeric(alder), sex = factor(gender), bmi_v = as.numeric(bmi), batch = factor(plasma_batch_effects))]
X <- t(as.matrix(filt[match(PANEL, ProteinAccessions), s1$unique_identifier, with = FALSE]))
colnames(X) <- names(PANEL)
dat <- cbind(s1[, .(patient_name, y, age, sex, bmi_v, batch)], as.data.table(X))
n_before <- nrow(dat)
dat <- dat[complete.cases(dat[, .(age, sex, bmi_v, batch)])]
cat(sprintf("S1 patients %d, complete covariates %d, F2+ %d, <F2 %d\n", n_before, nrow(dat), sum(dat$y), sum(dat$y == 0)))
panel_cols <- names(PANEL)
n_missing_panel <- colSums(is.na(dat[, ..panel_cols]))
for (g in names(PANEL)) dat[[g]] <- as.numeric(scale(dat[[g]]))   # unsupervised z-score

auc <- function(y, s) {
  r <- rank(s); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
forms <- list(
  covariates = y ~ age + sex + bmi_v + batch,
  panel      = y ~ C7 + ALDOB + ICAM1 + CTSD + IGF2,
  panel_plus_covariates = y ~ C7 + ALDOB + ICAM1 + CTSD + IGF2 + age + sex + bmi_v + batch)
impute <- function(tr, te) {                      # training-fold means for any missing panel value
  for (g in names(PANEL)) { m <- mean(tr[[g]], na.rm = TRUE); tr[[g]][is.na(tr[[g]])] <- m; te[[g]][is.na(te[[g]])] <- m }
  list(tr = tr, te = te)
}
oof <- lapply(forms, function(f) matrix(NA_real_, nrow(dat), N_REPEAT))
for (r in seq_len(N_REPEAT)) {
  set.seed(SEED + r)
  fold <- integer(nrow(dat))
  for (k in 0:1) { i <- which(dat$y == k); fold[i] <- sample(rep_len(seq_len(N_FOLD), length(i))) }
  for (f in seq_len(N_FOLD)) {
    tr <- dat[fold != f]; te <- dat[fold == f]
    ii <- impute(tr, te)
    for (nm in names(forms)) {
      fit <- suppressWarnings(glm(forms[[nm]], data = ii$tr, family = binomial()))
      oof[[nm]][fold == f, r] <- suppressWarnings(predict(fit, newdata = ii$te, type = "link"))
    }
  }
}
cv_auc <- function(P, idx = seq_len(nrow(P)), y = dat$y) mean(vapply(seq_len(ncol(P)), function(r) auc(y[idx], P[idx, r]), numeric(1)))
obs <- vapply(oof, cv_auc, numeric(1))
set.seed(SEED)
boot <- replicate(N_BOOT, {
  idx <- sample.int(nrow(dat), replace = TRUE)
  if (length(unique(dat$y[idx])) < 2) return(rep(NA_real_, length(forms)))
  vapply(oof, function(P) cv_auc(P, idx, dat$y), numeric(1))
})
rownames(boot) <- names(forms)
ci <- function(v) quantile(v, c(0.025, 0.975), na.rm = TRUE)
d_boot <- boot["panel_plus_covariates", ] - boot["covariates", ]
res <- data.table(model = names(forms), cv_auc = round(obs, 4),
                  ci_lo = round(apply(boot, 1, function(v) ci(v)[1]), 4), ci_hi = round(apply(boot, 1, function(v) ci(v)[2]), 4))
res <- rbind(res, data.table(model = "panel_plus_covariates minus covariates", cv_auc = round(obs[["panel_plus_covariates"]] - obs[["covariates"]], 4),
                             ci_lo = round(ci(d_boot)[1], 4), ci_hi = round(ci(d_boot)[2], 4)))
dat_c <- dat[complete.cases(dat[, ..panel_cols])]          # likelihood-ratio test on patients with all five proteins
fit0 <- glm(forms$covariates, data = dat_c, family = binomial()); fit1 <- glm(forms$panel_plus_covariates, data = dat_c, family = binomial())
lrt <- anova(fit0, fit1, test = "Chisq")
pass <- obs[["panel"]] > 0.7 && obs[["panel_plus_covariates"]] > obs[["covariates"]]
summ <- list(n_patients = nrow(dat), n_f2_plus = sum(dat$y), n_lt_f2 = sum(dat$y == 0), unit = "patient",
             n_dropped_incomplete_covariates = n_before - nrow(dat), panel_missing_values = as.list(n_missing_panel),
             repeats = N_REPEAT, folds = N_FOLD, bootstrap = N_BOOT, seed = SEED,
             lrt_n_patients = nrow(dat_c), lrt_panel_5df_chisq = lrt$Deviance[2], lrt_df = lrt$Df[2], lrt_p = lrt$`Pr(>Chi)`[2],
             positive_control_pass = pass, note = "panel derived in this cohort; not independent")
fwrite(res, file.path(out, "positive_control_cv_auc.tsv"), sep = "\t")
writeLines(jsonlite::toJSON(summ, auto_unbox = TRUE, pretty = TRUE, na = "null"), file.path(out, "positive_control_summary.json"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(res); str(summ)
