#!/usr/bin/env Rscript
# CONF1 C1c: inflamed-hepatocyte bulk score against ballooning grade in two cohorts not used for it
# (GSE267145, E-MTAB-9815). Scoring and model follow DIG1 job 58b; cohort loading follows H8.
# Prespecification: Analysis/Hypotheses/CONF1-20261109T221846Z/PRESPEC_C1c.md
# Usage: Rscript c1c_ballooning_heldout.R <new_out_dir>   (the directory must not exist)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1]]
if (dir.exists(out)) stop("Output directory exists: ", out)
dir.create(out, recursive = TRUE)

SEED <- 20261109L
N_NULL <- 2000L
N_PERM <- 10000L
MIN_DETECT_FRAC <- 0.20
DISCOVERY_EFFECT <- 0.24

repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
kam <- file.path(repo, "RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-20260818T020627Z/external/kamzolas_v1/data")
inputs <- list(
  gse267145_counts = list(
    path = file.path(repo, "Analysis/MASLD_Model_Benchmark/executions/gse267145-matrix-audit-21064333/raw/GSE269412_RNA.txt.gz"),
    sha = "fed9333cb63ecca817c8dbc6f318c8b0795852df2567685bf509b90e3e1690f3"),
  gse267145_meta = list(
    path = file.path(repo, "Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv"),
    sha = "60bc52b86d7c7547c1371927513ef848c9062a7d3b6e17c032c04862d9d30b91"),
  ucam_counts = list(
    path = file.path(kam, "merged_counts.csv"),
    sha = "7e7b4a4616565d00d06ec42cf4849b47fe321af264b0bfc580e0f159f50894cc"),
  ucam_meta = list(
    path = file.path(kam, "metadata.csv"),
    sha = "4d940394042a5f74fb0246861d096f0290e2777b9460fc6f942f97cf7e79dbb9"),
  gene_metadata = list(
    path = file.path(repo, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"),
    sha = NA_character_)
)
sha256 <- function(p) sub(" .*", "", system2("sha256sum", shQuote(p), stdout = TRUE))
input_log <- rbindlist(lapply(names(inputs), function(k) {
  obs <- sha256(inputs[[k]]$path)
  if (!is.na(inputs[[k]]$sha) && obs != inputs[[k]]$sha) stop("Checksum mismatch: ", k)
  data.table(input = k, path = inputs[[k]]$path, sha256 = obs)
}))
fwrite(input_log, file.path(out, "inputs_sha256.tsv"), sep = "\t")

# Symbols on the primary assembly only; a set gene must map to exactly one Ensembl gene.
gene_meta <- fread(inputs$gene_metadata$path, select = c("gene_name", "chromosome", "ensembl_base"))
gene_meta <- unique(gene_meta[grepl("^chr", chromosome) & !grepl("_PAR_Y", ensembl_base)], by = "ensembl_base")

SETS <- list(
  INFLAMED_HEP = c("FCAMR", "TRIM31", "CXCL10", "IL32", "CCL20"),
  HEPONLY = c("FCAMR", "TRIM31"),
  IMMUNE = c("PTPRC", "CD3E", "CD2", "CD68", "CD163"),
  HEPMASS = c("ALB", "APOA1", "CYP2E1", "HP")
)
set_symbols <- unique(unlist(SETS))
sym2ens <- gene_meta[gene_name %in% set_symbols, .(n = .N, ens = ensembl_base[1]), by = gene_name]
stopifnot(all(set_symbols %in% sym2ens$gene_name), all(sym2ens$n == 1L))
ens_of <- setNames(sym2ens$ens, sym2ens$gene_name)

strip_version <- function(m) {
  ids <- rownames(m)
  m <- m[!grepl("_PAR_Y", ids), , drop = FALSE]
  base <- sub("\\..*$", "", rownames(m))
  dup <- base[duplicated(base)]
  keep <- !(base %in% dup)
  m <- m[keep, , drop = FALSE]
  rownames(m) <- base[keep]
  m
}

# ------------------------------------------------------------------ cohorts (H8 loader)
cohorts <- list()
pj <- fread(inputs$gse267145_meta$path)
x269 <- fread(inputs$gse267145_counts$path)
m269 <- as.matrix(x269[, -1]); rownames(m269) <- x269[[1]]; rm(x269)
stopifnot(all(pj$participant_id %in% colnames(m269)), !anyDuplicated(pj$participant_id))
cohorts$GSE267145 <- list(
  counts = strip_version(m269[, pj$participant_id, drop = FALSE]),
  meta = pj[, .(person = participant_id, ballooning = as.integer(ballooning),
                steatosis = as.integer(steatosis), fibrosis = as.integer(fibrosis),
                inflammation = as.integer(lobular_inflammation), sex = sex)])
rm(m269)

um <- fread(inputs$ucam_meta$path)[Dataset == "UCAM"]
uc <- fread(inputs$ucam_counts$path)
mu <- as.matrix(uc[, -1]); rownames(mu) <- uc[[1]]; rm(uc)
stopifnot(all(um[["Sample name"]] %in% colnames(mu)), !anyDuplicated(um[["Sample name"]]))
cohorts[["E-MTAB-9815"]] <- list(
  counts = strip_version(mu[, um[["Sample name"]], drop = FALSE]),
  meta = um[, .(person = `Sample name`, ballooning = as.integer(Ballooning),
                steatosis = as.integer(Steatosis), fibrosis = as.integer(Fibrosis),
                inflammation = as.integer(Inflammation), sex = Sex)])
rm(mu); invisible(gc())

cohort_order <- c("GSE267145", "E-MTAB-9815")
for (k in cohort_order) {
  m <- cohorts[[k]]$meta
  ok <- complete.cases(m[, .(ballooning, steatosis, fibrosis, inflammation, sex)])
  if (!all(ok)) message(k, ": dropping ", sum(!ok), " persons with missing covariates")
  cohorts[[k]]$meta <- m[ok]
  cohorts[[k]]$counts <- cohorts[[k]]$counts[, m$person[ok], drop = FALSE]
  cohorts[[k]]$meta[, sex := factor(sex, levels = c("F", "M"))]
  stopifnot(!anyDuplicated(cohorts[[k]]$meta$person), !anyNA(cohorts[[k]]$meta$sex),
            identical(colnames(cohorts[[k]]$counts), cohorts[[k]]$meta$person))
}

# ------------------------------------------------------------------ helpers
fit_hc3 <- function(f, dat) {
  mod <- lm(f, data = dat); X <- model.matrix(mod); e <- residuals(mod); h <- hatvalues(mod)
  XtXi <- solve(crossprod(X)); V <- XtXi %*% crossprod(X * (e / (1 - h))) %*% XtXi
  b <- coef(mod); se <- sqrt(diag(V))
  data.frame(term = names(b), estimate = b, se = se, t = b / se,
             p = 2 * pt(-abs(b / se), df.residual(mod)), n = nobs(mod), row.names = NULL)
}
# HC3 t statistic of one column of X for outcome y (used inside the null loops).
hc3_coef <- function(X, y, j) {
  XtXi <- solve(crossprod(X)); b <- XtXi %*% crossprod(X, y)
  e <- as.vector(y - X %*% b); h <- rowSums((X %*% XtXi) * X)
  V <- XtXi %*% crossprod(X * (e / (1 - h))) %*% XtXi
  c(estimate = b[j], t = b[j] / sqrt(V[j, j]))
}
rint <- function(x) qnorm((rank(x) - 0.5) / length(x))

comp <- "~ steatosis + inflammation + ballooning + fibrosis"
adj <- " + IMMUNE + HEPMASS"

gate_rows <- list(); model_rows <- list(); null_rows <- list(); grade_rows <- list()
set.seed(SEED)

for (k in cohort_order) {
  cts <- cohorts[[k]]$counts; meta <- as.data.frame(cohorts[[k]]$meta)
  cts <- cts[rowSums(cts) > 0, , drop = FALSE]
  d <- calcNormFactors(DGEList(cts))
  lc <- cpm(d, log = TRUE, prior.count = 2)
  usable <- rowMeans(cts > 0) >= MIN_DETECT_FRAC
  set_ens <- lapply(SETS, function(s) { e <- ens_of[s]; e[e %in% rownames(lc)[usable]] })
  n_use <- sapply(set_ens, length)
  hep_ok <- any(c("FCAMR", "TRIM31") %in% names(set_ens$INFLAMED_HEP))
  testable <- n_use[["INFLAMED_HEP"]] >= 4 && hep_ok && n_use[["IMMUNE"]] >= 3 && n_use[["HEPMASS"]] >= 3
  gate_rows[[k]] <- data.table(
    cohort = k, n_persons = nrow(meta), n_genes_usable = sum(usable),
    balloon_0 = sum(meta$ballooning == 0), balloon_1 = sum(meta$ballooning == 1), balloon_2 = sum(meta$ballooning == 2),
    inflamed_hep_usable = paste(names(set_ens$INFLAMED_HEP), collapse = ","),
    immune_usable = paste(names(set_ens$IMMUNE), collapse = ","),
    hepmass_usable = paste(names(set_ens$HEPMASS), collapse = ","), testable = testable)
  if (!testable) next

  z <- t(scale(t(lc[usable, , drop = FALSE])))
  for (s in names(SETS)) meta[[s]] <- colMeans(z[set_ens[[s]], , drop = FALSE])
  for (g in names(set_ens$INFLAMED_HEP)) meta[[paste0("g_", g)]] <- z[set_ens$INFLAMED_HEP[[g]], ]
  meta$INFLAMED_RINT <- rint(meta$INFLAMED_HEP)

  specs <- list(
    list(outcome = "INFLAMED_HEP", f = paste("INFLAMED_HEP", comp, adj), model = "primary"),
    list(outcome = "INFLAMED_HEP", f = paste("INFLAMED_HEP", comp), model = "unadjusted"),
    list(outcome = "INFLAMED_HEP", f = paste("INFLAMED_HEP", comp, adj, "+ sex"), model = "primary_plus_sex"),
    list(outcome = "INFLAMED_RINT", f = paste("INFLAMED_RINT", comp, adj), model = "primary"),
    list(outcome = "IMMUNE", f = paste("IMMUNE", comp), model = "unadjusted"))
  if (n_use[["HEPONLY"]] == 2)
    specs[[length(specs) + 1]] <- list(outcome = "HEPONLY", f = paste("HEPONLY", comp, adj), model = "primary")
  for (g in names(set_ens$INFLAMED_HEP))
    specs[[length(specs) + 1]] <- list(outcome = paste0("g_", g), f = paste0("g_", g, " ", comp, adj), model = "primary")
  for (sp in specs) {
    r <- fit_hc3(as.formula(sp$f), meta); r$cohort <- k; r$outcome <- sp$outcome; r$model <- sp$model
    model_rows[[length(model_rows) + 1]] <- r[r$term != "(Intercept)", ]
  }

  # Design matrix of the primary model; the ballooning column is tested.
  X <- model.matrix(as.formula(paste(comp, adj)), meta)
  jb <- which(colnames(X) == "ballooning")
  obs <- hc3_coef(X, meta$INFLAMED_HEP, jb)

  # Freedman-Lane permutation: permute residuals of the model without ballooning.
  X0 <- X[, -jb, drop = FALSE]
  b0 <- solve(crossprod(X0), crossprod(X0, meta$INFLAMED_HEP))
  fit0 <- as.vector(X0 %*% b0); res0 <- meta$INFLAMED_HEP - fit0
  t_perm <- replicate(N_PERM, hc3_coef(X, fit0 + sample(res0), jb)[["t"]])
  p_fl <- (1 + sum(abs(t_perm) >= abs(obs[["t"]]))) / (1 + N_PERM)

  # Matched-set null: one gene per score gene from the same expression decile.
  ave <- rowMeans(lc[usable, , drop = FALSE])
  dec <- cut(ave, quantile(ave, 0:10 / 10), include.lowest = TRUE, labels = FALSE)
  names(dec) <- rownames(z)
  pool <- setdiff(rownames(z), unlist(set_ens))
  target_dec <- dec[set_ens$INFLAMED_HEP]
  pools <- lapply(target_dec, function(dd) pool[dec[pool] == dd])
  null_est <- replicate(N_NULL, {
    g <- vapply(pools, function(p) sample(p, 1L), character(1))
    hc3_coef(X, colMeans(z[g, , drop = FALSE]), jb)[["estimate"]]
  })
  null_rows[[k]] <- data.table(
    cohort = k, observed = obs[["estimate"]], observed_t = obs[["t"]], freedman_lane_p = p_fl,
    null_mean = mean(null_est), null_q025 = quantile(null_est, 0.025), null_q975 = quantile(null_est, 0.975),
    frac_null_ge_observed = mean(null_est >= obs[["estimate"]]),
    matched_p_two_sided = (1 + sum(abs(null_est - mean(null_est)) >= abs(obs[["estimate"]] - mean(null_est)))) / (1 + N_NULL),
    n_null = N_NULL, n_perm = N_PERM, seed = SEED)

  grade_rows[[k]] <- as.data.table(meta)[, .(n = .N, mean_INFLAMED_HEP = mean(INFLAMED_HEP),
                                              mean_IMMUNE = mean(IMMUNE)), by = ballooning][order(ballooning)][, cohort := k]
}

gate <- rbindlist(gate_rows); fwrite(gate, file.path(out, "cohort_gate.tsv"), sep = "\t")
if (!any(gate$testable)) {
  writeLines("No cohort passed the gate; untestable.", file.path(out, "STOPPED.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  quit(save = "no", status = 0)
}
models <- rbindlist(model_rows); fwrite(models, file.path(out, "models.tsv"), sep = "\t")
fwrite(rbindlist(null_rows), file.path(out, "matched_null.tsv"), sep = "\t")
fwrite(rbindlist(grade_rows), file.path(out, "mean_by_grade.tsv"), sep = "\t")

# ------------------------------------------------------------------ primary test
pr <- models[outcome == "INFLAMED_HEP" & model == "primary" & term == "ballooning"]
pr[, holm_p := p.adjust(p, method = "holm")]
w <- 1 / pr$se^2
est <- sum(w * pr$estimate) / sum(w); se <- sqrt(1 / sum(w)); zc <- est / se
mde <- 2.8 * se
is_confirmed <- est > 0 && 2 * pnorm(-abs(zc)) < 0.05 && all(pr$estimate > 0)
decision <- "not_confirmed"
if (mde > 2 * DISCOVERY_EFFECT) decision <- "uninformative"
if (is_confirmed) decision <- "confirmed"
combined <- data.table(
  cohorts = paste(pr$cohort, collapse = ","), n_cohorts = nrow(pr),
  per_cohort_estimate = paste(signif(pr$estimate, 4), collapse = ","),
  per_cohort_se = paste(signif(pr$se, 4), collapse = ","),
  per_cohort_p = paste(signif(pr$p, 4), collapse = ","),
  per_cohort_holm_p = paste(signif(pr$holm_p, 4), collapse = ","),
  combined_estimate = est, combined_se = se, combined_z = zc, combined_p = 2 * pnorm(-abs(zc)),
  ci_low = est - 1.96 * se, ci_high = est + 1.96 * se,
  min_detectable_effect = mde, discovery_estimate = DISCOVERY_EFFECT, decision = decision)
fwrite(combined, file.path(out, "combined.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))

options(width = 200)
print(gate); print(pr); print(combined); print(rbindlist(null_rows))
