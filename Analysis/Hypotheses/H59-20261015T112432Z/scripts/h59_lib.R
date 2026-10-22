# H59 shared library: inputs, universes, model fits, permutation machinery.
# Implements PRESPEC_H59.md (copied in the run dir; sha256 verified at start).
suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
  library(parallel)
})
setDTthreads(1)
options(stringsAsFactors = FALSE, warn = 1)

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- Sys.getenv("H59_RUN_DIR", "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H59-20261015T112432Z")
PRESPEC_SHA <- "9f04681ecd855b0e07eabb591ab679233ae70255cd36a5427da46bfebeac21ef"
SEED <- 20260923L
B_FULL <- 19999L
B <- as.integer(Sys.getenv("H59_B", B_FULL))
NCORES <- as.integer(Sys.getenv("H59_CORES", "16"))
OUT <- Sys.getenv("H59_OUT", file.path(RUN, "out"))
N_FAMILY <- 117L
STAGES <- c("Healthy", "Steatosis", "Steatohepatitis")
FIT_DATASETS <- c("GSE202379", "GSE244832")
DROP_DATASET <- "GSE189600"
DROP_DONORS <- c("SRR22677587", "SRR22677595")
NAMED <- data.table(cell_type = c("hepatocytes", "hepatocytes"), module = c(20L, 8L),
                    label = c("Ductular injury", "Stromal ECM"),
                    printed_beta = c(0.402086, 0.351525))

H <- "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot"
FST <- "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014"
INPUTS <- data.table(
  role = c("scores_primary", "scores_equal_run", "metadata", "registry", "primary_effects",
           "ambient_effects", "ambient_scores", "verifier_60", "producer_512", "producer_selection",
           "fig4c_pdf"),
  path = c(file.path(H, "donor_program_scores_primary.tsv"),
           file.path(H, "donor_program_scores_equal_run.tsv"),
           "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv",
           file.path(FST, "program_registry_v2.tsv"),
           file.path(FST, "primary_effects.tsv"),
           file.path(FST, "ambient_program_effects.tsv"),
           "Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/results/donor_program_scores.tsv.gz",
           "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R04-01-reproduce/refit_drop_gse189600_family.tsv",
           "Analysis/SingleCell/scripts/hotspot_modules/512_hotspot_v2_donor_refit.R",
           "scripts/analysis/bulk_refit_followup/refit_qualified_program_diagnosis.R",
           "figures/main/fig4_singlecell_programs/panels/fig4c_hotspot_stage_heatmap.pdf"),
  sha256 = c("97fdad8d6841420fad92ab85495ae33f0ee34b0bda7da17137539f38e4b28579",
             "a5d23bd1b8ea0e897af6fe09a40e221816af0c0162c1d976fa5cafd1b1b2f410",
             "f44f4daa5502983f7f424f87addb4dc6ee748f10a114b5767801b2c245cb6414",
             "d1eae3c0d02532acaaed14ae772580bafea4741bab76ca5efc668d1d84b50bfb",
             "d512433ecd9e74ece936f24b7cc1d4566b0316fb55be2ca08603db33e18b23c3",
             "d3cc0115864006b51e4d690f64eef861ea5aa84904c20786b88c964f1b112e75",
             "72c7a2bf3e80df912f11b8908597da9d92290367a912dc30a624dbc3437059d3",
             "c8d727feec9f307feedd831eb1b39d7c16fad74acd7b6d3547b04279e59eff91",
             "fa67d13e4d28f523b2f2d119ce9805dddea92c83df75422750434388e0e59d88",
             "bd59de5b485c17d4de3db2afdacb0a49f170b476b7ceb5b9f348903e7775d3f9",
             "4cf739c2a87ec4156271e12af7cf9e08a3bf7a87c347c9ded8025780acfabe31"))

msg <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush.console() }
inp <- function(role) file.path(ROOT, INPUTS$path[INPUTS$role == role])
wtsv <- function(x, name) {
  path <- file.path(OUT, name)
  if (file.exists(path)) stop("refusing to overwrite ", path)
  fwrite(x, path, sep = "\t")
  invisible(path)
}
sha256 <- function(path) sub(" .*$", "", system2("sha256sum", shQuote(path), stdout = TRUE))

check_hashes <- function() {
  ps <- sha256(file.path(RUN, "PRESPEC_H59.md"))
  if (ps != PRESPEC_SHA) stop("PRESPEC_H59.md sha256 mismatch: ", ps)
  obs <- vapply(INPUTS$role, function(r) sha256(inp(r)), "")
  res <- copy(INPUTS)[, observed := obs][, match := observed == sha256]
  sz <- file.size(inp("scores_primary"))
  if (!all(res$match) || sz != 1869675) {
    print(res[match == FALSE])
    stop("INPUT HASH MISMATCH; refusing to run")
  }
  msg("input hashes OK (", nrow(res), " files + prespec); primary size ", sz)
  res
}

bh117 <- function(p) {
  p <- ifelse(is.finite(p), p, 1)
  p.adjust(p, method = "BH", n = N_FAMILY)
}

# ---------------------------------------------------------------- data
load_data <- function() {
  meta <- fread(inp("metadata"))
  setnames(meta, "sample", "donor")
  stopifnot(!anyDuplicated(meta$donor))
  meta[, excl := tolower(as.character(exclude_stage_analysis)) %in% c("true", "t", "1")]
  elig <- meta[!excl & disease_stage_coarse %in% STAGES,
               .(donor, dataset, disease_stage_coarse)]
  elig[, stage := match(disease_stage_coarse, STAGES) - 1]
  rd <- function(role) {
    x <- fread(inp(role))
    x[, c("dataset", "disease_stage_coarse", "stage_ordinal", "exclude_stage_analysis") := NULL]
    x <- merge(x[, .(cell_type, module = as.integer(module), donor, score)], elig, by = "donor")
    x
  }
  p <- rd("scores_primary"); eq <- rd("scores_equal_run")
  stopifnot(uniqueN(p$donor) == 62L, !any(p$donor %in% c("SRR14602057", "SRR14602058")))
  reg <- fread(inp("registry"))
  programs <- unique(reg[, .(cell_type, module = as.integer(module), module_name,
                             stability_median)])
  setorder(programs, cell_type, module)
  stopifnot(nrow(programs) == N_FAMILY)
  donors62 <- sort(unique(p$donor), method = "radix")
  stopifnot(all(DROP_DONORS %in% donors62),
            setequal(p[dataset == DROP_DATASET, unique(donor)], DROP_DONORS))
  list(meta = meta, elig = elig, p = p, eq = eq, reg = reg, programs = programs,
       donors62 = donors62)
}

prog_data <- function(scores, ct, mod, datasets = NULL) {
  d <- scores[cell_type == ct & module == mod & is.finite(score)]
  if (!is.null(datasets)) d <- d[dataset %in% datasets]
  setorderv(d, c("dataset", "donor"))  # data.table ordering is C-locale
  d[]
}

# ------------------------------------------------------------ permutation
key_matrix <- function(nrow_draws = B) {
  set.seed(SEED)
  K <- matrix(runif(B_FULL * 62), nrow = B_FULL)
  K[seq_len(nrow_draws), , drop = FALSE]
}

# Permuted stage vectors for one program: B x n matrix. Donors in d are sorted
# by dataset then donor; within a dataset the donor with the j-th smallest key
# receives s[j] where s is the stage vector in sorted donor-ID order.
perm_stage <- function(d, K, donors62) {
  cols <- match(d$donor, donors62)
  stopifnot(!anyNA(cols))
  out <- matrix(NA_real_, nrow(K), nrow(d))
  for (ds in unique(d$dataset)) {
    idx <- which(d$dataset == ds)
    s <- d$stage[idx]                      # already sorted donor-ID order within dataset
    Kd <- K[, cols[idx], drop = FALSE]
    ord <- t(apply(Kd, 1, order))           # ord[k, j] = position of j-th smallest key
    if (length(idx) == 1L) ord <- matrix(1L, nrow(K), 1L)
    pm <- matrix(NA_real_, nrow(K), length(idx))
    pm[cbind(rep(seq_len(nrow(K)), length(idx)), as.vector(ord))] <- rep(s, each = nrow(K))
    out[, idx] <- pm
  }
  out
}

# ---------------------------------------------------------------- GLS fit
GLS_CTRL <- glsControl(maxIter = 200, msMaxIter = 200)
fit_gls <- function(y, stage, dataset) {
  df <- data.frame(y = y, stage = stage, dataset = factor(dataset))
  f <- tryCatch(gls(y ~ stage + dataset, data = df,
                    weights = varIdent(form = ~ 1 | dataset),
                    method = "REML", control = GLS_CTRL),
                error = function(e) NULL)
  if (is.null(f)) return(c(beta = NA, se = NA, t = NA, t3 = NA, s1 = NA, s2 = NA))
  beta <- unname(coef(f)["stage"]); se <- sqrt(vcov(f)["stage", "stage"])
  vf <- coef(f$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE)
  lev <- levels(df$dataset)
  sd_d <- f$sigma * vf[lev]
  # Display robustness: WLS with weights 1/sigma_d^2, HC3 sandwich.
  w <- 1 / sd_d[as.integer(df$dataset)]^2
  X <- model.matrix(~ stage + dataset, df)
  XtWX_inv <- solve(crossprod(X * w, X))
  bw <- XtWX_inv %*% crossprod(X * w, y)
  e <- as.vector(y - X %*% bw)
  h <- rowSums((X %*% XtWX_inv) * X) * w
  t3 <- NA_real_
  if (all(h < 1 - 1e-10)) {
    meat <- crossprod(X * (w * e / (1 - h)), X * (w * e / (1 - h)))
    V3 <- XtWX_inv %*% meat %*% XtWX_inv
    t3 <- unname(bw[2] / sqrt(V3[2, 2]))
  }
  c(beta = beta, se = se, t = beta / se, t3 = t3, s1 = unname(sd_d[1]), s2 = unname(sd_d[2]))
}

# ------------------------------------------------------- OLS closed form (FWL)
# y: n vector, X: matrix of stage vectors (rows = draws), g: dataset labels
ols_t_mat <- function(y, Xs, g) {
  g <- as.character(g); n <- length(y); k <- length(unique(g)); rdf <- n - k - 1
  yc <- y - ave(y, g)
  XC <- Xs
  for (gg in unique(g)) { idx <- which(g == gg); XC[, idx] <- Xs[, idx, drop = FALSE] - rowMeans(Xs[, idx, drop = FALSE]) }
  sxx <- rowSums(XC^2)
  beta <- as.vector(XC %*% yc) / sxx
  E <- matrix(yc, nrow(XC), n, byrow = TRUE) - beta * XC
  se <- sqrt(rowSums(E^2) / rdf / sxx)
  ng <- as.vector(table(g)[g])
  Hm <- matrix(1 / ng, nrow(XC), n, byrow = TRUE) + XC^2 / sxx
  W <- ifelse(Hm < 1 - 1e-10, E^2 / (1 - Hm)^2, 0)
  se3 <- sqrt(rowSums(XC^2 * W)) / sxx
  t <- beta / se
  t[!is.finite(sxx) | sxx <= sqrt(.Machine$double.eps)] <- NA
  list(beta = beta, se = se, t = t, t3 = beta / se3, rdf = rdf)
}

# --------------------------------------------- per-dataset diagnostics (D2, G2)
cell_centred <- function(d) d$score - ave(d$score, d$dataset, d$stage)
sigma_cells <- function(d) {
  r <- cell_centred(d)
  sapply(FIT_DATASETS, function(ds) {
    i <- d$dataset == ds; ncell <- uniqueN(d$stage[i])
    sqrt(sum(r[i]^2) / (sum(i) - ncell))
  })
}
sxx_by <- function(d) sapply(FIT_DATASETS, function(ds) { s <- d$stage[d$dataset == ds]; sum((s - mean(s))^2) })
brown_forsythe <- function(x, g) {
  g <- factor(g); z <- abs(x - ave(x, g, FUN = median))
  if (nlevels(g) < 2) return(NA_real_)
  a <- anova(lm(z ~ g)); a[["Pr(>F)"]][1]
}
within_slope <- function(d, ds) {
  x <- d[dataset == ds]
  f <- summary(lm(score ~ stage, x))$coefficients
  c(slope = f["stage", 1], se = f["stage", 2], p = f["stage", 4], n = nrow(x))
}
