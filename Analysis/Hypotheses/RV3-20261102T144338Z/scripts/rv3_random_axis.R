# RV3: random-axis null for the continuum program coefficients. Follows PRESPEC_RV3.md.
# Unit: genetic-identity cluster (one per person).
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
REL <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-corrected774-20261011T152106Z")
DGE <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/merged_dge.rds")
MNC <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/prespec/mnc-fib-v1-20260828T133158Z")
SCORES <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261011T152106Z/programs/hotspot/hotspot_participant_scores.tsv.gz")
BASIS <- file.path(ROOT, "figures/candidates/continuum-composition-774-20261012T012137Z/basis_774/continuum")
FRAC <- file.path(ROOT, "Analysis/Deconvolution/corrected_774_20261009T141331Z/results")
COHORTS <- c("GSE162694", "GSE213621")
AXES <- c("fixed_projection", "signature_pc1")
FOCAL <- c(hotspot_hepatocytes_f05c535ae5bbc0b9 = "Stromal ECM", hotspot_hepatocytes_48f39dd4d817a10e = "Ductular injury")
wr <- function(x, f) fwrite(x, file.path(OUT, f), sep = "\t", na = "NA")
base_id <- function(x) sub("\\.[0-9]+$", "", gsub('"', "", as.character(x)))
say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
checks <- list()
check <- function(name, value, ok) {
  checks[[length(checks) + 1L]] <<- data.table(check = name, value = value, pass = ok)
  wr(rbindlist(checks), "repro_checks.tsv")
  if (!ok && !grepl("signature_pc1", name)) stop("reproduction check failed: ", name, " = ", value)
}

# ---------------- inputs ----------------
model <- readRDS(file.path(REL, "reproduction/discovery_signature_model.rds"))
loading <- model$discovery_loading_oriented; names(loading) <- base_id(names(loading))
centre <- model$discovery_center; names(centre) <- base_id(names(centre))
M <- readRDS(file.path(REL, "reproduction/discovery_normalized_expression.rds")); rownames(M) <- base_id(rownames(M))
rep <- fread(file.path(REL, "reproduction/participant_scores.tsv"))
released_pc1 <- as.numeric(rep[match(colnames(M), sample_id), pc1_released_operational])
stopifnot(all(is.finite(released_pc1)))
reg <- fread(file.path(REL, "projection/score_registry.tsv"))

dge <- readRDS(DGE)
counts <- dge$counts; rownames(counts) <- base_id(rownames(counts))
stopifnot(!anyDuplicated(rownames(counts)))
hist <- fread(file.path(BASIS, "participant_histology.tsv"))

normalise_cohort <- function(cnt) {
  keep <- rowSums(cnt) > ncol(cnt)
  le <- log2(cnt[keep, , drop = FALSE] + 1)
  out <- normalizeQuantiles(le); dimnames(out) <- dimnames(le); out
}
X <- setNames(lapply(COHORTS, function(cid) normalise_cohort(counts[, hist[dataset == cid, sample_id], drop = FALSE])), COHORTS)

project <- function(mat, ld, ct) {
  g <- names(ld)[names(ld) %in% rownames(mat)]
  setNames(as.numeric(crossprod(ld[g], sweep(mat[g, , drop = FALSE], 1L, ct[g], FUN = "-"))), colnames(mat))
}
cohort_pc1 <- function(mat, ld) {
  g <- names(ld)[names(ld) %in% rownames(mat)]
  cp <- prcomp(t(mat[g, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 1L)
  rl <- cp$rotation[, 1L]
  cs <- sum(rl * ld[g]) / sqrt(sum(rl^2) * sum(ld[g]^2))
  if (!is.finite(cs) || cs == 0) return(NULL)
  setNames(as.numeric(cp$x[, 1L]) * sign(cs), colnames(mat))
}

# ---------------- check 1 and 3: stored axes ----------------
common <- names(loading)
obs_axis <- list()
for (cid in COHORTS) {
  r <- reg[gsub('"', "", dataset) == cid]; r[, sample_id := gsub('"', "", sample_id)]
  fp <- project(X[[cid]], loading, centre)
  d1 <- max(abs(fp[r$sample_id] - r$fixed_projection_raw))
  check(paste0("fixed_projection_raw:", cid), sprintf("%.3g", d1), d1 < 1e-6)
  sp <- cohort_pc1(X[[cid]], loading)
  d3 <- max(abs(sp[r$sample_id] - r$signature_pc1_raw))
  check(paste0("signature_pc1_raw:", cid), sprintf("%.3g", d3), d3 < 1e-6)
  obs_axis[[cid]] <- list(fixed_projection = setNames(r$fixed_projection_raw, r$sample_id),
                          signature_pc1 = setNames(r$signature_pc1_raw, r$sample_id))
}
axes_used <- if (all(rbindlist(checks)[grepl("signature_pc1", check), pass])) AXES else "fixed_projection"

# ---------------- fitting set and program scores ----------------
frac_ok <- unlist(lapply(COHORTS, function(cid) {
  f <- fread(file.path(FRAC, cid, paste0(cid, "_bayesprism_proportions.tsv")))
  f[[1]][complete.cases(f)]
}))
fitset <- hist[dataset %in% COHORTS & source_control == FALSE & sample_id %in% frac_ok & !is.na(inferred_sex) & !is.na(fibrosis_stage)]
check("fit_n_GSE162694", nrow(fitset[dataset == "GSE162694"]), nrow(fitset[dataset == "GSE162694"]) == 109)
check("fit_n_GSE213621", nrow(fitset[dataset == "GSE213621"]), nrow(fitset[dataset == "GSE213621"]) == 258)

sc <- fread(SCORES)
W <- dcast(sc[sample_id %in% fitset$sample_id], sample_id ~ program_uid, value.var = "outcome_z")
programs <- setdiff(names(W), "sample_id")

prep_cohort <- function(cid) {
  d <- fitset[dataset == cid]
  Y <- as.matrix(W[match(d$sample_id, sample_id), ..programs])
  keep <- colSums(!is.finite(Y)) == 0
  list(d = d, Y = Y[, keep, drop = FALSE])
}
PC <- setNames(lapply(COHORTS, prep_cohort), COHORTS)

ols_hc3 <- function(Y, X) {
  XtXi <- solve(crossprod(X)); Mx <- XtXi %*% t(X)
  B <- Mx %*% Y; E <- Y - X %*% B
  h <- rowSums((X %*% XtXi) * X)
  se <- sqrt(colSums((Mx[2, ]^2) * (E / (1 - h))^2))
  list(beta = B[2, ], se = se)
}
fit_axis <- function(cid, axis_vals) {
  p <- PC[[cid]]; d <- p$d
  a <- axis_vals[d$sample_id]; a <- (a - mean(a)) / sd(a)
  X <- cbind(1, a, model.matrix(~ factor(stage_group) + factor(inferred_sex), d)[, -1, drop = FALSE])
  ols_hc3(p$Y, X)
}
meta_fit <- function(axis_by_cohort) {
  f <- lapply(COHORTS, function(cid) fit_axis(cid, axis_by_cohort[[cid]]))
  pr <- Reduce(intersect, lapply(f, function(x) names(x$beta)))
  b <- sapply(f, function(x) x$beta[pr]); s <- sapply(f, function(x) x$se[pr])
  w <- 1 / s^2
  beta <- rowSums(w * b) / rowSums(w); se <- sqrt(1 / rowSums(w))
  p <- 2 * pnorm(-abs(beta / se))
  q <- p.adjust(p, "BH", n = 117)
  data.table(program_uid = pr, beta = beta, se = se, p = p, q = q)
}

# ---------------- check 2: observed coefficients ----------------
stored <- fread(file.path(BASIS, "continuum_meta_effects.tsv"))
OBS <- list()
for (ax in axes_used) {
  o <- meta_fit(lapply(obs_axis, `[[`, ax)); o[, axis := ax]; OBS[[ax]] <- o
  for (u in names(FOCAL)) {
    ref <- stored[feature_id == u & axis == ax & arm == "disease_only" & composition == FALSE & robust == TRUE, beta]
    dd <- abs(o[program_uid == u, beta] - ref)
    check(paste0("meta_beta:", FOCAL[[u]], ":", ax), sprintf("%.3g", dd), dd < 1e-6)
  }
}
OBS <- rbindlist(OBS)
wr(OBS, "observed.tsv")
say("observed calls: %s", paste(OBS[, sum(q < 0.05), by = axis][, paste(axis, V1)], collapse = "; "))

# ---------------- null draws ----------------
read_sets <- function(f) { tb <- fread(f); s <- split(tb$gene_id_base, tb$draw_index); s[as.character(sort(as.integer(names(s))))] }
NULLSETS <- list(A = read_sets(file.path(MNC, "null_gene_sets.tsv.gz")), B = read_sets(file.path(MNC, "null_gene_sets_caliper.tsv.gz")))

focal_rows <- list(); call_rows <- list(); prog_null <- list()
for (nl in names(NULLSETS)) {
  sets <- NULLSETS[[nl]]
  for (b in seq_along(sets)) {
    gs <- intersect(sets[[b]], rownames(M))
    pca <- prcomp(t(M[gs, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 1L)
    sgn <- sign(cor(pca$x[, 1L], released_pc1))
    if (sgn == 0) next
    ld <- setNames(pca$rotation[, 1L] * sgn, gs); ct <- pca$center[gs]
    for (ax in axes_used) {
      av <- lapply(COHORTS, function(cid) if (ax == "fixed_projection") project(X[[cid]], ld, ct) else cohort_pc1(X[[cid]], ld))
      if (any(vapply(av, is.null, logical(1)))) next
      names(av) <- COHORTS
      n_present <- min(vapply(COHORTS, function(cid) sum(gs %in% rownames(X[[cid]])), integer(1)))
      r_obs <- median(vapply(COHORTS, function(cid) { s <- PC[[cid]]$d$sample_id; cor(av[[cid]][s], obs_axis[[cid]][[ax]][s]) }, numeric(1)))
      m <- meta_fit(av)
      focal_rows[[length(focal_rows) + 1L]] <- data.table(null = nl, draw = b, axis = ax, program_uid = names(FOCAL),
                                                          beta = m[match(names(FOCAL), program_uid), beta], n_genes_min = n_present, r_with_observed_axis = r_obs)
      call_rows[[length(call_rows) + 1L]] <- data.table(null = nl, draw = b, axis = ax, n_called = sum(m$q < 0.05),
                                                        median_beta = median(m$beta), r_with_observed_axis = r_obs)
      if (nl == "A") prog_null[[length(prog_null) + 1L]] <- m[, .(draw = b, axis = ax, program_uid, beta)]
    }
  }
  say("null %s done", nl)
}
NF <- rbindlist(focal_rows); NC <- rbindlist(call_rows); NP <- rbindlist(prog_null)
wr(NF, "null_focal.tsv.gz"); wr(NC, "null_calls.tsv.gz")

# ---------------- primary ----------------
prim <- merge(NF, OBS[program_uid %in% names(FOCAL), .(program_uid, axis, observed = beta)], by = c("program_uid", "axis"))
prim <- prim[, .(observed = observed[1], n_draws = .N, null_median = median(beta), null_q95 = quantile(beta, 0.95),
                 perm_p = (1 + sum(beta >= observed[1])) / (1 + .N)), by = .(null, program_uid, axis)]
prim[, holm_p := p.adjust(perm_p, "holm"), by = null]
prim[, program := FOCAL[program_uid]]
wr(prim, "primary.tsv")

exc <- merge(NP, OBS[, .(program_uid, axis, observed = beta, q_obs = q)], by = c("program_uid", "axis"))
exc <- exc[, .(observed = observed[1], q_obs = q_obs[1], null_median = median(beta), null_q95 = quantile(beta, 0.95),
               pct = mean(beta < observed[1])), by = .(program_uid, axis)]
exc[, excess := observed - null_median]
exc[, rank_excess := frank(-excess), by = axis]
wr(exc, "program_excess.tsv")

sink(file.path(OUT, "summary.txt"))
print(rbindlist(checks))
cat("\nPrimary (null A Holm over", nrow(prim[null == "A"]), ")\n"); print(prim[, .(null, program, axis, observed, null_median, null_q95, n_draws, perm_p, holm_p)])
cat("\nPrograms called (BH q < 0.05 within 117)\n")
print(merge(OBS[, .(observed_called = sum(q < 0.05)), by = axis],
            NC[, .(null_median_called = median(n_called), null_q05 = quantile(n_called, 0.05), null_q95 = quantile(n_called, 0.95),
                   median_r_with_observed_axis = median(r_with_observed_axis)), by = .(null, axis)], by = "axis"))
cat("\nFocal programs: excess over null median, rank among programs\n"); print(exc[program_uid %in% names(FOCAL)])
cat("\nTop 10 programs by excess (fixed_projection)\n"); print(exc[axis == "fixed_projection"][order(-excess)][1:10])
cat(sprintf("\nPrograms with observed above null q95 (fixed_projection): %d of %d\n", exc[axis == "fixed_projection", sum(observed > null_q95)], exc[axis == "fixed_projection", .N]))
sink()
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("done\n")
