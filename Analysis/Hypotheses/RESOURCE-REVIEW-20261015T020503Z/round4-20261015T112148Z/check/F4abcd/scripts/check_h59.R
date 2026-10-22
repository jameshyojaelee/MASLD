# Blind re-derivation of H59 (F4abcd check). Reads only release inputs
# (donor scores, source-diagnosis metadata, registry, ambient scores).
# MODE = prespec : permutation keys built as the prespec text describes
#                  (set.seed(20260923); runif(19999*62) by donor-ID order)
# MODE = indep   : independent Monte Carlo stream (set.seed(20260923); rnorm keys)
# ARM  = primary : primary scores, 117 programs (+ OLS 62/60, PC4 for named)
# ARM  = ambient : ambient-corrected scores (E4 arm)
suppressPackageStartupMessages({ library(data.table); library(nlme); library(parallel) })
setDTthreads(1)
MODE <- Sys.getenv("MODE", "prespec"); ARM <- Sys.getenv("ARM", "primary")
B <- as.integer(Sys.getenv("BDRAW", "19999")); NC <- as.integer(Sys.getenv("NCORES", "8"))
OUT <- Sys.getenv("OUTDIR"); stopifnot(nzchar(OUT)); dir.create(OUT, showWarnings = FALSE)
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
STAGES <- c("Healthy", "Steatosis", "Steatohepatitis")
NMIN <- if (ARM == "ambient") 8L else 6L
cat("MODE", MODE, "ARM", ARM, "B", B, "\n")

meta <- fread(paste0(R, "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"))
meta[, excl := toupper(as.character(exclude_stage_analysis)) %in% c("TRUE", "T", "1")]
u62 <- meta[!excl & disease_stage_coarse %in% STAGES, .(donor = sample, dataset, stage = match(disease_stage_coarse, STAGES) - 1)]
cat("62-universe:", nrow(u62), "\n"); print(u62[, .N, by = .(dataset, stage)][order(dataset, stage)])
# 60-universe: drop datasets whose ambient correction failed (GSE189600 is the only one in u62)
u60 <- u62[dataset != "GSE189600"]
nst <- u60[, .(k = uniqueN(stage)), by = dataset]
fitds <- nst[k >= 2, dataset]; cat("fit datasets:", fitds, "\n")
ufit <- u60[dataset %in% fitds]

reg <- fread(paste0(R, "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv"))
prog <- unique(reg[, .(cell_type, module = as.integer(module), module_name, stability_median)])
setorder(prog, cell_type, module); stopifnot(nrow(prog) == 117)
prog[, pid := paste(cell_type, module, sep = "_")]

if (ARM == "primary") {
  sc <- fread(paste0(R, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv"))
  sc <- sc[, .(cell_type, module = as.integer(module), donor, score)]
} else {
  sc <- fread(paste0(R, "Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/results/donor_program_scores.tsv.gz"))
  sc <- sc[, .(cell_type, module = as.integer(module), donor, score = corrected)]
}
sc <- sc[is.finite(score)]

getd <- function(univ, ct, md) {
  d <- merge(sc[cell_type == ct & module == md], univ, by = "donor")
  d[order(dataset, donor)]   # C-locale order in data.table
}

gls_t <- function(y, st, ds) {
  f <- tryCatch(gls(y ~ st + ds, weights = varIdent(form = ~ 1 | ds), method = "REML",
                    control = glsControl(maxIter = 200, msMaxIter = 200)), error = function(e) NULL)
  if (is.null(f)) return(c(NA, NA))
  b <- unname(coef(f)[2]); se <- sqrt(vcov(f)[2, 2]); c(b, se)
}
ols_fit <- function(d) {
  X <- model.matrix(~ stage + factor(dataset), d); y <- d$score
  if (qr(X)$rank < ncol(X)) return(c(beta = NA, p = NA, p3 = NA))
  XtXi <- solve(crossprod(X)); b <- XtXi %*% crossprod(X, y); e <- as.vector(y - X %*% b)
  df <- nrow(X) - ncol(X); s2 <- sum(e^2) / df; se <- sqrt(s2 * XtXi[2, 2])
  h <- rowSums((X %*% XtXi) * X)
  if (any(h > 1 - 1e-10)) p3 <- NA else {
    V3 <- XtXi %*% crossprod(X * (e / (1 - h))) %*% XtXi; p3 <- 2 * pt(abs(b[2] / sqrt(V3[2, 2])), df, lower.tail = FALSE) }
  c(beta = b[2], p = 2 * pt(abs(b[2] / se), df, lower.tail = FALSE), p3 = unname(p3))
}

# ---------------- observed
obs <- rbindlist(lapply(seq_len(nrow(prog)), function(i) {
  d <- getd(ufit, prog$cell_type[i], prog$module[i])
  n <- nrow(d)
  g <- if (n >= NMIN && uniqueN(d$dataset) == 2) gls_t(d$score, d$stage, factor(d$dataset)) else c(NA, NA)
  sl <- lapply(fitds, function(ds) { x <- d[dataset == ds]; if (nrow(x) < 3 || uniqueN(x$stage) < 2) return(c(NA, NA))
    s <- summary(lm(score ~ stage, x))$coefficients; c(s[2, 1], s[2, 2]) })
  Q <- (sl[[1]][1] - sl[[2]][1])^2 / (sl[[1]][2]^2 + sl[[2]][2]^2)
  out <- data.table(pid = prog$pid[i], n_fit = n, n_ds1 = sum(d$dataset == fitds[1]), n_ds2 = sum(d$dataset == fitds[2]),
             beta = g[1], se = g[2], t = g[1] / g[2],
             ci_lo = g[1] - qt(.975, n - 3) * g[2], ci_hi = g[1] + qt(.975, n - 3) * g[2],
             slope1 = sl[[1]][1], se1 = sl[[1]][2], slope2 = sl[[2]][1], se2 = sl[[2]][2], Qp = pchisq(Q, 1, lower.tail = FALSE))
  if (ARM == "primary") {
    o62 <- ols_fit(getd(u62, prog$cell_type[i], prog$module[i])); o60 <- ols_fit(getd(u60, prog$cell_type[i], prog$module[i]))
    out[, `:=`(n62 = nrow(getd(u62, prog$cell_type[i], prog$module[i])), ols62_beta = o62[1], ols62_p = o62[2], ols62_p3 = o62[3],
               ols60_beta = o60[1], ols60_p = o60[2], ols60_p3 = o60[3])]
  }
  out
}))
obs <- cbind(prog, obs[, -1])

# ---------------- permutation keys
donors62 <- sort(u62$donor, method = "radix")
set.seed(20260923)
K <- if (MODE == "prespec") matrix(runif(19999 * 62), nrow = 19999) else matrix(rnorm(19999 * 62), nrow = 19999)
K <- K[seq_len(B), , drop = FALSE]
permmat <- function(d) {
  cols <- match(d$donor, donors62); P <- matrix(NA_real_, B, nrow(d))
  for (ds in unique(d$dataset)) {
    idx <- which(d$dataset == ds); s <- d$stage[idx]   # donor-ID order within dataset
    for (k in seq_len(B)) { o <- order(K[k, cols[idx]]); v <- numeric(length(idx)); v[o] <- s; P[k, idx] <- v }
  }
  P
}
nullT <- mclapply(seq_len(nrow(prog)), function(i) {
  d <- getd(ufit, prog$cell_type[i], prog$module[i])
  if (nrow(d) < NMIN || uniqueN(d$dataset) < 2) return(rep(NA_real_, B))
  P <- permmat(d); ds <- factor(d$dataset)
  vapply(seq_len(B), function(k) { g <- gls_t(d$score, P[k, ], ds); g[1] / g[2] }, 0)
}, mc.cores = NC)
Tn <- do.call(cbind, nullT)
fwrite(as.data.table(Tn), file.path(OUT, "null_t.tsv.gz"))
pp <- vapply(seq_len(nrow(prog)), function(j) { if (!is.finite(obs$t[j])) return(NA_real_)
  v <- abs(Tn[, j]); v <- v[is.finite(v)]; (1 + sum(v >= abs(obs$t[j]) * (1 - 1e-12))) / (1 + length(v)) }, 0)
obs[, perm_p := pp]
obs[, perm_q := p.adjust(ifelse(is.finite(perm_p), perm_p, 1), "BH", n = 117)]
obs[, n_fail := colSums(!is.finite(Tn))]
obs[, crit := apply(abs(Tn), 2, function(v) if (all(!is.finite(v))) NA_real_ else quantile(v, 1 - 0.05 / 117, na.rm = TRUE, type = 7))]
obs[, fw_lo := beta - crit * se][, fw_hi := beta + crit * se]
obs[, selected := perm_q < 0.05 & stability_median >= 0.5]
if (ARM == "primary") {
  obs[, `:=`(ols62_q = p.adjust(ifelse(is.finite(ols62_p), ols62_p, 1), "BH", n = 117),
             ols62_q3 = p.adjust(ifelse(is.finite(ols62_p3), ols62_p3, 1), "BH", n = 117),
             ols60_q = p.adjust(ifelse(is.finite(ols60_p), ols60_p, 1), "BH", n = 117))]
  # PC4 detection for the two named programs, rows 1..1000 of K
  pc4 <- rbindlist(lapply(c("hepatocytes_20", "hepatocytes_8"), function(id) {
    i <- match(id, prog$pid); d <- getd(ufit, prog$cell_type[i], prog$module[i]); P <- permmat(d)[1:min(1000, B), , drop = FALSE]
    ds <- factor(d$dataset)
    rbindlist(lapply(c(0.2, 0.3, if (id == "hepatocytes_20") 0.402086 else 0.351525), function(dl) {
      tt <- vapply(seq_len(nrow(P)), function(k) { g <- gls_t(d$score + dl * P[k, ], P[k, ], ds); g[1] / g[2] }, 0)
      data.table(pid = id, delta = dl, crit = obs[pid == id, crit], detection = mean(abs(tt) > obs[pid == id, crit], na.rm = TRUE))
    }))
  }))
  fwrite(pc4, file.path(OUT, "pc4_check.tsv"), sep = "\t")
}
fwrite(obs, file.path(OUT, paste0("check_", ARM, "_", MODE, ".tsv")), sep = "\t")
cat("selected:", sum(obs$selected, na.rm = TRUE), " min q:", min(obs$perm_q), obs$pid[which.min(obs$perm_q)], "\n")
print(obs[pid %in% c("hepatocytes_20", "hepatocytes_8")])
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
