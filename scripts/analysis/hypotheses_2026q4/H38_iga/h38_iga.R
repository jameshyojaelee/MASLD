#!/usr/bin/env Rscript
# H38: IgA fraction of heavy-chain constant-region reads in bulk liver RNA.
# Prespecification: Analysis/Hypotheses/H38-20261009T180100Z/PRESPEC_H38.md (written before any outcome was read).
# Usage: h38_iga.R <new_output_dir>     (H38_SMOKE=1 runs a reduced bootstrap and permutation count)
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(metafor); library(aod); library(parallel); library(jsonlite)
})
SEED <- 20260923L
set.seed(SEED)
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 1L) stop("Usage: h38_iga.R output_dir")
out <- argv[[1]]
if (dir.exists(out)) stop("Refusing to overwrite existing output dir: ", out)
dir.create(out, recursive = TRUE)
SMOKE <- nzchar(Sys.getenv("H38_SMOKE"))
B_BOOT <- if (SMOKE) 8L else 500L
N_PERM <- if (SMOKE) 20L else 1000L
NCORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
MIN_N <- 10; MIN_N_SENS <- 30
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
READS <- file.path(ROOT, "Analysis/Hypotheses/READS-20261009T141017Z")
CAND <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z")
MANIFEST5 <- file.path(ROOT, "figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/five_cohort_sample_manifest.tsv")
GENEMETA <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
ROSTER_COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
STAGE_COHORTS <- c("GSE130970", "GSE135251", "GSE162694", "GSE213621")
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n", sep = " ")
writeLines(c(capture.output(sessionInfo()), "", paste("smoke:", SMOKE), paste("cores:", NCORES)), file.path(out, "sessionInfo.txt"))

# ------------------------------------------------------------------ library-level heavy-chain counts
cls_of_set <- function(gs) {
  vapply(strsplit(gs, ",", fixed = TRUE), function(x) {
    cl <- c(A = any(x %in% c("IGHA1", "IGHA2")), G = any(x %in% paste0("IGHG", 1:4)), M = any(x == "IGHM"),
            D = any(x == "IGHD"), E = any(x == "IGHE"))
    if (sum(cl) == 1L) names(cl)[cl] else "X"
  }, character(1))
}
lib_dirs <- Sys.glob(file.path(READS, "libraries", "*", "*"))
log_msg("libraries found:", length(lib_dirs))
read_one <- function(d) {
  lib <- basename(d); coh <- basename(dirname(d))
  f <- fread(file.path(d, "ig_fragments.tsv"), colClasses = list(character = c("library", "gene_set", "mapping")))
  q <- fread(file.path(d, "library_qc.tsv"), select = c("library", "star_input_reads", "layout"))
  list(frag = f[, `:=`(cohort = coh)], qc = q[, cohort := coh])
}
res <- mclapply(lib_dirs, read_one, mc.cores = NCORES)
frag <- rbindlist(lapply(res, `[[`, "frag"))
qc <- rbindlist(lapply(res, `[[`, "qc"))
setnames(frag, "all_alignments_in_constant_exons", "inside")
frag[, cls := cls_of_set(gene_set)]
frag_in <- frag[inside == 1L]
cnt <- dcast(frag_in[, .(n = sum(n_fragments)), by = .(library, cls)], library ~ cls, value.var = "n", fill = 0)
for (cc in c("A", "G", "M", "D", "E", "X")) if (!cc %in% names(cnt)) cnt[, (cc) := 0]
setnames(cnt, c("A", "G", "M", "D", "E", "X"), c("n_A", "n_G", "n_M", "n_D", "n_E", "n_mixed"))
lib <- merge(qc, cnt, by = "library", all.x = TRUE)
for (cc in c("n_A", "n_G", "n_M", "n_D", "n_E", "n_mixed")) lib[is.na(get(cc)), (cc) := 0]
lib[, n_outside := frag[inside == 0L, sum(n_fragments), by = library][match(lib$library, library), V1]]
lib[is.na(n_outside), n_outside := 0]
lib[, `:=`(k = n_A, n = n_A + n_G + n_M)]
# proportional allocation of class-mixed sets (sensitivity 2)
mixed <- frag_in[cls == "X"]
alloc <- mixed[, {
  g <- strsplit(gene_set, ",", fixed = TRUE)[[1]]
  .(hasA = any(g %in% c("IGHA1", "IGHA2")), hasG = any(g %in% paste0("IGHG", 1:4)), hasM = any(g == "IGHM"), n_fragments = n_fragments)
}, by = .(library, gene_set)]
alloc <- merge(alloc, lib[, .(library, n_A, n_G, n_M)], by = "library")
alloc[, den := hasA * n_A + hasG * n_G + hasM * n_M]
alloc <- alloc[den > 0 & (hasA + hasG + hasM) >= 1]
alloc[, `:=`(add_A = n_fragments * hasA * n_A / den, add_all = n_fragments)]
al <- alloc[, .(add_A = sum(add_A), add_all = sum(add_all)), by = library]
lib <- merge(lib, al, by = "library", all.x = TRUE)
lib[is.na(add_A), `:=`(add_A = 0, add_all = 0)]
lib[, `:=`(k2 = k + add_A, n2 = n + add_all)]
lib[, ig_rate := log10((n + 1) / (star_input_reads / 1e6))]
fwrite(lib, file.path(out, "class_counts_by_library.tsv"), sep = "\t")

# ------------------------------------------------------------------ persons and covariates
roster <- fread(file.path(CAND, "primary_unit_roster.tsv"))
man5 <- fread(MANIFEST5)
stopifnot(nrow(roster) == 774L, all(roster$sample_id %in% man5$sample_id), anyDuplicated(roster$analysis_unit_id) == 0L)
per <- merge(man5[, .(sample_id, analysis_unit_id, cohort = dataset, group_binary, inferred_sex, fibrosis_stage, nas_score)],
             lib, by.x = c("sample_id", "cohort"), by.y = c("library", "cohort"), all.x = TRUE)
stopifnot(nrow(per) == 774L, !anyNA(per$n))
per[, sexM := as.integer(inferred_sex == "M")]
per[, disease := as.integer(group_binary == "Disease")]
per[, stage := suppressWarnings(as.numeric(fibrosis_stage))]
per[, nas := suppressWarnings(as.numeric(nas_score))]

# expression covariates from the corrected fragment counts (TMM log-CPM), z within cohort
dge <- readRDS(file.path(CAND, "primary_source_controls/merged_dge.rds"))
counts <- if (inherits(dge, "DGEList")) dge$counts else as.matrix(dge)
gm <- fread(GENEMETA)
if (any(grepl("^ENSG", rownames(counts)))) {
  pick <- function(sym) {
    g <- gm[gene_name == sym]
    g <- g[order(gene_biotype != "protein_coding")]
    stopifnot(nrow(g) >= 1L); g$gene_id[[1]]
  }
  ids <- vapply(c("CD79A", "MS4A1", "MZB1", "TNFRSF17"), pick, character(1))
  stopifnot(all(ids %in% rownames(counts)))
} else {
  ids <- setNames(c("CD79A", "MS4A1", "MZB1", "TNFRSF17"), c("CD79A", "MS4A1", "MZB1", "TNFRSF17"))
  stopifnot(all(ids %in% rownames(counts)))
}
y <- calcNormFactors(DGEList(counts[, per$sample_id]))
lcpm <- cpm(y, log = TRUE, prior.count = 1)[ids, , drop = FALSE]
zc <- function(v, g) ave(v, g, FUN = function(x) (x - mean(x)) / sd(x))
for (s in names(ids)) per[, (paste0("z_", s)) := zc(lcpm[ids[[s]], match(sample_id, colnames(lcpm))], cohort)]
per[, Bscore := (z_CD79A + z_MS4A1) / 2]
per[, Pscore := (z_MZB1 + z_TNFRSF17) / 2]
stopifnot(!anyNA(per$Bscore), !anyNA(per$Pscore))
fwrite(per, file.path(out, "person_table.tsv"), sep = "\t")
depth <- per[, .(persons = .N, n_ge10 = as.numeric(sum(n >= 10)), n_ge30 = as.numeric(sum(n >= 30)), n_ge100 = as.numeric(sum(n >= 100)),
                 median_n = as.numeric(median(n)), q25_n = as.numeric(quantile(n, .25)), q75_n = as.numeric(quantile(n, .75)), max_n = as.numeric(max(n))), by = cohort]
fwrite(depth, file.path(out, "depth_by_cohort.tsv"), sep = "\t")
print(depth)

# ------------------------------------------------------------------ beta-binomial machinery
lbb <- function(k, n, mu, rho) {
  s <- (1 - rho) / rho; a <- mu * s; b <- (1 - mu) * s
  lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + lbeta(k + a, n - k + b) - lbeta(a, b)
}
rho_of <- function(t) 1e-4 + (1 - 2e-4) * plogis(t)
bb_negll <- function(par, k, n, X) {
  p <- ncol(X); mu <- plogis(drop(X %*% par[1:p]))
  mu <- pmin(pmax(mu, 1e-9), 1 - 1e-9)
  -sum(lbb(k, n, mu, rho_of(par[p + 1])))
}
bb_fit <- function(k, n, X, start = NULL) {
  p <- ncol(X)
  if (is.null(start)) {
    b0 <- tryCatch(suppressWarnings(coef(glm.fit(X, cbind(k, n - k), family = binomial()))), error = function(e) rep(0, p))
    b0[!is.finite(b0)] <- 0
    start <- c(b0, qlogis((0.05 - 1e-4) / (1 - 2e-4)))
  }
  o <- tryCatch(optim(start, bb_negll, k = k, n = n, X = X, method = "BFGS", control = list(maxit = 500, reltol = 1e-10)),
                error = function(e) NULL)
  if (is.null(o)) return(list(ok = FALSE))
  H <- tryCatch(optimHess(o$par, bb_negll, k = k, n = n, X = X), error = function(e) NULL)
  V <- if (is.null(H)) NULL else tryCatch(solve(H), error = function(e) NULL)
  se <- if (is.null(V)) rep(NA_real_, length(o$par)) else sqrt(pmax(diag(V), 0))
  list(ok = o$convergence == 0 && all(is.finite(se[1:p])), par = o$par, se = se, ll = -o$value, p = p)
}
design <- function(d, xcol, covs) {
  X <- cbind(`(Intercept)` = 1, x = d[[xcol]])
  for (cv in covs) {
    v <- d[[cv]]
    if (length(unique(v)) > 1L) X <- cbind(X, structure(matrix(v, ncol = 1), dimnames = list(NULL, cv)))
  }
  X
}
pool_dl <- function(b, se) {
  ok <- is.finite(b) & is.finite(se) & se > 0
  b <- b[ok]; se <- se[ok]; k <- length(b)
  if (k < 2L) return(list(b = if (k == 1L) b else NA_real_, se = if (k == 1L) se else NA_real_, z = NA_real_, p = NA_real_, tau2 = NA_real_, I2 = NA_real_, Q = NA_real_, Qp = NA_real_, k = k))
  w <- 1 / se^2; bfe <- sum(w * b) / sum(w); Q <- sum(w * (b - bfe)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  ws <- 1 / (se^2 + tau2); bp <- sum(ws * b) / sum(ws); sp <- 1 / sqrt(sum(ws))
  list(b = bp, se = sp, z = bp / sp, p = 2 * pnorm(-abs(bp / sp)), tau2 = tau2,
       I2 = if (Q > 0) max(0, (Q - (k - 1)) / Q) else 0, Q = Q, Qp = pchisq(Q, k - 1, lower.tail = FALSE), k = k)
}
cohort_slopes <- function(d, xcol, covs, kcol = "k", ncol_ = "n", min_n = MIN_N, min_persons = 20L, cohorts = STAGE_COHORTS) {
  rows <- lapply(cohorts, function(cc) {
    x <- d[cohort == cc & get(ncol_) >= min_n & !is.na(get(xcol)) & !is.na(sexM)]
    if (nrow(x) < min_persons || length(unique(x[[xcol]])) < 2L) return(NULL)
    X <- design(x, xcol, covs)
    f <- bb_fit(x[[kcol]], x[[ncol_]], X)
    if (!f$ok) return(data.table(cohort = cc, persons = nrow(x), b = NA_real_, se = NA_real_, rho = NA_real_, ok = FALSE))
    data.table(cohort = cc, persons = nrow(x), b = f$par[2], se = f$se[2], rho = rho_of(f$par[f$p + 1]), ok = TRUE)
  })
  rbindlist(rows)
}
pool_table <- function(tab, label) {
  p <- pool_dl(tab$b, tab$se)
  data.table(analysis = label, cohorts = p$k, persons = sum(tab$persons[is.finite(tab$b)]), pooled_b = p$b, se = p$se, z = p$z, p = p$p,
             tau2 = p$tau2, I2 = p$I2, Q = p$Q, Q_p = p$Qp)
}

# ------------------------------------------------------------------ controls
per[, `:=`(fA = fifelse(n > 0, n_A / n, NA_real_), fG = fifelse(n > 0, n_G / n, NA_real_), fM = fifelse(n > 0, n_M / n, NA_real_))]
pc1 <- per[n > 0, max(abs(fA + fG + fM - 1))]
log_msg("PC1 max |sum-1|:", pc1)
cor_row <- function(sc, cc) {
  x <- per[cohort == cc]
  r <- suppressWarnings(cor(x$ig_rate, x[[sc]], method = "spearman"))
  data.table(score = sc, cohort = cc, persons = nrow(x), rho = r, z = atanh(r), vz = 1 / (nrow(x) - 3))
}
cors <- rbindlist(lapply(c("Bscore", "Pscore"), function(sc) rbindlist(lapply(ROSTER_COHORTS, function(cc) cor_row(sc, cc)))))
pc2 <- rbindlist(lapply(unique(cors$score), function(sc) {
  x <- cors[score == sc]
  m <- rma(yi = x$z, vi = x$vz, method = "DL")
  data.table(score = sc, pooled_rho = tanh(as.numeric(m$beta)), all_positive = all(x$rho > 0), pass = all(x$rho > 0) && tanh(as.numeric(m$beta)) >= 0.30)
}))
fwrite(cors, file.path(out, "pc2_ig_rate_correlations.tsv"), sep = "\t")
fwrite(pc2, file.path(out, "pc2_summary.tsv"), sep = "\t")
print(pc2)

# ------------------------------------------------------------------ T1: IgA-high subgroup (mixture)
log_msg("T1 mixture")
d1 <- per[n >= MIN_N]
d1[, ci := match(cohort, ROSTER_COHORTS)]
C <- length(ROSTER_COHORTS)
fit_m1 <- function(d) {
  fits <- lapply(seq_len(C), function(i) {
    x <- d[ci == i]; X <- matrix(1, nrow(x), 1)
    f <- bb_fit(x$k, x$n, X); list(mu = f$par[1], th = f$par[2], ll = f$ll, ok = f$ok)
  })
  list(mu = vapply(fits, `[[`, 0, "mu"), th = vapply(fits, `[[`, 0, "th"), ll = sum(vapply(fits, `[[`, 0, "ll")),
       ok = all(vapply(fits, `[[`, TRUE, "ok")))
}
negll_m2 <- function(par, k, n, ci) {
  m <- par[1:C][ci]; rho <- rho_of(par[(C + 1):(2 * C)])[ci]
  delta <- exp(par[2 * C + 1]); w <- 0.02 + 0.48 * plogis(par[2 * C + 2])
  mu1 <- pmin(pmax(plogis(m), 1e-9), 1 - 1e-9); mu2 <- pmin(pmax(plogis(m + delta), 1e-9), 1 - 1e-9)
  l1 <- log(1 - w) + lbb(k, n, mu1, rho); l2 <- log(w) + lbb(k, n, mu2, rho)
  mx <- pmax(l1, l2)
  -sum(mx + log(exp(l1 - mx) + exp(l2 - mx)))
}
fit_m2 <- function(d, m1) {
  starts <- list(c(1.5, 0.1), c(2.5, 0.1), c(1.5, 0.3), c(3, 0.05))
  best <- NULL
  for (s in starts) {
    st <- c(m1$mu, m1$th, log(s[1]), qlogis((s[2] - 0.02) / 0.48))
    o <- tryCatch(optim(st, negll_m2, k = d$k, n = d$n, ci = d$ci, method = "BFGS", control = list(maxit = 400, reltol = 1e-9)), error = function(e) NULL)
    if (!is.null(o) && is.finite(o$value) && (is.null(best) || o$value < best$value)) best <- o
  }
  best
}
m1 <- fit_m1(d1); m2 <- fit_m2(d1, m1)
lrt_obs <- max(0, 2 * (-m2$value - m1$ll))
boot_one <- function(r) {
  set.seed(SEED + 100000L + r)
  a_b <- {
    rho <- rho_of(m1$th)[d1$ci]; mu <- plogis(m1$mu)[d1$ci]; s <- (1 - rho) / rho
    pp <- rbeta(nrow(d1), mu * s, (1 - mu) * s); rbinom(nrow(d1), d1$n, pp)
  }
  db <- copy(d1)[, k := a_b]
  f1 <- fit_m1(db); f2 <- fit_m2(db, f1)
  if (is.null(f2) || !f1$ok) return(NA_real_)
  max(0, 2 * (-f2$value - f1$ll))
}
lrt_null <- unlist(mclapply(seq_len(B_BOOT), boot_one, mc.cores = NCORES))
boot_p <- (1 + sum(lrt_null >= lrt_obs, na.rm = TRUE)) / (1 + sum(is.finite(lrt_null)))
par2 <- m2$par
delta_hat <- exp(par2[2 * C + 1]); w_hat <- 0.02 + 0.48 * plogis(par2[2 * C + 2])
mix_rows <- data.table(cohort = ROSTER_COHORTS, base_logit_mean = par2[1:C], base_fraction = plogis(par2[1:C]),
                       high_fraction = plogis(par2[1:C] + delta_hat), rho = rho_of(par2[(C + 1):(2 * C)]),
                       persons = d1[, .N, by = ci][order(ci)]$N)
fwrite(mix_rows, file.path(out, "t1_mixture_components.tsv"), sep = "\t")
sub_supported <- is.finite(boot_p) && boot_p < 0.05 && w_hat >= 0.05 && w_hat <= 0.5 && delta_hat >= 1
t1 <- data.table(persons = nrow(d1), ll_m1 = m1$ll, ll_m2 = -m2$value, lrt = lrt_obs, boot_B = sum(is.finite(lrt_null)), boot_p = boot_p,
                 w = w_hat, delta = delta_hat, delta_odds_ratio = exp(delta_hat), supported = sub_supported,
                 lrt_null_q95 = quantile(lrt_null, .95, na.rm = TRUE), lrt_null_max = max(lrt_null, na.rm = TRUE))
fwrite(t1, file.path(out, "t1_mixture_test.tsv"), sep = "\t")
print(t1)
# posterior membership and, if supported, prevalence by stage group (CMH)
m <- par2
mu1 <- plogis(m[1:C][d1$ci]); mu2 <- plogis(m[1:C][d1$ci] + delta_hat); rho <- rho_of(m[(C + 1):(2 * C)])[d1$ci]
l1 <- log(1 - w_hat) + lbb(d1$k, d1$n, mu1, rho); l2 <- log(w_hat) + lbb(d1$k, d1$n, mu2, rho)
d1[, post_high := 1 / (1 + exp(l1 - l2))]
d1[, high := post_high > 0.5]
fwrite(d1[, .(sample_id, analysis_unit_id, cohort, k, n, post_high, high)], file.path(out, "t1_person_membership.tsv"), sep = "\t")
cmh <- NULL
if (sub_supported) {
  dd <- d1[!is.na(stage) & cohort %in% STAGE_COHORTS]
  dd[, late := factor(stage >= 2, levels = c(FALSE, TRUE))]
  dd[, hi := factor(high, levels = c(FALSE, TRUE))]
  arr <- table(dd$late, dd$hi, dd$cohort)
  cmh <- tryCatch(mantelhaen.test(arr), error = function(e) NULL)
  if (!is.null(cmh)) fwrite(data.table(statistic = unname(cmh$statistic), p = cmh$p.value, odds_ratio = unname(cmh$estimate)), file.path(out, "t1_cmh_high_by_late_stage.tsv"), sep = "\t")
}

# ------------------------------------------------------------------ T2: primary stage association, secondary, sensitivities
log_msg("T2")
COVS <- c("sexM", "Bscore")
st_tab <- cohort_slopes(per, "stage", COVS)
t2 <- pool_table(st_tab, "T2_primary_fibrosis_stage_per_step")
fwrite(st_tab, file.path(out, "t2_per_cohort.tsv"), sep = "\t")
stop_rule <- per[cohort %in% STAGE_COHORTS & !is.na(stage) & n >= MIN_N, .N]
log_msg("persons with stage and n>=10:", stop_rule)
print(st_tab); print(t2)
# cross-check against metafor
mchk <- tryCatch(rma(yi = st_tab$b[st_tab$ok], sei = st_tab$se[st_tab$ok], method = "DL"), error = function(e) NULL)
if (!is.null(mchk)) writeLines(sprintf("metafor DL pooled b=%.6f se=%.6f ; own b=%.6f se=%.6f", as.numeric(mchk$beta), mchk$se, t2$pooled_b, t2$se), file.path(out, "t2_metafor_check.txt"))

nas_tab <- cohort_slopes(per, "nas", COVS, cohorts = ROSTER_COHORTS)
s1 <- pool_table(nas_tab, "S1_NAS_per_point")
dis_tab <- cohort_slopes(per, "disease", COVS, cohorts = ROSTER_COHORTS)
s2 <- pool_table(dis_tab, "S2_disease_vs_control")
sec <- rbind(s1, s2)
sec[, q_BH := p.adjust(p, "BH")]
fwrite(nas_tab, file.path(out, "s1_per_cohort.tsv"), sep = "\t"); fwrite(dis_tab, file.path(out, "s2_per_cohort.tsv"), sep = "\t")
fwrite(sec, file.path(out, "secondary_family_BH2.tsv"), sep = "\t")
print(sec)

# sensitivities
sens <- list()
sens[["depth>=30"]] <- pool_table(cohort_slopes(per, "stage", COVS, min_n = MIN_N_SENS), "sens1_depth_ge30")
per_alloc <- copy(per)[, `:=`(k = k2, n = n2)]
sens[["alloc"]] <- pool_table(cohort_slopes(per_alloc, "stage", COVS), "sens2_proportional_allocation_of_mixed_sets")
sens[["plasma"]] <- pool_table(cohort_slopes(per, "stage", c("sexM", "Bscore", "Pscore", "ig_rate")), "sens3_plus_plasma_score_and_ig_rate")
per[, late := as.numeric(stage >= 2)]
sens[["binary"]] <- pool_table(cohort_slopes(per, "late", COVS), "sens6_F2plus_vs_F0F1")
for (lo in STAGE_COHORTS) sens[[paste0("loco_", lo)]] <- pool_table(st_tab[cohort != lo], paste0("sens5_leave_out_", lo))
# shared-slope pooled model with cohort intercepts, own optimiser versus aod::betabin
dp <- per[cohort %in% STAGE_COHORTS & n >= MIN_N & !is.na(stage) & !is.na(sexM)]
dp[, cohort := factor(cohort)]
Xp <- model.matrix(~ cohort + stage + sexM + Bscore, dp)
fp <- bb_fit(dp$k, dp$n, Xp)
ix <- which(colnames(Xp) == "stage")
aodfit <- tryCatch(betabin(cbind(k, n - k) ~ cohort + stage + sexM + Bscore, ~ 1, data = as.data.frame(dp)), error = function(e) NULL)
aodb <- if (!is.null(aodfit)) tryCatch(unname(aodfit@param["stage"]), error = function(e) NA_real_) else NA_real_
sens[["shared"]] <- data.table(analysis = "sens4_shared_slope_cohort_intercepts_own_optimiser", cohorts = 4L, persons = nrow(dp),
                               pooled_b = fp$par[ix], se = fp$se[ix], z = fp$par[ix] / fp$se[ix], p = 2 * pnorm(-abs(fp$par[ix] / fp$se[ix])),
                               tau2 = NA_real_, I2 = NA_real_, Q = NA_real_, Q_p = NA_real_)
writeLines(sprintf("own shared-slope b=%.6f ; aod::betabin stage coefficient=%.6f", fp$par[ix], aodb), file.path(out, "sens4_aod_check.txt"))
sens_tab <- rbindlist(sens, fill = TRUE)
fwrite(sens_tab, file.path(out, "sensitivities.tsv"), sep = "\t")
print(sens_tab[, .(analysis, cohorts, persons, pooled_b, se, p)])

# ------------------------------------------------------------------ negative control: stage permuted within cohort
log_msg("permutations", N_PERM)
obs_z <- t2$z
prep <- lapply(STAGE_COHORTS, function(cc) {
  x <- per[cohort == cc & n >= MIN_N & !is.na(stage) & !is.na(sexM)]
  if (nrow(x) < 20L) return(NULL)
  list(cc = cc, x = x, X = design(x, "stage", COVS))
})
prep <- Filter(Negate(is.null), prep)
perm_one <- function(r) {
  set.seed(SEED + r)
  bs <- vapply(prep, function(pp) {
    X <- pp$X; X[, 2] <- sample(X[, 2])
    f <- bb_fit(pp$x$k, pp$x$n, X)
    if (!f$ok) c(NA_real_, NA_real_) else c(f$par[2], f$se[2])
  }, numeric(2))
  pl <- pool_dl(bs[1, ], bs[2, ])
  c(z = pl$z, p = pl$p)
}
perm <- do.call(rbind, mclapply(seq_len(N_PERM), perm_one, mc.cores = NCORES))
perm <- perm[is.finite(perm[, "z"]), , drop = FALSE]
perm_p <- (1 + sum(abs(perm[, "z"]) >= abs(obs_z))) / (1 + nrow(perm))
nc <- data.table(permutations = nrow(perm), perm_p_T2 = perm_p, frac_nominal_p_lt_0.05 = mean(perm[, "p"] < 0.05),
                 null_passes_rule_0.03_to_0.07 = mean(perm[, "p"] < 0.05) >= 0.03 & mean(perm[, "p"] < 0.05) <= 0.07,
                 abs_z_q95 = quantile(abs(perm[, "z"]), .95), obs_abs_z = abs(obs_z))
fwrite(nc, file.path(out, "negative_control_permutation.tsv"), sep = "\t")
print(nc)

# ------------------------------------------------------------------ descriptive: other five cohorts (library level, no tests)
other <- lib[!cohort %in% ROSTER_COHORTS]
other_sum <- other[, .(libraries = .N, n_ge10 = sum(n >= 10), median_fraction = median((n_A / n)[n >= 10]),
                       q25 = quantile((n_A / n)[n >= 10], .25), q75 = quantile((n_A / n)[n >= 10], .75),
                       share_iga_gt_0.5 = mean((n_A / n)[n >= 10] > 0.5)), by = cohort]
fwrite(other_sum, file.path(out, "descriptive_other_cohorts_library_level.tsv"), sep = "\t")
roster_sum <- per[n >= MIN_N, .(persons = .N, median_fraction = median(n_A / n), q25 = quantile(n_A / n, .25), q75 = quantile(n_A / n, .75),
                                share_iga_gt_0.5 = mean(n_A / n > 0.5)), by = cohort]
fwrite(roster_sum, file.path(out, "iga_fraction_by_cohort_persons.tsv"), sep = "\t")
print(roster_sum); print(other_sum)

# ------------------------------------------------------------------ summary
dropped <- per[, .(persons = .N, below_min_depth = sum(n < MIN_N), na_sex = sum(is.na(sexM)), na_stage = sum(is.na(stage))), by = cohort]
fwrite(dropped, file.path(out, "person_exclusions_by_cohort.tsv"), sep = "\t")
summ <- list(smoke = SMOKE, seed = SEED, persons_total = nrow(per), persons_n_ge10 = per[n >= MIN_N, .N],
             stage_persons_n_ge10 = stop_rule, stop_rule_underpowered = stop_rule < 100L,
             pc1_max_abs_sum_minus_1 = pc1, pc2 = pc2, t1 = t1, t2 = t2, secondary = sec, negative_control = nc,
             bootstrap_B = B_BOOT, permutations = N_PERM, script = "h38_iga.R")
write_json(summ, file.path(out, "summary.json"), auto_unbox = TRUE, digits = 8, pretty = TRUE)
log_msg("done")
