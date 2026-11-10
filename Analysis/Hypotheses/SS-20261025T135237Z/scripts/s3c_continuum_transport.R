# SS job S3c: SteatoSITE as a third control-free transport cohort for the frozen 139-gene continuum loading.
# Fork of scripts/analysis/histology_anchored_continuum/control_free_transfer/01_project_control_free_cohorts.R
# (CFT-v1, 2026-08-28), which gave rho 0.510 in 93 GSE174478 samples and 0.710 in 66 GSE240729 samples.
# Unchanged from CFT-v1: frozen loading and discovery centre (hac-continuum-20260818T024923Z, same guards), gene
# universe = the 23,370-gene merged_dge substrate of that run (matched on unversioned Ensembl ID), within-cohort
# normalisation (row sum > n samples, log2(count + 1), limma normalizeQuantiles), projection without refit,
# statistic = Spearman rho of the projection with recorded fibrosis stage, sealed 1,000-draw matched nulls A
# (primary) and B (caliper), one-sided empirical p = (1 + #{null rho >= obs}) / (1 + draws retained), partial
# Spearman on log10 library size, the three competitors (log library size, hepatocyte fraction, 10-fold
# out-of-fold 6-lineage CLR; 10,000 bootstrap draws, seed 20260828) when fractions are given.
# Changed: input is the S2 layout (counts_all_genes.rds, units.tsv; stage joined from --stage); one cohort; the
# counting-convention audit and released-cohort reproduction of CFT-v1 are not repeated (the GSE174478 mock run is
# the reproduction check); strata are S1 all staged units (primary) and S_noF4 (stage < 4; explants are F4);
# added partial Spearman sensitivities on batch, batch + inferred_sex, and the prespecified quality covariates.
# SteatoSITE stage code 1 is F0-F1, so rho is on codes 1-4. Unit: one library per genetic individual (S2 include).
# Cross-sectional. Statistic and pass rule: PRESPEC_SS_addendum_S3c.md.
# Usage: Rscript s3c_continuum_transport.R --input <S2 dir> --stage <stage tsv> --out <new dir>
#        [--fractions <S3e fractions.tsv>] [--cohort <label>]
suppressPackageStartupMessages({library(data.table); library(limma)})
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/ss_common.R")
A <- ss_args(extra = c("fractions", "cohort")); OUT <- A$out
COHORT <- if (is.null(A$cohort)) "SteatoSITE" else A$cohort

H <- file.path(SS_ROOT, "RNA-seq/results/histology_anchored_continuum")
RELEASE <- file.path(H, "candidates/hac-continuum-20260818T024923Z")
PRESPEC <- file.path(H, "prespec/mnc-fib-v1-20260828T133158Z")
DGE <- file.path(SS_ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION/arms/F_five/results/integration/merged_dge.rds")
IN_SHA <- c(
  setNames("1370ace1951f46847fc08d15cd70075159ae582dc8bf950ab2c5c4a0ae061d57", file.path(RELEASE, "projection/fixed_projection_loadings.tsv")),
  setNames("0371dd3013cc18385786b9246a9720cbb9b91498412ec2481aaa4fab50cd2b7c", file.path(RELEASE, "reproduction/discovery_signature_model.rds")),
  setNames("5e32bbbacfc94f0bac2b2a6a8296b8a8cbb82749839854b48af7ed4e4cdc5206", file.path(RELEASE, "reproduction/discovery_normalized_expression.rds")),
  setNames("4ffb6133fdc44854bec37a6572d4f51b0c211f3de9aa5ade1f4eb2b92d35fed4", file.path(RELEASE, "reproduction/participant_scores.tsv")),
  setNames("1ac467d18e0bc45a8cdab900b178e38a1730fa634fead3c6b22e5010ff308b4d", file.path(PRESPEC, "null_gene_sets.tsv.gz")),
  setNames("91e88126efb3ce13fbac31705b9ece56943caddd6b843b896b142d65342a5576", file.path(PRESPEC, "null_gene_set_digests.tsv")),
  setNames("48375e7a42c8a9f0c63effcf5557444e11d2d136504cde33ca0d93b7e832fd29", file.path(PRESPEC, "null_gene_sets_caliper.tsv.gz")),
  setNames("1c8303499580e3596af484cc66fd40339ddb9e3e562132bb72135a900dd780ad", file.path(PRESPEC, "null_gene_set_digests_caliper.tsv")),
  setNames("51414c9676e59b5d1de755164aac02123d4d6f3e083d7b67845be41246a5789b", file.path(PRESPEC, "sampling_universe.tsv")),
  setNames("bc3ea19c8b064861339eba862abf58e82ddb5460c622b550dc45aa3b5fcd1f76", DGE))
ss_check_sha(IN_SHA)
ss_frozen(c("stage_continuum/comparison_summary.tsv"))  # frozen/ manifest intact

B_DRAWS <- 1000L; BOOT_B <- 10000L; BOOT_SEED <- 20260828L; CV_SEED <- 20260828L
MIN_DRAWS <- 800L; MIN_N <- 30L
LIN6 <- c("Hepatocytes", "Macrophages", "Mono+mono derived cells", "cDC2s", "T cells", "pDCs")
FLOOR_A <- c(lo = 0.465325, hi = 0.616087)
DISCOVERY <- data.table(cohort = c("GSE174478", "GSE240729"), n = c(93L, 66L), rho = c(0.510385240165136, 0.710030538764475),
                        p_emp_nullA = c(0.116883116883117, 0.031968031968032))

say <- function(...) { cat(sprintf(...), "\n", sep = ""); flush.console() }
abort <- function(...) stop(paste0("S3c-ABORT: ", sprintf(...)), call. = FALSE)
chk <- function(cond, ...) if (!isTRUE(cond)) abort(...)
base_id <- function(x) sub("\\.[0-9]+$", "", as.character(x))
sp_rho <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
close_to <- function(a, b, tol) is.finite(a) && is.finite(b) && abs(a - b) <= tol
guard <- list()
record <- function(check, expected, observed, ok, fatal = TRUE) {
  guard[[length(guard) + 1L]] <<- data.table(check = check, expected = as.character(expected), observed = as.character(observed), pass = ok, fatal = fatal)
  if (!ok && fatal) abort("guard failed: %s (expected %s, observed %s)", check, expected, observed)
}

## frozen axis (CFT-v1 tier 1, unchanged)
load_tab <- fread(file.path(RELEASE, "projection", "fixed_projection_loadings.tsv"))
load_tab[, gene_id_base := base_id(gsub('"', "", gene_id_base))]
sig139 <- load_tab$gene_id_base[load_tab$used_for_fixed_projection %in% c(TRUE, "TRUE")]
record("n_used_for_fixed_projection", 139L, length(sig139), length(sig139) == 139L)
sig_sha <- digest::digest(paste0(paste(sort(sig139), collapse = "\n"), "\n"), algo = "sha256", serialize = FALSE)
record("signature_sorted_139_sha256", "165b9168648651b16d6ffe59b4d1f1c18515308e2fd63701cc21d1c6a43ba8ff", sig_sha,
       identical(sig_sha, "165b9168648651b16d6ffe59b4d1f1c18515308e2fd63701cc21d1c6a43ba8ff"))
model <- readRDS(file.path(RELEASE, "reproduction", "discovery_signature_model.rds"))
loading <- model$discovery_loading_oriented; names(loading) <- base_id(names(loading))
centre <- model$discovery_center; names(centre) <- base_id(names(centre))
record("frozen_loading_sum", 9.774235120123, sprintf("%.12f", sum(loading)), close_to(sum(loading), 9.774235120123, 1e-9))
record("frozen_loading_sumsq", 1, sprintf("%.12f", sum(loading^2)), close_to(sum(loading^2), 1, 1e-9))
record("discovery_centre_sum", 771.590375354, sprintf("%.9f", sum(centre)), close_to(sum(centre), 771.590375354, 1e-8))
chk(setequal(names(loading), sig139), "frozen loading names are not the 139-gene set")
M <- readRDS(file.path(RELEASE, "reproduction", "discovery_normalized_expression.rds")); rownames(M) <- base_id(rownames(M))
rp <- fread(file.path(RELEASE, "reproduction", "participant_scores.tsv"))[match(colnames(M), sample_id)]
released_pc1 <- as.numeric(rp$pc1_released_operational); chk(all(is.finite(released_pc1)), "released PC1 non-finite")
project_set <- function(mat, g, ld, ct) as.numeric(crossprod(ld[g], sweep(mat[g, , drop = FALSE], 1L, ct[g], FUN = "-")))

substrate <- base_id(rownames(readRDS(DGE)$counts))
record("merged_dge_substrate_n_genes", 23370L, length(substrate), length(substrate) == 23370L && !anyDuplicated(substrate))

## sealed nulls (unchanged)
read_sets <- function(gzf, digf, coll_sha) {
  tb <- fread(gzf); chk(nrow(tb) == B_DRAWS * 139L, "%s: wrong row count", gzf)
  sets <- split(tb$gene_id_base, tb$draw_index)[as.character(seq_len(B_DRAWS))]
  per <- vapply(sets, function(s) digest::digest(paste(sort(s), collapse = "\n"), algo = "sha256", serialize = FALSE), character(1))
  chk(identical(unname(per), fread(digf)$set_sha256), "%s: per-draw digests differ", digf)
  chk(identical(digest::digest(paste(unname(per), collapse = "\n"), algo = "sha256", serialize = FALSE), coll_sha), "%s: collection digest", gzf)
  lapply(sets, as.character)
}
setsA <- read_sets(file.path(PRESPEC, "null_gene_sets.tsv.gz"), file.path(PRESPEC, "null_gene_set_digests.tsv"),
                   "1b6e2f040329e25faf8e9be85e0a8e5ed2920fdb2c232b14553143733d3f7e36")
setsB <- read_sets(file.path(PRESPEC, "null_gene_sets_caliper.tsv.gz"), file.path(PRESPEC, "null_gene_set_digests_caliper.tsv"),
                   "275454cce7a32de013b5efdb23a21c8ae3e1bba1d734908cc1518d5fd89c1fdb")
precompute <- function(sets) {
  LD <- CT <- matrix(NA_real_, 139L, B_DRAWS); GN <- matrix(NA_character_, 139L, B_DRAWS); bad <- integer(0)
  for (b in seq_len(B_DRAWS)) {
    pca <- prcomp(t(M[sets[[b]], , drop = FALSE]), center = TRUE, scale. = FALSE)
    sgn <- sign(cor(pca$x[, 1L], released_pc1)); if (sgn == 0) { bad <- c(bad, b); next }
    ld <- pca$rotation[, 1L] * sgn; GN[, b] <- names(ld); LD[, b] <- unname(ld); CT[, b] <- unname(pca$center[names(ld)])
  }
  list(LD = LD, CT = CT, GN = GN, bad = bad)
}
say("[null] precomputing discovery PC1 loadings for the sealed draws")
PRE <- list(A = precompute(setsA), B = precompute(setsB)); rm(M); gc()

## cohort: S2 counts -> substrate -> frozen normalisation -> projection
u <- ss_units(A$input, A$stage)
cm <- readRDS(file.path(A$input, "counts_all_genes.rds"))
if (inherits(cm, "DGEList")) cm <- cm$counts
cm <- as.matrix(cm); chk(all(u$sample_alias %in% colnames(cm)), "included units missing from counts_all_genes")
cm <- cm[, u$sample_alias, drop = FALSE]
rownames(cm) <- base_id(rownames(cm)); cm <- cm[!duplicated(rownames(cm)) & rownames(cm) %in% substrate, , drop = FALSE]
chk(nrow(cm) > 10000L, "only %d genes shared with the substrate", nrow(cm))
chk(!anyNA(cm) && all(cm >= 0) && all(cm == floor(cm)), "counts are not raw non-negative integers")
lib <- colSums(cm)
keep <- rowSums(cm) > ncol(cm); chk(sum(keep) >= 1000L, "fewer than 1,000 genes after the row-sum filter")
Xc <- normalizeQuantiles(log2(cm[keep, , drop = FALSE] + 1)); dimnames(Xc) <- dimnames(cm[keep, , drop = FALSE])
gg <- sig139[sig139 %in% rownames(Xc)]
u[, fixed_projection := project_set(Xc, gg, loading, centre)]
u[, log10_libsize := log10(as.numeric(lib[sample_alias]))]
NULL_UNIVERSE <- fread(file.path(PRESPEC, "sampling_universe.tsv"))$gene_id_base
cov <- data.table(cohort = COHORT, n_included = nrow(u), n_staged = sum(is.finite(u$stage_num)), n_genes_in_universe = nrow(cm),
                  n_genes_after_rowsum_filter = nrow(Xc), n_signature_present = length(gg),
                  signature_weight_present = sum(loading[gg]^2),
                  n_null_universe_genes_present = sum(NULL_UNIVERSE %in% rownames(Xc)), n_null_universe_total = length(NULL_UNIVERSE))
record("signature_coverage", 139L, length(gg), length(gg) == 139L, fatal = FALSE)
rm(cm); gc()

fr <- NULL
if (!is.null(A$fractions)) {
  fr <- fread(A$fractions, check.names = FALSE); chk(all(LIN6 %in% names(fr)), "fractions lack a 6-lineage column")
  fr <- fr[match(u$sample_alias, fr$sample_alias)]; for (L in LIN6) u[[L]] <- as.numeric(fr[[L]])
} else for (L in LIN6) u[[L]] <- NA_real_
u[, has_composition := rowSums(is.finite(as.matrix(.SD))) == length(LIN6), .SDcols = LIN6]
fwrite(u, file.path(OUT, "per_unit_scores.tsv"), sep = "\t")
fwrite(cov, file.path(OUT, "coverage.tsv"), sep = "\t")

## helpers (CFT-v1)
clr6 <- function(P) {
  P <- as.matrix(P); P[!is.finite(P)] <- 0; P <- P / rowSums(P)
  P <- P + 0.5 * min(P[P > 0]); P <- P / rowSums(P); L <- log(P); L - rowMeans(L)
}
partial_spearman <- function(x, y, Z) {
  rx <- rank(x); ry <- rank(y); Zr <- apply(as.matrix(Z), 2L, rank)
  suppressWarnings(cor(residuals(lm(rx ~ Zr)), residuals(lm(ry ~ Zr))))
}
tie_ceiling <- function(outcome, tiebreak) {
  star <- numeric(length(outcome))
  for (g in split(seq_along(outcome), outcome)) star[g] <- outcome[g[1]] + rank(tiebreak[g], ties.method = "first") / (length(g) + 1)
  sp_rho(star, outcome)
}
mde80 <- function(q95, n) tanh(atanh(q95) + qnorm(0.80) * 1.06 / sqrt(n - 3))
dummies <- function(f) { f <- droplevels(factor(f)); if (nlevels(f) < 2) NULL else model.matrix(~ f)[, -1, drop = FALSE] }
# batch-type covariates enter as dummies (rank of a 0/1 column is an affine map of it, so partial_spearman is unchanged)
partial_sets <- function(dd) {
  Zb <- dummies(dd$batch); Zs <- if (all(dd$inferred_sex %in% c("F", "M"))) dummies(dd$inferred_sex) else NULL
  q <- as.matrix(dd[, .(log10_assigned, unique_pct, intron_exon_ratio)]); if (any(!is.finite(q))) q <- NULL
  Filter(Negate(is.null), list(log_libsize = cbind(dd$log10_libsize), batch = Zb,
                               batch_sex = if (!is.null(Zb) || !is.null(Zs)) cbind(Zb, Zs) else NULL,
                               quality = if (!is.null(q)) cbind(q, Zb) else NULL))
}

## strata and matched nulls
STR <- list(S1_all_staged = u[is.finite(stage_num)], S_noF4 = u[is.finite(stage_num) & stage_num < 4])
nullres <- list(); partrows <- list(); comp_rows <- list()
for (st in names(STR)) {
  dd <- STR[[st]]; yk <- dd$stage_num
  if (nrow(dd) < 10L || length(unique(yk)) < 2L) next
  obs <- sp_rho(dd$fixed_projection, yk); obs_pl <- partial_spearman(dd$fixed_projection, yk, dd$log10_libsize)
  ceil <- tie_ceiling(yk, dd$fixed_projection); pos <- match(dd$sample_alias, colnames(Xc)); chk(!anyNA(pos), "unit not in matrix")
  for (nl in c("A", "B")) {
    P <- PRE[[nl]]; rho <- rhop <- rep(NA_real_, B_DRAWS); n_incomplete <- 0L
    for (b in seq_len(B_DRAWS)) {
      if (b %in% P$bad) next
      g <- P$GN[, b]; if (!all(g %in% rownames(Xc))) { n_incomplete <- n_incomplete + 1L; next }
      sv <- project_set(Xc, g, setNames(P$LD[, b], g), setNames(P$CT[, b], g))[pos]
      rho[b] <- sp_rho(sv, yk); rhop[b] <- partial_spearman(sv, yk, dd$log10_libsize)
    }
    rv <- rho[is.finite(rho)]; rpv <- rhop[is.finite(rhop)]
    q <- if (length(rv)) quantile(rv, c(0.05, 0.5, 0.95, 0.99), type = 7, names = FALSE) else rep(NA_real_, 4)
    nullres[[length(nullres) + 1L]] <- data.table(
      cohort = COHORT, stratum = st, outcome = "fibrosis_stage", null = nl, n = nrow(dd),
      n_draws_retained = length(rv), n_draws_incomplete_coverage = n_incomplete,
      testable = length(rv) >= MIN_DRAWS && nrow(dd) >= MIN_N,
      observed_rho = obs, null_median = q[2], null_q05 = q[1], null_q95 = q[3], null_q99 = q[4],
      null_mean = mean(rv), null_sd = sd(rv), margin_over_median = obs - q[2],
      p_emp_signed = (1 + sum(rv >= obs)) / (1 + length(rv)),
      observed_partial_rho = obs_pl, null_partial_q95 = if (length(rpv)) quantile(rpv, 0.95, type = 7, names = FALSE) else NA_real_,
      p_emp_partial = (1 + sum(rpv >= obs_pl)) / (1 + length(rpv)),
      tie_ceiling = ceil, obs_over_ceiling = obs / ceil, mde80_vs_null_q95 = mde80(q[3], nrow(dd)),
      inside_within_disease_floor_a = obs >= FLOOR_A[["lo"]] & obs <= FLOOR_A[["hi"]])
    fwrite(data.table(draw_index = seq_len(B_DRAWS), rho_signed = rho, rho_partial_libsize = rhop),
           file.path(OUT, sprintf("nullvec_%s_%s_null%s.tsv", COHORT, st, nl)), sep = "\t")
  }
  pr <- data.table(cohort = COHORT, stratum = st, n = nrow(dd), raw_spearman = obs)
  for (nm in names(ps <- partial_sets(dd))) pr[[paste0("partial_adj_", nm)]] <- partial_spearman(dd$fixed_projection, yk, ps[[nm]])
  if (all(dd$has_composition)) {
    pr$partial_adj_hepatocyte <- partial_spearman(dd$fixed_projection, yk, dd$Hepatocytes)
    pr$partial_adj_6clr_joint <- partial_spearman(dd$fixed_projection, yk, clr6(dd[, LIN6, with = FALSE]))
  }
  partrows[[st]] <- pr

  ## competitors (CFT-v1 floor c), only with fractions
  if (all(dd$has_composition) && nrow(dd) >= 20L) {
    Lc <- clr6(dd[, LIN6, with = FALSE]); colnames(Lc) <- paste0("CLR_", make.names(LIN6)); nn <- nrow(dd)
    set.seed(CV_SEED); fk <- sample(rep_len(1:10, nn)); dfcv <- data.frame(y = yk, Lc); oof <- numeric(nn)
    for (k in 1:10) { tr <- fk != k; oof[!tr] <- predict(lm(y ~ ., data = dfcv[tr, , drop = FALSE]), newdata = dfcv[!tr, , drop = FALSE]) }
    cands <- list(continuum = dd$fixed_projection, C_libsize = dd$log10_libsize, C_hepatocyte = dd$Hepatocytes, C_lineage6_oof = oof)
    set.seed(BOOT_SEED); IDX <- matrix(sample.int(nn, nn * BOOT_B, replace = TRUE), nrow = nn, ncol = BOOT_B)
    absr <- lapply(cands, function(v) vapply(seq_len(BOOT_B), function(j) {
      i <- IDX[, j]; s <- yk[i]; if (length(unique(s)) < 2L) NA_real_ else abs(sp_rho(v[i], s)) }, numeric(1)))
    rr <- rbindlist(lapply(c("C_libsize", "C_hepatocyte", "C_lineage6_oof"), function(nm) {
      del <- absr$continuum - absr[[nm]]; del <- del[is.finite(del)]; ci <- quantile(del, c(0.025, 0.975), type = 7, names = FALSE)
      data.table(cohort = COHORT, stratum = st, competitor = nm, n = nn, continuum_abs_rho = abs(obs),
                 competitor_rho = sp_rho(cands[[nm]], yk), diff_abs_rho = abs(obs) - abs(sp_rho(cands[[nm]], yk)),
                 boot_diff_median = median(del), boot_ci_lo = ci[1], boot_ci_hi = ci[2], win = ci[1] > 0,
                 boot_p = max(2 * min(mean(del <= 0), mean(del >= 0)), 1 / (length(del) + 1)), n_boot_used = length(del))
    }))
    rr[, boot_p_holm := p.adjust(boot_p, method = "holm")]; comp_rows[[st]] <- rr; rm(IDX, absr); gc()
  }
}

## verdict: primary = S1_all_staged, null A, one-sided alpha 0.05 (addendum S3c)
NR <- rbindlist(nullres, fill = TRUE); CMP <- rbindlist(comp_rows, fill = TRUE); PR <- rbindlist(partrows, fill = TRUE)
prim <- NR[stratum == "S1_all_staged" & null == "A"]
prim[, pass := testable & p_emp_signed < 0.05]
prim[, beats_null_q95 := observed_rho > null_q95]
if (nrow(CMP)) prim[, beats_all_competitors := CMP[stratum == "S1_all_staged", all(win & boot_p_holm <= 0.05)]]
three <- rbind(DISCOVERY[, .(cohort, n, observed_rho = rho, p_emp_signed = p_emp_nullA)], prim[, .(cohort, n, observed_rho, p_emp_signed)])
three[, p_holm_3cohort := p.adjust(p_emp_signed, method = "holm")]
fwrite(NR, file.path(OUT, "matched_null_results.tsv"), sep = "\t")
fwrite(PR, file.path(OUT, "partial_spearman_adjustments.tsv"), sep = "\t")
if (nrow(CMP)) fwrite(CMP, file.path(OUT, "competitor_comparisons.tsv"), sep = "\t")
fwrite(prim, file.path(OUT, "PRIMARY_transfer_verdict.tsv"), sep = "\t")
fwrite(three, file.path(OUT, "three_cohort_descriptive_holm.tsv"), sep = "\t")
fwrite(rbindlist(guard), file.path(OUT, "reproduction_guard.tsv"), sep = "\t")
writeLines(c(sprintf("S3c control-free transport of the frozen 139-gene loading: %s", COHORT), "",
             "COVERAGE:", capture.output(print(cov)), "", "PRIMARY (S1_all_staged, fibrosis, null A, one-sided alpha 0.05):",
             capture.output(print(prim)), "", "ALL NULL CELLS:",
             capture.output(print(NR[, .(stratum, null, n, n_draws_retained, observed_rho, null_median, null_q95, p_emp_signed, tie_ceiling)])),
             "", "PARTIALS:", capture.output(print(PR)), "", "COMPETITORS:", capture.output(print(CMP)),
             "", "WITH THE TWO 2026-08-28 COHORTS (descriptive Holm):", capture.output(print(three))),
           file.path(OUT, "VERDICT.txt"))
ss_finish(OUT)
