# SS job S3e (PRESPEC_SS.md section 5; exact statistics in PRESPEC_SS_addendum_S3e.md, written before any SteatoSITE
# count was read). Composition and two earlier hypotheses in the SteatoSITE cohort.
# 1. Cell-type fractions: the release estimator unchanged. counts_all_genes.rds -> GENCODE v49 symbols, rows sharing a
#    symbol summed (scripts/analysis/bulk_refit_followup/prepare_corrected_composition_counts.R), then
#    Analysis/Deconvolution/scripts/12_run_bayesprism.R (InstaPrism, label key celltype) with the release reference
#    reference_human_scalesc_sce.rds, in the music_deconv env, one run for the whole cohort. Stage is not read here.
# 2. 7-part lineage ILR within cohort, H51 coding (h51_common.R LINEAGE, ilr_basis, load_ilr: floor 1e-6). Written
#    per unit for the composition-adjusted fits of S3a/S3d; descriptive: share below floor, R2 of bin on the ILR.
# 3. BF4 P1 (PRESPEC_BF4.md): continuum position = within-cohort signature PC1 of the corrected-774 continuum
#    (02_build_unsupervised_axes.R rules: rowSums > n filter, log2(x + 1), quantile normalisation, PC1 of the 143
#    frozen signature genes present, sign by cosine with the frozen discovery loading), standardized. Frozen BF4
#    markers (30 per type), score = mean within-cohort gene z, standardized. A: y ~ factor(stage) + sex + batch;
#    B: A + bile-duct + fibroblast. Share = (RSS_A - RSS_B) / RSS_A; 2,000 bootstraps within batch; N1 = 1,000
#    expression-decile-matched random pairs (CG0 AveExpr), seed 20261024. S1, S2 and S3 as BF4.
# 4. H38 T2 (PRESPEC_H38.md), only with --ig: IgA / (IgA + IgG + IgM) class-unambiguous constant-exon fragments;
#    beta-binomial, logit(mu) = a + b stage + sex + B-cell score + batch; one-sided b > 0. NC: stage permuted within
#    batch x sex, 1,000 permutations, seed 20260923. PC2 and the H38 sensitivities, plus a fit without F4.
# Unit: one library per genetic individual (S2 include). Cross-sectional. Batch replaces dataset; one cohort, so no
# meta-analysis. Usage: Rscript s3e_composition.R --input <S2 dir> --stage <stage tsv> --out <new dir> [--ig <dir>]
#   --ig: <dir>/<sample_alias>/{ig_fragments.tsv,library_qc.tsv} from the unchanged extract_reads_one_library.py.
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/ss_common.R")
suppressPackageStartupMessages({library(edgeR); library(limma); library(jsonlite); library(parallel)})
A <- ss_args(extra = "ig"); out <- A$out
NCORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
ROOT <- SS_ROOT; COHORT <- "SteatoSITE"
DECONV <- file.path(ROOT, "Analysis/Deconvolution")
REF_SCE <- file.path(DECONV, "reference/reference_human_scalesc_sce.rds")
SCRIPT12 <- file.path(DECONV, "scripts/12_run_bayesprism.R")
MUSIC_ENV <- file.path(ROOT, "RNA-seq/.mamba/music_deconv")
GTF <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
HAC <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-corrected774-20261009T171017Z/reproduction")
BF4 <- file.path(ROOT, "Analysis/Hypotheses/BF4-20261024T010648Z/out-20261024T011020Z")
CG0G <- file.path(ROOT, "Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z/genes.tsv")
GENEMETA <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
HASH <- setNames(c("d7b7598cc42238938403d95f59ce2c69cbb95bd4abc27277d32a98105051f317",
                   "720e01674dac7a8b87175a7e77a2de46d991aa1e9f270cd3f7714a81dda48a79",
                   "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9",
                   "f1e602fdc8bc71b99813f32876c21004a04c8cafdaebffae9cae7eb4da42eb3d",
                   "232bf31cdf3b074ef34a0afea7723b47790e8a3fc6415d59fe9ea97d377bfb47",
                   "8552c0d44a09b7a40872865a3642d5bc4b32c6c7bf484bf119c1dc7cb12e8562",
                   "5ff4804243b9187f14c641196fbfd99d47362ba76fed3c60526c4ea75beeaf94",
                   "49dcaae77323ffb20f5fd8399bb2899cd553e3634755a37b60a133d62549c189"),
                 c(REF_SCE, SCRIPT12, GTF, file.path(HAC, "discovery_signature_model.rds"), file.path(BF4, "markers_main.tsv"),
                   file.path(BF4, "markers_s2.tsv"), CG0G, GENEMETA))
ss_check_sha(HASH)
L <- readLines(file.path(SS_FROZEN, "SHA256SUMS")); invisible(ss_frozen(sub("^[0-9a-f]{64} +", "", L)))
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
base_id <- function(x) sub("\\.[0-9]+$", "", as.character(x))
RES <- list(n_cores = NCORES)

u <- ss_units(A$input, A$stage)
dge <- readRDS(file.path(A$input, "dge_included.rds")); stopifnot(nrow(dge) == 24115L, setequal(colnames(dge), u$sample_alias))
u <- u[match(colnames(dge), sample_alias)]
call <- readRDS(file.path(A$input, "counts_all_genes.rds")); stopifnot(setequal(colnames(call), u$sample_alias))
call <- call[, u$sample_alias]
u[, sexM := as.integer(inferred_sex == "M")]
RES$n_units <- nrow(u); RES$n_staged <- u[!is.na(stage_num), .N]; RES$n_batches <- nlevels(droplevels(u$batch))
log_msg("units", nrow(u), "staged", RES$n_staged)
bterm <- if (nlevels(droplevels(u$batch)) > 1) " + batch" else ""   # reproduction mocks carry one batch level

# ---------------------------------------------------------------- 1. InstaPrism fractions (stage-blind)
gtf_lines <- fread(cmd = paste("zcat", GTF, "| grep -P '^chr.*\\tgene\\t'"), header = FALSE, sep = "\t")
ens2sym <- unique(data.table(ensembl = sub("\\..*$", "", gsub('.*gene_id "([^"]+)".*', "\\1", gtf_lines$V9)),
                             symbol = gsub('.*gene_name "([^"]+)".*', "\\1", gtf_lines$V9)), by = "ensembl")
stopifnot(is.numeric(call), all(call >= 0), all(call == round(call)))
sym <- ens2sym$symbol[match(sub("\\..*$", "", rownames(call)), ens2sym$ensembl)]; mapped <- !is.na(sym)
bt <- data.table(SYMBOL = sym[mapped], as.data.table(call[mapped, , drop = FALSE]))
bt <- bt[, lapply(.SD, sum), by = SYMBOL, .SDcols = setdiff(names(bt), "SYMBOL")]
for (col in setdiff(names(bt), "SYMBOL")) set(bt, j = col, value = as.integer(round(bt[[col]])))
RES$symbol_audit <- list(n_count_rows = nrow(call), n_mapped = sum(mapped), n_symbols = nrow(bt),
                         total_fragments = sum(call), mapped_fragments = sum(call[mapped, ]))
dir.create(file.path(out, "bulk")); cnt_tsv <- file.path(out, "bulk", paste0(COHORT, "_counts.tsv"))
fwrite(bt, cnt_tsv, sep = "\t", quote = FALSE); rm(bt, call); invisible(gc())
prism_dir <- file.path(out, "instaprism")
cmd <- paste0('eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)" && micromamba activate ', MUSIC_ENV,
              ' && export R_LIBS_USER="" R_LIBS_SITE=', MUSIC_ENV, '/lib/R/library LD_LIBRARY_PATH=', MUSIC_ENV, '/lib:${LD_LIBRARY_PATH:-}',
              ' && ', MUSIC_ENV, '/bin/R --no-save --no-restore --slave --file=', SCRIPT12, ' --args ', REF_SCE, ' ', cnt_tsv, ' ',
              prism_dir, ' ', COHORT, ' celltype')
log_msg("InstaPrism start"); rc <- system2("bash", c("-c", shQuote(cmd))); stopifnot(rc == 0); log_msg("InstaPrism done")
pf <- file.path(prism_dir, paste0(COHORT, "_bayesprism_proportions.tsv"))
fr <- read.delim(pf, check.names = FALSE)   # ragged header: row names are sample_alias
CT <- colnames(fr); stopifnot(setequal(rownames(fr), u$sample_alias), length(CT) == 16L)
fr <- as.matrix(fr[u$sample_alias, ])
FR <- data.table(sample_alias = u$sample_alias, analysis_unit_id = u$analysis_unit_id, batch = as.character(u$batch), fr, check.names = FALSE)
fwrite(FR, file.path(out, "fractions.tsv"), sep = "\t")

# ---------------------------------------------------------------- 2. 7-part ILR (H51 coding, copied from h51_common.R)
LINEAGE <- list(hep_chol = c("Hepatocytes", "Cholangiocytes"), endothelial = "Endothelial cells",
                fibroblasts = "Fibroblasts", mac_mono = c("Macrophages", "Mono+mono derived cells"),
                other_myeloid = c("cDC1s", "cDC2s", "pDCs", "Neutrophils", "Basophils"),
                t_nk = c("T cells", "Circulating NK/NKT", "Resident NK"), b_plasma = c("B cells", "Plasma cells"))
ilr_basis <- function(D = 7L) {
  V <- matrix(0, D, D - 1L)
  for (k in seq_len(D - 1L)) { V[seq_len(k), k] <- -1 / sqrt(k * (k + 1)); V[k + 1L, k] <- k / sqrt(k * (k + 1)) }
  V
}
stopifnot(setequal(unlist(LINEAGE), CT))
P <- sapply(LINEAGE, function(cols) rowSums(fr[, cols, drop = FALSE])); P <- P / rowSums(P)
below <- mean(P < 1e-6); P <- pmax(P, 1e-6); P <- P / rowSums(P)
Lg <- log(P); ilr <- (Lg - rowMeans(Lg)) %*% ilr_basis(7L); colnames(ilr) <- paste0("ilr", 1:6)
fwrite(data.table(sample_alias = u$sample_alias, P, ilr), file.path(out, "ilr7.tsv"), sep = "\t")
st_ok <- !is.na(u$stage_num)
RES$ilr <- list(share_below_floor = below, n_staged = sum(st_ok),
                r2_bin_on_ilr = summary(lm(as.integer(u$bin[st_ok]) ~ ilr[st_ok, ]))$r.squared,
                lineage_median = as.list(apply(P, 2, median)))
log_msg("ILR below floor", below)

# ---------------------------------------------------------------- 3. BF4 P1
set.seed(20261024)
cts <- dge$counts; rownames(cts) <- base_id(rownames(cts)); stopifnot(!anyDuplicated(rownames(cts)))
keep <- rowSums(cts) > ncol(cts); X <- normalizeQuantiles(log2(cts[keep, , drop = FALSE] + 1)); dimnames(X) <- dimnames(cts[keep, ])
model <- readRDS(file.path(HAC, "discovery_signature_model.rds"))
common <- unique(base_id(model$common_gene_ids)); sigall <- unique(base_id(model$signature_gene_ids_all))
ld <- model$discovery_loading_oriented; nm <- names(ld); ld <- as.numeric(ld)
names(ld) <- if (is.null(nm)) common else base_id(nm)
used <- common[common %in% rownames(X)]
RES$BF4_signature <- list(n_used = length(used), n_total = length(sigall), coverage = length(used) / length(sigall))
stopifnot(length(used) >= ceiling(0.90 * 145))
pc <- prcomp(t(X[used, ]), center = TRUE, scale. = FALSE, rank. = 2L)
sh <- intersect(used, names(ld)); cs <- sum(pc$rotation[sh, 1] * ld[sh]) / sqrt(sum(pc$rotation[sh, 1]^2) * sum(ld[sh]^2))
axis_raw <- pc$x[, 1] * ifelse(cs < 0, -1, 1); RES$BF4_signature$orientation_cosine <- abs(cs)
Xs <- X[apply(X, 1, sd) > 0, ]; Z <- t(scale(t(Xs)))
MK <- list(main = fread(file.path(BF4, "markers_main.tsv")), s2 = fread(file.path(BF4, "markers_s2.tsv")))
TYPES <- c(bile_duct = "Cholangiocytes", fibroblast = "Fibroblasts", endothelial = "Endothelial cells",
           macrophage = "Macrophages", t_cell = "T cells", hepatocyte = "Hepatocytes")
score <- function(ens) { e <- intersect(ens, rownames(Z)); as.numeric(scale(colMeans(Z[e, , drop = FALSE]))) }
D <- data.table(sample_alias = u$sample_alias, batch = u$batch, inferred_sex = u$inferred_sex, stage_num = u$stage_num,
                y = as.numeric(scale(axis_raw)))
mk_cov <- list()
for (k in names(TYPES)) for (s in names(MK)) {
  e <- MK[[s]][cell_type == TYPES[[k]]]$ensembl; cn <- if (s == "main") k else paste0(k, "_s2")
  D[, (cn) := score(e)]; mk_cov[[cn]] <- c(n_markers = length(e), n_present = sum(e %in% rownames(Z)))
}
RES$BF4_marker_coverage <- mk_cov
partial <- function(d, vars) {
  a <- lm(as.formula(paste0("y ~ factor(stage_num) + inferred_sex", bterm)), data = d)
  b <- lm(as.formula(paste0("y ~ factor(stage_num) + inferred_sex", bterm, " + ", paste(vars, collapse = " + "))), data = d)
  (sum(resid(a)^2) - sum(resid(b)^2)) / sum(resid(a)^2)
}
run_bf4 <- function(d, vars, nboot = 2000) {
  pt <- partial(d, vars)
  bs <- replicate(nboot, partial(d[, .SD[sample.int(.N, .N, replace = TRUE)], by = batch], vars))
  list(share = pt, ci = unname(quantile(bs, c(0.025, 0.975))), n = nrow(d))
}
staged <- D[!is.na(stage_num) & !is.na(inferred_sex)]
RES$BF4_P1 <- run_bf4(staged, c("bile_duct", "fibroblast"))
RES$BF4_P1$reading <- with(RES$BF4_P1, if (share >= 0.5) "mostly portal-tract content" else if (share >= 0.2) "partly" else "little")
log_msg("BF4 P1", RES$BF4_P1$share)
cg <- fread(CG0G, select = c("g", "AveExpr", "signature_gene"))
cg[, sig := as.character(signature_gene) %in% c("TRUE", "True", "true")]
cg[, dec := cut(rank(AveExpr, ties.method = "first"), 10, labels = FALSE)]
allmk <- unique(unlist(lapply(MK, `[[`, "ensembl")))
pool <- cg[g %in% rownames(Z) & !sig & !(g %in% allmk)]
draw_like <- function(ens) { need <- cg[g %in% ens, .N, by = dec]
  unlist(lapply(seq_len(nrow(need)), function(i) sample(pool[dec == need$dec[i]]$g, need$N[i]))) }
eb <- intersect(MK$main[cell_type == "Cholangiocytes"]$ensembl, rownames(Z)); ef <- intersect(MK$main[cell_type == "Fibroblasts"]$ensembl, rownames(Z))
null <- replicate(1000, {
  d <- copy(D); d[, `:=`(r1 = score(draw_like(eb)), r2 = score(draw_like(ef)))]
  partial(d[!is.na(stage_num) & !is.na(inferred_sex)], c("r1", "r2"))
})
fwrite(data.table(null_partial_r2 = null), file.path(out, "bf4_n1_null.tsv"), sep = "\t")
RES$BF4_N1 <- list(null_median = median(null), null_q95 = unname(quantile(null, 0.95)), null_max = max(null),
                   observed_percentile = mean(null < RES$BF4_P1$share))
RES$BF4_S1_all_six <- run_bf4(staged, names(TYPES))
RES$BF4_S2_strict <- run_bf4(staged, c("bile_duct_s2", "fibroblast_s2"))
RES$BF4_noF4 <- run_bf4(staged[stage_num < 4], c("bile_duct", "fibroblast"))
D[, c("P", "H") := { f <- lm(y ~ bile_duct + fibroblast); list(fitted(f), resid(f)) }]
RES$BF4_S3 <- D[!is.na(stage_num), .(n = .N, r_position = cor(y, stage_num, method = "spearman"),
                                     r_P = cor(P, stage_num, method = "spearman"), r_H = cor(H, stage_num, method = "spearman"))]
RES$BF4_S3_r2_position_on_P <- summary(lm(y ~ bile_duct + fibroblast, D))$r.squared
fwrite(D, file.path(out, "bf4_unit_scores.tsv"), sep = "\t")

# ---------------------------------------------------------------- 4. H38 T2 (only with --ig)
if (!is.null(A$ig)) {
  cls_of_set <- function(gs) vapply(strsplit(gs, ",", fixed = TRUE), function(x) {
    cl <- c(A = any(x %in% c("IGHA1", "IGHA2")), G = any(x %in% paste0("IGHG", 1:4)), M = any(x == "IGHM"),
            D = any(x == "IGHD"), E = any(x == "IGHE"))
    if (sum(cl) == 1L) names(cl)[cl] else "X" }, character(1))
  rd <- function(s) {
    d <- file.path(A$ig, s)
    list(f = fread(file.path(d, "ig_fragments.tsv"), colClasses = list(character = c("library", "gene_set", "mapping")))[, sample_alias := s],
         q = fread(file.path(d, "library_qc.tsv"), select = "star_input_reads")[, sample_alias := s])
  }
  rr <- lapply(u$sample_alias, rd)
  frag <- rbindlist(lapply(rr, `[[`, "f")); qc <- rbindlist(lapply(rr, `[[`, "q"))
  setnames(frag, "all_alignments_in_constant_exons", "inside"); frag[, cls := cls_of_set(gene_set)]
  fin <- frag[inside == 1L]
  cnt <- dcast(fin[, .(n = sum(n_fragments)), by = .(sample_alias, cls)], sample_alias ~ cls, value.var = "n", fill = 0)
  for (cc in c("A", "G", "M", "D", "E", "X")) if (!cc %in% names(cnt)) cnt[, (cc) := 0]
  setnames(cnt, c("A", "G", "M"), c("n_A", "n_G", "n_M"))
  per <- merge(u[, .(sample_alias, batch, inferred_sex, sexM, stage_num)], merge(qc, cnt[, .(sample_alias, n_A, n_G, n_M)], by = "sample_alias", all.x = TRUE), by = "sample_alias")
  for (cc in c("n_A", "n_G", "n_M")) per[is.na(get(cc)), (cc) := 0]
  per[, `:=`(k = n_A, n = n_A + n_G + n_M)]
  mixed <- fin[cls == "X"][, { g <- strsplit(gene_set, ",", fixed = TRUE)[[1]]
    .(hasA = any(g %in% c("IGHA1", "IGHA2")), hasG = any(g %in% paste0("IGHG", 1:4)), hasM = any(g == "IGHM"), n_fragments = n_fragments) },
    by = .(sample_alias, gene_set)]
  mixed <- merge(mixed, per[, .(sample_alias, n_A, n_G, n_M)], by = "sample_alias")
  mixed[, den := hasA * n_A + hasG * n_G + hasM * n_M]; mixed <- mixed[den > 0]
  al <- mixed[, .(add_A = sum(n_fragments * hasA * n_A / den), add_all = sum(n_fragments)), by = sample_alias]
  per <- merge(per, al, by = "sample_alias", all.x = TRUE); per[is.na(add_A), `:=`(add_A = 0, add_all = 0)]
  per[, `:=`(k2 = k + add_A, n2 = n + add_all, ig_rate = log10((n + 1) / (star_input_reads / 1e6)))]
  gm <- fread(GENEMETA)
  pick <- function(s) { g <- gm[gene_name == s][order(gene_biotype != "protein_coding")]; g$gene_id[[1]] }
  ids <- vapply(c("CD79A", "MS4A1", "MZB1", "TNFRSF17"), pick, character(1)); stopifnot(all(ids %in% rownames(dge)))
  lc <- cpm(calcNormFactors(DGEList(dge$counts)), log = TRUE, prior.count = 1)[ids, per$sample_alias]
  zc <- function(v) (v - mean(v)) / sd(v)
  per[, Bscore := (zc(lc[ids[["CD79A"]], ]) + zc(lc[ids[["MS4A1"]], ])) / 2]
  per[, Pscore := (zc(lc[ids[["MZB1"]], ]) + zc(lc[ids[["TNFRSF17"]], ])) / 2]
  lbb <- function(k, n, mu, rho) { s <- (1 - rho) / rho; a <- mu * s; b <- (1 - mu) * s
    lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + lbeta(k + a, n - k + b) - lbeta(a, b) }
  rho_of <- function(t) 1e-4 + (1 - 2e-4) * plogis(t)
  bb_negll <- function(par, k, n, X) { p <- ncol(X); mu <- pmin(pmax(plogis(drop(X %*% par[1:p])), 1e-9), 1 - 1e-9)
    -sum(lbb(k, n, mu, rho_of(par[p + 1]))) }
  bb_fit <- function(k, n, X) {   # h38_iga.R bb_fit unchanged
    p <- ncol(X)
    b0 <- tryCatch(suppressWarnings(coef(glm.fit(X, cbind(k, n - k), family = binomial()))), error = function(e) rep(0, p))
    b0[!is.finite(b0)] <- 0; start <- c(b0, qlogis((0.05 - 1e-4) / (1 - 2e-4)))
    o <- tryCatch(optim(start, bb_negll, k = k, n = n, X = X, method = "BFGS", control = list(maxit = 500, reltol = 1e-10)), error = function(e) NULL)
    if (is.null(o)) return(list(ok = FALSE))
    H <- tryCatch(optimHess(o$par, bb_negll, k = k, n = n, X = X), error = function(e) NULL)
    V <- if (is.null(H)) NULL else tryCatch(solve(H), error = function(e) NULL)
    se <- if (is.null(V)) rep(NA_real_, length(o$par)) else sqrt(pmax(diag(V), 0))
    list(ok = o$convergence == 0 && all(is.finite(se[1:p])), par = o$par, se = se, p = p)
  }
  design <- function(d, xcol, covs) {
    X <- cbind(`(Intercept)` = 1, x = d[[xcol]])
    for (cv in covs) if (length(unique(d[[cv]])) > 1L) X <- cbind(X, structure(matrix(d[[cv]], ncol = 1), dimnames = list(NULL, cv)))
    if (nlevels(droplevels(d$batch)) > 1) X <- cbind(X, model.matrix(~ batch, droplevels(d))[, -1, drop = FALSE])
    X
  }
  fit_x <- function(d, xcol, covs, kc = "k", nc = "n", min_n = 10) {
    x <- d[get(nc) >= min_n & !is.na(get(xcol)) & !is.na(sexM)]
    f <- bb_fit(x[[kc]], x[[nc]], design(x, xcol, covs))
    data.table(persons = nrow(x), b = if (f$ok) f$par[2] else NA_real_, se = if (f$ok) f$se[2] else NA_real_,
               rho = if (f$ok) rho_of(f$par[f$p + 1]) else NA_real_, ok = f$ok)
  }
  COVS <- c("sexM", "Bscore")
  t2 <- fit_x(per, "stage_num", COVS)
  t2[, `:=`(z = b / se, p_one_sided = pnorm(-b / se), p_two_sided = 2 * pnorm(-abs(b / se)))]
  stop_rule <- per[!is.na(stage_num) & n >= 10 & !is.na(sexM), .N]
  pc2 <- per[, .(rho_B = cor(ig_rate, Bscore, method = "spearman"), rho_P = cor(ig_rate, Pscore, method = "spearman"))]
  pc2[, `:=`(pass_B = rho_B >= 0.30, pass_P = rho_P >= 0.30)]
  per[, late := as.numeric(stage_num >= 2)]
  sens <- rbindlist(list(
    cbind(analysis = "depth_ge30", fit_x(per, "stage_num", COVS, min_n = 30)),
    cbind(analysis = "proportional_allocation", fit_x(per, "stage_num", COVS, kc = "k2", nc = "n2")),
    cbind(analysis = "plus_plasma_and_ig_rate", fit_x(per, "stage_num", c("sexM", "Bscore", "Pscore", "ig_rate"))),
    cbind(analysis = "F2plus_vs_F0F1", fit_x(per, "late", COVS)),
    cbind(analysis = "without_F4", fit_x(per[is.na(stage_num) | stage_num < 4], "stage_num", COVS))))
  sens[, p_one_sided := pnorm(-b / se)]
  x <- per[n >= 10 & !is.na(stage_num) & !is.na(sexM)]; X0 <- design(x, "stage_num", COVS)
  strata <- split(seq_len(nrow(x)), paste(x$batch, x$sexM))
  perm_one <- function(r) { set.seed(20260923 + r); X <- X0; v <- X[, 2]
    for (s in strata) v[s] <- v[s][sample.int(length(s))]; X[, 2] <- v
    f <- bb_fit(x$k, x$n, X); if (f$ok) f$par[2] / f$se[2] else NA_real_ }
  pz <- unlist(mclapply(seq_len(1000), perm_one, mc.cores = NCORES)); pz <- pz[is.finite(pz)]
  nc <- list(permutations = length(pz), perm_p_one_sided = (1 + sum(pz >= t2$z)) / (1 + length(pz)),
             frac_nominal_two_sided_p_lt_0.05 = mean(2 * pnorm(-abs(pz)) < 0.05))
  nc$null_passes_0.03_to_0.07 <- nc$frac_nominal_two_sided_p_lt_0.05 >= 0.03 && nc$frac_nominal_two_sided_p_lt_0.05 <= 0.07
  RES$H38 <- list(T2 = t2, stop_rule_persons = stop_rule, underpowered = stop_rule < 100, PC2 = pc2, NC = nc, sensitivities = sens,
                  median_iga_fraction = per[n >= 10, median(n_A / n)], depth_median = per[, median(n)])
  fwrite(per, file.path(out, "h38_person_table.tsv"), sep = "\t"); fwrite(sens, file.path(out, "h38_sensitivities.tsv"), sep = "\t")
  log_msg("H38 T2 b", t2$b, "se", t2$se)
}

write_json(RES, file.path(out, "results.json"), auto_unbox = TRUE, pretty = TRUE, digits = 8)
print(RES[setdiff(names(RES), "BF4_marker_coverage")])
ss_finish(out)
