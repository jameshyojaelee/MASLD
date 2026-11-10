# SS job S3a (prespecified; PRESPEC_SS.md section 5, families F-SS1, F-SS2, F-SS3; rules fixed in
# PRESPEC_SS_addendum_S3a.md). Held-out replication in SteatoSITE of frozen three-cohort fibrosis-stage effects.
# Forked from H51 (stage1_blinded.R, stage2_outcome.R); slope, null, bootstrap, MDE and per-gene state code unchanged.
# F-SS1 (2 tests, one-sided, Holm): b = OLS slope of G_SS log2FC on frozen D0 log2FC over D0 BH hits kept by the
#   SteatoSITE filterByExpr(group = bin); contrasts F3F4 and F2 vs F0F1. Null: 1,000 permutations of the stage code
#   within batch x sex strata (seed 20260923). Interval: 200-draw bootstrap within bin x sex cells (seed 20261007).
#   Testable if 2.80 x null SD <= 0.30. Controls PC, NC1, NC-perm, NC-quality. Beside: female-only, quality
#   covariates, F3-only (drops F4), composition (ILR file, optional).
# F-SS2 (3 tests, one-sided, Holm): same statistic for adjacent steps S2-S1, S3-S2, S4-S3 of `~ batch + sex + stage4`
#   against frozen Fig 3E F1_to_F2, F2_to_F3, F3_to_F4 hits. Bootstrap within stage4 x sex cells.
# F-SS3: hero trio stage_num (one-sided decrease, Holm 3); IGFBP7 S3-S2 (one-sided increase); 29 Fig 3F late-step
#   marks (one-sided in discovery sign, BH 29); stage-supported sign concordance vs permutation null (descriptive);
#   per-gene states for kept F3F4 hits (H51 definitions, frozen s_c).
# Models: limma voomWithQualityWeights, lmFit, eBayes; `~ factor(batch) + factor(inferred_sex) + <stage term>`;
#   filterByExpr(group = bin) on staged units, then RLE on kept genes (--norm stored keeps the input factors; mock only).
# Unit: one library per genetic individual (S2 include = TRUE), staged (codes 1-4), sex known. Cross-sectional.
# Usage: Rscript s3a_replication.R --input <S2 dir> --stage <stage tsv> --out <new dir> [--strata-column batch]
#        [--norm rle|stored] [--ilr <ilr tsv>] [--block-size 48] [--resume] [--debug] [--allow-unhashed-prespec]
# Checkpointing: all random draws are made first from the prespecified seeds; fits run in blocks saved to
# <out>/checkpoints/. --resume reuses finished blocks (same inputs, scripts, arguments and draws required), so a
# restarted run gives the same numbers as an uninterrupted one. --stop-after-blocks N exits after N new blocks (test only).
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3a_common.R")
t0 <- Sys.time()
A <- parse_args(commandArgs(trailingOnly = TRUE), c("input", "stage", "out"), c("debug", "allow-unhashed-prespec", "resume"),
                list("strata-column" = "batch", norm = "rle", ilr = NA_character_, "block-size" = "48", "stop-after-blocks" = "0"))
stopifnot(A$norm %in% c("rle", "stored"))
DEBUG <- isTRUE(A$debug)
B_PERM <- if (DEBUG) 8L else 1000L
B_BOOT <- if (DEBUG) 4L else 200L
BLOCK <- as.integer(A$`block-size`); STOP_AFTER <- as.integer(A$`stop-after-blocks`)
RUN_ARGS <- A[setdiff(sort(names(A)), c("resume", "stop-after-blocks"))]
if (isTRUE(A$resume) && dir.exists(A$out)) {
  OUT <- normalizePath(A$out)
  old <- jsonlite::read_json(file.path(OUT, "run_args.json"), simplifyVector = TRUE)
  as_txt <- function(L) vapply(L, function(x) if (is.null(x)) "NULL" else paste(x, collapse = ","), "")
  if (!identical(as_txt(old[names(RUN_ARGS)]), as_txt(RUN_ARGS))) stop("--resume: arguments differ from the original run")
  log_msg("RESUME in", OUT)
} else {
  OUT <- new_outdir(A$out)
  jsonlite::write_json(RUN_ARGS, file.path(OUT, "run_args.json"), auto_unbox = TRUE, pretty = TRUE)
}
CK <- file.path(OUT, "checkpoints"); dir.create(CK, showWarnings = FALSE)
notes <- character(); note <- function(...) { m <- paste0(...); notes <<- c(notes, m); log_msg("NOTE:", m) }
log_msg("S3a start. OUT =", OUT, " input =", A$input, " stage =", A$stage, " strata =", A$`strata-column`, " norm =", A$norm,
        " DEBUG =", DEBUG, " B_PERM =", B_PERM, " B_BOOT =", B_BOOT, " cores =", NCORES)
verify_prespec(OUT, isTRUE(A$`allow-unhashed-prespec`))
verify_frozen(OUT)
inp <- c(units = file.path(A$input, "units.tsv"), dge = file.path(A$input, "dge_included.rds"), stage = A$stage,
         ilr = if (!is.na(A$ilr) && file.exists(A$ilr)) A$ilr else NA)
sc <- list.files(file.path(SS, "scripts"), pattern = "^s3a_.*\\.R$", full.names = TRUE)
SHA_IN <- data.table(id = names(inp), path = inp, sha256 = vapply(inp, function(f) if (is.na(f)) NA_character_ else sha256(f), ""))
SHA_SC <- data.table(script = basename(sc), sha256 = vapply(sc, sha256, ""))
if (file.exists(file.path(OUT, "input_sha256.tsv"))) {  # resume: inputs and scripts must be unchanged
  if (!identical(fread(file.path(OUT, "input_sha256.tsv"))$sha256, SHA_IN$sha256)) stop("--resume: input sha256 changed")
  if (!identical(fread(file.path(OUT, "script_sha256.tsv"))$sha256, SHA_SC$sha256)) stop("--resume: script sha256 changed")
}
fwrite(SHA_IN, file.path(OUT, "input_sha256.tsv"), sep = "\t"); fwrite(SHA_SC, file.path(OUT, "script_sha256.tsv"), sep = "\t")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

## ---------- Units, labels, genes ----------
ann <- load_ann()
units <- fread(inp[["units"]])
dge <- readRDS(inp[["dge"]])
if (!"group" %in% names(dge$samples)) dge$samples <- cbind(group = factor(rep(1L, ncol(dge))), dge$samples)  # DGEList subsetting needs it
stopifnot(nrow(dge) == 24115L, setequal(rownames(dge), ann$gene_id_versioned), !anyDuplicated(units$sample_alias))
ann <- ann[match(rownames(dge), gene_id_versioned)]
stopifnot(setequal(colnames(dge), units[include == TRUE, sample_alias]), A$`strata-column` %in% names(units))
st <- read_stage(A$stage)
u <- merge(units[include == TRUE], st, by = "sample_alias", all.x = TRUE)
note("Included libraries: ", nrow(u), "; unstaged (no inferential test): ", sum(is.na(u$stage_code)),
     "; staged with unknown sex (excluded): ", sum(!is.na(u$stage_code) & is.na(u$inferred_sex)))
u <- u[!is.na(stage_code) & !is.na(inferred_sex)]
setorder(u, sample_alias)
u <- add_labels(u, u$stage_code)
u[, stratum := paste(get(A$`strata-column`), inferred_sex, sep = "_")]
u[, cell := paste(bin, inferred_sex, sep = "_")]
u[, cell4 := paste(stage4, inferred_sex, sep = "_")]
n_tab <- dcast(u[, .N, by = .(batch, inferred_sex, stage4)], batch + inferred_sex ~ stage4, value.var = "N", fill = 0)
fwrite(u, file.path(OUT, "units_used.tsv"), sep = "\t")
log_msg("Staged units:", nrow(u), "; bins", paste(names(table(u$bin)), table(u$bin), collapse = " "), "; stage4", paste(names(table(u$stage4)), table(u$stage4), collapse = " "))

dg_all <- dge[, u$sample_alias]
keep_l <- filterByExpr(dg_all, group = u$bin)
KEEP <- rownames(dg_all)[keep_l]
dg <- dg_all[keep_l, , keep.lib.sizes = TRUE]
if (A$norm == "rle") dg <- calcNormFactors(dg, method = "RLE") else note("norm = stored: input normalization factors used (mock reproduction only)")
log_msg("filterByExpr(group = bin) kept", length(KEEP), "of 24115")
iK <- function(gs) match(gs, KEEP)
gname <- function(gs) ann$gene_name[match(gs, ann$gene_id_versioned)]
gid <- function(nm) { g <- ann[gene_name == nm, gene_id_versioned]; stopifnot(length(g) == 1L); g }
base_to_id <- setNames(ann$gene_id_versioned, ann$gene_id_base)

## ---------- Frozen discovery ----------
D0 <- fread(file.path(FROZEN, "h51/D0_discovery_grouped_per_gene.tsv.gz"))
stopifnot(nrow(D0) == 2L * 24115L)
xD0 <- function(cf, gs) D0[contrast == cf][match(gs, gene), logFC]
HC <- lapply(CF, function(cf) D0[contrast == cf & FDR < 0.05, gene])
HCk <- lapply(HC, function(x) intersect(x, KEEP))
SE3 <- fread(file.path(FROZEN, "stage/stage_extension_all_gene_results.tsv"), select = c("gene_id_versioned", "axis", "contrast", "logFC", "FDR"))[axis == "fibrosis_adjacent"]
HS <- lapply(STEPS, function(s) SE3[contrast == s & FDR < 0.05, gene_id_versioned])
HSk <- lapply(HS, function(x) intersect(x, KEEP))
xSE <- function(s, gs) SE3[contrast == STEPS[[s]]][match(gs, gene_id_versioned), logFC]
SC <- unique(fread(file.path(FROZEN, "h51/per_gene_states.tsv.gz"), select = c("contrast", "s_c"))[!is.na(s_c)])
stopifnot(nrow(SC) == 2L); S_C <- setNames(SC$s_c, SC$contrast)
TSV <- fread(file.path(FROZEN, "stage_continuum/transcript_stage_vs_continuum.tsv.gz"))[stage_supported == TRUE]
stopifnot(nrow(TSV) == 11130L)
TSV[, gene := unname(base_to_id[gene_id_base])]
TSVk <- TSV[gene %in% KEEP]
F3F <- fread(file.path(FROZEN, "stage/fig3f.tsv"))[mark == "cell" & dot %in% c("filled", "open") & contrast != "F0_to_F1"]
stopifnot(nrow(F3F) == 29L)
F3F[, gene_id := unname(base_to_id[ensembl])]
F3F[, step := names(STEPS)[match(contrast, STEPS)]]
autos <- paste0("chr", 1:22)
nc1_genes <- HCk$F3F4[ann$chromosome[match(HCk$F3F4, ann$gene_id_versioned)] %in% autos]

# Gene sets carried through the null and the bootstrap: model, coefficient, rows of KEEP, discovery x
SETS <- c(
  lapply(setNames(names(CF), names(CF)), function(cn) list(model = "bin", cf = CF[[cn]], i = iK(HCk[[cn]]), x = xD0(CF[[cn]], HCk[[cn]]))),
  unlist(lapply(names(CF), function(cn) { gs <- HCk[[cn]]; x <- xD0(CF[[cn]], gs)
    setNames(list(list(model = "bin", cf = CF[[cn]], i = iK(gs[x > 0]), x = x[x > 0]), list(model = "bin", cf = CF[[cn]], i = iK(gs[x < 0]), x = x[x < 0])),
             paste0(cn, c("_up", "_down"))) }), recursive = FALSE),
  lapply(setNames(names(STEPS), names(STEPS)), function(s) list(model = "step", cf = s, i = iK(HSk[[s]]), x = xSE(s, HSk[[s]]))),
  list(stage_supported = list(model = "num", cf = "stage_num", i = iK(TSVk$gene), x = TSVk$beta)))
set_stats <- function(cm) unlist(lapply(names(SETS), function(n) { S <- SETS[[n]]; y <- cm[[S$model]][S$i, S$cf]
  setNames(c(ols_slope(y, S$x), conc(y, S$x)), paste0(c("slope.", "conc."), n)) }))

## ---------- Random draws (serial, fixed order) ----------
set.seed(SEED_PERM, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
strata <- sort(unique(u$stratum))
perm <- matrix(NA_integer_, B_PERM, nrow(u), dimnames = list(NULL, u$sample_alias))
for (k in seq_len(B_PERM)) {
  lab <- u$stage_code
  for (s in strata) { i <- which(u$stratum == s); lab[i] <- lab[i][sample.int(length(i))] }
  perm[k, ] <- lab
}
stopifnot(all(apply(perm, 1, function(l) all(table(u$stratum, factor(l, 1:4)) == table(u$stratum, factor(u$stage_code, 1:4))))))
set.seed(SEED_BOOT, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
bootA <- t(vapply(seq_len(B_BOOT), function(k) boot_draw(u$cell), integer(nrow(u))))
bootB <- t(vapply(seq_len(B_BOOT), function(k) boot_draw(u$cell4), integer(nrow(u))))
stopifnot(all(apply(bootA, 1, function(i) all(table(u$cell[i]) == table(u$cell)))), all(apply(bootB, 1, function(i) all(table(u$cell4[i]) == table(u$cell4)))))
fwrite(data.table(draw = seq_len(B_PERM), perm), file.path(OUT, "perm_stage_codes.tsv.gz"), sep = "\t")
fwrite(rbind(data.table(set = "bin_x_sex", draw = rep(seq_len(B_BOOT), each = nrow(u)), position = rep(seq_len(nrow(u)), B_BOOT), index = as.vector(t(bootA))),
             data.table(set = "stage4_x_sex", draw = rep(seq_len(B_BOOT), each = nrow(u)), position = rep(seq_len(nrow(u)), B_BOOT), index = as.vector(t(bootB)))),
       file.path(OUT, "boot_indices.tsv.gz"), sep = "\t")
log_msg("Draws written: strata", length(strata), "perms", B_PERM, "boot", B_BOOT, "x 2")

## ---------- Null fits ----------
base_info <- u[, .(batch, inferred_sex)]
RB <- function(tag, ids, fun, key) run_blocks(tag, ids, fun, key, CK, BLOCK, STOP_AFTER, NCORES)
null_res <- RB("null", seq_len(B_PERM), function(k) {
  info <- add_labels(base_info, perm[k, ])
  cm <- list(bin = fit_qw(dg, design_of(info, "bin"), FALSE)$coefficients,
             step = step_coefs(fit_qw(dg, design_of(info, "stage4"), FALSE)$coefficients),
             num = fit_qw(dg, design_of(info, "stage_num"), FALSE)$coefficients[, "stage_num", drop = FALSE])
  set_stats(cm)
}, function(ids) perm[ids, , drop = FALSE])
NULLS <- as.data.table(do.call(rbind, null_res)); rm(null_res); invisible(gc())
fwrite(cbind(draw = seq_len(B_PERM), NULLS), file.path(OUT, "null_set_statistics.tsv.gz"), sep = "\t")
log_msg("Null fits done", format(Sys.time() - t0))

## ---------- Observed, sensitivity and bootstrap fits ----------
ilr <- NULL
if (is.na(inp[["ilr"]])) note("Composition-adjusted fit skipped: no ILR file (", A$ilr, ")") else {
  il <- fread(inp[["ilr"]]); icol <- grep("^ilr[0-9]+$", names(il), value = TRUE)
  if (!"sample_alias" %in% names(il) || length(icol) == 0L || !all(u$sample_alias %in% il$sample_alias)) {
    note("Composition-adjusted fit skipped: ILR file lacks sample_alias/ilr columns or some staged units") } else {
    ilr <- as.data.table(il[match(u$sample_alias, sample_alias), ..icol]); note("Composition-adjusted fit uses ", length(icol), " ILR coordinates from ", inp[["ilr"]]) } }
qual_ok <- all(QUAL %in% names(u)) && !anyNA(u[, ..QUAL])
if (!qual_ok) note("Quality-covariate fit skipped: ", paste(QUAL, collapse = ", "), " missing or NA for some staged units")
n4 <- sum(u$stage_code == 4L)
if (n4 == 0L) note("No S4 (code 4) units: S4-S3 not estimable; F3-only fit equals the primary")
pack <- function(f, cfs) list(coef = f$coefficients[, cfs, drop = FALSE], tab = gene_tab(f, cfs))
CFS3 <- c("binF2", "binF3F4", if (length(unique(u$inferred_sex)) > 1L) SEXCF)
jobs <- c("bin", "step", "num", "female", "quality", "F3only", "comp", paste0("bootA", seq_len(B_BOOT)), paste0("bootB", seq_len(B_BOOT)))
fitfun <- function(j) {
  if (j == "bin") return(pack(fit_qw(dg, design_of(u, "bin")), CFS3))
  if (j == "step") { des <- design_of(u, "stage4"); v <- voomWithQualityWeights(dg, des, plot = FALSE)
    f <- eBayes(contrasts.fit(lmFit(v, des), step_matrix(des))); return(pack(f, colnames(f$coefficients))) }
  if (j == "num") return(pack(fit_qw(dg, design_of(u, "stage_num")), "stage_num"))
  if (j == "female") { w <- which(u$inferred_sex == "F"); return(pack(fit_qw(dg[, w], design_of(u[w], "bin")), c("binF2", "binF3F4"))) }
  if (j == "quality") { if (!qual_ok) return(NULL); return(pack(fit_qw(dg, design_of(u, "bin", QUAL)), c("binF2", "binF3F4"))) }
  if (j == "F3only") { w <- which(u$stage_code != 4L); return(pack(fit_qw(dg[, w], design_of(u[w], "bin")), c("binF2", "binF3F4"))) }
  if (j == "comp") { if (is.null(ilr)) return(NULL); info <- cbind(u[, .(batch, inferred_sex, bin)], ilr)
    return(pack(fit_qw(dg, design_of(info, "bin", names(ilr))), c("binF2", "binF3F4"))) }
  if (startsWith(j, "bootA")) { idx <- bootA[as.integer(sub("bootA", "", j)), ]; d <- dg[, idx]; colnames(d) <- make.unique(colnames(d))
    return(fit_qw(d, design_of(u[idx], "bin"), FALSE)$coefficients[, CFS3, drop = FALSE]) }
  idx <- bootB[as.integer(sub("bootB", "", j)), ]; d <- dg[, idx]; colnames(d) <- make.unique(colnames(d))
  step_coefs(fit_qw(d, design_of(u[idx], "stage4"), FALSE)$coefficients)
}
res <- c(RB("observed", jobs[1:7], fitfun, function(ids) ids),
         RB("bootA", jobs[grep("^bootA", jobs)], fitfun, function(ids) bootA[as.integer(sub("bootA", "", ids)), , drop = FALSE]),
         RB("bootB", jobs[grep("^bootB", jobs)], fitfun, function(ids) bootB[as.integer(sub("bootB", "", ids)), , drop = FALSE]))
stopifnot(length(res) == length(jobs))
names(res) <- jobs
FIT <- res[1:7]
BOOTA <- lapply(setNames(CFS3, CFS3), function(cf) sapply(res[grep("^bootA", names(res))], function(m) m[, cf]))
BOOTB <- lapply(setNames(names(STEPS), names(STEPS)), function(s) sapply(res[grep("^bootB", names(res))], function(m) m[, s]))
rm(res); invisible(gc())
saveRDS(FIT, file.path(OUT, "observed_fits.rds"))
for (nm in c("bin", "step", "num")) { tb <- FIT[[nm]]$tab; tb[, gene_name := gname(gene)]
  fwrite(tb, file.path(OUT, paste0("gene_", nm, ".tsv.gz")), sep = "\t") }
log_msg("Observed and bootstrap fits done", format(Sys.time() - t0))


## ---------- Slope blocks ----------
stepcm <- sapply(names(STEPS), function(s) if (s %in% colnames(FIT$step$coef)) FIT$step$coef[, s] else rep(NA_real_, length(KEEP)))
CM <- list(bin = FIT$bin$coef, step = stepcm, num = FIT$num$coef)
slope_block <- function(set, cm, boot) {
  S <- SETS[[set]]; y <- cm[[S$model]][S$i, S$cf]; b <- ols_slope(y, S$x)
  nl <- NULLS[[paste0("slope.", set)]]; nc <- NULLS[[paste0("conc.", set)]]
  bb <- if (is.null(boot)) NA_real_ else apply(boot[S$i, , drop = FALSE], 2, ols_slope, x = S$x)
  data.table(n_genes = length(S$i), slope = b, null_mean = mean(nl), null_sd = sd(nl), MDE_slope = MDE_SLOPE_K * sd(nl),
             p_perm_onesided = if (is.finite(b)) (1 + sum(nl >= b)) / (length(nl) + 1) else NA_real_,
             boot_se = sd(bb), ci_norm_lo = b - z975 * sd(bb), ci_norm_hi = b + z975 * sd(bb),
             ci_pct_lo = unname(quantile(bb, 0.025, na.rm = TRUE)), ci_pct_hi = unname(quantile(bb, 0.975, na.rm = TRUE)),
             conc_obs = conc(y, S$x), conc_null_mean = mean(nc), conc_null_sd = sd(nc))
}
point_slope <- function(set, fit) { S <- SETS[[set]]; if (is.null(fit)) NA_real_ else ols_slope(fit$coef[S$i, S$cf], S$x) }

## ---------- F-SS1 ----------
F1 <- rbindlist(lapply(names(CF), function(cn) cbind(data.table(contrast = cn, D0_hits = length(HC[[cn]])), slope_block(cn, CM, BOOTA[[CF[[cn]]]]))))
F1[, testable := is.finite(MDE_slope) & MDE_slope <= 0.30]
F1[, p_holm := p.adjust(ifelse(testable, p_perm_onesided, NA_real_), "holm", n = 2)]
F1[, holm_pass := testable & !is.na(p_holm) & p_holm < 0.05]
F1[, nc_perm_clean := abs(null_mean) <= 0.05]
updn <- rbindlist(lapply(c("F3F4_up", "F3F4_down", "F2_up", "F2_down"), function(s)
  cbind(data.table(contrast = sub("_.*", "", s), analysis = paste(sub(".*_", "", s), "hits")), slope_block(s, CM, BOOTA[[SETS[[s]]$cf]]))))
beside <- rbindlist(lapply(names(CF), function(cn) data.table(contrast = cn,
  analysis = c("female-only (~ batch + bin)", "quality covariates (+ log10 assigned, unique %, intron/exon)", "F3-only (drops code 4 units)", "composition-adjusted (+ ILR)"),
  n_units = c(sum(u$inferred_sex == "F"), if (qual_ok) nrow(u) else NA, sum(u$stage_code != 4L), if (!is.null(ilr)) nrow(u) else NA),
  slope = c(point_slope(cn, FIT$female), point_slope(cn, FIT$quality), point_slope(cn, FIT$F3only), point_slope(cn, FIT$comp)))))
beside <- beside[!(contrast == "F2" & grepl("^F3-only", analysis))]
ncq <- merge(F1[, .(contrast, slope)], beside[grepl("^quality", analysis), .(contrast, quality_slope = slope)], by = "contrast")
ncq[, sensitive_to_rna_quality := is.finite(quality_slope) & quality_slope < 0.5 * slope]

# PC: all kept genes have the discovery sign and >= 4 have nominal p < 0.05 (G_SS F3F4 vs F0F1)
pc <- rbindlist(lapply(names(PCG), function(nm) { g <- gid(nm); k <- g %in% KEEP; tb <- FIT$bin$tab[coef == "binF3F4"][gene == g]
  data.table(gene_name = nm, gene = g, expected_sign = PCG[[nm]], kept = k, D0_F3F4_log2FC = xD0("binF3F4", g),
             G_F3F4_log2FC = if (k) tb$logFC else NA_real_, G_F3F4_p = if (k) tb$P else NA_real_) }))
pc[, sign_ok := kept & sign(G_F3F4_log2FC) == sign(D0_F3F4_log2FC)]
PC_pass <- all(pc[kept == TRUE, sign_ok]) && sum(pc$kept & pc$G_F3F4_p < 0.05, na.rm = TRUE) >= 4L

# NC1: male-vs-female coefficient on D0 F3F4 log2FC over autosomal kept F3F4 hits; fails if the 95% interval excludes 0
if (SEXCF %in% CFS3) {
  ii <- iK(nc1_genes); xn <- xD0("binF3F4", nc1_genes); bn <- ols_slope(FIT$bin$coef[ii, SEXCF], xn)
  bbn <- apply(BOOTA[[SEXCF]][ii, , drop = FALSE], 2, ols_slope, x = xn)
  NC1 <- data.table(n_genes = length(ii), slope = bn, boot_se = sd(bbn), ci_norm_lo = bn - z975 * sd(bbn), ci_norm_hi = bn + z975 * sd(bbn))
  NC1[, fails := ci_norm_lo > 0 | ci_norm_hi < 0]
} else { NC1 <- data.table(n_genes = length(nc1_genes), slope = NA_real_, fails = NA); note("NC1 not estimable: one sex level") }

## ---------- F-SS2 ----------
F2s <- rbindlist(lapply(names(STEPS), function(s) cbind(data.table(step = s, label = STEP_LABEL[[s]], discovery = STEPS[[s]], discovery_hits = length(HS[[s]])),
  slope_block(s, CM, BOOTB[[s]]))))
F2s[, estimable := is.finite(slope)]
F2s[, testable := estimable & is.finite(MDE_slope) & MDE_slope <= 0.30]
F2s[, p_holm := p.adjust(ifelse(testable, p_perm_onesided, NA_real_), "holm", n = 3)]
F2s[, holm_pass := testable & !is.na(p_holm) & p_holm < 0.05]
F2s[, nc_perm_clean := abs(null_mean) <= 0.05]

## ---------- F-SS3 ----------
gstat <- function(tab, cf, g, sign) { r <- tab[coef == cf & gene == g]
  if (nrow(r) == 0L) return(data.table(kept = g %in% KEEP, logFC = NA_real_, SE = NA_real_, t = NA_real_, df = NA_real_, p_onesided = NA_real_))
  data.table(kept = TRUE, logFC = r$logFC, SE = r$SE, t = r$t, df = r$df, p_onesided = p_one(r$t, r$df, sign)) }
hero <- rbindlist(lapply(HERO, function(nm) cbind(data.table(gene_name = nm, test = "stage_num < 0"), gstat(FIT$num$tab, "stage_num", gid(nm), -1))))
hero[, p_holm := p.adjust(p_onesided, "holm", n = 3)]; hero[, pass := !is.na(p_holm) & p_holm < 0.05]
igf <- cbind(data.table(gene_name = "IGFBP7", test = "S3 - S2 > 0 (F2 -> F3)"), gstat(FIT$step$tab, "S3mS2", gid("IGFBP7"), 1))
igf[, pass := !is.na(p_onesided) & p_onesided < 0.05]
marks <- rbindlist(lapply(seq_len(nrow(F3F)), function(r) { m <- F3F[r]
  cbind(m[, .(gene, ensembl, block, discovery_step = contrast, discovery_logFC = logFC, discovery_FDR = FDR, dot)], step = m$step,
        gstat(FIT$step$tab, m$step, m$gene_id, sign(m$logFC))) }))
marks[, q_BH29 := p.adjust(p_onesided, "BH", n = 29)]; marks[, pass := !is.na(q_BH29) & q_BH29 < 0.05]
ssup <- slope_block("stage_supported", CM, NULL)[, .(n_genes, conc_obs, conc_null_mean, conc_null_sd)]
ssup[, p_perm_onesided := (1 + sum(NULLS$conc.stage_supported >= conc_obs)) / (B_PERM + 1)]
ssup[, n_stage_supported := nrow(TSV)]

# Per-gene states (H51 definitions; s_c frozen from H51 S5; BH within kept hits of the contrast)
pg <- rbindlist(lapply(names(CF), function(cn) { cf <- CF[[cn]]; gs <- HCk[[cn]]; x <- xD0(cf, gs)
  tb <- FIT$bin$tab[coef == cf][match(gs, gene)]; q <- p.adjust(tb$P, "BH"); mde <- MDE_GENE_K * tb$SE; e <- S_C[[cn]] * abs(x)
  st <- fifelse(q < 0.05 & sign(tb$logFC) == sign(x), "replicated", fifelse(q < 0.05, "opposite sign",
          fifelse(mde <= e, "not replicated at adequate power", "underpowered")))
  kept <- data.table(gene = gs, gene_name = gname(gs), contrast = cn, D0_log2FC = x, G_log2FC = tb$logFC,
    G_lo = tb$logFC - qt(0.975, tb$df) * tb$SE, G_hi = tb$logFC + qt(0.975, tb$df) * tb$SE, G_SE = tb$SE, G_p = tb$P, G_BH_q = q,
    MDE_g = mde, e_g = e, s_c = S_C[[cn]], state = st)
  miss <- setdiff(HC[[cn]], gs)
  rbind(kept, data.table(gene = miss, gene_name = gname(miss), contrast = cn, D0_log2FC = xD0(cf, miss), state = "not measurable in SteatoSITE"), fill = TRUE) }))
pg[, family_reported := contrast == "F3F4" & isTRUE(F1[contrast == "F3F4", holm_pass]) & PC_pass]
fwrite(pg, file.path(OUT, "per_gene_states.tsv.gz"), sep = "\t")
stc <- pg[, .N, by = .(contrast, state)][order(contrast, state)]

## ---------- Write tables, RESULT, json ----------
fwrite(F1, file.path(OUT, "F_SS1_primary.tsv"), sep = "\t"); fwrite(rbind(updn, beside, fill = TRUE), file.path(OUT, "F_SS1_beside.tsv"), sep = "\t")
fwrite(pc, file.path(OUT, "PC.tsv"), sep = "\t"); fwrite(NC1, file.path(OUT, "NC1.tsv"), sep = "\t"); fwrite(ncq, file.path(OUT, "NC_quality.tsv"), sep = "\t")
fwrite(F2s, file.path(OUT, "F_SS2_steps.tsv"), sep = "\t"); fwrite(rbind(hero, igf, fill = TRUE), file.path(OUT, "F_SS3_hero_igfbp7.tsv"), sep = "\t")
fwrite(marks, file.path(OUT, "F_SS3_fig3f_marks.tsv"), sep = "\t"); fwrite(ssup, file.path(OUT, "F_SS3_stage_supported.tsv"), sep = "\t")
r <- c("# RESULT_S3a (SteatoSITE held-out replication: F-SS1, F-SS2, F-SS3)", "",
  if (DEBUG) "**DEBUG RUN: 8 permutations, 4 bootstraps; numbers are not results.**\n" else NULL,
  paste0("Input: ", A$input, "; stage file: ", A$stage, "; strata column: ", A$`strata-column`, " x inferred_sex; norm: ", A$norm, "."),
  paste0("Unit: one library per genetic individual, staged, sex known. n = ", nrow(u), " (", paste(names(table(u$bin)), table(u$bin), collapse = ", "), "; ",
         paste(names(table(u$stage4)), table(u$stage4), collapse = ", "), "). Genes are the axis of the slope, not replicates. Cross-sectional."),
  paste0("Genes kept by filterByExpr(group = bin): ", length(KEEP), " of 24,115. Discovery D0: frozen H51 (371 disease units, 3 cohorts)."), "",
  "Units by batch x sex x stage:", "", tbl(n_tab), "",
  "## Controls", "", paste0("PC (must pass before F-SS1 is read): ", if (PC_pass) "PASS" else "FAIL"), "", tbl(fmt(pc)), "",
  "NC1 (male-vs-female coefficient on D0 F3F4 log2FC, autosomal kept F3F4 hits; fails if the 95% interval excludes 0):", "", tbl(fmt(NC1)), "",
  "## F-SS1 (2 one-sided permutation tests, Holm, FWER 0.05)", "", tbl(fmt(F1)), "",
  "Beside (descriptive):", "", tbl(fmt(rbind(updn, beside, fill = TRUE))), "", "NC-quality:", "", tbl(fmt(ncq)), "",
  "## F-SS2 (adjacent steps; 3 one-sided permutation tests, Holm)", "", tbl(fmt(F2s)), "",
  "## F-SS3", "", "Hero trio (stage_num, one-sided decrease, Holm 3) and IGFBP7 (S3 - S2, one-sided increase):", "", tbl(fmt(rbind(hero, igf, fill = TRUE))), "",
  paste0("Fig 3F late-step marks (one-sided in discovery sign, BH within 29): ", sum(marks$pass), " of 29 pass; ", sum(is.na(marks$p_onesided)), " not estimable or not kept."), "",
  tbl(fmt(marks[, .(gene, discovery_step, discovery_logFC, dot, kept, logFC, p_onesided, q_BH29, pass)])), "",
  "Stage-supported genes, sign concordance of stage_num with the frozen trend (descriptive; these genes include GSE213621):", "", tbl(fmt(ssup)), "",
  paste0("Per-gene states (kept D0 hits; BH within contrast; family reported only if F-SS1 F3F4 passes Holm and PC passes: ",
         isTRUE(F1[contrast == "F3F4", holm_pass]) && PC_pass, "):"), "", tbl(stc), "",
  "## Notes", "", paste0("- ", notes), "", paste0("Elapsed: ", format(Sys.time() - t0)))
writeLines(r, file.path(OUT, "RESULT_S3a.md"))
jsonlite::write_json(list(debug = DEBUG, n_units = nrow(u), n_kept_genes = length(KEEP), PC_pass = PC_pass, F_SS1 = F1, NC1 = NC1, F_SS2 = F2s,
  hero = hero, igfbp7 = igf, stage_supported = ssup, notes = notes, norm = A$norm, strata_column = A$`strata-column`),
  file.path(OUT, "result.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
log_msg("S3a DONE", format(Sys.time() - t0))
