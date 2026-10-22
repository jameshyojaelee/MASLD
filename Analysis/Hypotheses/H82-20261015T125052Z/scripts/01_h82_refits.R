#!/usr/bin/env Rscript
# H82: is the loss of Figure 3E stage DEG calls under library-complexity adjustment larger than adding
# three uninformative covariates causes, and is it confined to F1 versus F0?
# Prespecification: ../PRESPEC_H82.md (sha256 checked below; never edited).
#
# Assay: bulk liver RNA-seq fragment counts, five-cohort 774-unit release (merged_dge.rds, 24,115 genes,
# stored norm factors). Method native to counts: limma voomWithQualityWeights + lmFit + eBayes (defaults),
# exactly the producer's fit_contrast() (scripts/figures/build_pi_stage_extensions.R:142-226): units in
# C-locale sample_id order, design ~ dataset + inferred_sex [+ covariates] + group, BH over all 24,115
# genes per contrast, no fold-change filter. The only change is the added covariate columns.
# Model code follows RUN/jobs/1721/stage_complexity_refit.R (copy in ../reference/), which was not edited.
# Unit of replication: bulk analysis unit (genetic-identity cluster, not a verified person). Genes are
# not replicates; the retained fraction summarizes a unit-level refit over genes.
#
# Stages: A = fits with known reproduction targets (PC-1, stop rule S2); B = blind real arms (NAS
# covariate arms, D3, D4 joint arms, PC-3); C = permutation null (11 contrasts x B). Then summaries,
# T1/T2 tests, readings and RESULT_H82.md.
#
# Flags: --out DIR (must exist, empty of outputs) --seed 20260923 --workers 16 --B 200
#        [--resume EARLIER_RUN_DIR] [--stage all|gate] [--synthetic]
# --stage gate stops after the PC-1 reproduction gate (debug of the reproduction gate on real data).
# --synthetic runs every code path on permuted-label data with planted effects on a 2,000-gene subset
#   (debug only; never a result).

suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(parallel)})
options(digits = 10, scipen = 999, stringsAsFactors = FALSE)
invisible(Sys.setlocale("LC_COLLATE", "C"))

## ---------- arguments and paths ----------
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[[i + 1L]] }
OUT <- getarg("--out")
SEED <- as.integer(getarg("--seed", "20260923"))
WORKERS <- as.integer(getarg("--workers", "16"))
B <- as.integer(getarg("--B", "200"))
RESUME <- getarg("--resume", NA_character_)
STAGE <- getarg("--stage", "all")
SYNTH <- "--synthetic" %in% args
if (is.null(OUT) || !dir.exists(OUT)) stop("--out must name an existing directory")
OUT <- normalizePath(OUT)
if (file.exists(file.path(OUT, "diagnostics.tsv"))) stop("Refusing overwrite: outputs already exist in ", OUT)
stopifnot(STAGE %in% c("all", "gate"), B >= 1L, WORKERS >= 1L)

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
H82 <- dirname(dirname(script_path))
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
PRESPEC_SHA <- "4db60207123c649e29c74587d29dab6e9ad47f967f636d4bb033f97c73b433b8"
PRESPEC_ORIG <- file.path(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/PRESPEC_H82.md")
FITS <- file.path(OUT, "fits"); dir.create(FITS)
NULLDIR <- if (is.na(RESUME)) file.path(OUT, "null") else file.path(normalizePath(RESUME), "null")
if (is.na(RESUME)) dir.create(NULLDIR)

LOG <- file.path(OUT, "h82.log")
logm <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " ", paste0(..., collapse = ""))
  cat(msg, "\n", sep = ""); cat(msg, "\n", sep = "", file = LOG, append = TRUE)
}
DIAG <- file.path(OUT, "diagnostics.tsv")
diag_add <- function(step, item, value, expected = NA, pass = NA) {
  x <- data.table(step = step, item = item, value = as.character(value), expected = as.character(expected), pass = pass)
  fwrite(x, DIAG, sep = "\t", append = file.exists(DIAG), quote = FALSE, na = "NA")
  for (k in seq_len(nrow(x))) logm("DIAG ", step, " | ", x$item[k], " | value=", x$value[k], " | expected=", x$expected[k], " | pass=", x$pass[k])
}
wt <- function(x, name) {
  p <- file.path(OUT, name); if (file.exists(p)) stop("Refusing overwrite: ", p)
  fwrite(x, p, sep = "\t", quote = FALSE, na = "NA"); invisible(p)
}
sha <- function(p) strsplit(system2("sha256sum", shQuote(p), stdout = TRUE)[1], "[[:space:]]+")[[1]][1]
zc <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v)) / sd(v))
t0 <- Sys.time(); el <- function() round(as.numeric(difftime(Sys.time(), t0, units = "secs")))

logm("H82 start; OUT=", OUT, " seed=", SEED, " workers=", WORKERS, " B=", B, " stage=", STAGE,
     " synthetic=", SYNTH, " resume=", RESUME)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
scr <- list.files(file.path(H82, "scripts"), full.names = TRUE)
wt(data.table(script = scr, sha256 = vapply(scr, sha, "")), "script_sha256_at_start.tsv")

## ---------- result writers ----------
RESULT <- file.path(OUT, "RESULT_H82.md")
write_json <- function(x, name) jsonlite::write_json(x, file.path(OUT, name), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
finish_stopped <- function(rule, detail) {
  status <- if (rule == "GATE_ONLY") "debug_gate_only" else "stopped_by_rule"
  lines <- c("# RESULT H82", "",
             if (SYNTH) "SYNTHETIC DEBUG RUN: permuted labels with planted effects. Not a result." else NULL,
             paste0("Verdict: ", if (rule == "GATE_ONLY") "debug run ended after the PC-1 reproduction gate (no outcome computed)."
                    else paste0("STOPPED by rule ", rule, ". No H82 number enters the paper; review numbers stay review-only.")),
             "", paste0("Detail: ", detail), "",
             paste0("Run directory: ", OUT), paste0("Diagnostics: ", DIAG))
  writeLines(lines, RESULT)
  write_json(list(status = status, rule = rule, detail = detail, synthetic = SYNTH), "verdict.json")
  logm("FINISH ", status, " ", rule, ": ", detail)
  quit(save = "no", status = 0L)
}

## ---------- diagnostic 1: input hashes (S1) ----------
inp <- fread(file.path(H82, "scripts/inputs_sha256.tsv"))
pre_copy <- file.path(H82, "PRESPEC_H82.md")
pre_txt <- paste(readLines(pre_copy, warn = FALSE), collapse = "\n")
h_copy <- sha(pre_copy); h_orig <- sha(PRESPEC_ORIG)
diag_add("1_sha256", "PRESPEC_H82.md (run-dir copy)", h_copy, PRESPEC_SHA, h_copy == PRESPEC_SHA)
diag_add("1_sha256", "PRESPEC_H82.md (plans original)", h_orig, PRESPEC_SHA, h_orig == PRESPEC_SHA)
inp[, quoted_in_prespec := vapply(sha256, function(h) grepl(substr(h, 1, 8), pre_txt, fixed = TRUE), logical(1))]
inp[, observed := vapply(path, sha, "")]
inp[, pass := observed == sha256 & quoted_in_prespec]
diag_add("1_sha256", paste0(inp$role, " ", inp$path), inp$observed, inp$sha256, inp$pass)
wt(inp, "inputs_sha256_check.tsv")
if (h_copy != PRESPEC_SHA || h_orig != PRESPEC_SHA || !all(inp$pass)) {
  finish_stopped("S1", paste("sha256 mismatch:", paste(c(if (h_copy != PRESPEC_SHA) "prespec copy", if (h_orig != PRESPEC_SHA) "prespec original",
                                                      inp[pass == FALSE, role]), collapse = ", ")))
}
P <- setNames(inp$path, inp$role)

## ---------- load inputs ----------
d <- readRDS(P[["count_matrix"]])
stopifnot(inherits(d, "DGEList"), identical(dim(d), c(24115L, 774L)))
rel <- fread(P[["released_stage_results"]], select = c("gene_id_versioned", "gene_name", "logFC", "t", "P.Value", "FDR", "axis", "contrast"))
man <- fread(P[["unit_roster_per_contrast"]])
stg <- fread(P[["fibrosis_roster_1721"]])
cx <- fread(P[["complexity_covariates"]])
r1721 <- fread(P[["roster_1721_z"]])
logm("inputs loaded; ", el(), " s")

CT <- data.table(
  contrast = c("F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0", "F1_to_F2", "F2_to_F3", "F3_to_F4",
               "NAS1_2_vs_NAS0", "NAS3_4_vs_NAS0", "NAS5_8_vs_NAS0", "NAS5_8_vs_NAS3_4", "NAS3_4_vs_NAS1_2"),
  rel_axis = c(rep("fibrosis", 4), rep("fibrosis_adjacent", 3), rep("NAS", 3), rep("NAS_adjacent", 2)),
  roster = c(rep("fibrosis", 7), rep("NAS", 5)),
  ref = c("F0", "F0", "F0", "F0", "F1", "F2", "F3", "NAS0", "NAS0", "NAS0", "NAS3-4", "NAS1-2"),
  cmp = c("F1", "F2", "F3", "F4", "F2", "F3", "F4", "NAS1-2", "NAS3-4", "NAS5-8", "NAS5-8", "NAS3-4"),
  n_ref_exp = c(103L, 103L, 103L, 103L, 102L, 88L, 69L, 46L, 46L, 46L, 117L, 65L),
  n_cmp_exp = c(102L, 88L, 69L, 25L, 88L, 69L, 25L, 65L, 117L, 167L, 167L, 117L),
  released_exp = c(1338L, 4266L, 6190L, 7043L, 398L, 512L, 281L, 5338L, 9337L, 9841L, 4305L, 3L),
  null = c(rep(TRUE, 11), FALSE))
CT[, testable_family := fifelse(contrast == "F1_vs_F0", "T1", fifelse(null, "T2", "none (G2: 3 released calls)"))]
COVJ <- c("z_loglib", "z_gini", "z_nf")

## ---------- rosters, covariates, units ----------
mk_roster <- function(axes) {
  m <- unique(man[axis %in% axes, .(sample_id, analysis_unit_id, unit_evidence, dataset, inferred_sex, group)])
  if (anyDuplicated(m$sample_id)) stop("A unit has more than one group on axes ", paste(axes, collapse = ","))
  m <- merge(m, cx[, .(sample_id, lib_size, gini, norm_factor, group_binary, cx_sex = inferred_sex)], by = "sample_id", all.x = TRUE)
  if (anyNA(m$lib_size)) stop("Covariates missing for some units")
  m[, `:=`(z_loglib = zc(log2(lib_size), dataset), z_gini = zc(gini, dataset), z_nf = zc(norm_factor, dataset),
           r_loglib = log2(lib_size), r_gini = gini, r_nf = norm_factor,
           control_status = as.integer(group_binary == "Control"))]
  setorder(m, sample_id)
  m
}
ros <- list(fibrosis = mk_roster(c("fibrosis", "fibrosis_adjacent")), NAS = mk_roster(c("NAS", "NAS_adjacent")))

## ---------- diagnostic 2: unit counts (S3) ----------
n_exp <- rbind(
  data.table(roster = "fibrosis", dataset = "GSE130970", group = paste0("F", 0:4), n = c(23L, 27L, 9L, 10L, 2L)),
  data.table(roster = "fibrosis", dataset = "GSE135251", group = paste0("F", 0:4), n = c(45L, 45L, 52L, 51L, 14L)),
  data.table(roster = "fibrosis", dataset = "GSE162694", group = paste0("F", 0:4), n = c(35L, 30L, 27L, 8L, 9L)),
  data.table(roster = "NAS", dataset = "GSE130970", group = c("NAS0", "NAS1-2", "NAS3-4", "NAS5-8"), n = c(4L, 13L, 32L, 22L)),
  data.table(roster = "NAS", dataset = "GSE135251", group = c("NAS0", "NAS1-2", "NAS3-4", "NAS5-8"), n = c(10L, 31L, 61L, 105L)),
  data.table(roster = "NAS", dataset = "GSE162694", group = c("NAS0", "NAS1-2", "NAS3-4", "NAS5-8"), n = c(32L, 21L, 24L, 40L)))
n_obs <- rbindlist(lapply(names(ros), function(r) ros[[r]][, .(n = .N), by = .(dataset, group)][, roster := r]))
nc <- merge(n_exp, n_obs, by = c("roster", "dataset", "group"), all = TRUE, suffixes = c("_exp", "_obs"))
nc[, pass := !is.na(n_exp) & !is.na(n_obs) & n_exp == n_obs]
ctrl_exp <- data.table(roster = c("fibrosis", "fibrosis", "fibrosis", "fibrosis", "NAS", "NAS", "NAS", "NAS"),
                       group = c("F0", "F0", "F1", "F2", "NAS0", "NAS0", "NAS1-2", "NAS1-2"),
                       group_binary = c("Control", "Disease", "Control", "Control", "Control", "Disease", "Control", "Disease"),
                       n = c(14L, 89L, 1L, 1L, 45L, 1L, 2L, 63L))
ctrl_obs <- rbindlist(lapply(names(ros), function(r) ros[[r]][, .(n = .N), by = .(group, group_binary)][, roster := r]))
cc2 <- merge(ctrl_exp, ctrl_obs, by = c("roster", "group", "group_binary"), all.x = TRUE, suffixes = c("_exp", "_obs"))
cc2[, pass := !is.na(n_obs) & n_exp == n_obs]
units <- setNames(lapply(seq_len(nrow(CT)), function(i) {
  x <- ros[[CT$roster[i]]][group %in% c(CT$ref[i], CT$cmp[i])]; setorder(x, sample_id); x }), CT$contrast)
uc <- rbindlist(lapply(seq_len(nrow(CT)), function(i) {
  u <- units[[i]]; m <- man[axis == CT$rel_axis[i] & contrast == CT$contrast[i]]
  data.table(contrast = CT$contrast[i], n_ref = sum(u$group == CT$ref[i]), n_cmp = sum(u$group == CT$cmp[i]),
             n_ref_exp = CT$n_ref_exp[i], n_cmp_exp = CT$n_cmp_exp[i], n_cohorts = uniqueN(u$dataset),
             min_group_n = min(table(u$group)),
             manifest_match = setequal(paste(u$sample_id, u$group), paste(m$sample_id, m$group)),
             sex_match_cx = all(u$inferred_sex == u$cx_sex))
}))
uc[, pass := n_ref == n_ref_exp & n_cmp == n_cmp_exp & n_cohorts == 3L & min_group_n >= 20L & manifest_match & sex_match_cx]
fib_match <- setequal(paste(ros$fibrosis$sample_id, ros$fibrosis$group), paste(stg$sample_id, stg$fib_group))
diag_add("2_units", paste0(nc$roster, " ", nc$dataset, " ", nc$group), nc$n_obs, nc$n_exp, nc$pass)
diag_add("2_units", paste0(cc2$roster, " ", cc2$group, " ", cc2$group_binary), cc2$n_obs, cc2$n_exp, cc2$pass)
diag_add("2_units", paste0(uc$contrast, " ref/cmp, 3 cohorts, >=20 per group, manifest and sex agree"),
         paste0(uc$n_ref, "/", uc$n_cmp, " cohorts=", uc$n_cohorts, " min=", uc$min_group_n),
         paste0(uc$n_ref_exp, "/", uc$n_cmp_exp), uc$pass)
diag_add("2_units", "roster sizes fibrosis/NAS; fibrosis roster equals stage_sample_manifest.tsv",
         paste0(nrow(ros$fibrosis), "/", nrow(ros$NAS), " match=", fib_match), "387/395 match=TRUE",
         nrow(ros$fibrosis) == 387L && nrow(ros$NAS) == 395L && fib_match)
wt(rbind(nc[, .(check = "dataset_x_group", roster, dataset, group, group_binary = NA_character_, n_obs, n_exp, pass)],
         cc2[, .(check = "control_status", roster, dataset = NA_character_, group, group_binary, n_obs, n_exp, pass)]), "unit_counts.tsv")
wt(uc, "unit_counts_per_contrast.tsv")
if (!all(nc$pass) || !all(cc2$pass) || !all(uc$pass) || nrow(ros$fibrosis) != 387L || nrow(ros$NAS) != 395L || !fib_match) {
  finish_stopped("S3", "unit counts differ from the prespecified n tables (see unit_counts.tsv)")
}

## released calls (G2) with released rows aligned to the count-matrix gene order
REL <- setNames(lapply(seq_len(nrow(CT)), function(i) {
  x <- rel[axis == CT$rel_axis[i] & contrast == CT$contrast[i]][match(rownames(d), gene_id_versioned)]
  if (anyNA(x$gene_id_versioned) || nrow(x) != 24115L) stop("released table join failed for ", CT$contrast[i]); x }), CT$contrast)
rel_calls <- vapply(REL, function(x) sum(x$FDR < 0.05), 1L)
diag_add("2_units", paste0("G2 released BH calls ", CT$contrast), rel_calls, CT$released_exp, rel_calls == CT$released_exp)

## ---------- diagnostic 3: covariate recomputation (S5) ----------
cnt <- d$counts; storage.mode(cnt) <- "double"
lib <- colSums(cnt)
gini_lorenz <- function(x) { x <- sort(x); n <- length(x); L <- cumsum(x) / sum(x); 1 - (2 * sum(L) - 1) / n }
gl <- apply(cnt, 2, gini_lorenz)
cxm <- cx[match(colnames(d), sample_id)]
rec <- data.table(covariate = c("lib_size (colSums)", "lib_size (stored d$samples$lib.size vs colSums)", "gini (Lorenz, 4201 formula)", "norm_factor (stored)"),
                  max_abs_diff = c(max(abs(lib - cxm$lib_size)), max(abs(d$samples$lib.size - lib)),
                                   max(abs(gl - cxm$gini)), max(abs(d$samples$norm.factors - cxm$norm_factor))))
rec[, pass := max_abs_diff < 1e-10]
diag_add("3_covariate_recompute", rec$covariate, signif(rec$max_abs_diff, 4), "< 1e-10", rec$pass)
wt(rec, "covariate_recompute.tsv")
S5_cov <- !all(rec$pass)
rm(cnt); invisible(gc())
zchk <- merge(ros$fibrosis[, .(sample_id, z_loglib, z_gini, z_nf)], r1721[, .(sample_id, a = z_loglib, b = z_gini, c = z_nf)], by = "sample_id")
diag_add("3_covariate_recompute", "fibrosis-roster z vs 1721 roster z (max abs diff over 3 columns)",
         signif(max(abs(c(zchk$z_loglib - zchk$a, zchk$z_gini - zchk$b, zchk$z_nf - zchk$c))), 4), "< 1e-10 (not a stop rule)",
         nrow(zchk) == 387L && max(abs(c(zchk$z_loglib - zchk$a, zchk$z_gini - zchk$b, zchk$z_nf - zchk$c))) < 1e-10)

## ---------- synthetic debug data (never a result) ----------
if (SYNTH) {
  set.seed(1)
  named <- c("COL1A1", "COL1A2", "COL3A1", "LUM", "THBS2", "EFEMP1", "FBLN5", "LOX", "GNMT", "MAT1A", "CYP2C19")
  ids_named <- unique(REL[[1]][gene_name %in% named, gene_id_versioned])
  pool <- setdiff(rownames(d)[rowSums(d$counts) > 5000], ids_named)
  keep <- c(ids_named, sample(pool, 2000L - length(ids_named)))
  keep <- rownames(d)[rownames(d) %in% keep]
  d <- d[keep, , keep.lib.sizes = TRUE]
  REL <- lapply(REL, function(x) x[match(rownames(d), gene_id_versioned)])
  # permute count columns within dataset (breaks every label-expression link), then plant stage effects
  cm <- d$counts; ds <- d$samples$dataset
  if (is.null(ds)) ds <- cx[match(colnames(d), sample_id), dataset]
  for (g in unique(ds)) { j <- which(ds == g); cm[, j] <- cm[, j[sample.int(length(j))]] }
  lev_f <- setNames(as.integer(sub("F", "", ros$fibrosis$group)), ros$fibrosis$sample_id)
  lev_n <- setNames(match(ros$NAS$group, c("NAS0", "NAS1-2", "NAS3-4", "NAS5-8")) - 1L, ros$NAS$sample_id)
  gf <- sample(setdiff(keep, ids_named), 250); gn <- sample(setdiff(keep, c(ids_named, gf)), 250)
  lf <- lev_f[colnames(cm)]; lf[is.na(lf)] <- 0L; ln <- lev_n[colnames(cm)]; ln[is.na(ln)] <- 0L
  cm[gf, ] <- round(sweep(cm[gf, , drop = FALSE], 2, 1.25^lf, `*`)); cm[gn, ] <- round(sweep(cm[gn, , drop = FALSE], 2, 1.25^ln, `*`))
  d$counts <- cm
  logm("SYNTHETIC: ", nrow(d), " genes; permuted columns within dataset; planted 250 fibrosis and 250 NAS genes")
}
NG <- nrow(d)

## ---------- fitting helpers ----------
build_X <- function(u, ref, cmp, covdf) {
  info <- data.frame(dataset = droplevels(factor(u$dataset)), inferred_sex = droplevels(factor(u$inferred_sex)), row.names = u$sample_id)
  for (nm in names(covdf)) info[[nm]] <- covdf[[nm]]
  info$group <- factor(u$group, levels = c(ref, cmp))
  rhs <- c("dataset", "inferred_sex", names(covdf), "group")
  model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), info)
}
fit_X <- function(u, X, cmp) {
  dge <- d[, match(u$sample_id, colnames(d))]
  fit <- eBayes(lmFit(voomWithQualityWeights(dge, X, plot = FALSE), X))
  tt <- topTable(fit, coef = paste0("group", cmp), number = Inf, sort.by = "none")
  g <- data.table(gene_id_versioned = rownames(tt), logFC = tt$logFC, t = tt$t, P.Value = tt$P.Value,
                  FDR = p.adjust(tt$P.Value, method = "BH"))
  stopifnot(identical(g$gene_id_versioned, rownames(d)))
  covn <- setdiff(colnames(X), c("(Intercept)", grep("^dataset|^inferred_sex|^group", colnames(X), value = TRUE)))
  covbh <- vapply(covn, function(cc) sum(p.adjust(fit$p.value[, cc], "BH") < 0.05), 1L)
  list(g = g, covbh = covbh)
}
vs_rel <- function(g, rc) {
  sig <- g$FDR < 0.05; R <- rc$FDR < 0.05
  data.table(calls = sum(sig), up = sum(sig & g$logFC > 0), down = sum(sig & g$logFC < 0), released = sum(R),
             retained = sum(sig & R), retained_frac = if (sum(R) > 0) sum(sig & R) / sum(R) else NA_real_,
             new_calls = sum(sig & !R), median_logFC_ratio_on_released = if (sum(R) > 0) median(g$logFC[R] / rc$logFC[R]) else NA_real_,
             max_abs_logFC_diff_vs_released = max(abs(g$logFC - rc$logFC)))
}
CI <- function(ct) match(ct, CT$contrast)

## covariate data frames per arm
arm_cov <- function(arm, u) {
  switch(arm,
    base = data.frame(row.names = seq_len(nrow(u)))[, FALSE, drop = FALSE],
    lib = data.frame(z_loglib = u$z_loglib), gini = data.frame(z_gini = u$z_gini), nf = data.frame(z_nf = u$z_nf),
    joint = data.frame(z_loglib = u$z_loglib, z_gini = u$z_gini, z_nf = u$z_nf),
    joint_raw = data.frame(r_loglib = u$r_loglib, r_gini = u$r_gini, r_nf = u$r_nf),
    joint_zsub = data.frame(y_loglib = zc(log2(u$lib_size), u$dataset), y_gini = zc(u$gini, u$dataset), y_nf = zc(u$norm_factor, u$dataset)),
    base_diseaseF0 = data.frame(row.names = seq_len(nrow(u)))[, FALSE, drop = FALSE],
    base_ctrlterm = data.frame(control_status = u$control_status),
    joint_diseaseF0 = data.frame(z_loglib = u$z_loglib, z_gini = u$z_gini, z_nf = u$z_nf),
    joint_ctrlterm = data.frame(z_loglib = u$z_loglib, z_gini = u$z_gini, z_nf = u$z_nf, control_status = u$control_status),
    stop("unknown arm ", arm))
}
arm_units <- function(arm, ct) {
  u <- units[[ct]]
  if (arm %in% c("base_diseaseF0", "joint_diseaseF0")) u <- u[!(group == "F0" & group_binary == "Control")]
  u
}
make_task <- function(arm, ct, stage) list(id = paste0(arm, "__", ct), arm = arm, contrast = ct, stage = stage)

run_real <- function(tk) {
  i <- CI(tk$contrast); u <- arm_units(tk$arm, tk$contrast)
  covdf <- if (!is.null(tk$covdf)) tk$covdf else arm_cov(tk$arm, u)
  X <- build_X(u, CT$ref[i], CT$cmp[i], covdf)
  base_row <- data.table(arm = tk$arm, contrast = tk$contrast, stage = tk$stage, n_units = nrow(u),
                         n_ref = sum(u$group == CT$ref[i]), n_cmp = sum(u$group == CT$cmp[i]),
                         covariates = paste(names(covdf), collapse = "+"), design_cols = ncol(X), design_rank = qr(X)$rank,
                         design_kappa = kappa(X, exact = TRUE))
  if (qr(X)$rank < ncol(X)) return(data.table(base_row, status = "not_estimable_rank_deficient"))
  ts <- Sys.time(); f <- fit_X(u, X, CT$cmp[i])
  fwrite(f$g, file.path(FITS, paste0(tk$id, ".genes.tsv.gz")), sep = "\t")
  s <- vs_rel(f$g, REL[[tk$contrast]])
  row <- data.table(base_row, status = "fitted", s,
                    cov_bh_sig = paste0(names(f$covbh), ":", f$covbh, collapse = ";"),
                    fit_s = round(as.numeric(difftime(Sys.time(), ts, units = "secs"))))
  fwrite(row, file.path(FITS, paste0(tk$id, ".summary.tsv")), sep = "\t")
  row
}
run_batch <- function(tasks, label) {
  ord <- order(-vapply(tasks, function(tk) nrow(units[[tk$contrast]]), 1L))
  logm("stage ", label, ": ", length(tasks), " fits on ", WORKERS, " workers")
  res <- mclapply(tasks[ord], function(tk) tryCatch(run_real(tk), error = function(e)
    data.table(arm = tk$arm, contrast = tk$contrast, stage = tk$stage, status = paste("ERROR", conditionMessage(e)))),
    mc.cores = WORKERS, mc.preschedule = FALSE)
  bad <- vapply(res, function(r) !is.data.table(r) || grepl("^ERROR", r$status[1]), logical(1))
  if (any(bad)) { print(res[bad]); stop("stage ", label, ": ", sum(bad), " fits failed") }
  out <- rbindlist(res, fill = TRUE)
  logm("stage ", label, " done; ", el(), " s")
  out
}

## ---------- synthetic: released table := synthetic base fits ----------
if (SYNTH) {
  sb <- run_batch(lapply(CT$contrast, function(ct) make_task("base", ct, "synthetic_release")), "synthetic_release")
  REL <- setNames(lapply(CT$contrast, function(ct) {
    g <- fread(file.path(FITS, paste0("base__", ct, ".genes.tsv.gz")))
    g[, gene_name := REL[[ct]][match(g$gene_id_versioned, gene_id_versioned), gene_name]]; g }), CT$contrast)
  for (ct in CT$contrast) file.rename(file.path(FITS, paste0("base__", ct, c(".genes.tsv.gz", ".summary.tsv"))),
                                     file.path(FITS, paste0("synthrelease__", ct, c(".genes.tsv.gz", ".summary.tsv"))))
  rel_calls <- vapply(REL, function(x) sum(x$FDR < 0.05), 1L)
  logm("SYNTHETIC released calls: ", paste(names(rel_calls), rel_calls, collapse = "; "))
}

## ---------- permutations and PC-3 draws (seeded once, fixed order) ----------
PERMF <- file.path(if (is.na(RESUME)) OUT else normalizePath(RESUME), "perm_index.rds")
if (is.na(RESUME)) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(SEED)
  perms <- list()
  for (ct in CT[null == TRUE, contrast]) {
    u <- units[[ct]]; dss <- sort(unique(u$dataset), method = "radix")
    Pm <- matrix(NA_integer_, nrow(u), B)
    for (b in seq_len(B)) {
      p <- seq_len(nrow(u))
      for (g in dss) { idx <- which(u$dataset == g); stopifnot(length(idx) > 1L); p[idx] <- sample(idx) }
      Pm[, b] <- p
    }
    perms[[ct]] <- Pm
  }
  plant_e <- matrix(NA_real_, nrow(units[["F1_vs_F0"]]), 20L)
  for (k in 1:20) plant_e[, k] <- rnorm(nrow(units[["F1_vs_F0"]]))
  saveRDS(list(seed = SEED, B = B, rngkind = RNGkind(), perms = perms, plant_e = plant_e,
               unit_order = lapply(units, `[[`, "sample_id")), PERMF)
  writeLines(sha(PERMF), file.path(OUT, "perm_index.sha256"))
} else {
  saved <- readLines(file.path(normalizePath(RESUME), "perm_index.sha256"))
  if (sha(PERMF) != saved) stop("--resume: perm_index.rds sha256 differs from its saved value")
  pi_ <- readRDS(PERMF); if (pi_$B != B || pi_$seed != SEED) stop("--resume: B or seed differs")
  if (!identical(pi_$unit_order, lapply(units, `[[`, "sample_id"))) stop("--resume: unit order differs")
  perms <- pi_$perms; plant_e <- pi_$plant_e
  file.copy(PERMF, file.path(OUT, "perm_index_resumed_from.rds"))
}
diag_add("seed", "perm_index.rds sha256 (written before any fit)", sha(PERMF), NA, NA)

## ---------- diagnostic 4: covariate balance ----------
bal <- rbindlist(lapply(seq_len(nrow(CT)), function(i) {
  u <- units[[i]]
  rbindlist(lapply(COVJ, function(v) rbindlist(lapply(c(sort(unique(u$dataset)), "all (within-cohort z)"), function(g) {
    q <- if (g == "all (within-cohort z)") u else u[dataset == g]
    a <- q[group == CT$ref[i]][[v]]; b2 <- q[group == CT$cmp[i]][[v]]
    data.table(contrast = CT$contrast[i], covariate = v, dataset = g, n_ref = length(a), n_cmp = length(b2),
               smd = (mean(b2) - mean(a)) / sqrt((var(a) + var(b2)) / 2), wilcox_p = suppressWarnings(wilcox.test(b2, a)$p.value))
  }))))
}))
r2 <- rbindlist(lapply(seq_len(nrow(CT)), function(i) {
  u <- copy(units[[i]]); u[, gi := as.integer(group == CT$cmp[i])]
  data.table(contrast = CT$contrast[i], r2_dataset_plus_3z = summary(lm(gi ~ dataset + z_loglib + z_gini + z_nf, u))$r.squared,
             r2_dataset_only = summary(lm(gi ~ dataset, u))$r.squared)
}))
wt(bal, "covariate_balance.tsv"); wt(r2, "covariate_balance_r2.tsv")
diag_add("4_balance", paste0(bal$contrast, " ", bal$covariate, " ", bal$dataset, " SMD | Wilcoxon p"),
         paste0(round(bal$smd, 4), " | ", signif(bal$wilcox_p, 3)), "descriptive", NA)
diag_add("4_balance", paste0(r2$contrast, " R2 group ~ dataset + 3 z (dataset only)"),
         paste0(round(r2$r2_dataset_plus_3z, 4), " (", round(r2$r2_dataset_only, 4), ")"), "descriptive", NA)

## ---------- diagnostic 5: design rank (real arms and null replicates) ----------
ARMS_MAIN <- c("base", "lib", "gini", "nf", "joint")
real_specs <- c(lapply(CT$contrast, function(ct) lapply(c(ARMS_MAIN, "joint_raw", "joint_zsub"), function(a) c(a, ct))),
                list(list(c("base_diseaseF0", "F1_vs_F0"), c("base_ctrlterm", "F1_vs_F0"),
                          c("joint_diseaseF0", "F1_vs_F0"), c("joint_ctrlterm", "F1_vs_F0"))))
real_specs <- unlist(real_specs, recursive = FALSE)
dr <- rbindlist(lapply(real_specs, function(s) {
  i <- CI(s[2]); u <- arm_units(s[1], s[2]); X <- build_X(u, CT$ref[i], CT$cmp[i], arm_cov(s[1], u))
  data.table(arm = s[1], contrast = s[2], n_units = nrow(u), cols = ncol(X), rank = qr(X)$rank, kappa = kappa(X, exact = TRUE))
}))
null_rank <- rbindlist(lapply(CT[null == TRUE, contrast], function(ct) {
  i <- CI(ct); u <- units[[ct]]
  rk <- vapply(seq_len(B), function(b) { X <- build_X(u, CT$ref[i], CT$cmp[i], as.data.frame(u[perms[[ct]][, b], ..COVJ])); qr(X)$rank == ncol(X) }, logical(1))
  data.table(contrast = ct, B = B, rank_deficient = sum(!rk))
}))
wt(dr, "design_rank_real_arms.tsv"); wt(null_rank, "design_rank_null.tsv")
diag_add("5_design_rank", paste0(dr$arm, "__", dr$contrast, " rank/cols (kappa)"), paste0(dr$rank, "/", dr$cols, " (", signif(dr$kappa, 4), ")"),
         "full rank", dr$rank == dr$cols)
diag_add("5_design_rank", paste0("null ", null_rank$contrast, " rank-deficient replicates of ", B), null_rank$rank_deficient, "<= 10 (S4)",
         null_rank$rank_deficient <= 10L)

## ---------- stage A: fits with known reproduction targets ----------
fib7 <- CT$contrast[1:7]
tasksA <- c(lapply(CT$contrast, function(ct) make_task("base", ct, "A")),
            unlist(lapply(c("lib", "gini", "nf", "joint"), function(a) lapply(fib7, function(ct) make_task(a, ct, "A"))), recursive = FALSE),
            list(make_task("joint_raw", "F1_vs_F0", "A"), make_task("joint_zsub", "F1_vs_F0", "A"),
                 make_task("base_diseaseF0", "F1_vs_F0", "A"), make_task("base_ctrlterm", "F1_vs_F0", "A")))
resA <- run_batch(tasksA, "A")

## ---------- diagnostic 6: PC-1 reproduction (S2) ----------
tg <- rbindlist(list(
  data.table(arm = "base", contrast = CT$contrast, metric = "calls", target = CT$released_exp, source = "release (prespec n table)"),
  data.table(arm = "joint", contrast = fib7, metric = "calls", target = c(104, 2386, 5564, 6025, 769, 873, 166), source = "1721"),
  data.table(arm = "lib", contrast = c("F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0", "F3_to_F4"), metric = "calls", target = c(1350, 3931, 6064, 7155, 643), source = "1721"),
  data.table(arm = "gini", contrast = c("F1_vs_F0", "F2_vs_F0", "F4_vs_F0", "F3_to_F4"), metric = "calls", target = c(485, 3728, 7361, 115), source = "1721"),
  data.table(arm = "nf", contrast = c("F1_vs_F0", "F2_vs_F0", "F4_vs_F0", "F3_to_F4"), metric = "calls", target = c(280, 3358, 7129, 117), source = "1721"),
  data.table(arm = c("joint", "lib", "gini", "nf"), contrast = "F1_vs_F0", metric = "retained", target = c(99, 1272, 459, 273), source = "1721/4201"),
  data.table(arm = c("joint_raw", "joint_zsub"), contrast = "F1_vs_F0", metric = "calls", target = c(188, 123), source = "4201"),
  data.table(arm = c("joint_raw", "joint_zsub"), contrast = "F1_vs_F0", metric = "retained", target = c(180, 118), source = "4201"),
  data.table(arm = c("base_diseaseF0", "base_ctrlterm"), contrast = "F1_vs_F0", metric = "calls", target = c(783, 814), source = "RESULTS:1120"),
  data.table(arm = "joint", contrast = c("F2_vs_F0", "F3_vs_F0", "F4_vs_F0", "F1_to_F2", "F2_to_F3", "F3_to_F4"), metric = "retained",
             target = c(2113, 4797, 5307, 356, 497, 114), source = "A1 refit_vs_released_counts.tsv"),
  data.table(arm = "lib", contrast = c("F2_vs_F0", "F3_vs_F0", "F4_vs_F0", "F3_to_F4"), metric = "retained", target = c(3637, 5688, 6188, 267), source = "A1"),
  data.table(arm = "gini", contrast = c("F2_vs_F0", "F4_vs_F0", "F3_to_F4"), metric = "retained", target = c(3694, 6088, 107), source = "A1"),
  data.table(arm = "nf", contrast = c("F2_vs_F0", "F4_vs_F0", "F3_to_F4"), metric = "retained", target = c(3322, 6014, 110), source = "A1"),
  data.table(arm = "lib", contrast = c("F1_to_F2", "F2_to_F3"), metric = "calls", target = c(846, 546), source = "check 9050"),
  data.table(arm = "gini", contrast = c("F1_to_F2", "F2_to_F3", "F3_vs_F0"), metric = "calls", target = c(807, 696, 6090), source = "check 9050"),
  data.table(arm = "nf", contrast = c("F1_to_F2", "F2_to_F3", "F3_vs_F0"), metric = "calls", target = c(745, 722, 5922), source = "check 9050")))
tg[, observed := mapply(function(a, ct, m) { r <- resA[arm == a & contrast == ct]; if (nrow(r) != 1L) NA_real_ else as.numeric(r[[m]]) }, arm, contrast, metric)]
tg[, pass := !is.na(observed) & observed == target]
lfc <- resA[arm == "base", .(arm, contrast, metric = "max_abs_logFC_diff_vs_released", target = 1e-6, source = "release",
                              observed = max_abs_logFC_diff_vs_released, pass = max_abs_logFC_diff_vs_released < 1e-6)]
pc1 <- rbind(tg, lfc)
f1j <- resA[arm == "joint" & contrast == "F1_vs_F0"]
pc1_extra <- data.table(arm = "joint", contrast = "F1_vs_F0", metric = "median_logFC_ratio_on_released (not a stop rule)",
                        target = 0.7726, source = "4201", observed = round(f1j$median_logFC_ratio_on_released, 4),
                        pass = round(f1j$median_logFC_ratio_on_released, 4) == 0.7726)
wt(rbind(pc1, pc1_extra), "pc1_reproduction.tsv")
diag_add("6_PC1", paste0(pc1$arm, "__", pc1$contrast, " ", pc1$metric, " [", pc1$source, "]"), pc1$observed, pc1$target, pc1$pass)
diag_add("6_PC1", paste0(pc1_extra$arm, "__", pc1_extra$contrast, " ", pc1_extra$metric), pc1_extra$observed, pc1_extra$target, pc1_extra$pass)
PC1_pass <- all(pc1$pass)
logm("PC-1 reproduction: ", sum(pc1$pass), "/", nrow(pc1), " pass -> ", if (PC1_pass) "PASS" else "FAIL")

## D5 fibrosis (reproduction targets) and D6 F1 targets, computed on stage-A fits
d5_axis <- function(rn) {
  r <- ros[[rn]]
  lc <- cpm(d[, r$sample_id], log = TRUE, prior.count = 2)
  slopes <- sapply(COVJ, function(v) {
    X <- model.matrix(as.formula(paste("~ dataset + inferred_sex + group +", v)), r)
    lmFit(lc, X)$coefficients[, v] })
  rownames(slopes) <- rownames(lc)
  rbindlist(lapply(CT[roster == rn, contrast], function(ct) {
    rc <- REL[[ct]]; sig <- rc$FDR < 0.05; e <- rc$logFC[sig]; a <- slopes[rc$gene_id_versioned[sig], , drop = FALSE]
    data.table(axis_roster = rn, n_roster = nrow(r), contrast = ct, n_released_calls = sum(sig),
               rho_lib = if (sum(sig) > 2) cor(e, a[, "z_loglib"], method = "spearman") else NA_real_,
               rho_gini = if (sum(sig) > 2) cor(e, a[, "z_gini"], method = "spearman") else NA_real_,
               rho_nf = if (sum(sig) > 2) cor(e, a[, "z_nf"], method = "spearman") else NA_real_)
  }))
}
d5f <- d5_axis("fibrosis")
d5_tg <- data.table(contrast = c("F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0"),
                    t_nf = c(0.689, 0.738, 0.647, 0.550), t_gini = c(-0.675, -0.721, -0.639, -0.537), t_lib = c(0.406, 0.135, 0.061, 0.206))
d5c <- merge(d5f, d5_tg, by = "contrast")
d5c[, pass := abs(rho_nf - t_nf) <= 5e-4 + 1e-9 & abs(rho_gini - t_gini) <= 5e-4 + 1e-9 & abs(rho_lib - t_lib) <= 5e-4 + 1e-9]
diag_add("D5_reproduction", paste0(d5c$contrast, " rho nf/gini/lib"), paste(round(d5c$rho_nf, 4), round(d5c$rho_gini, 4), round(d5c$rho_lib, 4), sep = "/"),
         paste(d5c$t_nf, d5c$t_gini, d5c$t_lib, sep = "/"), d5c$pass)

if (!PC1_pass && !SYNTH) finish_stopped("S2", paste0("PC-1 reproduction failed for ", sum(!pc1$pass), " of ", nrow(pc1), " targets (pc1_reproduction.tsv)"))
if (STAGE == "gate") finish_stopped("GATE_ONLY", paste0("PC-1 ", if (PC1_pass) "passed" else "failed", " (", sum(pc1$pass), "/", nrow(pc1), "); stage B and the null were not run"))

## ---------- stage B: blind real arms and PC-3 ----------
nas5 <- CT$contrast[8:12]
uF1 <- units[["F1_vs_F0"]]
g_ind <- zc(as.numeric(uF1$group == "F1"), uF1$dataset)
pc3_z <- sapply(1:20, function(k) zc(0.6 * g_ind + 0.8 * plant_e[, k], uF1$dataset))
tasksB <- c(unlist(lapply(c("lib", "gini", "nf", "joint"), function(a) lapply(nas5, function(ct) make_task(a, ct, "B"))), recursive = FALSE),
            unlist(lapply(c("joint_raw", "joint_zsub"), function(a) lapply(setdiff(CT$contrast, "F1_vs_F0"), function(ct) make_task(a, ct, "B"))), recursive = FALSE),
            list(make_task("joint_diseaseF0", "F1_vs_F0", "B"), make_task("joint_ctrlterm", "F1_vs_F0", "B")),
            lapply(1:20, function(k) { tk <- make_task(sprintf("pc3_%02d", k), "F1_vs_F0", "B"); tk$covdf <- data.frame(z_plant = pc3_z[, k]); tk }))
resB <- run_batch(tasksB, "B")
REAL <- rbind(resA, resB, fill = TRUE)
wt(REAL, "real_arm_summaries.tsv")

## ---------- stage C: permutation null ----------
null_task <- function(ct, b) {
  dirc <- file.path(NULLDIR, ct); f <- file.path(dirc, paste0(b, ".tsv"))
  if (file.exists(f)) return("exists")
  i <- CI(ct); u <- units[[ct]]; Z <- as.data.frame(u[perms[[ct]][, b], ..COVJ])
  X <- build_X(u, CT$ref[i], CT$cmp[i], Z)
  if (qr(X)$rank < ncol(X)) {
    row <- data.table(contrast = ct, b = b, calls = NA_integer_, retained = NA_integer_, released = rel_calls[[ct]], rank_ok = FALSE,
                      cov_z_loglib = NA_integer_, cov_z_gini = NA_integer_, cov_z_nf = NA_integer_, fit_s = 0)
  } else {
    ts <- Sys.time(); f2 <- fit_X(u, X, CT$cmp[i]); s <- vs_rel(f2$g, REL[[ct]])
    row <- data.table(contrast = ct, b = b, calls = s$calls, retained = s$retained, released = s$released, rank_ok = TRUE,
                      cov_z_loglib = f2$covbh[["z_loglib"]], cov_z_gini = f2$covbh[["z_gini"]], cov_z_nf = f2$covbh[["z_nf"]],
                      fit_s = round(as.numeric(difftime(Sys.time(), ts, units = "secs"))))
  }
  tmp <- paste0(f, ".part.", Sys.getpid()); fwrite(row, tmp, sep = "\t")
  if (file.exists(f)) { unlink(tmp); return("exists") }
  file.rename(tmp, f); "written"
}
for (ct in CT[null == TRUE, contrast]) dir.create(file.path(NULLDIR, ct), showWarnings = FALSE)
null_specs <- unlist(lapply(seq_len(B), function(b) lapply(CT[null == TRUE, contrast], function(ct) list(ct = ct, b = b))), recursive = FALSE)
logm("stage C (null): ", length(null_specs), " fits on ", WORKERS, " workers")
progress_every <- max(1L, length(CT[null == TRUE, contrast]) * 10L)
resC <- mclapply(seq_along(null_specs), function(k) {
  s <- null_specs[[k]]
  r <- tryCatch(null_task(s$ct, s$b), error = function(e) paste("ERROR", conditionMessage(e)))
  if (k %% progress_every == 0L) cat(format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), "null progress: task", k, "of", length(null_specs), "\n")
  r }, mc.cores = WORKERS, mc.preschedule = FALSE)
nerr <- sum(grepl("^ERROR", unlist(resC)))
if (nerr > 0L) { print(head(unlist(resC)[grepl("^ERROR", unlist(resC))])); stop(nerr, " null fits failed") }
NUL <- rbindlist(lapply(CT[null == TRUE, contrast], function(ct) rbindlist(lapply(seq_len(B), function(b) fread(file.path(NULLDIR, ct, paste0(b, ".tsv")))))))
stopifnot(nrow(NUL) == B * sum(CT$null))
wt(NUL, "null_replicates.tsv")
logm("stage C done; ", el(), " s")

## ---------- null summaries (computed before diagnostics 7 to 9 are written) ----------
nsum <- NUL[, {
  R <- retained[rank_ok] / released[1]
  list(B = .N, rank_deficient = sum(!rank_ok), mean_R = mean(R), median_R = median(R),
    q025_R = unname(quantile(R, 0.025)), q975_R = unname(quantile(R, 0.975)), min_R = min(R), max_R = max(R),
    mean_calls = mean(calls[rank_ok]), median_calls = median(calls[rank_ok]), q025_calls = unname(quantile(calls[rank_ok], 0.025)),
    q975_calls = unname(quantile(calls[rank_ok], 0.975)), min_calls = min(calls[rank_ok]), max_calls = max(calls[rank_ok]))
}, by = contrast]

## ---------- diagnostic 7: PC-3 planted confounder ----------
pc3 <- REAL[grepl("^pc3_", arm), .(draw = arm, calls, retained, released, retained_frac)]
pc3[, k := as.integer(sub("^pc3_", "", draw))]; setorder(pc3, k)
stopifnot(identical(pc3$k, 1:20))
pc3[, cor_zplant_group := vapply(k, function(kk) cor(pc3_z[, kk], as.numeric(uF1$group == "F1")), 0)]
q025_F1 <- nsum[contrast == "F1_vs_F0", q025_R]
pc3[, below_half := retained_frac < 0.5]; pc3[, below_null_q025 := retained_frac < q025_F1]
PC3_pass <- sum(pc3$below_half) >= 18L && all(pc3$below_null_q025)
if (SYNTH && nzchar(Sys.getenv("H82_SYNTH_FORCE_PC3_PASS"))) PC3_pass <- TRUE  # synthetic debug only: exercises the non-S5 reading code
wt(pc3, "pc3_planted.tsv")
diag_add("7_PC3", paste0(pc3$draw, " retained fraction (cor with group)"), paste0(round(pc3$retained_frac, 4), " (", round(pc3$cor_zplant_group, 3), ")"),
         "< 0.5 and < F1 null 2.5th pct", pc3$below_half & pc3$below_null_q025)
diag_add("7_PC3", "PC-3 overall: draws below 0.5; draws below F1 null 2.5th percentile",
         paste0(sum(pc3$below_half), "/20; ", sum(pc3$below_null_q025), "/20 (q025=", round(q025_F1, 4), ")"),
         ">= 18/20; 20/20", PC3_pass)

## ---------- diagnostic 8: null summaries ----------
wt(nsum, "null_summary.tsv")
diag_add("8_null", paste0(nsum$contrast, " null R mean/median/q025/q975/min/max"),
         paste(round(nsum$mean_R, 4), round(nsum$median_R, 4), round(nsum$q025_R, 4), round(nsum$q975_R, 4), round(nsum$min_R, 4), round(nsum$max_R, 4), sep = "/"),
         "NC-1: median >= 0.80 (F1 0.85-0.97)", nsum$median_R >= 0.80)
diag_add("8_null", paste0(nsum$contrast, " null calls mean/median/q025/q975/min/max"),
         paste(round(nsum$mean_calls, 1), nsum$median_calls, nsum$q025_calls, nsum$q975_calls, nsum$min_calls, nsum$max_calls, sep = "/"), "descriptive", NA)

## ---------- diagnostic 9: covariate-coefficient BH counts ----------
d9_real <- REAL[status == "fitted" & covariates != "", .(contrast, arm, cov_bh_sig)]
d9_null <- NUL[rank_ok == TRUE, .(mean_z_loglib = mean(cov_z_loglib), max_z_loglib = max(cov_z_loglib), mean_z_gini = mean(cov_z_gini),
                                  max_z_gini = max(cov_z_gini), mean_z_nf = mean(cov_z_nf), max_z_nf = max(cov_z_nf)), by = contrast]
wt(d9_real, "d9_covariate_bh_real.tsv"); wt(d9_null, "d9_covariate_bh_null.tsv")
diag_add("9_cov_bh", paste0(d9_real$arm, "__", d9_real$contrast), d9_real$cov_bh_sig, "descriptive", NA)
diag_add("9_cov_bh", paste0("null ", d9_null$contrast, " mean/max BH<0.05 genes on z_loglib, z_gini, z_nf"),
         paste0(round(d9_null$mean_z_loglib, 1), "/", d9_null$max_z_loglib, ", ", round(d9_null$mean_z_gini, 1), "/", d9_null$max_z_gini, ", ",
                round(d9_null$mean_z_nf, 1), "/", d9_null$max_z_nf), "near zero (prediction)", NA)

## ---------- primary statistic, T1 and T2 ----------
S5 <- S5_cov || !PC3_pass
S5_reason <- c(if (S5_cov) "covariate recomputation differs by >= 1e-10", if (!PC3_pass) "PC-3 failed")
prim <- rbindlist(lapply(CT[null == TRUE, contrast], function(ct) {
  rj <- REAL[arm == "joint" & contrast == ct]; nn <- NUL[contrast == ct]; ns <- nsum[contrast == ct]
  k_rd <- sum(!nn$rank_ok)
  ge <- sum(nn$retained[nn$rank_ok] <= rj$retained)
  data.table(contrast = ct, family = CT[contrast == ct, testable_family], released = rj$released, joint_calls = rj$calls,
             joint_retained = rj$retained, R = rj$retained_frac, null_median = ns$median_R, null_q025 = ns$q025_R, null_q975 = ns$q975_R,
             Delta = ns$median_R - rj$retained_frac, n_null_le_R = ge, rank_deficient = k_rd, untestable_S4 = k_rd > 10L,
             p = (1 + ge + k_rd) / (B + 1))
}))
prim[, p_BH := NA_real_]
t2 <- prim[family == "T2" & untestable_S4 == FALSE, .(contrast, p)]
t2[, p_BH := p.adjust(p, method = "BH")]
prim[t2, p_BH := i.p_BH, on = "contrast"]
prim[, rejected := fifelse(untestable_S4, NA, fifelse(family == "T1", p < 0.05, p_BH <= 0.05))]
prim[, nc1_label := fifelse(null_median < 0.50, "not specific to complexity", fifelse(null_median < 0.80, "weak null", "null ok"))]
cls_of <- function(rej, delta, nm) {
  if (is.na(rej)) return("untestable (S4)")
  if (nm < 0.50) return("not specific to complexity")
  if (!rej || delta < 0.10) "robust" else if (delta < 0.40) "moderate" else "severe"
}
prim[, class := mapply(cls_of, rejected, Delta, null_median)]
if (S5) prim[, class := paste0("untestable (S5): ", class)]
wt(prim, "primary_tests.tsv")

## ---------- descriptive outputs ----------
D1 <- REAL[arm %in% c(ARMS_MAIN, "joint_raw", "joint_zsub"),
           .(contrast, arm, n_units, covariates, calls, up, down, released, retained, retained_frac, new_calls,
             median_logFC_ratio_on_released, max_abs_logFC_diff_vs_released, cov_bh_sig, design_rank, design_cols)]
D1[, contrast := factor(contrast, levels = CT$contrast)]; D1[, arm := factor(arm, levels = c(ARMS_MAIN, "joint_raw", "joint_zsub"))]
setorder(D1, contrast, arm); D1[, `:=`(contrast = as.character(contrast), arm = as.character(arm))]
wt(D1, "d1_arm_table.tsv")
d2 <- rbindlist(lapply(CT$contrast, function(ct) {
  rc <- REL[[ct]]; g <- fread(file.path(FITS, paste0("joint__", ct, ".genes.tsv.gz")))
  R <- rc$FDR < 0.05; at <- abs(rc$t)
  st <- cut(at, c(-Inf, 4, 5, Inf), right = FALSE, labels = c("[threshold,4)", "[4,5)", "[5,inf)"))
  rbindlist(lapply(levels(st), function(s) { k <- R & st == s
    data.table(contrast = ct, stratum = s, released = sum(k), retained = sum(k & g$FDR < 0.05),
               retained_frac = if (sum(k) > 0) sum(k & g$FDR < 0.05) / sum(k) else NA_real_) }))
}))
wt(d2, "d2_t_strata.tsv")
d4 <- REAL[contrast == "F1_vs_F0" & arm %in% c("base", "joint", "base_diseaseF0", "base_ctrlterm", "joint_diseaseF0", "joint_ctrlterm"),
           .(arm, n_units, calls, retained, released, retained_frac)]
wt(d4, "d4_f0_reference.tsv")
d5n <- d5_axis("NAS")
d5 <- rbind(d5f, d5n)
wt(d5, "d5_alignment.tsv")
named <- c("COL1A1", "COL1A2", "COL3A1", "LUM", "THBS2", "EFEMP1", "FBLN5", "LOX")
gene_rows <- function(ct, arm, genes) {
  g <- fread(file.path(FITS, paste0(arm, "__", ct, ".genes.tsv.gz"))); rc <- REL[[ct]]
  g[, gene_name := rc$gene_name[match(g$gene_id_versioned, rc$gene_id_versioned)]]
  out <- g[gene_name %in% genes, .(contrast = ct, arm, gene_name, gene_id_versioned, logFC, FDR)]
  out[, `:=`(released_logFC = rc$logFC[match(gene_id_versioned, rc$gene_id_versioned)], released_FDR = rc$FDR[match(gene_id_versioned, rc$gene_id_versioned)])]
  out
}
pc2 <- rbind(gene_rows("F4_vs_F0", "joint", named)[, pc := "PC-2a"], gene_rows("NAS5_8_vs_NAS0", "joint", named)[, pc := "PC-2b"])
pc2[, survives := FDR < 0.05 & logFC > 0]
wt(pc2, "pc2_markers.tsv")
PC2a_n <- pc2[pc == "PC-2a", sum(survives)]; PC2b_n <- pc2[pc == "PC-2b", sum(survives)]
PC2a_pass <- PC2a_n == 8L && pc2[pc == "PC-2a", .N] == 8L
d6 <- rbind(gene_rows("F1_vs_F0", "joint", c("GNMT", "MAT1A", "CYP2C19")), gene_rows("F4_vs_F0", "joint", c("GNMT", "MAT1A", "CYP2C19")))
d6_tg <- data.table(gene_name = c("GNMT", "MAT1A", "CYP2C19"), t_logFC = c(-0.225, -0.156, -0.184), t_FDR = c(0.36, 0.048, 0.75))
d6 <- merge(d6, d6_tg, by = "gene_name", all.x = TRUE)
d6[contrast == "F1_vs_F0", reproduces_A1 := round(logFC, 3) == t_logFC & signif(FDR, 2) == t_FDR]
wt(d6, "d6_named_genes.tsv")
diag_add("D6_reproduction", paste0(d6[contrast == "F1_vs_F0", gene_name], " F1 joint logFC/FDR"),
         paste0(round(d6[contrast == "F1_vs_F0", logFC], 4), "/", signif(d6[contrast == "F1_vs_F0", FDR], 3)),
         paste0(d6[contrast == "F1_vs_F0", t_logFC], "/", d6[contrast == "F1_vs_F0", t_FDR]), d6[contrast == "F1_vs_F0", reproduces_A1])

## ---------- readings (prespec "Reading" section, coded) ----------
F1 <- prim[contrast == "F1_vs_F0"]; T2 <- prim[family == "T2"]
t2_testable <- T2[untestable_S4 == FALSE]
cls_plain <- function(x) sub("^untestable \\(S5\\): ", "", x)
F1cls <- cls_plain(F1$class); T2cls <- cls_plain(t2_testable$class)
strata <- c("[threshold,4)", "[4,5)", "[5,inf)")
d2_f1_below <- vapply(strata, function(s) {
  f1 <- d2[contrast == "F1_vs_F0" & stratum == s, retained_frac]
  oth <- d2[contrast %in% t2_testable$contrast & stratum == s & released > 0, retained_frac]
  length(f1) == 1L && !is.na(f1) && length(oth) > 0L && all(f1 < oth) }, logical(1))
D2_cond <- sum(d2_f1_below) >= 2L
f1_alone <- F1cls == "severe" && length(T2cls) > 0L && all(T2cls == "robust") && D2_cond
graded <- F1cls %in% c("severe", "moderate") && any(T2cls %in% c("moderate", "severe"))
r_joint_F1 <- d4[arm == "joint", retained_frac]; r_dis <- d4[arm == "joint_diseaseF0", retained_frac]
D4_flag <- (r_dis - r_joint_F1) >= 0.25
PC2b_over <- PC2b_n < 6L

if (S5) {
  reading <- paste0("UNTESTABLE (S5: ", paste(S5_reason, collapse = "; "), "). T1 and T2 are untestable; D1 counts are descriptive only. No H82 number enters the paper.")
  status <- "untestable_S5"
} else {
  f1_read <- if (F1$null_median < 0.50 || F1$p >= 0.05) {
    "Not specific to complexity: the F1 count is fragile to any added covariates (F1 null median < 0.50 or T1 p >= 0.05); no complexity-specific wording; the F1 bar is still qualified with the null retention."
  } else if (F1$null_median < 0.80) {
    "T1 rejected but the F1 null is weak (median < 0.80): the F1 count is fragile to any three added columns; T1 reported with the label 'weak null'; complexity-specific wording not allowed."
  } else if (!PC2a_pass) {
    "T1 rejected with null median >= 0.80, but PC-2a failed: the complexity-specific F1 reading requires PC-2a, so that wording is not allowed."
  } else {
    "Complexity-specific F1 loss: T1 p < 0.05, F1 null median >= 0.80, PC-1, PC-2a and PC-3 pass."
  }
  pat <- if (f1_alone) {
    "Affects F1 alone: F1 severe, every T2 contrast robust, and D2 F1 retention below every T2 contrast in >= 2 of 3 |t| strata ('the sensitivity is confined to the F0-to-F1 comparison')."
  } else if (graded) {
    paste0("Graded: F1 ", F1cls, " and at least one T2 contrast moderate or severe (",
           paste(t2_testable[cls_plain(class) %in% c("moderate", "severe"), paste0(contrast, " ", cls_plain(class))], collapse = ", "),
           "). The legend lists retained fractions for every bar, names these contrasts, and drops any wording that F1 is uniquely fragile.")
  } else {
    paste0("Neither 'F1 alone' nor 'graded': F1 class ", F1cls, "; T2 classes ", paste(paste0(t2_testable$contrast, " ", T2cls), collapse = ", "),
           "; D2 condition ", D2_cond, ".")
  }
  reading <- paste(f1_read, pat)
  status <- "completed"
}
extra <- c(if (!S5 && PC2b_over) paste0("PC-2b failed (", PC2b_n, " of 8 ECM genes survive at NAS5_8_vs_NAS0 joint): NAS adjusted counts are reported as possibly over-adjusted and not used to qualify the NAS bars.") else NULL,
           if (!S5 && D4_flag) paste0("D4: disease-only F0 + joint retains ", round(r_dis, 4), " vs joint ", round(r_joint_F1, 4), " (difference >= 0.25): the text names the F0 controls as part of the mechanism.") else NULL)

## ---------- RESULT_H82.md ----------
fmt <- function(x, k = 4) formatC(x, digits = k, format = "f")
pct <- function(x) paste0(formatC(100 * x, digits = 1, format = "f"), "%")
ptab <- prim[, sprintf("| %s | %s | %d / %d | %s | %s | %s [%s, %s] | %s | %s | %s | %s | %s |", contrast, family, joint_retained, released, fmt(R),
                        fmt(Delta), fmt(null_median), fmt(null_q025), fmt(null_q975), formatC(p, digits = 6, format = "f"),
                        ifelse(is.na(p_BH), "-", formatC(p_BH, digits = 5, format = "f")), ifelse(is.na(rejected), "-", as.character(rejected)), nc1_label, class)]
nas3 <- REAL[arm == "joint" & contrast == "NAS3_4_vs_NAS1_2"]
d1tab <- dcast(D1[arm %in% ARMS_MAIN], contrast ~ arm, value.var = "retained_frac")[match(CT$contrast, contrast)]
d1calls <- dcast(D1[arm %in% ARMS_MAIN], contrast ~ arm, value.var = "calls")[match(CT$contrast, contrast)]
d1lines <- vapply(seq_len(nrow(d1tab)), function(k) sprintf("| %s | %d | %s | %s | %s | %s |", d1tab$contrast[k], d1calls$base[k],
  paste0(d1calls$lib[k], " (", fmt(d1tab$lib[k], 3), ")"), paste0(d1calls$gini[k], " (", fmt(d1tab$gini[k], 3), ")"),
  paste0(d1calls$nf[k], " (", fmt(d1tab$nf[k], 3), ")"), paste0(d1calls$joint[k], " (", fmt(d1tab$joint[k], 3), ")")), "")
d2w <- dcast(d2, contrast ~ stratum, value.var = c("released", "retained_frac"))[match(CT$contrast, contrast)]
d2lines <- vapply(seq_len(nrow(d2w)), function(k) sprintf("| %s | %d (%s) | %d (%s) | %d (%s) |", d2w$contrast[k],
  d2w[["released_[threshold,4)"]][k], fmt(d2w[["retained_frac_[threshold,4)"]][k], 3), d2w[["released_[4,5)"]][k], fmt(d2w[["retained_frac_[4,5)"]][k], 3),
  d2w[["released_[5,inf)"]][k], fmt(d2w[["retained_frac_[5,inf)"]][k], 3)), "")
legend_new <- if (!S5) sprintf(paste0("Adjusting for within-cohort library size, Gini coefficient and normalization factor retains %s of the\n",
  "1,338 F1 genes (null with permuted covariates: median %s); retention for F2-F4 is %s / %s / %s. The full\n",
  "fibrosis/NAS, composition, F0-reference and library-complexity sensitivities are shown in Figure S3."),
  format(F1$joint_retained, big.mark = ","), pct(F1$null_median), pct(prim[contrast == "F2_vs_F0", R]), pct(prim[contrast == "F3_vs_F0", R]),
  pct(prim[contrast == "F4_vs_F0", R])) else NA_character_
res_lines <- c(
  "# RESULT H82", "",
  if (SYNTH) "**SYNTHETIC DEBUG RUN: permuted labels with planted effects on a 2,000-gene subset. Not a result.**" else NULL,
  if (SYNTH) "" else NULL,
  paste0("Status: ", status, ". Prespec sha256 ", PRESPEC_SHA, " (verified). Run: ", OUT, "."), "",
  "## Verdict (computed by the prespec rule)", "", reading, "", if (length(extra)) c(paste0("- ", extra), "") else NULL,
  "## Primary statistic", "",
  paste0("R_c = released BH calls of contrast c that are also BH < 0.05 in the joint fit (within-cohort z of log2 library size, Gini, stored norm factor), divided by released calls. Null: the three z columns permuted together within cohort, B = ",
  B, ". One-sided p_c = (1 + #{b: R_c^(b) <= R_c} + rank-deficient replicates) / ", B + 1, " (floor ", fmt(1 / (B + 1), 6), "). Delta_c = null median - R_c."), "",
  "T1: family of 1 (F1_vs_F0), alpha 0.05, reject if p < 0.05. T2: family of the testable contrasts among 10, BH at q 0.05 (reject if BH-adjusted p <= 0.05). Within each refit, BH over all 24,115 genes defines calls (producer's family).", "",
  "| Contrast | Family | Joint retained / released | R_c | Delta_c | Null median [2.5%, 97.5%] | p_c | p_BH (T2) | Rejected | NC-1 label | Class |",
  "|---|---|---|---|---|---|---|---|---|---|---|", ptab, "",
  sprintf("NAS3_4_vs_NAS1_2 (G2: %d released calls; not tested): joint %d calls, %d of %d retained.", nas3$released, nas3$calls, nas3$retained, nas3$released), "",
  "## Controls", "",
  sprintf("- PC-1 reproduction (S2): %d of %d targets reproduced exactly (pc1_reproduction.tsv); base max |logFC diff| vs release %s.", sum(pc1$pass), nrow(pc1),
          signif(max(resA[arm == "base", max_abs_logFC_diff_vs_released]), 3)),
  sprintf("- PC-2a (F4_vs_F0 joint, 8 ECM genes BH < 0.05 and logFC > 0): %d of 8 -> %s.", PC2a_n, if (PC2a_pass) "pass" else "FAIL"),
  sprintf("- PC-2b (NAS5_8_vs_NAS0 joint; predicted >= 7 of 8, over-adjusted if < 6): %d of 8.", PC2b_n),
  sprintf("- PC-3 (planted covariate, corr with group about 0.6; 20 draws): %d of 20 below 0.5, %d of 20 below the F1 null 2.5th percentile %s; median retained %s -> %s.",
          sum(pc3$below_half), sum(pc3$below_null_q025), fmt(q025_F1), fmt(median(pc3$retained_frac)), if (PC3_pass) "pass" else "FAIL (S5)"),
  sprintf("- NC-1 (null itself; predicted F1 median 0.85 to 0.97, all 11 >= 0.80): F1 median %s; %d of 11 contrasts have null median >= 0.80 (min %s, %s).",
          fmt(F1$null_median), sum(nsum$median_R >= 0.80), fmt(min(nsum$median_R)), nsum[which.min(median_R), contrast]),
  sprintf("- Covariate recomputation (S5): max abs diff %s.", paste(paste0(rec$covariate, " ", signif(rec$max_abs_diff, 3)), collapse = "; ")),
  sprintf("- Rank-deficient null replicates: %s (S4 threshold > 10).", paste(paste0(null_rank$contrast, " ", null_rank$rank_deficient), collapse = ", ")),
  sprintf("- D5 alignment reproduction (F1-F4 vs F0): %d of 4 contrasts within 5e-4; D6 F1 joint lookups reproduce A1: %d of 3.", sum(d5c$pass), sum(d6[contrast == "F1_vs_F0", reproduces_A1], na.rm = TRUE)), "",
  "## D1: calls (retained fraction of released) per arm", "",
  "| Contrast | Released (base) | lib | gini | nf | joint |", "|---|---|---|---|---|---|", d1lines, "",
  "## D2: joint retention by released |t| stratum, released calls (retained fraction)", "",
  "| Contrast | [threshold,4) | [4,5) | [5,inf) |", "|---|---|---|---|", d2lines, "",
  sprintf("D2 condition for 'F1 alone' (F1 below every testable T2 contrast in >= 2 of 3 strata): strata met %s -> %s. The prespec gives no numeric rule for 'similar within-stratum retention'; that sentence is a judgement on the table above.",
          if (any(d2_f1_below)) paste(names(d2_f1_below)[d2_f1_below], collapse = ", ") else "none", D2_cond), "",
  "## D3, D4, D6", "",
  sprintf("- D3 F1_vs_F0: joint raw %d calls (%d retained); joint z within the 205 units %d (%d). All 12 contrasts in d1_arm_table.tsv.",
          REAL[arm == "joint_raw" & contrast == "F1_vs_F0", calls], REAL[arm == "joint_raw" & contrast == "F1_vs_F0", retained],
          REAL[arm == "joint_zsub" & contrast == "F1_vs_F0", calls], REAL[arm == "joint_zsub" & contrast == "F1_vs_F0", retained]),
  sprintf("- D4 F1_vs_F0 retained fraction: joint %s; joint + disease-only F0 (%d units) %s; joint + control-status term %s; base disease-only F0 %d calls; base control-status %d calls.",
          fmt(r_joint_F1), d4[arm == "joint_diseaseF0", n_units], fmt(r_dis), fmt(d4[arm == "joint_ctrlterm", retained_frac]),
          d4[arm == "base_diseaseF0", calls], d4[arm == "base_ctrlterm", calls]),
  paste0("- D6 joint logFC / FDR: ", paste(d6[, sprintf("%s %s %s / %s", gene_name, contrast, fmt(logFC, 3), signif(FDR, 2))], collapse = "; "), "."),
  "- D7 (quoted from 1721, not rerun): F1_vs_F0 leave-one-cohort-out calls 14 / 929 / 651 without GSE130970 / GSE135251 / GSE162694.", "",
  "## Proposed legend text (returned as data, not installed)", "",
  "File docs/manuscript/05_figure_legends.md, lines 145-146.", "",
  "old_text: \"The full fibrosis/NAS, composition and F0-reference sensitivities are\\nshown in Figure S3.\"", "",
  paste0("new_text: ", if (is.na(legend_new)) "none (untestable)" else paste0("\"", gsub("\n", "\\\\n", legend_new), "\"")), "",
  "This template is the prespec's. Its wording applies only under the reading above; under 'graded' the legend must also name the moderate or severe contrasts and drop any F1-unique wording. Released counts stay printed; adjusted counts are a sensitivity (lower bound), because Gini and norm factor can lie downstream of disease.", "",
  "## Units and scope", "",
  "Unit: bulk analysis unit (genetic-identity cluster, not a verified person); fibrosis roster 387 units, NAS roster 395 units, three cohorts. Genes are not replicates. Covariates in every model: dataset, inferred sex, plus the listed complexity columns. Cross-sectional stage groups from unrelated biopsies. Human bulk RNA-seq only; no Cas13 data, no MR, GSE289173 not opened, Hotspot programs not used.")
writeLines(res_lines, RESULT)
write_json(list(status = status, synthetic = SYNTH, reading = reading, extra = extra, S5 = S5, S5_reason = S5_reason,
                PC1_pass = PC1_pass, PC2a = PC2a_n, PC2b = PC2b_n, PC3_pass = PC3_pass, D2_condition = D2_cond, D4_flag = D4_flag,
                f1_alone = f1_alone, graded = graded, primary = prim, legend_new = legend_new), "verdict.json")
logm("FINISH status=", status, "; ", reading)
logm("done; ", el(), " s")
