#!/usr/bin/env Rscript
# H79: do the 1,854 canonical disease calls of the adopted 774-unit bulk fit hold when
# genetic ancestry PCs from RNA genotypes enter the model?
# Prespecification: <run dir>/prespec/PRESPEC_H79.md (sha256 checked below; never edited).
# Model steps follow scripts/v1_774_gc_covariate.template.R (copy of the V1-774 script,
# which copies fit_disease_model() of refit_corrected_controls.R).
#
# Environment: H79_OUT (existing output dir), H79_CORES, H79_SEED (must be 20260923),
# H79_DEBUG ("1" = code-path test: every 10th gene, tiny B, and the "observed" PCs and EUR
# set replaced by one extra within-cell permutation and one extra subset, so the real
# outcome is never computed), H79_RESUME_DRAWS (optional: draws/ dir of an earlier run with
# identical draw indices; finished draw files are copied, not refit).

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(parallel)
})
fail <- function(...) stop(..., call. = FALSE)
now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
say <- function(...) { cat(sprintf("[%s] ", now()), ..., "\n", sep = ""); flush(stdout()) }

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT <- Sys.getenv("H79_OUT")
if (!nzchar(OUT) || !dir.exists(OUT)) fail("H79_OUT must be an existing directory")
H79 <- dirname(dirname(OUT))
NCORES <- as.integer(Sys.getenv("H79_CORES", "1"))
SEED <- as.integer(Sys.getenv("H79_SEED", NA))
DEBUG <- identical(Sys.getenv("H79_DEBUG", "0"), "1")
RESUME <- Sys.getenv("H79_RESUME_DRAWS", "")
if (!identical(SEED, 20260923L)) fail("H79_SEED must be 20260923")
B_PERM <- if (DEBUG) 4L else 200L
B_SUB <- if (DEBUG) 3L else 100L
COEF <- "group_binaryDisease"
PCN <- c("PC1", "PC2", "PC3", "PC4")
say("H79 start; OUT=", OUT, " cores=", NCORES, " debug=", DEBUG, " B_PERM=", B_PERM, " B_SUB=", B_SUB)

w <- function(x, name) {
  p <- file.path(OUT, name)
  if (file.exists(p)) fail("Refusing to overwrite ", p)
  fwrite(x, p, sep = "\t")
}
wlines <- function(x, name) {
  p <- file.path(OUT, name)
  if (file.exists(p)) fail("Refusing to overwrite ", p)
  writeLines(x, p)
}

# Stop rules S1-S7: write RESULT_H79.md naming the rule and quit. No number is used.
stop_by_rule <- function(rule, detail) {
  say("STOP ", rule, ": ", detail)
  txt <- c("# RESULT_H79", "",
           if (DEBUG) "DEBUG RUN (code-path test; not the analysis)." else NULL,
           sprintf("Status: stopped_by_rule (%s)", rule), "",
           sprintf("Detail: %s", detail), "",
           "Reading (prespec): stopped by S1 to S7, so no number from the run is used; the bulk",
           "claims stay as adopted and the review records the analysis as not completed.",
           "", sprintf("Written %s.", now()))
  writeLines(txt, file.path(OUT, "RESULT_H79.md"))
  writeLines(rule, file.path(OUT, "STOPPED_BY_RULE"))
  writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------ S1: input sha256
FIT_REL <- "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
R4 <- "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z"
ANC_REL <- "RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/obs-20260926T155651Z/ancestry/full-20260929T142807Z/library_ancestry_pcs.tsv"
GENO_REL <- "Analysis/Hypotheses/H3-20261009T180331Z/results/geno-20261009T180507Z/person_genotypes.tsv"
CONC_REL <- "Analysis/Hypotheses/V1-774-20261012T012943Z/inputs/catalog_r8_concordant_genes.tsv"
inputs <- data.table(
  role = c("prespec (run-dir copy)", "prespec (plans original)", "dge", "meta", "deg_results",
           "lvqw_weights", "ancestry_pcs", "pnpla3_calls", "v1_template", "gate_count_R",
           "gate_count_out", "prespec_count_R", "prespec_count_out", "d5_genes"),
  path = c(file.path(H79, "prespec/PRESPEC_H79.md"),
           file.path(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/PRESPEC_H79.md"),
           file.path(ROOT, FIT_REL, c("merged_dge.rds", "meta_matched.rds", "deg_results.csv", "lvqw_sample_weights.tsv")),
           file.path(ROOT, c(ANC_REL, GENO_REL, "Analysis/Hypotheses/V1-774-20261012T012943Z/v1_774_gc_covariate.R",
                             file.path(R4, "gate/NG18/count.R"), file.path(R4, "gate/NG18/count.out"),
                             file.path(R4, "write/H79/count_n.R"), file.path(R4, "write/H79/count_n.out"), CONC_REL))),
  bytes = c(NA, NA, 33106829, 9082, 5049936, 22417, 485660, 143454, 15934, 2512, 3803, 3502, 3328, 1209),
  sha256 = c("86d21b240928759268c3116d6f8d12c1e3efdaa3b0c3569626e0e7d5785ec632",
             "86d21b240928759268c3116d6f8d12c1e3efdaa3b0c3569626e0e7d5785ec632",
             "73f3c53fcb9c847e98b0e4f0918c52a04a6b0b21744e5b448cb45e5cbcfcad3b",
             "429e27376ad421f5f5fae19c75f5e760030e727fab65c5768ef562a9939edfc3",
             "0491ed10024388c3497dccb9f99bf669683edd5903c9d21c5308850d7c10bde6",
             "4e99563bd8d59dcba7efc81ab59ecb7f6f6fd5c4c1e8453226d93e435f0ec615",
             "3fbf184947b47ad6614dbcab238c7a78d117dcc9ac0152951ebc2659354f0d8e",
             "cab013902f1df18ade819dcb9227fd557d0cb6e312475af6cd78c2176559c517",
             "4482ea63d9880f4f3f8a7ce667995ea0841a401c0dcf816190df11761bb16b35",
             "f14d427eee01bfe8bb893ba4a8446004eb7663ccf42317789c7b6686816ba9d1",
             "a3204dc12f3efb3f1c54362929abb6e3231ff83f1c98f46dff86c4b2f3da9b8c",
             "868b37202b4721b77f5adf8b5acdcbd1a71fbcb407efffccccc18304efd1fe73",
             "ea7680f2af688c9832b03e90231499de5f6b26ec766d9ec3e7dd4f8024608eef",
             "7cc1307f64a94074894d983fba17466c1aff93fd323595a75e97dff2428c0d53"))
sha <- function(p) if (file.exists(p)) sub(" .*", "", system2("sha256sum", shQuote(p), stdout = TRUE)) else NA_character_
inputs[, observed_sha256 := vapply(path, sha, character(1))]
inputs[, observed_bytes := file.size(path)]
inputs[, ok := !is.na(observed_sha256) & observed_sha256 == sha256 & (is.na(bytes) | observed_bytes == bytes)]
w(inputs, "s1_input_sha256.tsv")
if (!all(inputs$ok)) stop_by_rule("S1", paste("sha256 or size mismatch:", paste(inputs[ok == FALSE, role], collapse = ", ")))
say("S1 PASSED: ", nrow(inputs), " input sha256 values and sizes match the prespec table")

# ------------------------------------------------------------------ inputs
FIT <- file.path(ROOT, FIT_REL)
dge <- readRDS(file.path(FIT, "merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(FIT, "meta_matched.rds")))
stored <- fread(file.path(FIT, "deg_results.csv"))
pcs <- fread(file.path(ROOT, ANC_REL))
lvqw <- fread(file.path(FIT, "lvqw_sample_weights.tsv"))
meta <- meta[match(colnames(dge), sample_id)]
if (!identical(meta$sample_id, colnames(dge))) fail("meta not aligned to DGEList columns")
if (!identical(stored$gene, rownames(dge))) fail("deg_results rows not aligned to DGEList rows")
if (anyNA(meta$inferred_sex)) fail("NA inferred_sex")
N <- ncol(dge)

# ------------------------------------------------------------------ S3: join and tables
idx_pc <- match(colnames(dge), pcs$run)                 # never merge(): it reorders rows
pcu <- pcs[idx_pc]
s3 <- list()
s3$n_units <- N
s3$n_unmatched <- sum(is.na(idx_pc))
s3$aligned <- isTRUE(identical(pcu$run, colnames(dge)))
s3$na_pc <- sum(!complete.cases(as.data.frame(pcu[, ..PCN])))
s3$na_call <- sum(is.na(pcu$ancestry_call))
s3$na_sex <- sum(is.na(meta$inferred_sex))
s3$fail_min_sites <- sum(!(as.character(pcu$passes_min_sites) %in% c("TRUE", "True", "true")))
s3$dataset_eq_cohort <- isTRUE(all(pcu$cohort == as.character(dge$samples$dataset)))
info <- data.frame(
  group_binary = factor(dge$samples$group_binary, levels = c("Control", "Disease")),
  dataset = factor(dge$samples$dataset),
  inferred_sex = factor(meta$inferred_sex),
  row.names = colnames(dge))
if (!identical(as.character(info$group_binary), meta$group_binary)) fail("group label mismatch")
info$cell <- interaction(info$dataset, info$group_binary, drop = TRUE, lex.order = TRUE)
info$ancestry_call <- pcu$ancestry_call
info$anc3 <- ifelse(info$ancestry_call %in% c("AFR", "EUR"), info$ancestry_call, "other")
info$possibly_mixed <- as.character(pcu$possibly_mixed) %in% c("TRUE", "True", "true")
PCM <- as.matrix(pcu[, ..PCN]); rownames(PCM) <- colnames(dge)
expected_tab <- data.table(
  dataset = rep(c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"), each = 2),
  group_binary = rep(c("Control", "Disease"), 5),
  AFR = c(0, 0, 0, 10, 6, 2, 11, 7, 14, 18),
  EUR = c(22, 6, 5, 48, 2, 171, 18, 77, 50, 215),
  other = c(1, 1, 1, 7, 2, 24, 2, 25, 4, 25))
obs_tab <- as.data.table(info)[, .(AFR = sum(anc3 == "AFR"), EUR = sum(anc3 == "EUR"), other = sum(anc3 == "other")),
                               by = .(dataset = as.character(dataset), group_binary = as.character(group_binary))]
setkey(obs_tab, dataset, group_binary); setkey(expected_tab, dataset, group_binary)
s3$table_matches <- isTRUE(all.equal(as.data.frame(obs_tab), as.data.frame(expected_tab), check.attributes = FALSE))
s3$n_cells <- nlevels(info$cell)
w(obs_tab, "s3_ancestry_by_cohort_group.tsv")
s3dt <- data.table(check = names(s3), value = vapply(s3, function(x) as.character(x), character(1)))
w(s3dt, "s3_join_checks.tsv")
if (s3$n_units != 774 || s3$n_unmatched != 0 || !s3$aligned || s3$na_pc != 0 || s3$na_call != 0 ||
    s3$na_sex != 0 || s3$fail_min_sites != 0 || !s3$dataset_eq_cohort || !s3$table_matches || s3$n_cells != 10) {
  stop_by_rule("S3", paste(sprintf("%s=%s", s3dt$check, s3dt$value), collapse = "; "))
}
say("S3 PASSED: 774 of 774 joined by match(), aligned, ancestry-by-cohort-by-group table equals the prespec")

# Debug: every 10th gene (deterministic) so the code path runs in minutes.
if (DEBUG) {
  keep_g <- seq(1, nrow(dge), by = 10)
  dge <- dge[keep_g, , keep.lib.sizes = TRUE]
  stored <- stored[keep_g]
  say("DEBUG: gene universe reduced to ", nrow(dge))
}
G <- nrow(dge)

# ------------------------------------------------------------------ fit function
F_M0 <- ~ dataset + inferred_sex + group_binary
F_M1 <- ~ dataset + inferred_sex + PC1 + PC2 + PC3 + PC4 + group_binary
F_D3 <- ~ dataset + inferred_sex + PC1 + group_binary
F_AUX <- ~ 0 + cell + inferred_sex + PC1 + PC2 + PC3 + PC4
full_rank <- function(f, d) { X <- model.matrix(f, data = d); qr(X)$rank == ncol(X) }
is_call <- function(lfc, padj) !is.na(padj) & padj < 0.05 & abs(lfc) > 0.5

fit_one <- function(dg, d, f, aux = FALSE, keep_E = FALSE) {
  design <- model.matrix(f, data = d)
  if (qr(design)$rank != ncol(design)) fail("design is rank deficient")
  v <- voomWithQualityWeights(dg, design)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = COEF, number = Inf, sort.by = "none")
  if (!identical(rownames(tt), rownames(dg))) fail("topTable row order")
  out <- list(logFC = tt$logFC, padj = tt$adj.P.Val, sw = v$targets$sample.weights)
  pcn <- intersect(PCN, colnames(design))
  if (length(pcn)) {
    tf <- topTable(fit, coef = pcn, number = Inf, sort.by = "none")
    out$pcF_q <- tf$adj.P.Val
    out$pc1_coef <- unname(fit$coefficients[, "PC1"])
    out$pc1_t <- unname(fit$t[, "PC1"])
  }
  if (aux) {
    da <- model.matrix(F_AUX, data = d)
    if (qr(da)$rank != ncol(da)) fail("auxiliary design is rank deficient")
    fa <- eBayes(lmFit(v, da))                    # reuses this draw's voom weights
    ta <- topTable(fa, coef = PCN, number = Inf, sort.by = "none")
    out$aux_q <- ta$adj.P.Val
  }
  if (keep_E) out$E <- v$E
  out
}

# ------------------------------------------------------------------ draws (seed set once; fixed order; before any fit)
eur_true <- which(info$ancestry_call == "EUR")
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
cells <- split(seq_len(N), info$cell)                 # fixed level order (lex.order)
perm <- lapply(seq_len(B_PERM), function(k) {
  p <- seq_len(N)
  for (ci in cells) p[ci] <- ci[sample.int(length(ci))]
  p
})
eur_n <- table(info$cell[eur_true])[names(cells)]
draw_sub <- function() sort(unlist(lapply(names(cells), function(cn) {
  ci <- cells[[cn]]; ci[sample.int(length(ci), eur_n[[cn]])]
}), use.names = FALSE))
subs <- lapply(seq_len(B_SUB), function(k) draw_sub())
if (DEBUG) {
  fake_perm <- seq_len(N); for (ci in cells) fake_perm[ci] <- ci[sample.int(length(ci))]
  fake_eur <- draw_sub()
}
draw_file <- file.path(OUT, "draw_indices.rds")
if (file.exists(draw_file)) fail("Refusing to overwrite ", draw_file)
saveRDS(list(perm = perm, subs = subs, seed = SEED, rngkind = RNGkind(), cell_levels = names(cells),
             units = colnames(dge), eur_n = eur_n), draw_file)
say("draw indices saved: ", B_PERM, " within-cell PC permutations, then ", B_SUB, " within-cell subsets")

# ------------------------------------------------------------------ S2: PC-0 reproduction
say("PC-0: refitting the adopted model M0 on the stored DGEList")
pc0 <- fit_one(dge, info, F_M0, keep_E = TRUE)
pc0_call <- is_call(pc0$logFC, pc0$padj)
stored_call <- is_call(stored$logFC, stored$padj)
repro <- data.table(
  n_units = N, n_control = sum(info$group_binary == "Control"), n_disease = sum(info$group_binary == "Disease"),
  n_genes = G, calls = sum(pc0_call), up = sum(pc0_call & pc0$logFC > 0), down = sum(pc0_call & pc0$logFC < 0),
  max_abs_logFC_delta_vs_stored = max(abs(pc0$logFC - stored$logFC)),
  max_abs_padj_delta_vs_stored = max(abs(pc0$padj - stored$padj)),
  stored_calls = sum(stored_call), call_mismatch_vs_stored = sum(pc0_call != stored_call),
  cor_M0_weights_vs_adopted = cor(pc0$sw, lvqw$sample_weight[match(colnames(dge), lvqw$sample_id)]))
w(repro, "pc0_reproduction.tsv")
if (!DEBUG && (repro$calls != 1854L || repro$up != 1333L || repro$down != 521L ||
               repro$max_abs_logFC_delta_vs_stored > 1e-6)) {
  stop_by_rule("S2", sprintf("calls %d up %d down %d max|dlogFC| %.3g", repro$calls, repro$up, repro$down,
                             repro$max_abs_logFC_delta_vs_stored))
}
say(if (DEBUG) "PC-0 (debug, thresholds not applied): " else "S2 PASSED (PC-0): ",
    sprintf("calls %d (up %d, down %d), max |log2FC - stored| %.3g, call mismatches vs stored %d",
            repro$calls, repro$up, repro$down, repro$max_abs_logFC_delta_vs_stored, repro$call_mismatch_vs_stored))
# The adopted call set is the stored release (deg_results.csv). In debug it is the debug M0 fit.
adopted <- if (DEBUG) pc0_call else stored_call
base_lfc <- if (DEBUG) pc0$logFC else stored$logFC
low_lfc <- adopted & abs(base_lfc) < 0.6
gstm1_row <- match("GSTM1", stored$symbol)
gstm1_E <- if (!is.na(gstm1_row)) pc0$E[gstm1_row, ] else NULL
m0_sw <- pc0$sw
pc0$E <- NULL

# ------------------------------------------------------------------ S6: PC-A (rs738409)
geno <- fread(file.path(ROOT, GENO_REL))
gl <- geno[, .(lib = unlist(strsplit(libraries, "[,;|]"))), by = .(person_id, rs738409_call)]
z <- merge(data.table(sample_id = colnames(dge), ancestry_call = info$ancestry_call,
                      group_binary = as.character(info$group_binary)),
           gl, by.x = "sample_id", by.y = "lib")
pca <- list(rows_joined = nrow(z), units_joined = uniqueN(z$sample_id))
z <- z[!is.na(rs738409_call)]
z[, Gc := c(RR = 0, RA = 1, AA = 2)[rs738409_call]]
pca$callable <- nrow(z)
qa <- z[ancestry_call %in% c("AFR", "EUR")]
al <- rbind(qa[, .(allele = "G", wt = Gc, anc = ancestry_call, grp = group_binary)],
            qa[, .(allele = "C", wt = 2 - Gc, anc = ancestry_call, grp = group_binary)])
t4 <- xtabs(wt ~ factor(allele, levels = c("G", "C")) + factor(anc, levels = c("AFR", "EUR")) + grp, al)
mh <- mantelhaen.test(t4, exact = TRUE)
tb <- rbind(G = tapply(qa$Gc, factor(qa$ancestry_call, levels = c("AFR", "EUR")), sum),
            C = tapply(2 - qa$Gc, factor(qa$ancestry_call, levels = c("AFR", "EUR")), sum))
qa[, carrier := Gc > 0]
t5 <- table(factor(qa$carrier, levels = c(TRUE, FALSE)), factor(qa$ancestry_call, levels = c("AFR", "EUR")), qa$group_binary)
mh3 <- mantelhaen.test(t5, exact = TRUE)
pca_freq <- rbind(
  z[, .(stratum = "all", n = .N, G_freq = mean(Gc) / 2), by = ancestry_call],
  qa[, .(n = .N, G_freq = mean(Gc) / 2), by = .(ancestry_call, stratum = group_binary)], use.names = TRUE)
w(pca_freq[order(stratum, ancestry_call)], "pcA_rs738409_freq_by_call.tsv")
pca_tests <- data.table(
  test = c("allele_CMH_exact_by_group (decision)", "allele_Fisher_unstratified", "carrier_CMH_exact_by_group"),
  odds_ratio_AFR_vs_EUR = c(unname(mh$estimate), unname(fisher.test(tb)$estimate), unname(mh3$estimate)),
  ci95_low = c(mh$conf.int[1], fisher.test(tb)$conf.int[1], mh3$conf.int[1]),
  ci95_high = c(mh$conf.int[2], fisher.test(tb)$conf.int[2], mh3$conf.int[2]),
  p = c(mh$p.value, fisher.test(tb)$p.value, mh3$p.value),
  predicted = c("OR 0.451 (0.27 to 0.73), p 6.3e-04", "p 8.0e-06", "OR 0.359, p 3.2e-04"))
w(pca_tests, "pcA_rs738409_tests.tsv")
pca$cmh_or <- unname(mh$estimate); pca$cmh_p <- mh$p.value
if (!(pca$cmh_or < 1 && pca$cmh_p < 0.05)) stop_by_rule("S6", sprintf("allele CMH OR %.3f p %.3g", pca$cmh_or, pca$cmh_p))
say(sprintf("S6 PASSED (PC-A): group-stratified exact CMH allele OR %.3f (%.3f to %.3f), p %.3g; callable %d",
            pca$cmh_or, mh$conf.int[1], mh$conf.int[2], pca$cmh_p, pca$callable))

# ------------------------------------------------------------------ S4 (true designs) and design diagnostics
design_diag <- data.table(
  model = c("M0", "M1", "ME"),
  rank = c(qr(model.matrix(F_M0, info))$rank, qr(model.matrix(F_M1, cbind(info, PCM)))$rank,
           qr(model.matrix(F_M0, info[eur_true, ]))$rank),
  ncol = c(ncol(model.matrix(F_M0, info)), ncol(model.matrix(F_M1, cbind(info, PCM))),
           ncol(model.matrix(F_M0, info[eur_true, ]))))
dd <- cbind(info, PCM)
r20 <- summary(lm(as.numeric(group_binary == "Disease") ~ dataset + inferred_sex, dd))$r.squared
r21 <- summary(lm(as.numeric(group_binary == "Disease") ~ dataset + inferred_sex + PC1 + PC2 + PC3 + PC4, dd))$r.squared
design_diag[, `:=`(r2_disease_on_dataset_sex = r20, r2_disease_on_dataset_sex_PCs = r21,
                   variance_inflation_disease_coef = (1 - r20) / (1 - r21))]
w(design_diag, "design_ranks.tsv")
if (any(design_diag$rank != design_diag$ncol)) stop_by_rule("S4", "M0, M1 or ME rank deficient")
pc_shift <- rbindlist(lapply(PCN, function(p) {
  m <- summary(lm(as.formula(paste(p, "~ dataset + group_binary")), dd))$coefficients[COEF, ]
  data.table(pc = p, disease_minus_control_within_cohort = m[["Estimate"]], t = m[["t value"]])
}))
w(pc_shift, "pc_shifts_within_cohort.tsv")

# Observed PCs and EUR set (debug: replaced by one extra null permutation and one extra subset)
PC_obs <- if (DEBUG) PCM[fake_perm, , drop = FALSE] else PCM
eur_idx <- if (DEBUG) fake_eur else eur_true
make_d <- function(P) { d <- info; d[, PCN] <- unname(P); d }

# S7: permutation integrity
cm <- function(P) apply(P, 2, function(x) tapply(x, info$cell, mean))
cm_true <- cm(PCM)
s7_dev <- vapply(perm, function(p) max(abs(cm(PCM[p, , drop = FALSE]) - cm_true)), numeric(1))
s7_cor <- vapply(perm, function(p) cor(PCM[, "PC1"], PCM[p, "PC1"]), numeric(1))
if (max(s7_dev) > 1e-12) stop_by_rule("S7", sprintf("max cell-mean deviation %.3g", max(s7_dev)))
say(sprintf("S7 PASSED: max |cell-mean difference| over %d permutations %.3g", B_PERM, max(s7_dev)))

# S4: every permuted, auxiliary, subset and sensitivity design
rk <- c(
  vapply(perm, function(p) full_rank(F_M1, make_d(PCM[p, , drop = FALSE])) &&
                           full_rank(F_AUX, make_d(PCM[p, , drop = FALSE])), logical(1)),
  full_rank(F_M1, make_d(PC_obs)), full_rank(F_AUX, make_d(PC_obs)),
  vapply(subs, function(s) full_rank(F_M0, info[s, ]), logical(1)),
  full_rank(F_M0, info[eur_idx, ]),
  full_rank(F_M1, make_d(PC_obs)[!info$possibly_mixed, ]),
  full_rank(F_D3, make_d(PC_obs)))
if (!all(rk)) stop_by_rule("S4", sprintf("%d of %d designs rank deficient", sum(!rk), length(rk)))
sub_counts_ok <- all(vapply(subs, function(s) identical(as.vector(table(info$cell[s])[names(cells)]), as.vector(eur_n)), logical(1)))
say("S4 PASSED: ", length(rk), " designs full rank")

# ------------------------------------------------------------------ fits
DRAWDIR <- file.path(OUT, "draws"); dir.create(DRAWDIR)
if (nzchar(RESUME)) {
  old <- readRDS(file.path(dirname(RESUME), "draw_indices.rds"))
  if (!identical(old$perm, perm) || !identical(old$subs, subs)) fail("resume draw indices differ")
  fs <- list.files(RESUME, pattern = "[.]rds$", full.names = TRUE)
  file.copy(fs, DRAWDIR)
  say("resume: copied ", length(fs), " finished draw files from ", RESUME)
}
compare <- function(lfc, padj) {
  cl <- is_call(lfc, padj)
  lost <- adopted & !cl; gained <- !adopted & cl
  ca <- adopted
  data.table(calls = sum(cl), up = sum(cl & lfc > 0), down = sum(cl & lfc < 0), kept = sum(adopted & cl),
    lost = sum(lost), gained = sum(gained), C = sum(lost) + sum(gained),
    lost_up = sum(lost & base_lfc > 0), lost_down = sum(lost & base_lfc < 0),
    lost_padj_only = sum(lost & padj >= 0.05 & abs(lfc) > 0.5),
    lost_lfc_only = sum(lost & padj < 0.05 & abs(lfc) <= 0.5),
    lost_both = sum(lost & padj >= 0.05 & abs(lfc) <= 0.5),
    lost_sign_flipped = sum(lost & sign(lfc) != sign(base_lfc)),
    lost_low_lfc = sum(lost & low_lfc),
    pearson_logFC_all = cor(base_lfc, lfc), spearman_logFC_all = cor(base_lfc, lfc, method = "spearman"),
    pearson_logFC_calls = cor(base_lfc[ca], lfc[ca]),
    spearman_logFC_calls = cor(base_lfc[ca], lfc[ca], method = "spearman"),
    slope_logFC_calls_no_intercept = unname(coef(lm(lfc[ca] ~ 0 + base_lfc[ca]))),
    n_na_padj = sum(is.na(padj)))
}
jobs <- c(list(list(id = "t1_true", kind = "t1"), list(id = "t2_true", kind = "sub"),
               list(id = "d2_no_mixed", kind = "d2"), list(id = "d3_pc1_only", kind = "d3")),
          lapply(seq_len(B_PERM), function(k) list(id = sprintf("perm_%03d", k), kind = "t1", k = k)),
          lapply(seq_len(B_SUB), function(k) list(id = sprintf("sub_%03d", k), kind = "sub", k = k)))
true_ids <- c("t1_true", "t2_true", "d2_no_mixed", "d3_pc1_only")
run_job <- function(j) {
  fn <- file.path(DRAWDIR, paste0(j$id, ".rds"))
  if (file.exists(fn)) return("skipped")
  t0 <- Sys.time()
  if (j$kind == "t1") {
    P <- if (is.null(j[["k"]])) PC_obs else PCM[perm[[j[["k"]]]], , drop = FALSE]
    f <- fit_one(dge, make_d(P), F_M1, aux = TRUE)
    units <- seq_len(N)
  } else if (j$kind == "sub") {
    units <- if (is.null(j[["k"]])) eur_idx else subs[[j[["k"]]]]
    f <- fit_one(dge[, units, keep.lib.sizes = TRUE], info[units, ], F_M0)
  } else if (j$kind == "d2") {
    units <- which(!info$possibly_mixed)
    f <- fit_one(dge[, units, keep.lib.sizes = TRUE], make_d(PC_obs)[units, ], F_M1)
  } else if (j$kind == "d3") {
    units <- seq_len(N)
    f <- fit_one(dge, make_d(PC_obs), F_D3)
  }
  cl <- is_call(f$logFC, f$padj)
  r <- list(id = j$id, kind = j$kind, k = if (is.null(j[["k"]])) NA_integer_ else j[["k"]], n_units = length(units),
            summary = compare(f$logFC, f$padj), call = cl, changed = cl != adopted,
            aux_bh05 = if (!is.null(f$aux_q)) sum(f$aux_q < 0.05, na.rm = TRUE) else NA_integer_,
            m1F_bh05 = if (!is.null(f$pcF_q)) sum(f$pcF_q < 0.05, na.rm = TRUE) else NA_integer_,
            seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  if (j$id %in% true_ids) r$fit <- f[setdiff(names(f), "E")]
  if (j$id == "t1_true") r$sw_units <- colnames(dge)[units]
  tmp <- paste0(fn, ".part")
  saveRDS(r, tmp); file.rename(tmp, fn)
  "done"
}
run_safe <- function(j) tryCatch(run_job(j), error = function(e) paste("ERROR:", conditionMessage(e)))
say("fitting ", length(jobs), " models on ", NCORES, " cores")
res <- mclapply(jobs, run_safe, mc.cores = NCORES, mc.preschedule = FALSE)
missing <- which(!file.exists(file.path(DRAWDIR, paste0(vapply(jobs, `[[`, "", "id"), ".rds"))))
if (length(missing)) {
  say("S5: ", length(missing), " fits failed; first error: ", as.character(res[[missing[1]]]), "; rerunning once")
  chk <- readRDS(draw_file)
  if (!identical(chk$perm, perm) || !identical(chk$subs, subs)) fail("stored draw indices differ from memory")
  res2 <- mclapply(jobs[missing], run_safe, mc.cores = NCORES, mc.preschedule = FALSE)
  still <- missing[!file.exists(file.path(DRAWDIR, paste0(vapply(jobs[missing], `[[`, "", "id"), ".rds")))]
  if (length(still)) stop_by_rule("S5", sprintf("%d fits failed twice (%s); first error: %s", length(still),
                                   paste(head(vapply(jobs[still], `[[`, "", "id"), 5), collapse = ","),
                                   as.character(res2[[match(still[1], missing)]])))
}
say("all fits finished")

rd <- function(id) readRDS(file.path(DRAWDIR, paste0(id, ".rds")))
perm_r <- lapply(sprintf("perm_%03d", seq_len(B_PERM)), rd)
sub_r <- lapply(sprintf("sub_%03d", seq_len(B_SUB)), rd)
t1 <- rd("t1_true"); t2 <- rd("t2_true"); d2 <- rd("d2_no_mixed"); d3 <- rd("d3_pc1_only")
perm_sum <- rbindlist(lapply(perm_r, function(r) data.table(draw = r$k, r$summary, aux_bh05 = r$aux_bh05,
                                                       m1F_bh05 = r$m1F_bh05, seconds = r$seconds)))
sub_sum <- rbindlist(lapply(sub_r, function(r) data.table(draw = r$k, r$summary, seconds = r$seconds)))
w(perm_sum, "t1_reference_draws.tsv")
w(sub_sum, "t2_reference_draws.tsv")

# ------------------------------------------------------------------ diagnostics (written before C1, CE or any p)
q <- function(x, p) unname(quantile(x, p))
aw <- lvqw$sample_weight[match(colnames(dge), lvqw$sample_id)]
m1_sw <- t1$fit$sw
sub_afr_ctrl <- vapply(subs, function(s) mean(info$ancestry_call[s][info$group_binary[s] == "Control"] == "AFR"), numeric(1))
me_afr_ctrl <- mean(info$ancestry_call[eur_idx][info$group_binary[eur_idx] == "Control"] == "AFR")
diag <- rbindlist(list(
  data.table(n = 1, diagnostic = "input sha256", metric = inputs$role, value = as.character(inputs$ok)),
  data.table(n = 2, diagnostic = "join", metric = s3dt$check, value = s3dt$value),
  data.table(n = 2, diagnostic = "ancestry by cohort by group",
             metric = paste(obs_tab$dataset, obs_tab$group_binary, sep = "_"),
             value = sprintf("AFR %d / EUR %d / other %d", obs_tab$AFR, obs_tab$EUR, obs_tab$other)),
  data.table(n = 3, diagnostic = "PC-0 reproduction", metric = names(repro), value = vapply(repro, function(x) format(x, digits = 6), "")),
  data.table(n = 4, diagnostic = "PC-A rs738409", metric = c(paste0("freq_", pca_freq$stratum, "_", pca_freq$ancestry_call), pca_tests$test),
             value = c(sprintf("G %.4f (n %d)", pca_freq$G_freq, pca_freq$n),
                       sprintf("OR %.3f (%.3f to %.3f) p %.3g", pca_tests$odds_ratio_AFR_vs_EUR, pca_tests$ci95_low, pca_tests$ci95_high, pca_tests$p))),
  data.table(n = 5, diagnostic = "design", metric = c(paste0("rank_", design_diag$model), "r2_disease_on_dataset_sex",
                                                     "r2_with_PCs", "variance_inflation"),
             value = c(sprintf("%d of %d", design_diag$rank, design_diag$ncol), format(r20, digits = 4),
                       format(r21, digits = 4), format((1 - r20) / (1 - r21), digits = 4))),
  data.table(n = 6, diagnostic = "within-cohort PC shift (disease minus control)", metric = pc_shift$pc,
             value = sprintf("%.4g (t %.2f)", pc_shift$disease_minus_control_within_cohort, pc_shift$t)),
  data.table(n = 7, diagnostic = "S7 permutation integrity",
             metric = c("max_abs_cell_mean_diff", "median_abs_cor_true_vs_perm_PC1", "max_abs_cor_true_vs_perm_PC1"),
             value = format(c(max(s7_dev), median(abs(s7_cor)), max(abs(s7_cor))), digits = 4)),
  data.table(n = 8, diagnostic = "quality weights",
             metric = c("cor_M1_vs_adopted", "cor_M0_vs_adopted",
                        paste0("median_M0_", c("AFR", "EUR", "other")), paste0("median_M1_", c("AFR", "EUR", "other"))),
             value = format(c(cor(m1_sw, aw), cor(m0_sw, aw),
                              vapply(c("AFR", "EUR", "other"), function(a) median(m0_sw[info$anc3 == a]), 0),
                              vapply(c("AFR", "EUR", "other"), function(a) median(m1_sw[info$anc3 == a]), 0)), digits = 4)),
  data.table(n = 9, diagnostic = "PC-B / NC-2 (auxiliary cell-means fits, BH q < 0.05 on PC F-test)",
             metric = c("aux_true", "aux_perm_median", "aux_perm_q95", "aux_perm_max",
                        "M1F_true", "M1F_perm_median", "M1F_perm_q95", "M1F_perm_max"),
             value = format(c(t1$aux_bh05, median(perm_sum$aux_bh05), q(perm_sum$aux_bh05, .95), max(perm_sum$aux_bh05),
                              t1$m1F_bh05, median(perm_sum$m1F_bh05), q(perm_sum$m1F_bh05, .95), max(perm_sum$m1F_bh05)), digits = 6)),
  data.table(n = 10, diagnostic = "T2 subset composition",
             metric = c("all_subsets_have_EUR_cell_counts", "mean_AFR_fraction_subset_controls", "AFR_fraction_ME_controls"),
             value = c(as.character(sub_counts_ok), format(mean(sub_afr_ctrl), digits = 4), format(me_afr_ctrl, digits = 4)))
), use.names = TRUE)
w(diag, "diagnostics.tsv")
say("DIAGNOSTICS (prespec order 1-10):")
for (i in seq_len(nrow(diag))) say(sprintf("  D%02d | %s | %s | %s", diag$n[i], diag$diagnostic[i], diag$metric[i], diag$value[i]))
if (!sub_counts_ok) fail("a T2 subset does not have the EUR cell counts")

# ------------------------------------------------------------------ primary F1
emp_p <- function(obs, ref) (1 + sum(ref >= obs)) / (length(ref) + 1)
min_reject <- function(ref, alpha) {
  k <- floor(alpha * (length(ref) + 1) + 1e-9) - 1
  if (k < 0) return(NA_real_)
  if (k + 1 > length(ref)) return(0)
  sort(ref, decreasing = TRUE)[k + 1] + 1
}
C1 <- t1$summary$C; CE <- t2$summary$C
p1 <- emp_p(C1, perm_sum$C); p2 <- emp_p(CE, sub_sum$C)
holm <- p.adjust(c(p1, p2), method = "holm")
prim <- data.table(
  test = c("T1 (M1 with true PCs vs within-cell permuted PCs)", "T2 (ME, EUR-only vs within-cell random subsets)"),
  observed_C = c(C1, CE), lost = c(t1$summary$lost, t2$summary$lost), gained = c(t1$summary$gained, t2$summary$gained),
  B = c(B_PERM, B_SUB), ref_median = c(median(perm_sum$C), median(sub_sum$C)),
  ref_q025 = c(q(perm_sum$C, .025), q(sub_sum$C, .025)), ref_q95 = c(q(perm_sum$C, .95), q(sub_sum$C, .95)),
  ref_q975 = c(q(perm_sum$C, .975), q(sub_sum$C, .975)), ref_max = c(max(perm_sum$C), max(sub_sum$C)),
  n_ref_ge_obs = c(sum(perm_sum$C >= C1), sum(sub_sum$C >= CE)),
  p_empirical_one_sided = c(p1, p2), p_floor = 1 / (c(B_PERM, B_SUB) + 1), p_holm = holm,
  reject_holm_005 = holm < 0.05,
  min_C_reject_holm_step1 = c(min_reject(perm_sum$C, 0.025), min_reject(sub_sum$C, 0.025)),
  min_C_reject_holm_step2 = c(min_reject(perm_sum$C, 0.05), min_reject(sub_sum$C, 0.05)))
w(prim, "primary_F1.tsv")
say(sprintf("T1: C1 %d (lost %d, gained %d); permuted median %.1f, q95 %.1f; p %.4g; Holm p %.4g",
            C1, t1$summary$lost, t1$summary$gained, prim$ref_median[1], prim$ref_q95[1], p1, holm[1]))
say(sprintf("T2: CE %d (lost %d, gained %d); subset median %.1f, q95 %.1f; p %.4g; Holm p %.4g",
            CE, t2$summary$lost, t2$summary$gained, prim$ref_median[2], prim$ref_q95[2], p2, holm[2]))

beside <- function(obs, ref, label) {
  mets <- c("calls", "kept", "lost", "gained", "C", "lost_padj_only", "lost_lfc_only", "lost_both", "lost_sign_flipped",
            "lost_low_lfc", "pearson_logFC_all", "spearman_logFC_all", "pearson_logFC_calls", "spearman_logFC_calls",
            "slope_logFC_calls_no_intercept")
  rbindlist(lapply(mets, function(m) data.table(arm = label, metric = m, observed = obs[[m]],
    ref_median = median(ref[[m]]), ref_q025 = q(ref[[m]], .025), ref_q975 = q(ref[[m]], .975))))
}
w(rbind(beside(t1$summary, perm_sum, "T1"), beside(t2$summary, sub_sum, "T2")), "beside_tests.tsv")

# ------------------------------------------------------------------ controls
aux_q95 <- q(perm_sum$aux_bh05, .95)
pcb_pass <- t1$aux_bh05 > aux_q95
nc2_med <- median(perm_sum$aux_bh05)
nc2_pass <- nc2_med <= 50
nc1_med <- median(perm_sum$C)
g_lfc_m1 <- if (!is.na(gstm1_row)) t1$fit$logFC[gstm1_row] else NA_real_
g_pc1 <- if (!is.na(gstm1_row)) t1$fit$pc1_coef[gstm1_row] else NA_real_
g_afr <- if (!is.null(gstm1_E)) mean(gstm1_E[info$group_binary == "Control" & info$ancestry_call == "AFR"]) else NA_real_
g_eur <- if (!is.null(gstm1_E)) mean(gstm1_E[info$group_binary == "Control" & info$ancestry_call == "EUR"]) else NA_real_
pcc_pass <- isTRUE(g_pc1 < 0 && abs(g_lfc_m1) < abs(base_lfc[gstm1_row]))
controls <- data.table(
  control = c("PC-0 reproduction (S2)", "PC-A rs738409 (S6)", "PC-B within-cell PC signal (S8)",
              "PC-C GSTM1 (not a stop rule)", "PC-C GSTM1 AFR vs EUR controls log2 CPM", "NC-1 permuted C median",
              "NC-2 permuted auxiliary count median (S8)"),
  prediction = c("1,854 calls (1,333 up, 521 down); max |dlog2FC| < 1e-6",
                 "allele CMH OR 0.451 (0.27 to 0.73), p 6.3e-04; pass if OR < 1 and p < 0.05",
                 "true count exceeds every permuted count; pass if above permuted 95th percentile",
                 "PC1 coefficient negative and |log2FC| below 1.440 in M1",
                 "AFR controls higher than EUR controls",
                 "[Estimate] median C between 20 and 200",
                 "[Estimate] median at most 5; fails if above 50"),
  observed = c(sprintf("%d calls (%d up, %d down); max |dlog2FC| %.3g", repro$calls, repro$up, repro$down, repro$max_abs_logFC_delta_vs_stored),
               sprintf("OR %.3f (%.3f to %.3f), p %.3g", pca$cmh_or, mh$conf.int[1], mh$conf.int[2], pca$cmh_p),
               sprintf("true %d; permuted median %.1f, q95 %.1f, max %d", t1$aux_bh05, nc2_med, aux_q95, max(perm_sum$aux_bh05)),
               sprintf("PC1 coefficient %.4g; M1 log2FC %.4f (adopted %.4f)", g_pc1, g_lfc_m1, base_lfc[gstm1_row]),
               sprintf("AFR %.3f, EUR %.3f", g_afr, g_eur),
               sprintf("median %.1f", nc1_med),
               sprintf("median %.1f", nc2_med)),
  outcome = c(if (DEBUG) "debug (not applied)" else "pass",
              "pass",
              if (pcb_pass) (if (t1$aux_bh05 > max(perm_sum$aux_bh05)) "pass (prediction met: above every permuted count)" else "pass (above q95, not above every permuted count)") else "fail",
              if (pcc_pass) "pass" else "fail",
              if (isTRUE(g_afr > g_eur)) "prediction met" else "prediction not met",
              if (nc1_med >= 20 && nc1_med <= 200) "prediction met" else "prediction not met",
              if (nc2_pass) (if (nc2_med <= 5) "pass (prediction met)" else "pass (median above predicted 5, not above 50)") else "fail"))
w(controls, "controls.tsv")

# ------------------------------------------------------------------ D1 per-gene table and flags
perm_changed <- vapply(perm_r, `[[`, logical(G), "changed")
sub_changed <- vapply(sub_r, `[[`, logical(G), "changed")
n_perm_ch <- rowSums(perm_changed); n_sub_ch <- rowSums(sub_changed)
k_perm <- 9L; k_sub <- 4L                                   # under 5% of 200 and of 100
if (DEBUG) { k_perm <- 0L; k_sub <- 0L }
m1_call <- t1$call; me_call <- t2$call
flag_m1 <- t1$changed & n_perm_ch <= k_perm
flag_me <- t2$changed & n_sub_ch <= k_sub
loo_perm <- vapply(seq_len(B_PERM), function(j) sum(perm_changed[, j] & (n_perm_ch - perm_changed[, j]) <= k_perm), numeric(1))
loo_sub <- vapply(seq_len(B_SUB), function(j) sum(sub_changed[, j] & (n_sub_ch - sub_changed[, j]) <= k_sub), numeric(1))
status <- function(cl) ifelse(adopted & cl, "kept", ifelse(adopted, "lost", ifelse(cl, "gained", "not_called")))
gt <- data.table(gene = rownames(dge), symbol = stored$symbol, AveExpr = stored$AveExpr,
  adopted_logFC = base_lfc, adopted_padj = if (DEBUG) pc0$padj else stored$padj, adopted_call = adopted,
  M1_logFC = t1$fit$logFC, M1_padj = t1$fit$padj, M1_status = status(m1_call),
  ME_logFC = t2$fit$logFC, ME_padj = t2$fit$padj, ME_status = status(me_call),
  M1_PC_F_q = t1$fit$pcF_q, aux_PC_F_q = t1$fit$aux_q, M1_PC1_coef = t1$fit$pc1_coef, M1_PC1_t = t1$fit$pc1_t,
  n_perm_changed = n_perm_ch, n_perm = B_PERM, n_sub_changed = n_sub_ch, n_sub = B_SUB,
  D1_rule_M1 = flag_m1, D1_rule_ME = flag_me)
fwrite(gt, file.path(OUT, "D1_gene_table.tsv.gz"), sep = "\t")
d1 <- data.table(arm = c("M1", "ME"), rule = c(sprintf("changed in M1 and in <= %d of %d permuted fits", k_perm, B_PERM),
                                               sprintf("changed in ME and in <= %d of %d subsets", k_sub, B_SUB)),
  observed_flags = c(sum(flag_m1), sum(flag_me)),
  expected_false_flags_median = c(median(loo_perm), median(loo_sub)),
  expected_false_flags_q95 = c(q(loo_perm, .95), q(loo_sub, .95)))
w(d1, "D1_flag_counts.tsv")
w(gt[D1_rule_M1 | D1_rule_ME][order(-D1_rule_M1, n_perm_changed)], "D1_flagged_genes.tsv")

# ------------------------------------------------------------------ D2, D3
w(rbind(data.table(fit = "D2 M1 without 16 possibly_mixed units", n_units = d2$n_units, d2$summary),
        data.table(fit = "D3 M1 with PC1 only", n_units = d3$n_units, d3$summary),
        data.table(fit = "M1 true (T1)", n_units = t1$n_units, t1$summary),
        data.table(fit = "ME (T2)", n_units = t2$n_units, t2$summary)), "D2_D3_sensitivity.tsv")

# ------------------------------------------------------------------ D4
d4a <- as.data.table(info)[, .(n_control = sum(group_binary == "Control"), n_disease = sum(group_binary == "Disease"),
  AFR_frac_control = mean(ancestry_call[group_binary == "Control"] == "AFR"),
  AFR_frac_disease = mean(ancestry_call[group_binary == "Disease"] == "AFR")), by = .(dataset)][order(dataset)]
w(d4a, "D4_afr_fraction_by_cohort.tsv")
top100 <- order(-abs(t1$fit$pc1_t))[1:min(100, G)]
d4b <- data.table(top_n_by_abs_PC1_t = length(top100), adopted_calls = sum(adopted[top100]),
                  M1_lost = sum(adopted[top100] & !m1_call[top100]), M1_gained = sum(!adopted[top100] & m1_call[top100]),
                  M1_lost_all_genes = sum(adopted & !m1_call))
w(d4b, "D4_top100_PC1_genes.tsv")

# ------------------------------------------------------------------ D5
conc <- fread(file.path(ROOT, CONC_REL))
base_id <- sub("[.][0-9]+$", "", rownames(dge))
i26 <- match(conc$ensembl_id, base_id)
d5 <- data.table(ensembl_id = conc$ensembl_id, symbol = conc$symbol, in_universe = !is.na(i26),
                 adopted_call = adopted[i26], M1_call = m1_call[i26], ME_call = me_call[i26])
w(d5, "D5_concordant26.tsv")
d5s <- d5[, .(n = .N, in_universe = sum(in_universe), adopted_called = sum(adopted_call, na.rm = TRUE),
              kept_M1 = sum(adopted_call & M1_call, na.rm = TRUE), kept_ME = sum(adopted_call & ME_call, na.rm = TRUE))]

# ------------------------------------------------------------------ D6
n_low <- sum(low_lfc)
d6 <- data.table(n_low_lfc_calls = n_low, M1_lost = t1$summary$lost_low_lfc,
                 M1_frac_lost = t1$summary$lost_low_lfc / n_low,
                 perm_median_frac_lost = median(perm_sum$lost_low_lfc) / n_low,
                 perm_q025_frac_lost = q(perm_sum$lost_low_lfc, .025) / n_low,
                 perm_q975_frac_lost = q(perm_sum$lost_low_lfc, .975) / n_low)
w(d6, "D6_low_lfc_losses.tsv")

# ------------------------------------------------------------------ reading (verdict computed from the prespec rule)
t1_rej <- holm[1] < 0.05; t2_rej <- holm[2] < 0.05
fmt <- function(x) formatC(x, format = "g", digits = 3)
m1v <- prim$ref_median[1]; q1v <- prim$ref_q95[1]; m2v <- prim$ref_median[2]; q2v <- prim$ref_q95[2]
if (t1_rej && pcb_pass && nc2_pass) {
  t1_case <- "T1 rejects (PC-B and NC-2 pass)"
  t1_text <- sprintf(paste("Adding four genetic ancestry PCs to the disease-versus-control model changed %d of the 1,854 calls",
    "(lost %d, gained %d), more than within-cohort-and-group permutations of the same PCs (median %s, 95th percentile %s;",
    "Holm p %s). %d genes met the ancestry-sensitive rule (expected false flags from the reference: median %s, 95th percentile %s)."),
    C1, t1$summary$lost, t1$summary$gained, fmt(m1v), fmt(q1v), fmt(holm[1]), sum(flag_m1), fmt(d1$expected_false_flags_median[1]),
    fmt(d1$expected_false_flags_q95[1]))
  flags_apply <- TRUE
} else if (t1_rej) {
  t1_case <- "T1 rejects but PC-B or NC-2 fails"
  t1_text <- sprintf(paste("Adding four genetic ancestry PCs changed %d of the 1,854 calls, more than permuted PCs (Holm p %s),",
    "but %s, so the excess change is not attributed to ancestry. No genes are flagged as ancestry-sensitive and the call set",
    "stands as adopted."), C1, fmt(holm[1]),
    if (!nc2_pass) "the check of whether the PCs carry within-cell expression signal could not be read" else
      "the PCs showed no within-cell expression signal beyond permuted PCs")
  flags_apply <- FALSE
} else if (pcb_pass && nc2_pass) {
  t1_case <- "T1 does not reject, PC-B passes, NC-2 passes"
  t1_text <- sprintf(paste("Adjusting for four genetic ancestry PCs changed %d of the 1,854 calls (lost %d, gained %d), within the range",
    "of fits in which the same PCs were permuted within cohort and disease group (median %s, 95th percentile %s; Holm p %s).",
    "A change smaller than %s calls (95th percentile minus median) is not excluded. The established-state call set stands as adopted."),
    C1, t1$summary$lost, t1$summary$gained, fmt(m1v), fmt(q1v), fmt(holm[1]), fmt(q1v - m1v))
  flags_apply <- FALSE
} else if (nc2_pass) {
  t1_case <- "T1 does not reject and PC-B fails (NC-2 passes)"
  t1_text <- paste("Ancestry PCs from RNA genotypes carried no detectable within-cell expression signal in these units,",
                   "so this check is uninformative about ancestry effects on the calls.")
  flags_apply <- FALSE
} else {
  t1_case <- "T1 does not reject and NC-2 fails"
  t1_text <- paste("The test of whether the ancestry PCs carry expression signal was not calibrated in these units,",
                   "so this check is uninformative about ancestry effects on the calls.")
  flags_apply <- FALSE
}
if (t2_rej) {
  t2_case <- "T2 rejects"
  t2_text <- sprintf(paste("Restricting to EUR-called units changed %d calls (lost %d, gained %d), more than random subsets of the same",
    "cohort and group sizes (median %s, 95th percentile %s; Holm p %s)."), CE, t2$summary$lost, t2$summary$gained,
    fmt(m2v), fmt(q2v), fmt(holm[2]))
  if (!t1_rej) t2_text <- paste(t2_text, "Read with T1: the change follows from removing units, not from ancestry itself",
                                "being in the model; the EUR-only arm differs while PC adjustment does not.")
} else {
  t2_case <- "T2 does not reject"
  t2_text <- sprintf(paste("Restricting to %d EUR-called units (%d control, %d disease) changed %d calls (lost %d, gained %d),",
    "within the range of random subsets of the same composition (median %s, 95th percentile %s; Holm p %s). A change smaller than",
    "%s calls (95th percentile minus median) is not excluded."), length(eur_idx), sum(info$group_binary[eur_idx] == "Control"),
    sum(info$group_binary[eur_idx] == "Disease"), CE, t2$summary$lost, t2$summary$gained, fmt(m2v), fmt(q2v), fmt(holm[2]), fmt(q2v - m2v))
}
verdict <- sprintf("%s; %s", t1_case, t2_case)
say("VERDICT: ", verdict)

# ------------------------------------------------------------------ D7 panels
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })
FIGDIR <- file.path(OUT, "figures"); dir.create(FIGDIR)
FS <- 6
th <- theme_classic(base_size = FS, base_family = "Helvetica") +
  theme(text = element_text(size = FS, colour = "black", face = "plain"),
        axis.text = element_text(size = FS, colour = "black"), axis.title = element_text(size = FS, colour = "black"),
        legend.text = element_text(size = FS, colour = "black"), legend.title = element_text(size = FS, colour = "black"),
        axis.line = element_line(linewidth = 0.3, colour = "black"), axis.ticks = element_line(linewidth = 0.3, colour = "black"),
        legend.key.size = unit(0.25, "cm"))
save_pdf <- function(p, name, wd, ht) {
  fn <- file.path(FIGDIR, name)
  if (file.exists(fn)) fail("Refusing to overwrite ", fn)
  pdf(fn, width = wd, height = ht, useDingbats = FALSE, family = "Helvetica")
  print(p); invisible(dev.off())
}
hist_panel <- function(ref, obs, xlab, ylab) {
  bw <- max(1, ceiling(diff(range(c(ref, obs))) / 30))
  ggplot(data.frame(C = ref), aes(C)) +
    geom_histogram(binwidth = bw, fill = "#9E9E9E", colour = NA, boundary = 0) +
    geom_vline(xintercept = obs, colour = "black", linewidth = 0.4) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = xlab, y = ylab) + th
}
save_pdf(hist_panel(perm_sum$C, C1, "Calls changed vs the adopted 1,854 (lost + gained)", "Fits with permuted PCs"),
         "D7a_T1_changed_calls_vs_permuted.pdf", 2.6, 2.0)
save_pdf(hist_panel(sub_sum$C, CE, "Calls changed vs the adopted 1,854 (lost + gained)", "Random-subset fits"),
         "D7b_T2_changed_calls_vs_subsets.pdf", 2.6, 2.0)
sc <- data.table(adopted = base_lfc, m1 = t1$fit$logFC, status = status(m1_call), symbol = gt$symbol, flag = flag_m1)
sc[, status := factor(status, levels = c("not_called", "kept", "lost", "gained"))]
setorder(sc, status)
cols <- c(not_called = "#BDBDBD", kept = "#424242", lost = "#1565C0", gained = "#C9265E")
pt <- geom_point(size = 0.25, stroke = 0, alpha = 0.8)
if (requireNamespace("ggrastr", quietly = TRUE)) pt <- ggrastr::rasterise(pt, dpi = 300)
p7c <- ggplot(sc, aes(adopted, m1, colour = status)) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, colour = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, colour = "black") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.2, colour = "#9E9E9E") + pt +
  scale_colour_manual(values = cols, name = "Call status in M1") +
  geom_text_repel(data = sc[flag == TRUE], aes(label = symbol), colour = "black", size = FS / .pt,
                  seed = 20260923, max.overlaps = Inf, segment.size = 0.15, min.segment.length = 0) +
  labs(x = "Adopted log2 fold change (disease vs control)", y = "log2 fold change with ancestry PC1-PC4") + th
save_pdf(p7c, "D7c_adopted_vs_M1_log2FC.pdf", 3.6, 3.0)
captions <- c(
  sprintf("D7a. Number of the adopted 1,854 calls that change status (lost plus gained) in %d refits of M1 with PC1-PC4 rows permuted within the 10 cohort-by-group cells (gray). Black line: the fit with each unit's own PCs (C1 = %d).", B_PERM, C1),
  sprintf("D7b. The same count for %d M0 fits on random subsets with the EUR cohort-by-group counts (gray). Black line: the EUR-only fit (CE = %d).", B_SUB, CE),
  sprintf("D7c. Disease log2 fold change in the adopted fit against the fit with PC1-PC4, all %d genes. Dashed lines: |log2FC| = 0.5. Labels: genes meeting the D1 rule (changed in M1 and in at most %d of %d permuted fits); %s.", G, k_perm, B_PERM,
          if (flags_apply) "the reading permits the ancestry-sensitive label" else "the reading does not permit the ancestry-sensitive label; labels are descriptive"))
writeLines(captions, file.path(FIGDIR, "captions.txt"))
for (cp in captions) message(cp)

# ------------------------------------------------------------------ RESULT_H79.md
fmtp <- function(x) formatC(x, format = "g", digits = 4)
res_md <- c(
  "# RESULT_H79", "",
  if (DEBUG) "**DEBUG RUN: code-path test on every 10th gene, with the observed PCs and EUR set replaced by null draws. Not the analysis.**" else NULL,
  if (DEBUG) "" else NULL,
  sprintf("Written %s by `h79_ancestry_refit.R` (prespec sha256 86d21b24...c632). Seed 20260923, L'Ecuyer-CMRG.", now()), "",
  "Status: completed", "",
  sprintf("Verdict (prespec Reading rule, computed in code): %s.", verdict), "",
  "## Reading", "", t1_text, "", t2_text, "",
  "Whatever the result, GSE135251 per-cohort membership in Figure 3C cannot be ancestry matched (only 2 EUR controls).", "",
  "## Family F1 = {T1, T2}, Holm at family-wise alpha 0.05", "",
  "| Test | C | lost | gained | B | ref median | ref q95 | #ref >= obs | p (one-sided) | Holm p | reject | min C to reject (step 1 / step 2) |",
  "|---|---|---|---|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %d | %d | %d | %d | %s | %s | %d | %s | %s | %s | %s / %s |", c("T1", "T2"), prim$observed_C, prim$lost, prim$gained,
          prim$B, fmt(prim$ref_median), fmt(prim$ref_q95), prim$n_ref_ge_obs, fmtp(prim$p_empirical_one_sided), fmtp(prim$p_holm),
          prim$reject_holm_005, prim$min_C_reject_holm_step1, prim$min_C_reject_holm_step2), "",
  "Reference 2.5 to 97.5 percentiles of C: " ,
  sprintf("- T1: %s to %s; T2: %s to %s.", fmt(prim$ref_q025[1]), fmt(prim$ref_q975[1]), fmt(prim$ref_q025[2]), fmt(prim$ref_q975[2])), "",
  "## Controls", "",
  "| Control | Prediction | Observed | Outcome |", "|---|---|---|---|",
  sprintf("| %s | %s | %s | %s |", controls$control, controls$prediction, controls$observed, controls$outcome), "",
  sprintf("S8 (PC-B and NC-2): PC-B %s, NC-2 %s.", if (pcb_pass) "passes" else "fails", if (nc2_pass) "passes" else "fails"), "",
  "## Descriptive outputs (nominal; not tests)", "",
  sprintf("- D1: M1 rule met by %d genes (expected false flags median %s, 95th percentile %s); ME rule met by %d genes (expected false flags median %s, 95th percentile %s). %s",
          d1$observed_flags[1], fmt(d1$expected_false_flags_median[1]), fmt(d1$expected_false_flags_q95[1]),
          d1$observed_flags[2], fmt(d1$expected_false_flags_median[2]), fmt(d1$expected_false_flags_q95[2]),
          if (flags_apply) "The reading permits the ancestry-sensitive label." else "The reading does not permit the ancestry-sensitive label; these counts are descriptive only."),
  sprintf("- D2 (M1 without 16 possibly_mixed units, n %d): C %d (lost %d, gained %d).", d2$n_units, d2$summary$C, d2$summary$lost, d2$summary$gained),
  sprintf("- D3 (M1 with PC1 only): C %d (lost %d, gained %d).", d3$summary$C, d3$summary$lost, d3$summary$gained),
  sprintf("- D4: among the top 100 genes by |t| of the PC1 coefficient in M1, %d are adopted calls and %d are M1-lost calls.", d4b$adopted_calls, d4b$M1_lost),
  sprintf("- D5: of the 26 shared genes, %d in the universe, %d adopted calls, %d kept in M1, %d kept in ME.", d5s$in_universe, d5s$adopted_called, d5s$kept_M1, d5s$kept_ME),
  sprintf("- D6: of %d calls with adopted |log2FC| < 0.6, M1 lost %d (%.3f); permuted median fraction %.3f (2.5 to 97.5 percentile %.3f to %.3f).",
          n_low, d6$M1_lost, d6$M1_frac_lost, d6$perm_median_frac_lost, d6$perm_q025_frac_lost, d6$perm_q975_frac_lost),
  sprintf("- T1 log2FC agreement over the adopted calls: Pearson %.4f, slope %.4f; T2: Pearson %.4f, slope %.4f.",
          t1$summary$pearson_logFC_calls, t1$summary$slope_logFC_calls_no_intercept, t2$summary$pearson_logFC_calls, t2$summary$slope_logFC_calls_no_intercept),
  "- D7 panels: `figures/D7a_T1_changed_calls_vs_permuted.pdf`, `figures/D7b_T2_changed_calls_vs_subsets.pdf`, `figures/D7c_adopted_vs_M1_log2FC.pdf`.", "",
  "Unit of inference: 774 genetic-identity clusters (138 control, 636 disease) in five cohorts; genes are not replicates.",
  "Covariates: dataset, inferred_sex, group (M0); plus PC1-PC4 (M1). Per-gene BH over the fixed gene universe inside each fit is the call rule, not a test of H79.",
  "", "Files: diagnostics.tsv, primary_F1.tsv, beside_tests.tsv, controls.tsv, D1_gene_table.tsv.gz, D1_flag_counts.tsv, D1_flagged_genes.tsv, D2_D3_sensitivity.tsv, D4_*.tsv, D5_concordant26.tsv, D6_low_lfc_losses.tsv, t1_reference_draws.tsv, t2_reference_draws.tsv, draw_indices.rds, draws/.")
wlines(res_md, "RESULT_H79.md")
wlines(capture.output(sessionInfo()), "sessionInfo.txt")
if (!DEBUG) {
  root_res <- file.path(H79, "RESULT_H79.md")
  if (!file.exists(root_res)) file.copy(file.path(OUT, "RESULT_H79.md"), root_res)
}
wlines("COMPLETE", "COMPLETE")
say("H79 complete")
