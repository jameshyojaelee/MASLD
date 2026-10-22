#!/usr/bin/env Rscript
# H56 stage 2: fixed-weight composition-conditional refits, controls, labels, verdict.
# Implements PRESPEC_H56.md (sha256 cd786c6d...). Reads stage-1 outputs written by h56_prep.py.
# H56_MODE=debug replaces every real coordinate set by a within-cohort row permutation
# (null coordinates), uses tiny B, and writes to a debug directory.
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table); library(parallel)})

RUN <- Sys.getenv("H56_RUN"); OUT <- Sys.getenv("H56_OUT"); MODE <- Sys.getenv("H56_MODE", "final")
DEBUG <- MODE == "debug"
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
P <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
PRESPEC_SHA <- "cd786c6d17f60f2f88a7e3e11389820264cfa09e51e478d610516f38a847cf6d"
SEED <- 20260923L
NCORE <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
B_BOOT <- as.integer(Sys.getenv("H56_B", if (DEBUG) "4" else "1000"))
N_NULL <- as.integer(Sys.getenv("H56_NNULL", if (DEBUG) "2" else "10"))
SKIP_VQW <- DEBUG && Sys.getenv("H56_DEBUG_SKIP_VQW", "0") == "1"
FLAG <- c("A16_rebuilt", "A16_earlier", "C7_rebuilt", "C7_earlier", "D1_rebuilt", "D1_earlier",
          "DRmain_rebuilt", "DRmain_earlier")
GRP <- "group_binaryDisease"

logf <- file(file.path(OUT, "fit_log.txt"), open = "wt")
say <- function(...) { s <- paste0(...); cat(s, "\n"); writeLines(s, logf); flush(logf) }
sayt <- function(x) { s <- capture.output(print(x, row.names = FALSE)); cat(s, sep = "\n"); writeLines(s, logf); flush(logf) }
state <- list(mode = MODE, stops = list())
write_state <- function() writeLines(jsonlite::toJSON(state, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null"),
                                     file.path(OUT, "state.json"))
say("H56 fit stage; mode ", MODE, "; cores ", NCORE, "; B ", B_BOOT, "; null draws ", N_NULL, "; ", format(Sys.time(), tz = "UTC"))

# ------------------------------------------------------------------ S1 (recheck) and prespec hash
sha <- function(p) strsplit(system2("sha256sum", shQuote(p), stdout = TRUE), " ")[[1]][1]
if (sha(file.path(RUN, "PRESPEC_H56.md")) != PRESPEC_SHA) stop("prespec hash mismatch")
ht <- fread(file.path(RUN, "input_sha256.tsv"), header = FALSE, col.names = c("sha", "bytes", "path"))
ht[, got := vapply(file.path(ROOT, path), sha, "")]
if (nrow(ht) != 31L || any(ht$got != ht$sha)) {
  state$stops$S1 <- ht[got != sha, path]; write_state(); stop("S1 hash mismatch")
}
say("S1 recheck in R: 31 of 31 input hashes OK")
if (file.exists(file.path(OUT, "STOP_PREP.json"))) stop("stage 1 stopped; see STOP_PREP.json")

# ------------------------------------------------------------------ inputs
d <- readRDS(file.path(P, "merged_dge.rds"))
dz <- fread(file.path(P, "model_design.tsv"))
units <- fread(file.path(OUT, "units.tsv"))
cells <- fread(file.path(OUT, "cells.tsv"))
cls <- fread(file.path(OUT, "call_classes.tsv"))
prop <- fread(file.path(OUT, "propensity.tsv.gz"))
deg <- fread(file.path(P, "deg_results.csv"))
sw_stored <- fread(file.path(P, "lvqw_sample_weights.tsv"))
cells[, `:=`(estimable = as.logical(estimable), flag_cell = as.logical(flag_cell))]
units[, `:=`(in_DRmain = as.logical(in_DRmain), in_DRsep = as.logical(in_DRsep))]
ids <- colnames(d)
stopifnot(setequal(ids, dz$sample_id), setequal(ids, units$sample_id), nrow(dz) == 774L)
dz <- dz[match(ids, sample_id)]; units <- units[match(ids, sample_id)]; prop <- prop[match(ids, sample_id)]
sw_stored <- sw_stored[match(ids, sample_id)]
X774 <- as.matrix(dz[, -1]); rownames(X774) <- ids
meta <- as.data.table(readRDS(file.path(P, "meta_matched.rds")))[match(ids, sample_id)]
mm <- model.matrix(~ dataset + inferred_sex + group_binary,
                   data.frame(dataset = factor(meta$dataset), inferred_sex = factor(meta$inferred_sex),
                              group_binary = factor(meta$group_binary, levels = c("Control", "Disease"))))
design_same <- max(abs(unclass(mm)[, seq_len(ncol(mm))] - X774)) == 0 && identical(colnames(mm), colnames(X774))
say("model_design.tsv equals model.matrix(~ dataset + inferred_sex + group_binary) from meta_matched: ", design_same)
stopifnot(design_same)
coords <- list()
for (cl in cells$cell) {
  z <- fread(file.path(OUT, sprintf("coords_%s.tsv.gz", cl)))
  z <- as.matrix(z[match(ids, sample_id), -1]); rownames(z) <- ids
  coords[[cl]] <- z
}
dset <- units$dataset; grp <- units$group
if (DEBUG) {
  set.seed(777L)  # debug only: null coordinates
  for (cl in names(coords)) {
    perm <- seq_along(ids)
    for (g in unique(dset)) { i <- which(dset == g); perm[i] <- i[sample.int(length(i))] }
    z <- coords[[cl]][perm, , drop = FALSE]; rownames(z) <- ids; coords[[cl]] <- z
  }
  say("DEBUG: coordinates replaced by within-cohort row permutations")
}
calls <- deg[padj < 0.05 & abs(logFC) > 0.5]
stopifnot(nrow(calls) == 1854L)
calls <- merge(calls[, .(gene, symbol, logFC_stored = logFC, SE_stored = SE)], cls[, .(gene, sc_class, sc_class_detail, stromal_ecm)],
               by = "gene", all.x = TRUE, sort = FALSE)
calls[, s := sign(logFC_stored)]

# ------------------------------------------------------------------ S2: fixed-weight refit reproduces the stored fit
cache <- file.path(RUN, "tmp", "voom_refit.rds")
if (DEBUG && file.exists(cache)) {
  v <- readRDS(cache); say("DEBUG: reusing cached refit")
} else {
  t0 <- Sys.time()
  v <- voomWithQualityWeights(d, X774)
  say("voomWithQualityWeights refit: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  if (DEBUG) saveRDS(v, cache)
}
fit0 <- eBayes(lmFit(v, X774))
j0 <- match(GRP, colnames(X774))
ref <- data.table(gene = rownames(fit0), logFC = fit0$coefficients[, j0],
                  SE = fit0$stdev.unscaled[, j0] * sqrt(fit0$s2.post), P = fit0$p.value[, j0])
ref[, padj := p.adjust(P, "BH")]
cmp <- merge(ref, deg[, .(gene, logFC_s = logFC, SE_s = SE, padj_s = padj)], by = "gene")
swd <- max(abs(v$targets$sample.weights - sw_stored$sample_weight))
refcalls <- ref[padj < 0.05 & abs(logFC) > 0.5, gene]
s2 <- list(n_genes = nrow(cmp), max_abs_logFC = max(abs(cmp$logFC - cmp$logFC_s)), max_abs_SE = max(abs(cmp$SE - cmp$SE_s)),
           max_abs_padj = max(abs(cmp$padj - cmp$padj_s)), max_abs_sample_weight = swd,
           same_calls = setequal(refcalls, calls$gene), n_calls_refit = length(refcalls))
state$S2 <- s2
say("Diagnostic 1 / S2: ", paste(names(s2), unlist(s2), sep = "=", collapse = "; "))
S2_pass <- s2$n_genes == 24115 && s2$max_abs_logFC < 1e-6 && s2$max_abs_SE < 1e-6 && s2$max_abs_padj < 1e-6 &&
  swd < 1e-8 && s2$same_calls
state$S2_pass <- S2_pass; write_state()
if (!S2_pass) {
  state$verdict <- "untestable"; state$stops$S2 <- "refit does not reproduce the stored fit"; write_state()
  say("STOP S2"); quit(status = 4)
}
E <- v$E; W <- v$weights; dimnames(W) <- dimnames(E)
stopifnot(identical(colnames(E), ids))
gene_idx_calls <- match(calls$gene, rownames(E))
rm(fit0)

# ------------------------------------------------------------------ design helpers
cov_of <- function(cl) strsplit(cells[cell == cl, covariate_cols], ";")[[1]]
units_of <- function(arm) switch(arm, A16 = , C7 = , D1 = rep(TRUE, 774), DRmain = units$in_DRmain, DRsep = units$in_DRsep,
                                 U774 = rep(TRUE, 774), U537 = units$in_DRmain, U237 = units$in_DRsep)
comp_of <- function(cl) cells[cell == cl, comparator]
comp_cov <- list(U774 = cov_of("A16_rebuilt"), U537 = cov_of("DRmain_rebuilt"), U237 = cov_of("DRsep_rebuilt"))
design_u <- function(cov, rows) cbind(X774[rows, cov, drop = FALSE], group_binaryDisease = X774[rows, GRP])
design_a <- function(cov, Z, rows) cbind(X774[rows, cov, drop = FALSE], Z[rows, , drop = FALSE], group_binaryDisease = X774[rows, GRP])
fullrank <- function(X) qr(X)$rank == ncol(X)
fitg <- function(rows, X, genes = NULL, eb = TRUE) {
  Ei <- E[, rows, drop = FALSE]; Wi <- W[, rows, drop = FALSE]
  if (!is.null(genes)) { Ei <- Ei[genes, , drop = FALSE]; Wi <- Wi[genes, , drop = FALSE] }
  f <- lmFit(Ei, X, weights = Wi)
  j <- match(GRP, colnames(X))
  if (!eb) return(list(b = f$coefficients[, j]))
  e <- eBayes(f)
  list(b = e$coefficients[, j], u = e$stdev.unscaled[, j], s2 = e$s2.post, df = e$df.total, p = e$p.value[, j])
}
infl <- function(C0, Z, yv, wv) {
  rv <- function(X) { r <- lm.wfit(X, yv, wv)$residuals; sum(wv * r^2) }
  sqrt(rv(C0) / rv(cbind(C0, Z)))
}

# ------------------------------------------------------------------ S4, S5
cells[, fitted := estimable]
for (cl in cells$cell) {
  rows <- units_of(cells[cell == cl, arm])
  if (!fullrank(design_a(cov_of(cl), coords[[cl]], rows))) cells[cell == cl, fitted := FALSE]
}
flag_ok <- cells[cell %in% FLAG & fitted == TRUE, cell]
state$flag_cells_estimable <- flag_ok
say("S4/S5: estimable flag cells ", length(flag_ok), " of 8: ", paste(flag_ok, collapse = ", "))
flag_untestable <- length(flag_ok) < 6
if (flag_untestable) say("S5: fewer than six flag cells estimable; flag untestable")
write_state()

# ------------------------------------------------------------------ diagnostic 4: top-100 share on coordinates
cnt <- d$counts
t100 <- apply(cnt, 2, function(x) sum(sort(x, decreasing = TRUE)[1:100])) / colSums(cnt)
d4 <- rbindlist(lapply(cells$cell, function(cl) {
  rows <- units_of(cells[cell == cl, arm]); ds <- factor(dset[rows])
  f0 <- lm(t100[rows] ~ ds); f1 <- lm(t100[rows] ~ ds + coords[[cl]][rows, ])
  data.table(cell = cl, partial_R2_top100_on_coords = 1 - sum(resid(f1)^2) / sum(resid(f0)^2))
}))
say("Diagnostic 4: partial R2 of top-100 share on coordinates given dataset"); sayt(d4)
fwrite(d4, file.path(OUT, "diag4_complexity_r2.tsv"), sep = "\t")

# ------------------------------------------------------------------ diagnostic 6: single-cell benchmark (NC2 inputs)
source(file.path(ROOT, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
sc_res <- local({
  F <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk")
  hp <- fread(cmd = paste("zcat", shQuote(file.path(F, "pseudobulk_Hepatocytes.tsv.gz"))))
  m <- as.matrix(hp[, -1]); rownames(m) <- hp[[1]]
  m <- rowsum(m, rownames(m), reorder = FALSE)
  m <- collapse_counts_to_donor(m, ROOT)
  meta_sc <- fread(file.path(F, "sample_celltype_metadata.tsv"))[cell_type == "Hepatocytes"]
  map <- build_srr_to_donor_map(ROOT)
  meta_sc[, donor := ifelse(sample %in% names(map), map[sample], sample)]
  nc <- meta_sc[, .(n_cells = sum(n_cells)), by = .(donor, dataset)]
  ros <- fread(file.path(ROOT, "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/current_source_diagnosis_donor_roster.tsv"))
  ros <- ros[as.character(exclude) %in% c("FALSE", "False")]
  r <- merge(ros, nc, by = c("donor", "dataset"))
  r <- r[donor %in% colnames(m) & n_cells >= 50]
  r[, disease := as.integer(disease_stage_coarse != "Healthy")]
  n_h <- sum(r$disease == 0); n_d <- sum(r$disease == 1)
  out <- list(n_donors = nrow(r), n_healthy = n_h, n_disease = n_d,
              by_dataset = r[, .(healthy = sum(disease == 0), disease = sum(disease == 1)), by = dataset])
  if (nrow(r) != 70 || n_h != 12 || n_d != 58) { out$testable <- FALSE; out$reason <- "S8 donor count"; return(out) }
  y <- DGEList(counts = m[, r$donor])
  des <- model.matrix(~ factor(r$dataset) + r$disease); colnames(des)[ncol(des)] <- "disease"
  keep <- filterByExpr(y, group = r$disease)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  fv <- eBayes(lmFit(voom(y, des), des))
  tt <- data.table(symbol = rownames(fv), logFC_hep = fv$coefficients[, "disease"], p = fv$p.value[, "disease"])
  tt[, q := p.adjust(p, "BH")]
  hepc <- calls[sc_class == "hepatocyte"]
  j <- merge(hepc[, .(gene, symbol, logFC_stored)], tt, by = "symbol")
  hplus <- j[sign(logFC_hep) == sign(logFC_stored) & q < 0.10, gene]
  c(out, list(testable = length(hplus) >= 10, genes_tested = nrow(tt), genes_q10 = tt[q < 0.10, .N],
              hep_calls_tested = nrow(j), n_Hplus = length(hplus), Hplus = list(hplus), tab = list(tt)))
})
say("Diagnostic 6: single-cell hepatocyte benchmark donors ", sc_res$n_donors, " (healthy ", sc_res$n_healthy,
    ", disease ", sc_res$n_disease, "); genes tested ", sc_res$genes_tested, "; genes q<0.10 ", sc_res$genes_q10,
    "; hepatocyte-specific calls tested ", sc_res$hep_calls_tested, "; H+ ", sc_res$n_Hplus, "; NC2 testable ", sc_res$testable)
sayt(sc_res$by_dataset)
if (!is.null(sc_res$tab)) fwrite(sc_res$tab[[1]], file.path(OUT, "nc2_hepatocyte_sc_de.tsv.gz"), sep = "\t")
Hplus <- if (isTRUE(sc_res$testable)) sc_res$Hplus[[1]] else character(0)
state$NC2_inputs <- sc_res[c("n_donors", "n_healthy", "n_disease", "testable", "genes_tested", "genes_q10", "hep_calls_tested", "n_Hplus")]
write_state()

# ------------------------------------------------------------------ comparator fits (full gene set)
comp <- list()
for (cn in names(comp_cov)) {
  rows <- units_of(cn); comp[[cn]] <- fitg(which(rows), design_u(comp_cov[[cn]], rows))
}
say("U774 refit vs stored logFC max abs diff ", signif(max(abs(comp$U774$b[gene_idx_calls] - calls$logFC_stored)), 3))

# ------------------------------------------------------------------ diagnostic 9: calibration null (S11)
set.seed(SEED)
perms <- list()
for (cl in flag_ok) {
  rows <- which(units_of(cells[cell == cl, arm])); dd <- dset[rows]
  perms[[cl]] <- lapply(seq_len(N_NULL), function(b) {
    p <- seq_along(rows); for (g in unique(dd)) { i <- which(dd == g); p[i] <- i[sample.int(length(i))] }; p })
}
jobs <- expand.grid(cell = flag_ok, draw = seq_len(N_NULL), stringsAsFactors = FALSE)
nullres <- mclapply(seq_len(nrow(jobs)), function(k) {
  cl <- jobs$cell[k]; b <- jobs$draw[k]
  rlog <- units_of(cells[cell == cl, arm]); rows <- which(rlog)
  C <- design_u(cov_of(cl), rlog); Z <- coords[[cl]][rows, , drop = FALSE]
  wv <- sw_stored$sample_weight[rows]
  lf <- lm.wfit(C, Z, wv); Fv <- as.matrix(lf$fitted.values); R <- as.matrix(lf$residuals)
  Zp <- Fv + R[perms[[cl]][[b]], , drop = FALSE]
  colnames(Zp) <- colnames(Z); rownames(Zp) <- ids[rows]
  Xa <- cbind(C[, colnames(C) != GRP, drop = FALSE], Zp, group_binaryDisease = C[, GRP])
  fa <- fitg(rows, Xa); fu <- comp[[comp_of(cl)]]
  D <- fa$b - fu$b; se <- sqrt(fa$s2 * (fa$u^2 - fu$u^2)); tt <- D / se
  pp <- 2 * pt(-abs(tt), fa$df)
  data.table(cell = cl, draw = b, share_p05 = mean(pp < 0.05), lambda = median(tt^2) / 0.455,
             se_inflation = infl(C[, colnames(C) != GRP, drop = FALSE], Zp, C[, GRP], wv))
}, mc.cores = NCORE, mc.preschedule = FALSE)
if (any(vapply(nullres, inherits, TRUE, "try-error"))) stop("calibration null failed: ", paste(nullres[vapply(nullres, inherits, TRUE, "try-error")], collapse = " | "))
nullres <- rbindlist(nullres)
fwrite(nullres, file.path(OUT, "diag9_calibration_null_draws.tsv"), sep = "\t")
d9 <- nullres[, .(mean_share_p05 = mean(share_p05), mean_lambda = mean(lambda), mean_se_inflation = mean(se_inflation)), by = cell]
d9[, S11_stop := mean_share_p05 > 0.075 | mean_lambda > 1.2]
d9[, conservative := mean_share_p05 < 0.025 | mean_lambda < 0.8]
d9[, within_prediction := mean_share_p05 >= 0.04 & mean_share_p05 <= 0.06 & mean_lambda >= 0.9 & mean_lambda <= 1.1]
say("Diagnostic 9: calibration null (", N_NULL, " draws per cell)"); sayt(d9)
fwrite(d9, file.path(OUT, "diag9_calibration_null.tsv"), sep = "\t")
S11_stop <- any(d9$S11_stop)
state$S11_stop <- S11_stop; state$diag9 <- d9; write_state()

# ------------------------------------------------------------------ real cell fits (per-gene adjusted values read from here on)
say("Diagnostics 1 to 9 printed above. Fitting the real cells.")
real <- list()
for (cl in cells[fitted == TRUE, cell]) {
  rlog <- units_of(cells[cell == cl, arm]); rows <- which(rlog)
  real[[cl]] <- fitg(rows, design_a(cov_of(cl), coords[[cl]], rlog))
}
# diagnostic 5: variance identity on 20 seeded genes per cell
set.seed(SEED)
d5 <- rbindlist(lapply(names(real), function(cl) {
  rlog <- units_of(cells[cell == cl, arm]); rows <- which(rlog)
  Xa <- design_a(cov_of(cl), coords[[cl]], rlog); Xu <- design_u(comp_cov[[comp_of(cl)]], rlog)
  gs <- sample.int(nrow(E), 20)
  rbindlist(lapply(gs, function(g) {
    w <- W[g, rows]; yv <- E[g, rows]
    Aa <- solve(crossprod(Xa, w * Xa), t(w * Xa))[GRP, ]; Au <- solve(crossprod(Xu, w * Xu), t(w * Xu))[GRP, ]
    cvec <- Aa - Au; vd <- sum(cvec^2 / w)
    fa <- real[[cl]]; fu <- comp[[comp_of(cl)]]
    data.table(cell = cl, gene = rownames(E)[g], rel_diff_var = abs(vd - (fa$u[g]^2 - fu$u[g]^2)) / vd,
               abs_diff_bu = abs(sum(Au * yv) - fu$b[g]), abs_diff_ba = abs(sum(Aa * yv) - fa$b[g]))
  }))
}))
d5s <- d5[, .(max_rel_diff_var = max(rel_diff_var), max_abs_diff_bu = max(abs_diff_bu), max_abs_diff_ba = max(abs_diff_ba)), by = cell]
d5s[, pass := max_rel_diff_var < 1e-8 & max_abs_diff_bu < 1e-8]
say("Diagnostic 5: variance identity"); sayt(d5s)
fwrite(d5s, file.path(OUT, "diag5_variance_identity.tsv"), sep = "\t")

# per-gene tables
pergene <- rbindlist(lapply(names(real), function(cl) {
  fa <- real[[cl]]; cn <- comp_of(cl); fu <- comp[[cn]]
  bu_ref <- if (cn == "U774") { x <- fu$b; x[gene_idx_calls] <- calls$logFC_stored; x } else fu$b
  D <- fa$b - fu$b; seD <- sqrt(fa$s2 * (fa$u^2 - fu$u^2))
  data.table(cell = cl, gene = rownames(E), logFC_adj = fa$b, SE_adj = fa$u * sqrt(fa$s2), P_adj = fa$p,
             padj_adj = p.adjust(fa$p, "BH"), logFC_comparator = fu$b, logFC_unadj_reported = bu_ref,
             ratio = fa$b / bu_ref, D = D, SE_D = seD, t_D = D / seD, p_D = 2 * pt(-abs(D / seD), fa$df), df = fa$df,
             u_a = fa$u, u_u = fu$u, s2_a = fa$s2)
}))
pergene[, is_call := gene %in% calls$gene]
pergene[is_call == TRUE, q_D_calls := p.adjust(p_D, "BH"), by = cell]
fwrite(pergene[, !c("u_a", "u_u", "s2_a")], file.path(OUT, "per_gene_all_cells.tsv.gz"), sep = "\t")
adjcalls <- pergene[, .(adjusted_calls = sum(padj_adj < 0.05 & abs(logFC_adj) > 0.5),
                        of_stored_calls_kept = sum(is_call & padj_adj < 0.05 & abs(logFC_adj) > 0.5),
                        D_q05_among_calls = sum(is_call & q_D_calls < 0.05, na.rm = TRUE),
                        median_ratio_calls = median(ratio[is_call])), by = cell]
say("Per-cell descriptive counts"); sayt(adjcalls)

# ------------------------------------------------------------------ labels (IUT over the flag cells; BH over 1,854 per family)
pc <- pergene[is_call == TRUE & cell %in% flag_ok]
pc <- merge(pc, calls[, .(gene, s)], by = "gene")
pc[, tkeep := s * (D + 0.5 * logFC_unadj_reported) / SE_D]
pc[, p_keep := pt(tkeep, df, lower.tail = FALSE)]
pc[, p_lose := pt(tkeep, df, lower.tail = TRUE)]
pc[s * logFC_unadj_reported <= 0, `:=`(p_keep = 1, p_lose = 1)]
gl <- pc[, .(p_keep_iut = max(p_keep), p_lose_iut = max(p_lose), n_cells = .N,
             logFC_adj_min = min(logFC_adj), logFC_adj_max = max(logFC_adj),
             ratio_min = min(ratio), ratio_max = max(ratio)), by = gene]
gl <- merge(calls, gl, by = "gene", all.x = TRUE)
stopifnot(nrow(gl) == 1854L)
gl[, q_keep := p.adjust(p_keep_iut, "BH")]
gl[, q_lose := p.adjust(p_lose_iut, "BH")]
labels_computed <- (!S11_stop || (DEBUG && Sys.getenv("H56_DEBUG_FORCE_LABELS", "0") == "1")) && !flag_untestable
if (labels_computed) {
  gl[, label := fifelse(q_keep < 0.05, "at_least_half_retained_all_arms",
                        fifelse(q_lose < 0.05, "more_than_half_removed_all_arms", "depends_on_arm_or_not_resolved"))]
} else gl[, label := NA_character_]
wide <- dcast(pergene[is_call == TRUE, .(gene, cell, logFC_adj, ratio)], gene ~ cell, value.var = c("logFC_adj", "ratio"))
gl <- merge(gl, wide, by = "gene")
fwrite(gl, file.path(OUT, "call_labels_composition_sensitivity.tsv"), sep = "\t")
lab_counts <- if (labels_computed) gl[, .N, by = label] else data.table(label = "not computed", N = NA)
say("Labels computed: ", labels_computed); sayt(lab_counts)
n_dep <- if (labels_computed) gl[label == "depends_on_arm_or_not_resolved", .N] else NA_integer_
state$label_counts <- lab_counts; write_state()

# ------------------------------------------------------------------ controls
mr <- function(gs, cl) median(pergene[cell == cl & gene %in% gs, ratio])
ecm <- calls[stromal_ecm == TRUE, gene]; strom <- calls[sc_class == "stromal", gene]; unif <- calls[sc_class == "uniform", gene]
ctl <- rbindlist(lapply(flag_ok, function(cl) data.table(
  cell = cl,
  PC1_median_ratio = mr(ecm, cl), PC1_all_below_1 = all(pergene[cell == cl & gene %in% ecm, ratio] < 1),
  PC2_median_ratio = mr(strom, cl), NC1_median_ratio = mr(unif, cl),
  NC2_median_ratio = if (length(Hplus)) mr(Hplus, cl) else NA_real_)))
a16 <- merge(pergene[cell == "A16_rebuilt" & gene %in% ecm, .(gene, r_new = ratio)],
             pergene[cell == "A16_earlier" & gene %in% ecm, .(gene, r_old = ratio)], by = "gene")
pc1_more_att <- if (nrow(a16)) a16[, sum(r_new < r_old)] else NA_integer_
PC1_median_pass <- all(ctl$PC1_median_ratio < 0.5)
PC2_median_pass <- all(ctl$PC2_median_ratio < 0.5)
S6_fail <- !(PC1_median_pass && PC2_median_pass)
PC1_all_below_1 <- all(ctl$PC1_all_below_1)
PC1_a16_pass <- isTRUE(pc1_more_att >= 9)
PC2_retained_share <- if (labels_computed) gl[gene %in% strom, mean(label == "at_least_half_retained_all_arms")] else NA_real_
PC2_share_pass <- if (labels_computed) PC2_retained_share <= 0.40 else NA
ctl[, NC1_pass := NC1_median_ratio >= 0.7]
NC1_pass <- all(ctl$NC1_pass)
NC2_status <- if (!isTRUE(sc_res$testable)) "not testable" else {
  med_ok <- all(ctl$NC2_median_ratio >= 0.7)
  ret <- if (labels_computed) gl[gene %in% Hplus, mean(label == "at_least_half_retained_all_arms")] else NA_real_
  state$NC2_retained_share <- ret
  if (!med_ok) "fail" else if (is.na(ret)) "median passes; label part not evaluable (labels not computed)" else if (ret >= 0.5) "pass" else "fail"
}
say("Controls per flag cell"); sayt(ctl)
say("PC1: medians < 0.5 in every flag cell ", PC1_median_pass, "; all 12 ratios < 1 in every cell ", PC1_all_below_1,
    "; A16 rebuilt more attenuated than earlier ", pc1_more_att, " of ", nrow(a16), " (prediction >= 9: ", PC1_a16_pass, ")")
say("PC2: medians < 0.5 in every flag cell ", PC2_median_pass, "; share retained ", signif(PC2_retained_share, 3),
    " (prediction <= 0.40: ", PC2_share_pass, ")")
say("S6 (PC1 or PC2 median ratio >= 0.5 in any flag cell): ", S6_fail)
say("NC1 (uniform calls, median ratio >= 0.7 in every flag cell): ", NC1_pass, "; over-adjusting cells: ",
    paste(ctl[NC1_pass == FALSE, cell], collapse = ", "))
say("NC2: ", NC2_status)
fwrite(ctl, file.path(OUT, "controls_per_cell.tsv"), sep = "\t")
state$controls <- list(PC1_median_pass = PC1_median_pass, PC1_all_below_1 = PC1_all_below_1, PC1_a16_more_attenuated = pc1_more_att,
                       PC1_a16_pass = PC1_a16_pass, PC2_median_pass = PC2_median_pass, PC2_retained_share = PC2_retained_share,
                       PC2_share_pass = PC2_share_pass, S6_fail = S6_fail, NC1_pass = NC1_pass, NC2_status = NC2_status)
write_state()

# class medians
classes <- c("stromal", "hepatocyte", "immune", "cholangiocyte", "uniform", "other", "not_joined")
cm <- rbindlist(lapply(names(real), function(cl) rbindlist(lapply(classes, function(k) {
  gs <- calls[sc_class == k, gene]; data.table(cell = cl, class = k, n = length(gs), median_ratio = mr(gs, cl)) }))))
fwrite(dcast(cm, class + n ~ cell, value.var = "median_ratio"), file.path(OUT, "median_ratio_by_class.tsv"), sep = "\t")
say("Median ratio by class"); sayt(dcast(cm, class + n ~ cell, value.var = "median_ratio"))

# ------------------------------------------------------------------ bootstrap (S7, diagnostic 8)
strata <- paste(dset, grp)
say("Bootstrap strata (774): ", length(unique(strata)), "; minimum stratum size ", min(table(strata)))
set.seed(SEED)
boot_idx <- lapply(seq_len(B_BOOT), function(b) {
  unlist(lapply(sort(unique(strata)), function(s) { i <- which(strata == s); i[sample.int(length(i), length(i), replace = TRUE)] }))
})
boot_cells <- cells[fitted == TRUE, cell]
bres <- mclapply(seq_len(B_BOOT), function(b) {
  idx <- boot_idx[[b]]
  bu <- list()
  for (cn in names(comp_cov)) {
    ii <- idx[units_of(cn)[idx]]; X <- design_u(comp_cov[[cn]], ii)
    bu[[cn]] <- if (fullrank(X)) fitg(ii, X, gene_idx_calls, eb = FALSE)$b else NULL
  }
  ba <- matrix(NA_real_, length(gene_idx_calls), length(boot_cells), dimnames = list(NULL, boot_cells))
  med <- setNames(rep(NA_real_, length(boot_cells)), boot_cells)
  for (cl in boot_cells) {
    ii <- idx[units_of(cells[cell == cl, arm])[idx]]; X <- design_a(cov_of(cl), coords[[cl]], ii)
    cn <- comp_of(cl)
    if (!fullrank(X) || is.null(bu[[cn]])) next
    ba[, cl] <- fitg(ii, X, gene_idx_calls, eb = FALSE)$b
    med[cl] <- median(ba[, cl] / bu[[cn]])
  }
  list(med = med, ba = ba)
}, mc.cores = NCORE, mc.preschedule = TRUE)
berr <- vapply(bres, inherits, TRUE, "try-error")
medmat <- do.call(rbind, lapply(bres[!berr], `[[`, "med"))
failed <- colSums(is.na(medmat)) + sum(berr)
bt <- rbindlist(lapply(boot_cells, function(cl) {
  m <- medmat[, cl]; m <- m[!is.na(m)]; fr <- failed[cl] / B_BOOT
  data.table(cell = cl, median_ratio = adjcalls[cell == cl, median_ratio_calls], draws_ok = length(m), draws_failed = failed[cl],
             failed_share = fr, ci_low = if (fr <= 0.05) unname(quantile(m, 0.025)) else NA_real_,
             ci_high = if (fr <= 0.05) unname(quantile(m, 0.975)) else NA_real_)
}))
say("Diagnostic 8 / S7: bootstrap of the median ratio over the 1,854 calls"); sayt(bt)
fwrite(bt, file.path(OUT, "bootstrap_median_ratio.tsv"), sep = "\t")
# rebuilt minus earlier adjusted log2FC per coding, per call
diffci <- rbindlist(lapply(c("A16", "C7", "D1", "DRmain", "DRsep"), function(a) {
  cr <- paste0(a, "_rebuilt"); ce <- paste0(a, "_earlier")
  if (!all(c(cr, ce) %in% boot_cells)) return(NULL)
  dm <- do.call(cbind, lapply(bres[!berr], function(r) r$ba[, cr] - r$ba[, ce]))
  obs <- pergene[cell == cr & is_call == TRUE][match(calls$gene, gene), logFC_adj] -
    pergene[cell == ce & is_call == TRUE][match(calls$gene, gene), logFC_adj]
  data.table(coding = a, gene = calls$gene, symbol = calls$symbol, diff_rebuilt_minus_earlier = obs,
             ci_low = apply(dm, 1, quantile, 0.025, na.rm = TRUE), ci_high = apply(dm, 1, quantile, 0.975, na.rm = TRUE))
}))
fwrite(diffci, file.path(OUT, "rebuilt_minus_earlier_per_call.tsv.gz"), sep = "\t")
say("Rebuilt minus earlier adjusted log2FC: calls whose 95% interval excludes 0, by coding")
sayt(diffci[, .(n = .N, excl0 = sum(ci_low > 0 | ci_high < 0), median_abs_diff = median(abs(diff_rebuilt_minus_earlier))), by = coding])
rm(bres); invisible(gc())

# ------------------------------------------------------------------ DR-sep (reported apart)
say("DR-sep: median ratio over calls rebuilt ", signif(adjcalls[cell == "DRsep_rebuilt", median_ratio_calls], 4),
    "; earlier ", signif(adjcalls[cell == "DRsep_earlier", median_ratio_calls], 4))

# ------------------------------------------------------------------ library complexity arms
nested <- function(rows, Xs, Xb) {  # Xb nested in Xs; per-gene D and exact SE over calls
  fs <- fitg(rows, Xs); fb <- fitg(rows, Xb)
  D <- fs$b - fb$b; se <- sqrt(fs$s2 * (fs$u^2 - fb$u^2))
  list(ratio = fs$b[gene_idx_calls] / fb$b[gene_idx_calls], D = D[gene_idx_calls], se = se[gene_idx_calls])
}
cx <- list()
addX <- function(X, rows) { X <- cbind(X, top100_share = t100[rows]); X[, c(setdiff(colnames(X), GRP), GRP)] }
r774 <- rep(TRUE, 774); r537 <- units$in_DRmain
XU774 <- design_u(comp_cov$U774, r774); XX774 <- addX(XU774, r774)
XU537 <- design_u(comp_cov$U537, r537); XX537 <- addX(XU537, r537)
cx[["X774_vs_U774"]] <- nested(which(r774), XX774, XU774)
cx[["X537_vs_U537"]] <- nested(which(r537), XX537, XU537)
for (cl in intersect(c("A16_rebuilt", "A16_earlier"), names(real)))
  cx[[paste0(cl, "+X_vs_X774")]] <- nested(which(r774), addX(design_a(cov_of(cl), coords[[cl]], r774), r774), XX774)
for (cl in intersect(c("DRmain_rebuilt", "DRmain_earlier"), names(real)))
  cx[[paste0(cl, "+X_vs_X537")]] <- nested(which(r537), addX(design_a(cov_of(cl), coords[[cl]], r537), r537), XX537)
cxt <- rbindlist(lapply(names(cx), function(k) data.table(comparison = k, median_ratio = median(cx[[k]]$ratio),
  median_abs_D = median(abs(cx[[k]]$D)), n_D_q05 = sum(p.adjust(2 * pnorm(-abs(cx[[k]]$D / cx[[k]]$se)), "BH") < 0.05))))
say("Library complexity arms (median ratio over the 1,854 calls)"); sayt(cxt)
fwrite(cxt, file.path(OUT, "library_complexity_arms.tsv"), sep = "\t")
fwrite(rbindlist(lapply(names(cx), function(k) data.table(comparison = k, gene = calls$gene, ratio = cx[[k]]$ratio, D = cx[[k]]$D, SE_D = cx[[k]]$se))),
       file.path(OUT, "library_complexity_per_call.tsv.gz"), sep = "\t")

# ------------------------------------------------------------------ common support
csup <- rbindlist(lapply(intersect(c("A16_rebuilt", "A16_earlier", "DRmain_rebuilt", "DRmain_earlier"), names(real)), function(cl) {
  arm <- cells[cell == cl, arm]; rlog <- units_of(arm) & prop[[paste0("keep_", cl)]] %in% 1
  cn <- comp_of(cl); Xu <- design_u(comp_cov[[cn]], rlog); Xa <- design_a(cov_of(cl), coords[[cl]], rlog)
  if (!fullrank(Xu) || !fullrank(Xa)) return(data.table(cell = cl, n_kept = sum(rlog), median_ratio = NA_real_, spearman_vs_untrimmed = NA_real_))
  rows <- which(rlog)
  r <- fitg(rows, Xa, gene_idx_calls, eb = FALSE)$b / fitg(rows, Xu, gene_idx_calls, eb = FALSE)$b
  r0 <- pergene[cell == cl][match(calls$gene, gene), ratio]
  data.table(cell = cl, n_kept = sum(rlog), n_dropped = sum(units_of(arm)) - sum(rlog), median_ratio = median(r),
             spearman_vs_untrimmed = cor(r, r0, method = "spearman"))
}), fill = TRUE)
say("Common support (trimmed to overlapping propensity)"); sayt(csup)
fwrite(csup, file.path(OUT, "common_support.tsv"), sep = "\t")

# ------------------------------------------------------------------ weight re-estimation check
if (!SKIP_VQW && "A16_rebuilt" %in% names(real)) {
  t0 <- Sys.time()
  XA <- design_a(cov_of("A16_rebuilt"), coords[["A16_rebuilt"]], r774)
  vA <- voomWithQualityWeights(d, XA)
  fA <- eBayes(lmFit(vA, XA)); bA <- fA$coefficients[gene_idx_calls, GRP]
  bF <- real[["A16_rebuilt"]]$b[gene_idx_calls]
  wchk <- data.table(pearson = cor(bA, bF), median_abs_diff = median(abs(bA - bF)), max_abs_diff = max(abs(bA - bF)),
                     minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  say("Weight re-estimation check (A16 rebuilt, 1,854 calls)"); sayt(wchk)
  fwrite(wchk, file.path(OUT, "weight_reestimation_check.tsv"), sep = "\t")
  rm(vA, fA); invisible(gc())
}

# ------------------------------------------------------------------ F4 vs F0 (THRB, PPARA)
f4 <- local({
  SM <- file.path(ROOT, "figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions")
  sm <- fread(file.path(SM, "stage_extension_sample_manifest.tsv"))[contrast == "F4_vs_F0"]
  setorder(sm, sample_id)
  st <- fread(cmd = paste("grep -P 'F4_vs_F0|^gene_id_versioned'", shQuote(file.path(SM, "stage_extension_all_gene_results.tsv"))))
  st <- st[contrast == "F4_vs_F0" & gene_name %in% c("THRB", "PPARA")]
  dge <- d[, match(sm$sample_id, colnames(d))]
  info <- data.frame(dataset = droplevels(factor(sm$dataset)), inferred_sex = droplevels(factor(sm$inferred_sex)),
                     group = factor(sm$group, levels = c("F0", "F4")), row.names = sm$sample_id)
  des <- model.matrix(~ dataset + inferred_sex + group, info)
  v4 <- voomWithQualityWeights(dge, des, plot = FALSE)
  f <- eBayes(lmFit(v4, des))
  g <- st$gene_id_versioned; j <- match("groupF4", colnames(des))
  gi <- match(g, rownames(f$coefficients))
  b <- f$coefficients[gi, j]; se <- f$stdev.unscaled[gi, j] * sqrt(f$s2.post[gi]); cr <- qt(0.975, f$df.total[gi])
  s9 <- max(abs(b - st$logFC))
  res <- data.table(gene = st$gene_name, model = "unadjusted", logFC = b, ci_low = b - cr * se, ci_high = b + cr * se, stored = st$logFC)
  if (s9 < 1e-4) {
    for (cl in c("A16_rebuilt", "A16_earlier", "C7_rebuilt", "C7_earlier", "D1_rebuilt", "D1_earlier")) {
      Xa <- cbind(des[, colnames(des) != "groupF4"], coords[[cl]][sm$sample_id, ], groupF4 = des[, "groupF4"])
      if (!fullrank(Xa)) next
      fa <- eBayes(lmFit(v4$E, Xa, weights = v4$weights)); ja <- ncol(Xa)
      ba <- fa$coefficients[match(g, rownames(v4$E)), ja]; sea <- fa$stdev.unscaled[match(g, rownames(v4$E)), ja] * sqrt(fa$s2.post[match(g, rownames(v4$E))])
      cra <- qt(0.975, fa$df.total[match(g, rownames(v4$E))])
      res <- rbind(res, data.table(gene = st$gene_name, model = cl, logFC = ba, ci_low = ba - cra * sea, ci_high = ba + cra * sea, stored = NA_real_))
    }
  }
  list(s9 = s9, res = res, n = nrow(sm), n_f4 = sum(sm$group == "F4"))
})
say("F4 vs F0: units ", f4$n, " (F4 ", f4$n_f4, "); S9 max |refit - stored| THRB/PPARA ", signif(f4$s9, 3), "; extension run ", f4$s9 < 1e-4)
sayt(f4$res)
fwrite(f4$res, file.path(OUT, "f4_vs_f0_thrb_ppara.tsv"), sep = "\t")
state$S9 <- f4$s9

# ------------------------------------------------------------------ verdict (PRESPEC "Reading")
pc_fail <- S6_fail
nc2_fail <- identical(NC2_status, "fail")
if (flag_untestable) {
  verdict <- "untestable"; reason <- "fewer than six flag cells estimable (S5)"
} else if (pc_fail) {
  verdict <- "positive control failed"; reason <- "S6: PC1 or PC2 median ratio >= 0.5 in a flag cell"
} else {
  why <- c(if (S11_stop) "S11 calibration null stopped the labels",
           if (!is.na(n_dep) && n_dep > 927) sprintf("%d of 1,854 calls depends_on_arm_or_not_resolved (> 927)", n_dep),
           if (!NC1_pass) paste("NC1 failed in", paste(ctl[NC1_pass == FALSE, cell], collapse = ", ")),
           if (nc2_fail) "NC2 failed")
  if (length(why)) { verdict <- "null"; reason <- paste(why, collapse = "; ") } else {
    verdict <- "labels published (composition-conditional)"; reason <- "PC1, PC2, NC1 pass; NC2 pass or not testable; S11 pass; <= 927 unresolved"
  }
}
state$verdict <- verdict; state$reason <- reason; write_state()
say("VERDICT: ", verdict, " (", reason, ")")
saveRDS(list(ctl = ctl, bt = bt, gl = gl, d9 = d9, adjcalls = adjcalls, cm = cm, cxt = cxt, csup = csup, f4 = f4,
             sc = sc_res[setdiff(names(sc_res), c("tab", "Hplus"))], Hplus = Hplus, lab_counts = lab_counts, n_dep = n_dep,
             d2 = fread(file.path(OUT, "diag2_design.tsv")), d3 = fread(file.path(OUT, "diag3_overlap.tsv")), d4 = d4, d5s = d5s,
             s2 = s2, flag_ok = flag_ok, pc1_more_att = pc1_more_att, verdict = verdict, reason = reason,
             labels_computed = labels_computed, NC2_status = NC2_status, PC2_retained_share = PC2_retained_share,
             wchk = if (exists("wchk")) wchk else NULL),
        file.path(OUT, "summary.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
say("FIT_DONE ", format(Sys.time(), tz = "UTC"))
