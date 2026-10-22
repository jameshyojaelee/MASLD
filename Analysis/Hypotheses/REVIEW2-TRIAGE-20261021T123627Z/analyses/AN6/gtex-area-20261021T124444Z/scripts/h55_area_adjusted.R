# AN6 / PLAN A9: H55 GTEx fat concordance with slide tissue area added as a covariate.
# Post hoc sensitivity requested by review triage (REVIEW2-TRIAGE-20261021T123627Z). Not prespecified.
# Data loading, filter, TMM, family, voom weighting, Freedman-Lane machinery and the controls are copied
# from the producer: Analysis/Hypotheses/H55-20261015T112558Z/run-20261015T114811Z/scripts_used/h55_stats.R
# (sections 1, 3, 4, 5, 7). Permutation and bootstrap indices are the producer's saved matrices
# (perm_index_primary.rds, 10,000 draws; boot_index_Q2.rds, 2,000 stratified donor resamples).
# Covariate versions follow the producer's DEVIATIONS.md item 3: voomWithQualityWeights is rerun with the
# covariate-only design that includes the added covariate, and the primary permutation matrix is reused.
#
# Usage: Rscript h55_area_adjusted.R <out_dir>
suppressMessages({library(edgeR); library(limma); library(data.table)})
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
SEED <- 20260923L; set.seed(SEED)
OUTD <- commandArgs(trailingOnly = TRUE)[1]
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
H55 <- paste0(BASE, "Analysis/Hypotheses/H55-20261015T112558Z/")
H55RUN <- paste0(H55, "run-20261015T114811Z/")
RR <- paste0(BASE, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/")
PL <- paste0(RR, "scratch/plan-H55/")
LOG <- file(file.path(OUTD, "h55_area_adjusted.log"), open = "wt")
say <- function(...) { s <- paste0(...); cat(s, "\n"); cat(s, "\n", file = LOG); flush(LOG) }
say("start ", format(Sys.time(), tz = "UTC"))

## 0. Input hashes (producer's lists); refuse on mismatch
for (f in c("prespec.sha256", "inputs.sha256")) {
  rc <- system(paste0("sha256sum --quiet -c ", H55, "scripts/", f), intern = FALSE)
  if (rc != 0) { say("HASH MISMATCH in ", f, "; refusing"); quit(status = 2) }
}
say("hash check OK")

## 1. Donor set (producer section 1)
sa <- fread(paste0(BASE, "data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
            select = c("SAMPID", "SMTSD", "SMAFRZE", "SMRIN", "SMTSISCH", "SMCENTER", "SMEXNCRT"))
subj <- fread(paste0(BASE, "data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"))
lv <- sa[SMTSD == "Liver" & SMAFRZE == "RNASEQ"]; lv[, donor := sub("^(GTEX-[^-]+)-.*$", "\\1", SAMPID)]
sl <- fread(paste0(BASE, "Analysis/MASLD_Model_Benchmark/executions/he-molecular-c1-20260924T004124Z/liver_slides_matched.tsv"))
both <- intersect(lv$donor, sl$PatientID)
m <- merge(lv[donor %in% both], subj, by.x = "donor", by.y = "SUBJID")
m <- m[!is.na(DTHHRDY) & !is.na(SMCENTER) & SMCENTER != ""]
dm <- fread(paste0(PL, "design_meta_counts_only.tsv")); dm <- dm[!is.na(SMCENTER) & SMCENTER != ""]
setkey(m, donor); setkey(dm, donor)
stopifnot(nrow(m) == 222, identical(m$donor, dm$donor), identical(m$SAMPID, dm$SAMPID))
ng <- fread(paste0(RR, "scratch/gate-L02-01/note_grades.tsv"))
meta <- merge(dm[, .(donor, SAMPID, H, SMCENTER, AGE, SEX, SMRIN, SMTSISCH, SMEXNCRT, grade, pct, fib_or_nash)],
              ng[, .(donor, ord)], by = "donor", all.x = TRUE, sort = TRUE)
meta[, ageg := factor(fcase(AGE %in% c("20-29", "30-39"), "20-39", AGE == "40-49", "40-49", AGE == "50-59", "50-59",
                            AGE %in% c("60-69", "70-79"), "60-79"), levels = c("20-39", "40-49", "50-59", "60-79"))]
meta[, hg := fcase(H == 0, "H0", H %in% 1:2, "H12", H %in% 3:4, "H34")]
meta[, stratum := paste(SMCENTER, hg)]

## 2. Exposure and tissue area from the producer's final morphometry
mf <- fread(paste0(H55RUN, "morph/morph_fat.tsv"))
meta <- merge(meta, mf[status == "ok", .(donor, fat_frac, tissue_px)], by = "donor", all.x = TRUE, sort = TRUE)
fi <- which(is.finite(meta$fat_frac)); nf <- length(fi); stopifnot(nf == 222)
rk <- function(x) (rank(x, ties.method = "average") - 0.5) / length(x)
meta[, fat_rank := NA_real_]; meta[fi, fat_rank := rk(fat_frac)]
meta[, tissue_rank := NA_real_]; meta[fi, tissue_rank := rk(tissue_px)]
meta[, log10_tissue := log10(tissue_px)]
dex <- fread(paste0(H55RUN, "outputs/donor_exposure_table.tsv"))
stopifnot(identical(dex$donor, meta$donor), isTRUE(all.equal(dex$fat_rank, meta$fat_rank)), isTRUE(all.equal(dex$tissue_rank, meta$tissue_rank)))
rho_fat_area_donor <- cor(meta$fat_rank[fi], meta$tissue_rank[fi], method = "spearman")
say(sprintf("donors %d; donor-level Spearman(fat rank, tissue-area rank) %.4f", nf, rho_fat_area_donor))

## 3. Counts, filter, TMM, family (producer section 3)
cnt <- fread(paste0(PL, "inputs/gene_reads_2017-06-05_v8_liver.gct.gz"), skip = 2)
cnt <- cnt[!grepl("_PAR_Y$", Name)]
ens <- sub("\\..*$", "", cnt$Name); stopifnot(!anyDuplicated(ens), nrow(cnt) == 56156)
M <- as.matrix(cnt[, meta$SAMPID, with = FALSE]); rownames(M) <- ens; colnames(M) <- meta$donor
rm(cnt); invisible(gc())
y0 <- DGEList(M); keep <- filterByExpr(y0)
y <- calcNormFactors(y0[keep, , keep.lib.sizes = FALSE], method = "TMM"); rm(y0, M); invisible(gc())
say("filterByExpr pass: ", sum(keep))
fam <- fread(paste0(PL, "primary_family_genes.tsv")); stopifnot(nrow(fam) == 1002)
fg <- fam$ens; lfc <- fam$bulk_logFC

## 4. FL machinery (verbatim from producer)
fl_prep <- function(E, W, Zn, x, genes) {
  n <- ncol(E); G <- length(genes)
  e <- matrix(0, n, G); xr <- matrix(0, n, G); Ql <- vector("list", G)
  for (j in seq_len(G)) {
    sw <- sqrt(W[genes[j], ]); q <- qr(Zn * sw); Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
    yw <- E[genes[j], ] * sw; xw <- x * sw
    e[, j] <- yw - Q %*% crossprod(Q, yw); xr[, j] <- xw - Q %*% crossprod(Q, xw); Ql[[j]] <- Q
  }
  sxx <- colSums(xr^2)
  list(e = e, xr = xr, Q = Ql, sxx = sxx, slope = colSums(xr * e) / sxx, genes = genes, n = n, p = ncol(Zn))
}
null_slopes <- function(P, IDX) {
  out <- matrix(NA_real_, nrow(IDX), ncol(P$e))
  for (b in seq_len(nrow(IDX))) out[b, ] <- colSums(P$xr * P$e[IDX[b, ], , drop = FALSE]) / P$sxx
  out
}
rank_rows <- function(Mx) t(apply(Mx, 1, rank))
null_rho <- function(NS, ref) { rr <- rank(ref); R <- rank_rows(NS); as.vector(cor(t(R), rr)) }
flp <- function(obs, null) {
  Bn <- length(null); pu <- (1 + sum(null >= obs)) / (Bn + 1); pl <- (1 + sum(null <= obs)) / (Bn + 1)
  c(p_upper = pu, p_lower = pl, p_two = min(1, 2 * min(pu, pl)), null_q025 = unname(quantile(null, 0.025)), null_q975 = unname(quantile(null, 0.975)))
}
sp <- function(a, b) cor(a, b, method = "spearman")
perm_matrix <- function(strata, Bn) {
  set.seed(SEED); n <- length(strata); us <- sort(unique(strata))
  IDX <- matrix(0L, Bn, n)
  for (b in seq_len(Bn)) { p <- seq_len(n); for (s in us) { i <- which(strata == s); p[i] <- i[sample.int(length(i))] }; IDX[b, ] <- p }
  IDX
}
IDX <- readRDS(paste0(H55RUN, "outputs/perm_index_primary.rds"))
BIDX <- readRDS(paste0(H55RUN, "outputs/boot_index_Q2.rds"))
mf_ <- meta[fi]
stopifnot(identical(IDX, perm_matrix(mf_$stratum, nrow(IDX))))
say(sprintf("permutation matrix %d x %d reproduced from seed; bootstrap matrix %d x %d", nrow(IDX), ncol(IDX), nrow(BIDX), ncol(BIDX)))

## 5. Models
covf <- ~ ageg + factor(SEX) + factor(H) + SMTSISCH + SMRIN + SMCENTER         # producer primary nuisance design
cova <- update(covf, ~ . + tissue_rank)                                     # + slide tissue-area rank
covl <- update(covf, ~ . + log10_tissue)                                    # + log10 tissue area (descriptive)
voom_w <- function(form) {
  v <- voomWithQualityWeights(y, design = model.matrix(form, data = meta), plot = FALSE)
  E <- v$E[, fi]; W <- v$weights[, fi]; rownames(W) <- rownames(E); list(E = E, W = W)
}
V0 <- voom_w(covf); VA <- voom_w(cova); VL <- voom_w(covl)
Z0 <- model.matrix(covf, data = mf_); ZA <- model.matrix(cova, data = mf_); ZL <- model.matrix(covl, data = mf_)
say(sprintf("1 - R2 of fat rank on primary covariates %.3f; with tissue-area rank %.3f",
            1 - summary(lm(mf_$fat_rank ~ Z0[, -1]))$r.squared, 1 - summary(lm(mf_$fat_rank ~ ZA[, -1]))$r.squared))

pcg <- fread(paste0(PL, "pc_composite_genes.tsv")); stopifnot(nrow(pcg) == 197)
plin <- "ENSG00000147872"
run_version <- function(name, V, Zn, do_boot = FALSE) {
  P <- fl_prep(V$E, V$W, Zn, mf_$fat_rank, fg)
  r <- sp(lfc, P$slope); nul <- null_rho(null_slopes(P, IDX), lfc); p <- flp(r, nul)
  X <- cbind(Zn[, 1, drop = FALSE], fat_rank = mf_$fat_rank, Zn[, -1, drop = FALSE])
  eb <- eBayes(lmFit(V$E, X, weights = V$W))
  stopifnot(max(abs(P$slope - eb$coefficients[fg, "fat_rank"])) < 1e-6)
  Ppc <- fl_prep(V$E, V$W, Zn, mf_$fat_rank, pcg$ens)
  rpc <- sp(pcg$h8_mu, Ppc$slope); ppc <- flp(rpc, null_rho(null_slopes(Ppc, IDX), pcg$h8_mu))
  bt <- rep(NA_real_, nrow(BIDX)); n_drop <- 0L
  if (do_boot) {
    SW <- sqrt(V$W[fg, , drop = FALSE]); Y <- V$E[fg, , drop = FALSE]
    for (b in seq_len(nrow(BIDX))) {
      ib <- BIDX[b, ]; Xb <- X[ib, , drop = FALSE]
      kc <- c(TRUE, apply(Xb[, -1, drop = FALSE], 2, function(z) length(unique(z)) > 1))
      if (!all(kc)) n_drop <- n_drop + 1L
      Xb <- Xb[, kc, drop = FALSE]; jf <- which(colnames(Xb) == "fat_rank")
      bf <- numeric(length(fg))
      for (j in seq_along(fg)) { w <- SW[j, ib]; bf[j] <- qr.coef(qr(Xb * w), Y[j, ib] * w)[jf] }
      bt[b] <- sp(lfc, bf)
    }
  }
  out <- data.table(version = name, n_donors = nf, n_genes = length(fg), design_cols = ncol(X), resid_df = nf - ncol(X),
                    rho = r, fl_p_two = unname(p["p_two"]), fl_p_upper = unname(p["p_upper"]),
                    null_q025 = unname(p["null_q025"]), null_q975 = unname(p["null_q975"]),
                    boot_q025 = if (do_boot) unname(quantile(bt, 0.025)) else NA_real_,
                    boot_q975 = if (do_boot) unname(quantile(bt, 0.975)) else NA_real_,
                    boot_n_col_dropped = if (do_boot) n_drop else NA_integer_,
                    plin2_slope = eb$coefficients[plin, "fat_rank"], plin2_mod_p = eb$p.value[plin, "fat_rank"],
                    pc3_rho = rpc, pc3_fl_p_two = unname(ppc["p_two"]))
  say(sprintf("%-34s rho %.4f FL p %.4g null [%.3f, %.3f] boot [%s, %s] PLIN2 %.4f p %.3g PC3 %.4f p %.4g",
              name, r, p["p_two"], p["null_q025"], p["null_q975"], format(out$boot_q025, digits = 4), format(out$boot_q975, digits = 4),
              out$plin2_slope, out$plin2_mod_p, rpc, ppc["p_two"]))
  list(row = out, slope = P$slope, boot = bt)
}
R0 <- run_version("primary_reproduced", V0, Z0, do_boot = TRUE)
stopifnot(abs(R0$row$rho - 0.63800644940459) < 1e-10)             # producer h55_results.rds primary$rho
RA <- run_version("tissue_area_rank_covariate", VA, ZA, do_boot = TRUE)
RW <- run_version("tissue_area_rank_covariate_primary_weights", V0, ZA)
RL <- run_version("log10_tissue_area_covariate", VL, ZL)

# paired bootstrap difference (area-adjusted minus primary), same resamples
dlt <- RA$boot - R0$boot
dsum <- data.table(quantity = "rho_area_adjusted_minus_primary", observed = RA$row$rho - R0$row$rho,
                   boot_q025 = unname(quantile(dlt, 0.025)), boot_q975 = unname(quantile(dlt, 0.975)), n_boot = length(dlt))
say(sprintf("paired difference rho_area - rho_primary %.4f, bootstrap 95%% [%.4f, %.4f]", dsum$observed, dsum$boot_q025, dsum$boot_q975))

## 6. NC1 reproduction (producer section 7), as a check of the tissue-area exposure
Pnc <- fl_prep(V0$E, V0$W, cbind(Z0, fat_rank = mf_$fat_rank), mf_$tissue_rank, fg)
rho_nc1 <- sp(lfc, Pnc$slope)
say(sprintf("NC1 reproduced rho %.4f (producer 0.3537)", rho_nc1))
stopifnot(abs(rho_nc1 - 0.353721770316298) < 1e-10)

summ <- rbindlist(list(R0$row, RA$row, RW$row, RL$row))
summ[, donor_spearman_fat_vs_area := rho_fat_area_donor]
fwrite(summ, file.path(OUTD, "area_adjusted_summary.tsv"), sep = "\t")
fwrite(dsum, file.path(OUTD, "area_minus_primary_bootstrap.tsv"), sep = "\t")
fwrite(data.table(gene = fg, symbol = fam$symbol, bulk_logFC = lfc, slope_primary = R0$slope, slope_area_rank = RA$slope,
                  slope_area_rank_primary_weights = RW$slope, slope_log10_area = RL$slope, nc1_slope = Pnc$slope),
       file.path(OUTD, "family_slopes_area_adjusted.tsv"), sep = "\t")
fwrite(data.table(b = seq_along(R0$boot), rho_primary = R0$boot, rho_area_rank = RA$boot), file.path(OUTD, "bootstrap_rho.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUTD, "sessionInfo.txt"))
say("done ", format(Sys.time(), tz = "UTC"))
