#!/usr/bin/env Rscript
# A4 step 2: donor x lineage peak-gene links (Q1, Q2). See PRESPEC.md.
#   Units: donor x lineage pseudobulk (H24 build), >= 20 nuclei, lineage tested if >= 15 donors.
#   RNA filterByExpr + TMM logCPM; ATAC count >= 5 in >= 50% units + TMM logCPM.
#   Residuals on intercept + well + log(nuclei) + K own-modality PCs; Pearson r across donors.
#   Background z: all tested peaks of the pair's peak stratum (GC quintile x mean logCPM
#   decile x width tertile, merged to >= 250 peaks) on chromosomes other than the gene's,
#   computed exactly from per-stratum sums and Gram matrices.
#   K chosen on 2,000 own-promoter and 2,000 trans control pairs.
#   Q1 family: all cis pairs (+/- 500 kb); Q2 family: cis pairs of variant peaks
#   (credible-set noncoding PIP sum >= 0.1; sensitivity any member). BH within family.
#   Stability: N_BOOT seeded donor bootstrap replicates (K, features, strata fixed).
#   Split halves: 20 seeded splits, used only if residual df >= 8 in both halves.
# Seeds: 20261024 + offsets (controls 10k+1/10k+2, bootstrap 1000+100k+b, splits 5000+100k+s;
#   k = 1 hepatocyte, 2 fibroblast, 3 macrophage).
# Usage: 02_links.R <A4_DIR> <timing|full> <lineages,comma> <N_BOOT> <N_CORES>
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(jsonlite); library(Matrix)
  library(GenomicRanges); library(parallel)
})
args <- commandArgs(TRUE)
A4 <- args[1]; MODE <- args[2]; LINS <- strsplit(args[3], ",")[[1]]
N_BOOT <- as.integer(args[4]); NCORE <- as.integer(args[5])
SEED <- 20261024L; MIN_NUC <- 20L; MIN_DONORS <- 15L; N_CTRL <- 2000L; WIN <- 5e5; PROM <- 1000L
K_GRID <- c(0L, 2L, 4L, 6L, 8L, 10L); MIN_STRATUM <- 250L; N_SPLIT <- 20L; MIN_DF_SPLIT <- 8L
MIN_DF_BOOT <- 3L; ATAC_MIN_COUNT <- 5; ATAC_MIN_FRAC <- 0.5; Q_LINK <- 0.05; STABLE <- 0.8
PIP_MIN <- 0.1; LIN_INDEX <- c(hepatocyte = 1L, fibroblast = 2L, macrophage = 3L)
PRIMARY <- paste0("chr", c(1:22, "X"))
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
H24 <- file.path(ROOT, "Analysis/Hypotheses/H24-20261008T211750Z/gse296875")
W <- file.path(A4, "work"); RES <- file.path(A4, "results")
TAG <- if (MODE %in% c("timing", "posthoc")) paste0(MODE, "_") else ""
dir.create(RES, showWarnings = FALSE)
tic <- function(...) message(format(Sys.time(), "%H:%M:%S"), " ", ...)

read_f32 <- function(path, dims) matrix(readBin(path, "numeric", n = prod(dims), size = 4, endian = "little"), nrow = dims[1])
units <- fread(file.path(H24, "units.tsv"))
dims <- fromJSON(file.path(H24, "dims.json"))
rna <- read_f32(file.path(H24, "rna_gene_by_unit.f32"), dims$rna)
atac <- read_f32(file.path(H24, "atac_peak_by_unit.f32"), dims$atac)
genes <- fread(file.path(W, "genes_tss.tsv"))
pk <- fread(file.path(W, "peaks_annot.tsv"))
vpc <- fread(file.path(W, "variant_peak_cs.tsv"))
stopifnot(nrow(genes) == nrow(rna), all(genes$gene_index == seq_len(nrow(rna)) - 1L), nrow(pk) == nrow(atac))
tic("inputs loaded")

norm_logcpm <- function(M) cpm(calcNormFactors(DGEList(M)), log = TRUE, prior.count = 2)
design0 <- function(u) model.matrix(~ factor(well_id) + log(nuclei), data = u)
resid_K <- function(L, X0, K) {
  R0 <- t(qr.resid(qr(X0), t(L)))
  if (K > 0L) {
    S <- R0 / sqrt(rowSums(R0^2))
    U <- eigen(crossprod(S), symmetric = TRUE)$vectors[, seq_len(K), drop = FALSE]
    R0 <- R0 - (R0 %*% U) %*% t(U)
  }
  nr <- sqrt(rowSums(R0^2)); bad <- nr < 1e-10
  nr[bad] <- 1; R0 <- R0 / nr; R0[bad, ] <- 0
  R0
}
# Background mean and sd of r(g, q) over peaks q of stratum s on chromosomes other than g's.
bg_stats <- function(R, A, gchr, pchr, strat, S) {
  Ind <- sparseMatrix(i = seq_along(strat), j = strat, x = 1, dims = c(length(strat), S))
  MU <- R %*% t(as.matrix(crossprod(Ind, A)))
  NB <- matrix(tabulate(strat, S), nrow(R), S, byrow = TRUE)
  SS <- matrix(0, nrow(R), S)
  sidx <- split(seq_along(strat), factor(strat, levels = seq_len(S)))
  for (s in seq_len(S)) SS[, s] <- rowSums((R %*% crossprod(A[sidx[[s]], , drop = FALSE])) * R)
  for (cc in unique(gchr)) {
    gi <- which(gchr == cc); pj <- which(pchr == cc)
    if (!length(pj) || !length(gi)) next
    Rg <- R[gi, , drop = FALSE]
    MU[gi, ] <- MU[gi, ] - Rg %*% t(as.matrix(crossprod(Ind[pj, , drop = FALSE], A[pj, , drop = FALSE])))
    NB[gi, ] <- NB[gi, ] - matrix(tabulate(strat[pj], S), length(gi), S, byrow = TRUE)
    sc <- split(pj, strat[pj])
    for (s in names(sc)) {
      j <- as.integer(s)
      SS[gi, j] <- SS[gi, j] - rowSums((Rg %*% crossprod(A[sc[[s]], , drop = FALSE])) * Rg)
    }
  }
  mu <- MU / NB
  list(mu = mu, sd = sqrt(pmax((SS - NB * mu^2) / (NB - 1), 0)), nb_min = min(NB))
}
pair_stats <- function(R, A, g, p, bg, strat) {
  r <- numeric(length(g))
  for (ix in split(seq_along(g), ceiling(seq_along(g) / 5e5)))
    r[ix] <- rowSums(R[g[ix], , drop = FALSE] * A[p[ix], , drop = FALSE])
  m <- bg$mu[cbind(g, strat[p])]; s <- bg$sd[cbind(g, strat[p])]
  z <- (r - m) / s
  list(r = r, z = z, p = 2 * pnorm(-abs(z)))
}
mk_strata <- function(gc, acc, wid) {
  q <- function(x, k) cut(x, unique(quantile(x, 0:k / k)), include.lowest = TRUE, labels = FALSE)
  d <- data.table(i = seq_along(gc), g = q(gc, 5), a = q(acc, 10), w = q(wid, 3))
  d[, grp := {
    cnt <- tabulate(a, 10); lab <- integer(10); cur <- 1L; acc_n <- 0L
    for (j in 1:10) { lab[j] <- cur; acc_n <- acc_n + cnt[j]; if (acc_n >= MIN_STRATUM && j < 10) { cur <- cur + 1L; acc_n <- 0L } }
    if (acc_n < MIN_STRATUM && cur > 1L) lab[lab == cur] <- cur - 1L
    lab[a]
  }, by = .(g, w)]
  d[, s := .GRP, by = .(g, w, grp)]
  list(strat = d$s, S = max(d$s), gc_q = d$g, acc_d = d$a, wid_t = d$w)
}
bh_call <- function(st, idx = NULL) {
  if (is.null(idx)) idx <- seq_along(st$p)
  q <- p.adjust(st$p[idx], "BH")
  list(q = q, call = q < Q_LINK & st$z[idx] > 0)
}

summ <- list(); calib_all <- list(); hist_all <- list(); split_all <- list(); bootrep_all <- list()
for (lin in LINS) {
  k <- LIN_INDEX[[lin]]; t0 <- Sys.time()
  u <- units[lineage_id == lin & nuclei >= MIN_NUC]; n <- nrow(u)
  if (n < MIN_DONORS) { summ[[lin]] <- data.table(lineage = lin, donors = n, status = "untestable_donors"); next }
  cols <- u$unit_index + 1L
  kg <- which(filterByExpr(DGEList(rna[, cols]), design = matrix(1, n, 1)) & genes$has_tss)
  kp <- which(pk$primary & rowSums(atac[, cols] >= ATAC_MIN_COUNT) >= ATAC_MIN_FRAC * n)
  gt <- genes[kg]; pt <- pk[kp]
  gchr <- match(gt$chr, PRIMARY); pchr <- match(pt$chr, PRIMARY)
  X0 <- design0(u); rk0 <- qr(X0)$rank
  lr_full <- norm_logcpm(rna[kg, cols]); la_full <- norm_logcpm(atac[kp, cols])
  stt <- mk_strata(pt$gc, rowMeans(la_full), pt$width)
  strat <- stt$strat; S <- stt$S
  tic(lin, ": donors ", n, ", genes ", length(kg), ", peaks ", length(kp), ", strata ", S)

  # cis pairs
  ov <- findOverlaps(GRanges(pt$chr, IRanges(pt$centre, width = 1L)),
                     GRanges(gt$chr, IRanges(gt$tss, width = 1L)), maxgap = WIN)
  pr <- data.table(p = queryHits(ov), g = subjectHits(ov))
  pr[, dist := pt$centre[p] - gt$tss[g]]
  pr <- pr[abs(dist) <= WIN]
  vp1 <- unique(vpc[pip_sum >= PIP_MIN, peak_index]); vpa <- unique(vpc$peak_index)
  pr[, `:=`(in_q2 = pt$peak_index[p] %in% vp1, in_q2any = pt$peak_index[p] %in% vpa)]
  i2 <- which(pr$in_q2); i2a <- which(pr$in_q2any)
  tic(lin, ": cis pairs ", nrow(pr), ", Q2 pairs ", length(i2), ", Q2-any pairs ", length(i2a))

  # controls
  set.seed(SEED + 10000L * k + 1L)
  po <- findOverlaps(GRanges(gt$chr, IRanges(gt$tss - PROM, gt$tss + PROM)), GRanges(pt$chr, IRanges(pt$start, pt$end)))
  prom <- data.table(g = queryHits(po), p = subjectHits(po))
  prom[, d := abs(pt$centre[p] - gt$tss[g])]
  prom <- prom[order(g, d)][, .SD[1], by = g]
  prom <- prom[sort(sample.int(nrow(prom), min(N_CTRL, nrow(prom))))]
  set.seed(SEED + 10000L * k + 2L)
  tp <- sample.int(nrow(pt), N_CTRL)
  tg <- vapply(tp, function(p) { cand <- which(gchr != pchr[p]); cand[sample.int(length(cand), 1L)] }, 1L)
  trans <- data.table(g = tg, p = tp)
  cg <- sort(unique(c(prom$g, trans$g)))
  calib <- rbindlist(lapply(K_GRID, function(K) {
    R <- resid_K(lr_full, X0, K)[cg, , drop = FALSE]; A <- resid_K(la_full, X0, K)
    bg <- bg_stats(R, A, gchr[cg], pchr, strat, S)
    gm <- match(prom$g, cg); tm <- match(trans$g, cg)
    sp <- pair_stats(R, A, gm, prom$p, bg, strat); sx <- pair_stats(R, A, tm, trans$p, bg, strat)
    nu <- n - rk0 - K
    pp <- function(r) 2 * stats::pt(-abs(r * sqrt((nu - 1) / pmax(1 - r^2, 1e-12))), nu - 1)
    data.table(lineage = lin, K = K, resid_df = nu, bg_min_size = bg$nb_min,
               promoter_pairs = nrow(prom), trans_pairs = nrow(trans),
               promoter_rate = mean(sp$p < 0.05 & sp$r > 0), trans_rate = mean(sx$p < 0.05 & sx$r > 0),
               promoter_median_r = median(sp$r), trans_median_r = median(sx$r),
               promoter_rate_param = mean(pp(sp$r) < 0.05 & sp$r > 0),
               trans_rate_param = mean(pp(sx$r) < 0.05 & sx$r > 0))
  }))
  calib[, eligible := trans_rate <= 0.05]
  Kc <- if (any(calib$eligible)) calib[eligible == TRUE][which.max(promoter_rate), K] else NA_integer_
  calib[, chosen := !is.na(Kc) & K == Kc]
  informative <- !is.na(Kc) && calib[K == Kc, promoter_rate >= 2 * trans_rate]
  calib_all[[lin]] <- calib
  tic(lin, ": chosen K ", Kc, ", informative ", informative)
  Kuse <- if (is.na(Kc)) calib[which.min(trans_rate), K] else Kc

  # full data
  R <- resid_K(lr_full, X0, Kuse); A <- resid_K(la_full, X0, Kuse)
  bg <- bg_stats(R, A, gchr, pchr, strat, S)
  st <- pair_stats(R, A, pr$g, pr$p, bg, strat)
  c1 <- bh_call(st); c2 <- bh_call(st, i2); c2a <- bh_call(st, i2a)
  tic(lin, ": full data Q1 links ", sum(c1$call), ", Q2 links ", sum(c2$call), ", Q2-any links ", sum(c2a$call))
  t_full <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # bootstrap
  lk_idx <- sort(unique(c(which(c1$call), i2[c2$call], i2a[c2a$call])))
  boot_one <- function(b) {
    set.seed(SEED + 1000L + 100L * k + b)
    ib <- sample.int(n, n, replace = TRUE); ub <- u[ib]; Xb <- design0(ub)
    dfb <- uniqueN(ub$donor_id) - qr(Xb)$rank - Kuse
    if (dfb < MIN_DF_BOOT) return(list(b = b, ok = FALSE, df = dfb))
    Rb <- resid_K(norm_logcpm(rna[kg, cols[ib]]), Xb, Kuse)
    Ab <- resid_K(norm_logcpm(atac[kp, cols[ib]]), Xb, Kuse)
    bgb <- bg_stats(Rb, Ab, gchr, pchr, strat, S)
    sb <- pair_stats(Rb, Ab, pr$g, pr$p, bgb, strat)
    list(b = b, ok = TRUE, df = dfb, unique_donors = uniqueN(ub$donor_id),
         c1 = which(bh_call(sb)$call), c2 = which(bh_call(sb, i2)$call), c2a = which(bh_call(sb, i2a)$call),
         zl = sb$z[lk_idx], rl = sb$r[lk_idx])
  }
  rep1 <- integer(nrow(pr)); rep2 <- integer(length(i2)); rep2a <- integer(length(i2a)); nb_ok <- 0L
  t_b0 <- Sys.time()
  if (N_BOOT > 0L && (informative || MODE == "timing")) {
    bl <- mclapply(seq_len(N_BOOT), boot_one, mc.cores = NCORE, mc.preschedule = FALSE)
    if (any(vapply(bl, inherits, TRUE, "try-error"))) stop("bootstrap replicate failed: ", bl[vapply(bl, inherits, TRUE, "try-error")][[1]])
    for (x in bl) if (x$ok) { rep1[x$c1] <- rep1[x$c1] + 1L; rep2[x$c2] <- rep2[x$c2] + 1L; rep2a[x$c2a] <- rep2a[x$c2a] + 1L }
    nb_ok <- sum(vapply(bl, `[[`, TRUE, "ok"))
    bootrep_all[[lin]] <- rbindlist(lapply(bl, function(x) data.table(lineage = lin, b = x$b, ok = x$ok, resid_df = x$df,
      unique_donors = if (x$ok) x$unique_donors else NA_integer_,
      q1_calls = if (x$ok) length(x$c1) else NA_integer_, q2_calls = if (x$ok) length(x$c2) else NA_integer_)))
  }
  t_boot <- as.numeric(difftime(Sys.time(), t_b0, units = "secs"))
  # post hoc (amendment 2026-10-25): per-replicate z and r of the full-data links
  if (MODE == "posthoc" && N_BOOT > 0L && length(lk_idx))
    fwrite(rbindlist(lapply(bl, function(x) if (x$ok) data.table(lineage = lin, b = x$b, gene = gt$symbol[pr$g[lk_idx]],
      peak_id = pt$peak_id[pr$p[lk_idx]], z_full = st$z[lk_idx], r_full = st$r[lk_idx], z_boot = x$zl, r_boot = x$rl))),
      file.path(RES, paste0(TAG, "link_boot_z_", lin, ".tsv.gz")), sep = "\t")
  tic(lin, ": ", N_BOOT, " bootstrap replicates in ", round(t_boot), " s (", nb_ok, " non-degenerate)")
  # degenerate replicates count as not reproduced (denominator N_BOOT)
  f1 <- rep1 / max(N_BOOT, 1L); f2 <- rep2 / max(N_BOOT, 1L); f2a <- rep2a / max(N_BOOT, 1L)

  # split halves
  sp_rows <- list()
  if (MODE == "full" && informative) for (s in seq_len(N_SPLIT)) {
    set.seed(SEED + 5000L + 100L * k + s)
    o <- sample.int(n); hA <- sort(o[seq_len(ceiling(n / 2))]); hB <- sort(setdiff(o, hA))
    XA <- design0(u[hA]); XB <- design0(u[hB])
    dfA <- length(hA) - qr(XA)$rank - Kuse; dfB <- length(hB) - qr(XB)$rank - Kuse
    row <- data.table(lineage = lin, split = s, n_A = length(hA), n_B = length(hB), df_A = dfA, df_B = dfB,
                      used = dfA >= MIN_DF_SPLIT && dfB >= MIN_DF_SPLIT, links_A = NA_integer_, replication_rate = NA_real_)
    if (row$used) {
      RA <- resid_K(norm_logcpm(rna[kg, cols[hA]]), XA, Kuse); AA <- resid_K(norm_logcpm(atac[kp, cols[hA]]), XA, Kuse)
      sA <- pair_stats(RA, AA, pr$g, pr$p, bg_stats(RA, AA, gchr, pchr, strat, S), strat)
      la <- which(bh_call(sA)$call)
      RB <- resid_K(norm_logcpm(rna[kg, cols[hB]]), XB, Kuse); AB <- resid_K(norm_logcpm(atac[kp, cols[hB]]), XB, Kuse)
      sB <- pair_stats(RB, AB, pr$g[la], pr$p[la], bg_stats(RB, AB, gchr, pchr, strat, S), strat)
      row[, `:=`(links_A = length(la), replication_rate = if (length(la)) mean(sB$r > 0 & sB$p < 0.05) else NA_real_)]
    }
    sp_rows[[s]] <- row
  }
  split_all[[lin]] <- rbindlist(sp_rows)

  # outputs
  out <- data.table(lineage = lin, gene_index = gt$gene_index[pr$g], ensembl = gt$ensembl[pr$g], symbol = gt$symbol[pr$g],
                    gene_type = gt$gene_type[pr$g], peak_index = pt$peak_index[pr$p], peak_id = pt$peak_id[pr$p],
                    dist = pr$dist, r = st$r, z = st$z, p = st$p, q1 = c1$q, link_q1 = c1$call, repro_q1 = f1,
                    in_q2 = pr$in_q2, in_q2any = pr$in_q2any)
  out[i2, `:=`(q2 = c2$q, link_q2 = c2$call, repro_q2 = f2)]
  out[i2a, `:=`(q2any = c2a$q, link_q2any = c2a$call, repro_q2any = f2a)]
  fwrite(out[link_q1 == TRUE], file.path(RES, paste0(TAG, "links_q1_", lin, ".tsv.gz")), sep = "\t")
  fwrite(out[in_q2any == TRUE], file.path(RES, paste0(TAG, "q2_pairs_", lin, ".tsv.gz")), sep = "\t")
  fwrite(gt[, .(gene_index, ensembl, symbol, gene_type, chr, tss)], file.path(W, paste0(TAG, "tested_genes_", lin, ".tsv")), sep = "\t")
  fwrite(data.table(peak_index = pt$peak_index, stratum = strat, gc_q = stt$gc_q, acc_d = stt$acc_d, wid_t = stt$wid_t,
                    mean_logcpm = rowMeans(la_full)), file.path(W, paste0(TAG, "tested_peaks_", lin, ".tsv")), sep = "\t")
  brk <- seq(0, 1, 0.1)
  hist_all[[lin]] <- rbind(
    data.table(lineage = lin, family = "Q1", bin_lo = head(brk, -1), n = tabulate(findInterval(f1[c1$call], brk, rightmost.closed = TRUE), 10)),
    data.table(lineage = lin, family = "Q2", bin_lo = head(brk, -1), n = tabulate(findInterval(f2[c2$call], brk, rightmost.closed = TRUE), 10)))
  summ[[lin]] <- data.table(lineage = lin, status = "tested", donors = n, wells = uniqueN(u$well_id),
    genes = length(kg), peaks = length(kp), strata = S, cis_pairs = nrow(pr), K = Kc, K_used = Kuse,
    informative = informative,
    promoter_rate = calib[K == Kuse, promoter_rate], trans_rate = calib[K == Kuse, trans_rate],
    q1_links = sum(c1$call), q1_stable = sum(c1$call & f1 >= STABLE),
    q1_frac_stable = if (sum(c1$call)) mean(f1[c1$call] >= STABLE) else NA_real_,
    q1_median_repro = if (sum(c1$call)) median(f1[c1$call]) else NA_real_,
    q2_variant_peaks = uniqueN(pr$p[i2]), q2_pairs = length(i2), q2_links = sum(c2$call), q2_stable = sum(c2$call & f2 >= STABLE),
    q2any_pairs = length(i2a), q2any_links = sum(c2a$call), q2any_stable = sum(c2a$call & f2a >= STABLE),
    n_boot = N_BOOT, n_boot_nondegenerate = nb_ok,
    splits_used = if (nrow(split_all[[lin]])) sum(split_all[[lin]]$used) else 0L,
    split_median_replication = if (nrow(split_all[[lin]]) && any(split_all[[lin]]$used)) median(split_all[[lin]]$replication_rate, na.rm = TRUE) else NA_real_,
    seconds_full_data = round(t_full), seconds_per_boot_replicate = if (N_BOOT) round(t_boot / N_BOOT * min(NCORE, N_BOOT)) else NA_real_)
}
fwrite(rbindlist(summ, fill = TRUE), file.path(RES, paste0(TAG, "link_summary.tsv")), sep = "\t")
fwrite(rbindlist(calib_all), file.path(RES, paste0(TAG, "link_calibration.tsv")), sep = "\t")
fwrite(rbindlist(hist_all), file.path(RES, paste0(TAG, "link_reproduction_hist.tsv")), sep = "\t")
fwrite(rbindlist(split_all, fill = TRUE), file.path(RES, paste0(TAG, "link_split_half.tsv")), sep = "\t")
fwrite(rbindlist(bootrep_all), file.path(RES, paste0(TAG, "link_bootstrap_replicates.tsv")), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(A4, paste0("sessionInfo_02", if (MODE != "full") paste0("_", MODE), ".txt")))
tic("done")
