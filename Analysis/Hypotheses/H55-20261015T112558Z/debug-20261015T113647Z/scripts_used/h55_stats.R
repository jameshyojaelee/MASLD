# H55: GTEx v8 postmortem liver; do canonical bulk disease-state calls change with
# image-measured histologic fat in the same direction as disease-vs-control logFC?
# Implements PRESPEC_H55.md (sha256 293114cd...). Do not change logic without DEVIATIONS.md.
#
# Usage: Rscript h55_stats.R <run_dir> <morph_fat.tsv> <mode: final|debug> [B] [Bboot] [Bcal]
# mode=debug: gates are logged but do not stop; used only with a permuted (fake) exposure.
suppressMessages({library(edgeR); library(limma); library(data.table)})
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
a <- commandArgs(trailingOnly = TRUE)
RUN <- a[1]; MORPH <- a[2]; MODE <- a[3]
B <- if (length(a) >= 4) as.integer(a[4]) else 10000L
BBOOT <- if (length(a) >= 5) as.integer(a[5]) else 2000L
BCAL <- if (length(a) >= 6) as.integer(a[6]) else 1000L
SEED <- 20260923L
stopifnot(MODE %in% c("final", "debug"))
H55 <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H55-20261015T112558Z/"
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
RR <- paste0(BASE, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/")
PL <- paste0(RR, "scratch/plan-H55/")
DG <- file.path(RUN, "diagnostics"); dir.create(DG, showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(RUN, "outputs"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file(file.path(RUN, "h55_stats.log"), open = "wt")
say <- function(...) { s <- paste0(...); cat(s, "\n"); cat(s, "\n", file = LOG); flush(LOG) }
res <- list(mode = MODE, B = B, BBOOT = BBOOT, BCAL = BCAL)
save_res <- function() saveRDS(res, file.path(OUT, "h55_results.rds"))
STOP <- function(rule, msg) {
  say("STOP RULE ", rule, " FIRED: ", msg)
  res$stop_rule <<- rule; res$stop_msg <<- msg
  if (MODE == "final") { res$verdict <<- paste0("untestable (", rule, ")"); save_res(); quit(status = 0) }
}
say("H55 stats, mode ", MODE, ", B ", B, ", Bboot ", BBOOT, ", Bcal ", BCAL, ", start ", format(Sys.time(), tz = "UTC"))

## 0. Hashes: prespec and every listed input; refuse on mismatch
for (f in c("prespec.sha256", "inputs.sha256")) {
  rc <- system(paste0("sha256sum --quiet -c ", H55, "scripts/", f), intern = FALSE)
  if (rc != 0) { say("HASH MISMATCH in ", f, "; refusing to run"); quit(status = 2) }
}
say("hash check: prespec and ", length(readLines(paste0(H55, "scripts/inputs.sha256"))), " inputs OK")

## 1. Donor set, re-derived from the raw GTEx annotations
sa <- fread(paste0(BASE, "data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
            select = c("SAMPID", "SMTSD", "SMAFRZE", "SMRIN", "SMTSISCH", "SMCENTER", "SMEXNCRT"))
subj <- fread(paste0(BASE, "data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"))
lv <- sa[SMTSD == "Liver" & SMAFRZE == "RNASEQ"]; lv[, donor := sub("^(GTEX-[^-]+)-.*$", "\\1", SAMPID)]
sl <- fread(paste0(BASE, "Analysis/MASLD_Model_Benchmark/executions/he-molecular-c1-20260924T004124Z/liver_slides_matched.tsv"))
both <- intersect(lv$donor, sl$PatientID)
m <- merge(lv[donor %in% both], subj, by.x = "donor", by.y = "SUBJID")
say("liver RNASEQ donors ", uniqueN(lv$donor), "; slide donors ", uniqueN(sl$PatientID), "; both ", length(both),
    "; DTHHRDY present ", sum(!is.na(m$DTHHRDY)), "; SMCENTER present ", sum(m$SMCENTER != "" & !is.na(m$SMCENTER)))
m <- m[!is.na(DTHHRDY) & !is.na(SMCENTER) & SMCENTER != ""]
dm <- fread(paste0(PL, "design_meta_counts_only.tsv")); dm <- dm[!is.na(SMCENTER) & SMCENTER != ""]
setkey(m, donor); setkey(dm, donor)
chk <- nrow(m) == 222 && identical(m$donor, dm$donor) && identical(m$SAMPID, dm$SAMPID) &&
  all(m$DTHHRDY == dm$H) && all(m$SMCENTER == dm$SMCENTER) && all(m$AGE == dm$AGE) && all(m$SEX == dm$SEX) &&
  isTRUE(all.equal(m$SMRIN, dm$SMRIN)) && isTRUE(all.equal(m$SMTSISCH, dm$SMTSISCH)) && isTRUE(all.equal(m$SMEXNCRT, dm$SMEXNCRT))
if (!chk) { say("model donor table does not match design_meta_counts_only.tsv; refusing"); quit(status = 2) }
ng <- fread(paste0(RR, "scratch/gate-L02-01/note_grades.tsv"))
meta <- merge(dm[, .(donor, SAMPID, H, SMCENTER, AGE, SEX, SMRIN, SMTSISCH, SMEXNCRT, grade, pct, fib_or_nash)],
              ng[, .(donor, ord)], by = "donor", all.x = TRUE, sort = TRUE)
meta[, ageg := factor(fcase(AGE %in% c("20-29", "30-39"), "20-39", AGE == "40-49", "40-49", AGE == "50-59", "50-59",
                            AGE %in% c("60-69", "70-79"), "60-79"), levels = c("20-39", "40-49", "50-59", "60-79"))]
meta[, hg := fcase(H == 0, "H0", H %in% 1:2, "H12", H %in% 3:4, "H34")]
meta[, stratum := paste(SMCENTER, hg)]
meta[, fib := as.integer(fib_or_nash %in% c(TRUE, "True", "TRUE"))]
say("model donors ", nrow(meta), "; strata ", paste(names(table(meta$stratum)), table(meta$stratum), collapse = ", "),
    "; graded (ord present) ", sum(!is.na(meta$ord)), "; no_mention ", sum(meta$grade == "no_mention"), "; fib/NASH notes ", sum(meta$fib))

## 2. Morphometry gates S1, S2, S3, S3b and diagnostics 1-5
mf <- fread(MORPH)
order_tsv <- fread(paste0(RR, "jobs/3111/order.tsv"))
stopifnot(nrow(order_tsv) == 224, all(order_tsv$uid %in% sl$SeriesInstanceUID))
n_err <- sum(startsWith(mf$status, "error")); n_ok <- sum(mf$status == "ok")
st_center <- merge(mf[, .(donor, status = ifelse(startsWith(status, "error"), "error", status))],
                   meta[, .(donor, SMCENTER)], by = "donor", all.x = TRUE)
fwrite(st_center[, .N, by = .(SMCENTER, status)], file.path(DG, "d1_morph_status_by_center.tsv"), sep = "\t")
fwrite(mf[startsWith(status, "error")], file.path(DG, "d1_morph_errors.tsv"), sep = "\t")
meta <- merge(meta, mf[status == "ok", .(donor, fat_frac, tissue_px, white_frac, n_droplets)], by = "donor", all.x = TRUE, sort = TRUE)
n_fat <- sum(is.finite(meta$fat_frac))
say("S1: slides ", nrow(mf), ", ok ", n_ok, ", error ", n_err, "; model donors with fat fraction ", n_fat)
res$S1 <- list(slides = nrow(mf), ok = n_ok, error = n_err, model_with_fat = n_fat)
if (nrow(mf) != 224) STOP("S1", paste("morphometry table has", nrow(mf), "rows, not 224"))
if (n_err > 10 || n_fat < 210) STOP("S1", sprintf("errors %d (>10?) or model donors with fat %d (<210?)", n_err, n_fat))
pil <- fread(paste0(RR, "jobs/3111/morph_fat.tsv"))[status == "ok"]
cmp <- merge(pil[, .(uid, fat_pilot = fat_frac)], mf[, .(uid, status, fat_frac)], by = "uid", all.x = TRUE)
cmp[, rel := abs(fat_frac - fat_pilot) / pmax(abs(fat_pilot), .Machine$double.xmin)]
max_rel <- if (any(is.na(cmp$fat_frac))) Inf else max(cmp$rel)
say("S2: pilot slides ", nrow(cmp), "; reproduced ok ", sum(cmp$status == "ok", na.rm = TRUE), "; max relative difference ", signif(max_rel, 3))
fwrite(cmp, file.path(DG, "d2_pilot_reproduction.tsv"), sep = "\t")
res$S2 <- list(n = nrow(cmp), max_rel = max_rel)
if (nrow(cmp) != 119 || !is.finite(max_rel) || max_rel > 1e-9) STOP("S2", sprintf("max relative difference %.3g over %d pilot slides", max_rel, nrow(cmp)))
gr <- meta[!is.na(ord) & is.finite(fat_frac)]
s3 <- cor(gr$fat_frac, gr$ord, method = "spearman")
s3_pct <- with(meta[!is.na(pct) & is.finite(fat_frac)], c(n = length(pct), rho = cor(fat_frac, pct, method = "spearman")))
s3_c <- gr[, .(n = .N, rho = cor(fat_frac, ord, method = "spearman"), median_fat = median(fat_frac)), by = SMCENTER]
say(sprintf("S3: Spearman fat vs note ordinal, graded model donors n %d rho %.3f; stated pct n %d rho %.3f", nrow(gr), s3, s3_pct[1], s3_pct[2]))
res$S3 <- list(n = nrow(gr), rho = s3, pct = s3_pct, by_center = s3_c)
d3 <- rbind(meta[is.finite(fat_frac), .(what = "grade", level = grade, n = .N, median = median(fat_frac)), by = grade][, grade := NULL],
            meta[is.finite(fat_frac), .(what = "hardy_group", level = hg, n = .N, median = median(fat_frac)), by = hg][, hg := NULL],
            meta[is.finite(fat_frac), .(what = "center", level = SMCENTER, n = .N, median = median(fat_frac)), by = SMCENTER][, SMCENTER := NULL])
fwrite(d3, file.path(DG, "d3_fat_by_grade_hardy_center.tsv"), sep = "\t")
fwrite(s3_c, file.path(DG, "d3_spearman_within_center.tsv"), sep = "\t")
if (s3 < 0.5) STOP("S3", sprintf("Spearman %.3f < 0.5", s3))
nm <- meta[grade == "no_mention" & is.finite(fat_frac)]
s3b <- median(nm$fat_frac)
d4 <- rbind(meta[grade %in% c("no_mention", "unspecified") & is.finite(fat_frac),
                 .(n = .N, median = median(fat_frac), q25 = quantile(fat_frac, 0.25), q75 = quantile(fat_frac, 0.75)), by = .(grade, SMCENTER)],
            meta[grade %in% c("no_mention", "unspecified") & is.finite(fat_frac),
                 .(SMCENTER = "all", n = .N, median = median(fat_frac), q25 = quantile(fat_frac, 0.25), q75 = quantile(fat_frac, 0.75)), by = grade])
fwrite(d4, file.path(DG, "d4_no_mention_unspecified_fat.tsv"), sep = "\t")
say(sprintf("S3b: median fat of %d no-mention model donors %.5f (stop if > 0.0335; prediction <= 0.0056)", nrow(nm), s3b))
res$S3b <- list(n = nrow(nm), median = s3b, prediction_met = s3b <= 0.0056)
if (s3b > 0.0335) STOP("S3b", sprintf("median %.5f > 0.0335", s3b))

## exposure (ranks among donors with a fat fraction)
fi <- which(is.finite(meta$fat_frac)); nf <- length(fi)
rk <- function(x) (rank(x, ties.method = "average") - 0.5) / length(x)
meta[, fat_rank := NA_real_]; meta[fi, fat_rank := rk(fat_frac)]
meta[, tissue_rank := NA_real_]; meta[fi, tissue_rank := rk(tissue_px)]
meta[, isch_rank := NA_real_]; meta[fi, isch_rank := rk(SMTSISCH)]
meta[, fat_rank_center := NA_real_]; meta[fi, fat_rank_center := rk(fat_frac), by = SMCENTER]
meta[, logfat := log10(fat_frac + 1e-4)]
agenum <- as.integer(meta$ageg)
d5 <- data.table(var = c("SMTSISCH", "SMRIN", "age_group", "sex", "Hardy", "center_C1", "SMEXNCRT", "tissue_px"),
  rho = sapply(list(meta$SMTSISCH, meta$SMRIN, agenum, meta$SEX, meta$H, as.integer(meta$SMCENTER == "C1"), meta$SMEXNCRT, meta$tissue_px),
               function(z) cor(meta$fat_rank[fi], z[fi], method = "spearman")))
fwrite(d5, file.path(DG, "d5_fat_rank_vs_covariates.tsv"), sep = "\t")

## 3. Counts, filter, TMM, family (S4)
cnt <- fread(paste0(PL, "inputs/gene_reads_2017-06-05_v8_liver.gct.gz"), skip = 2)
cnt <- cnt[!grepl("_PAR_Y$", Name)]
ens <- sub("\\..*$", "", cnt$Name); stopifnot(!anyDuplicated(ens), nrow(cnt) == 56156)
sym <- setNames(cnt$Description, ens)
M <- as.matrix(cnt[, meta$SAMPID, with = FALSE]); rownames(M) <- ens; colnames(M) <- meta$donor
rm(cnt); invisible(gc())
y0 <- DGEList(M); keep <- filterByExpr(y0)
say("filterByExpr pass: ", sum(keep))
y <- calcNormFactors(y0[keep, , keep.lib.sizes = FALSE], method = "TMM"); rm(y0, M); invisible(gc())
deg <- fread(paste0(BASE, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"))
deg[, ens := sub("\\..*$", "", gene)]
can <- deg[padj < 0.05 & abs(logFC) > 0.5]
fam_now <- can[can$ens %in% rownames(y)]
fam <- fread(paste0(PL, "primary_family_genes.tsv"))
s4_ok <- nrow(can) == 1854 && setequal(fam_now$ens, fam$ens) && nrow(fam) == 1002 &&
  isTRUE(all.equal(fam_now[match(fam$ens, ens)]$logFC, fam$bulk_logFC))
say("S4: canonical ", nrow(can), "; in GTEx ", sum(can$ens %in% ens), "; family recomputed ", nrow(fam_now), "; equal to frozen ", s4_ok)
res$S4 <- list(canonical = nrow(can), in_gtex = sum(can$ens %in% ens), family = nrow(fam_now), equal = s4_ok)
if (!s4_ok) STOP("S4", "recomputed family differs from primary_family_genes.tsv")
fg <- fam$ens; lfc <- fam$bulk_logFC

## 4. Exposure-free voom weights (222 donors) and primary fit
covf <- ~ ageg + factor(SEX) + factor(H) + SMTSISCH + SMRIN + SMCENTER
Z0all <- model.matrix(covf, data = meta)
v <- voomWithQualityWeights(y, design = Z0all, plot = FALSE, save.plot = TRUE)
sw_q <- v$targets$sample.weights
say(sprintf("quality weights range %.3f to %.3f; donors with weight < 0.3: %d", min(sw_q), max(sw_q), sum(sw_q < 0.3)))
pdf(file.path(DG, "d7_voom_trend.pdf"), width = 3.5, height = 3, useDingbats = FALSE, pointsize = 6)
par(mar = c(3, 3, 1, 1), mgp = c(1.6, 0.4, 0), cex = 1)
plot(v$voom.xy, pch = 16, cex = 0.15, col = "#00000033", xlab = "log2(count size + 0.5)", ylab = "sqrt(standard deviation)")
lines(v$voom.line, col = "#C9265E"); dev.off()
fwrite(data.table(donor = meta$donor, sample_weight = sw_q, center = meta$SMCENTER, H = meta$H), file.path(DG, "d7_quality_weights.tsv"), sep = "\t")
res$d7 <- list(range = range(sw_q), n_below_0.3 = sum(sw_q < 0.3))

mf_ <- meta[fi]; Ef <- v$E[, fi, drop = FALSE]; Wf <- v$weights[, fi, drop = FALSE]
rownames(Wf) <- rownames(Ef)
Z0 <- model.matrix(covf, data = mf_)
Xfull <- model.matrix(~ fat_rank + ageg + factor(SEX) + factor(H) + SMTSISCH + SMRIN + SMCENTER, data = mf_)
d6 <- 1 - summary(lm(mf_$fat_rank ~ Z0[, -1]))$r.squared  # centred R^2 (Z0 column 1 is the intercept)
say(sprintf("regression donors %d; design columns %d; residual df %d; 1 - R2 of fat rank on covariates %.3f", nf, ncol(Xfull), nf - ncol(Xfull), d6))
res$d6 <- d6; res$n <- nf; res$df_res <- nf - ncol(Xfull)
fit <- lmFit(Ef, Xfull, weights = Wf); eb <- eBayes(fit)
res$eb_prior <- c(df.prior = eb$df.prior, s2.prior = eb$s2.prior)

## FL machinery
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
# moderated t under permutation with fixed prior (Family C calibration, statistic variant)
null_modt <- function(P, IDX, d0, s20) {
  Bn <- nrow(IDX); G <- ncol(P$e); dfr <- P$n - P$p - 1; out <- matrix(NA_real_, Bn, G); tI <- t(IDX)
  for (j in seq_len(G)) {
    ej <- P$e[, j]; Epi <- matrix(ej[tI], P$n, Bn)
    num <- drop(crossprod(P$xr[, j], Epi)); proj <- colSums(crossprod(P$Q[[j]], Epi)^2)
    rss <- sum(ej^2) - proj - num^2 / P$sxx[j]; s2 <- rss / dfr
    s2p <- if (is.finite(d0)) (d0 * s20 + dfr * s2) / (d0 + dfr) else rep(s20, Bn)
    out[, j] <- (num / P$sxx[j]) / sqrt(s2p / P$sxx[j])
  }
  out
}

## 5. Primary permutation matrix (six strata) and identity check S5
IDX <- perm_matrix(mf_$stratum, B)
saveRDS(IDX, file.path(OUT, "perm_index_primary.rds"))
P <- fl_prep(Ef, Wf, Z0, mf_$fat_rank, fg)
id_diff <- max(abs(P$slope - fit$coefficients[fg, "fat_rank"]))
say(sprintf("S5 identity check: max |FL slope - lmFit coefficient| = %.3g over %d family genes", id_diff, length(fg)))
res$S5 <- id_diff
if (!is.finite(id_diff) || id_diff > 1e-6) STOP("S5", sprintf("identity difference %.3g > 1e-6", id_diff))
# moderated t identity (diagnostic)
d0 <- eb$df.prior; s20 <- eb$s2.prior
t_chk <- null_modt(P, matrix(seq_len(nf), 1), d0, s20)[1, ]
res$modt_identity <- max(abs(t_chk - eb$t[fg, "fat_rank"]))
say(sprintf("moderated-t identity (diagnostic): max diff %.3g", res$modt_identity))

## 6. Primary statistic
rho_obs <- sp(lfc, P$slope)
NS <- null_slopes(P, IDX)
null_p <- null_rho(NS, lfc)
pp <- flp(rho_obs, null_p)
res$primary <- list(rho = rho_obs, p = pp, n_genes = length(fg), n_donors = nf,
                    null = c(mean = mean(null_p), sd = sd(null_p), q025 = unname(quantile(null_p, 0.025)), q975 = unname(quantile(null_p, 0.975))))
fwrite(data.table(rho_null = null_p), file.path(OUT, "null_rho_primary.tsv.gz"), sep = "\t")
say(sprintf("PRIMARY: rho %.4f, FL p two-sided %.4g (upper %.4g, lower %.4g); null mean %.4f sd %.4f [%.4f, %.4f]",
            rho_obs, pp["p_two"], pp["p_upper"], pp["p_lower"], mean(null_p), sd(null_p), quantile(null_p, 0.025), quantile(null_p, 0.975)))
res$d9_flag <- sd(null_p) < 0.15 || sd(null_p) > 0.30

## 7. Positive and negative controls
plin <- "ENSG00000147872"
res$PC2 <- list(slope = eb$coefficients[plin, "fat_rank"], p = eb$p.value[plin, "fat_rank"],
                pass = eb$coefficients[plin, "fat_rank"] > 0 && eb$p.value[plin, "fat_rank"] < 0.05)
say(sprintf("PC2 PLIN2: slope %.4f, moderated p %.3g, pass %s", res$PC2$slope, res$PC2$p, res$PC2$pass))
pcg <- fread(paste0(PL, "pc_composite_genes.tsv")); stopifnot(nrow(pcg) == 197, all(pcg$ens %in% rownames(Ef)))
Ppc <- fl_prep(Ef, Wf, Z0, mf_$fat_rank, pcg$ens)
rho_pc3 <- sp(pcg$h8_mu, Ppc$slope); null_pc3 <- null_rho(null_slopes(Ppc, IDX), pcg$h8_mu); p_pc3 <- flp(rho_pc3, null_pc3)
res$PC3 <- list(rho = rho_pc3, p = p_pc3, pass = rho_pc3 > 0 && p_pc3["p_two"] < 0.05, in_predicted_range = rho_pc3 >= 0.3 && rho_pc3 <= 0.8)
say(sprintf("PC3: rho %.4f, FL p %.4g, S6 fires %s", rho_pc3, p_pc3["p_two"], !res$PC3$pass))
Znc <- cbind(Z0, fat_rank = mf_$fat_rank)
Pnc <- fl_prep(Ef, Wf, Znc, mf_$tissue_rank, fg)
rho_nc1 <- sp(lfc, Pnc$slope); null_nc1 <- null_rho(null_slopes(Pnc, IDX), lfc); p_nc1 <- flp(rho_nc1, null_nc1)
nc1_fail <- p_nc1["p_two"] < 0.05 && abs(rho_nc1) >= 0.5 * abs(rho_obs)
res$NC1 <- list(rho = rho_nc1, p = p_nc1, fail = unname(nc1_fail), prediction_met = abs(rho_nc1) < 0.2 && p_nc1["p_two"] > 0.05)
say(sprintf("NC1 tissue area: rho %.4f, FL p %.4g, S7 fires %s", rho_nc1, p_nc1["p_two"], nc1_fail))

## 8. Family C: per-gene slopes, moderated t, calibration
NT_cal <- null_modt(P, IDX[seq_len(min(BCAL, B)), , drop = FALSE], d0, s20)
dft <- eb$df.total[match(fg, rownames(fit$coefficients))]
any_q <- apply(NT_cal, 1, function(t) any(p.adjust(2 * pt(-abs(t), dft), "BH") < 0.05))
cal_frac <- mean(any_q)
use_emp <- cal_frac > 0.07
p_lim <- eb$p.value[fg, "fat_rank"]
p_emp <- sapply(seq_along(fg), function(j) (1 + sum(abs(NS[, j]) >= abs(P$slope[j]))) / (B + 1))
p_col <- if (use_emp) p_emp else p_lim
q_col <- p.adjust(p_col, "BH")
res$famC <- list(cal_frac = cal_frac, use_empirical = use_emp, n_q05 = sum(q_col < 0.05),
                 agree = if (any(q_col < 0.05)) mean(sign(P$slope[q_col < 0.05]) == sign(lfc[q_col < 0.05])) else NA)
say(sprintf("Family C calibration: fraction of %d null draws with any q < 0.05 = %.3f; use empirical p %s; genes q<0.05 %d",
            nrow(NT_cal), cal_frac, use_emp, res$famC$n_q05))

## 9. Family B (Q1, Q2) and NC2
# Q1: Hardy 1-4, fixed filter/TMM/family, own voomWQW, global fat ranks
q1i <- which(meta$H != 0 & is.finite(meta$fat_frac)); mq1 <- meta[q1i]; mq1[, Hf := factor(H, levels = 1:4)]
Z0q1 <- model.matrix(~ ageg + factor(SEX) + Hf + SMTSISCH + SMRIN + SMCENTER, data = mq1)
vq1 <- voomWithQualityWeights(y[, q1i], design = Z0q1, plot = FALSE); rownames(vq1$weights) <- rownames(vq1$E)
Xq1 <- cbind(Z0q1[, 1, drop = FALSE], fat_rank = mq1$fat_rank, Z0q1[, -1])
fq1 <- lmFit(vq1$E, Xq1, weights = vq1$weights)
Pq1 <- fl_prep(vq1$E, vq1$weights, Z0q1, mq1$fat_rank, fg)
stopifnot(max(abs(Pq1$slope - fq1$coefficients[fg, "fat_rank"])) < 1e-6)
IDXq1 <- perm_matrix(mq1$stratum, B); saveRDS(IDXq1, file.path(OUT, "perm_index_Q1.rds"))
rho_q1 <- sp(lfc, Pq1$slope); p_q1 <- flp(rho_q1, null_rho(null_slopes(Pq1, IDXq1), lfc))
say(sprintf("Q1 Hardy 1-4: n %d, cols %d, df %d, rho %.4f, FL p two-sided %.4g", nrow(mq1), ncol(Xq1), nrow(mq1) - ncol(Xq1), rho_q1, p_q1["p_two"]))
# Q2: joint fat + ischemic-time ranks; covariates without SMTSISCH
covq2 <- ~ ageg + factor(SEX) + factor(H) + SMRIN + SMCENTER
vq2 <- voomWithQualityWeights(y, design = model.matrix(covq2, data = meta), plot = FALSE)
Eq2 <- vq2$E[, fi]; Wq2 <- vq2$weights[, fi]; rownames(Wq2) <- rownames(Eq2)
Z0q2 <- model.matrix(covq2, data = mf_)
Xq2 <- cbind(Z0q2[, 1, drop = FALSE], fat_rank = mf_$fat_rank, isch_rank = mf_$isch_rank, Z0q2[, -1])
fq2 <- lmFit(Eq2, Xq2, weights = Wq2)
sf <- fq2$coefficients[fg, "fat_rank"]; si <- fq2$coefficients[fg, "isch_rank"]
s_sign <- if (rho_obs < 0) -1 else 1
rho_f2 <- sp(lfc, sf); rho_i2 <- sp(lfc, si); D_obs <- s_sign * (rho_f2 - rho_i2)
# paired donor bootstrap within six strata
set.seed(SEED); us <- sort(unique(mf_$stratum))
BIDX <- matrix(0L, BBOOT, nf)
for (b in seq_len(BBOOT)) { p <- integer(nf); for (s in us) { i <- which(mf_$stratum == s); p[i] <- i[sample.int(length(i), length(i), replace = TRUE)] }; BIDX[b, ] <- p }
saveRDS(BIDX, file.path(OUT, "boot_index_Q2.rds"))
SWq2 <- sqrt(Wq2[fg, , drop = FALSE]); Yq2 <- Eq2[fg, , drop = FALSE]
Dstar <- numeric(BBOOT); n_drop <- 0L
for (b in seq_len(BBOOT)) {
  ib <- BIDX[b, ]; Xb <- Xq2[ib, , drop = FALSE]
  kc <- c(TRUE, apply(Xb[, -1, drop = FALSE], 2, function(z) length(unique(z)) > 1))
  if (!all(kc)) n_drop <- n_drop + 1L
  Xb <- Xb[, kc, drop = FALSE]; jf <- which(colnames(Xb) == "fat_rank"); ji <- which(colnames(Xb) == "isch_rank")
  bf <- numeric(length(fg)); bi <- numeric(length(fg))
  for (j in seq_along(fg)) { w <- SWq2[j, ib]; cf <- qr.coef(qr(Xb * w), Yq2[j, ib] * w); bf[j] <- cf[jf]; bi[j] <- cf[ji] }
  Dstar[b] <- s_sign * (sp(lfc, bf) - sp(lfc, bi))
}
p_q2 <- (1 + sum(Dstar <= 0)) / (BBOOT + 1)
# joint-permutation p for D (descriptive): FL on covariates without both exposures
Pj_e <- matrix(0, nf, length(fg)); A <- vector("list", length(fg))
for (j in seq_along(fg)) {
  sw <- sqrt(Wq2[fg[j], ]); Q <- qr.Q(qr(Z0q2 * sw)); yw <- Eq2[fg[j], ] * sw
  Xe <- cbind(mf_$fat_rank, mf_$isch_rank) * sw; Xr <- Xe - Q %*% crossprod(Q, Xe)
  Pj_e[, j] <- yw - Q %*% crossprod(Q, yw); A[[j]] <- solve(crossprod(Xr), t(Xr))
}
stopifnot(max(abs(sapply(seq_along(fg), function(j) (A[[j]] %*% Pj_e[, j])[1]) - sf)) < 1e-6)
tI <- t(IDX); NF2 <- matrix(0, B, length(fg)); NI2 <- matrix(0, B, length(fg))
for (j in seq_along(fg)) { Epi <- matrix(Pj_e[, j][tI], nf, B); s2 <- A[[j]] %*% Epi; NF2[, j] <- s2[1, ]; NI2[, j] <- s2[2, ] }
Dnull <- s_sign * (null_rho(NF2, lfc) - null_rho(NI2, lfc)); p_Dperm <- flp(D_obs, Dnull)
rm(NF2, NI2); invisible(gc())
say(sprintf("Q2: rho_fat %.4f rho_isch %.4f D %.4f; bootstrap one-sided p %.4g (resamples with a dropped column %d); joint-permutation p %.4g",
            rho_f2, rho_i2, D_obs, p_q2, n_drop, p_Dperm["p_two"]))
nc2_met <- rho_i2 >= -0.2 && rho_i2 <= 0.3 && rho_i2 < rho_f2
famA <- pp["p_two"] < 0.05
holm <- p.adjust(c(Q1 = unname(p_q1["p_two"]), Q2 = p_q2), "holm")
q1_sup <- famA && holm["Q1"] < 0.05 && sign(rho_q1) == sign(rho_obs)
q2_sup <- famA && holm["Q2"] < 0.05 && D_obs > 0
res$Q1 <- list(n = nrow(mq1), rho = rho_q1, p = p_q1, p_holm = holm["Q1"], supports = unname(q1_sup), tested = unname(famA))
res$Q2 <- list(rho_fat = rho_f2, rho_isch = rho_i2, D = D_obs, D_boot_ci = quantile(Dstar, c(0.025, 0.975)), p_boot = p_q2,
               p_holm = holm["Q2"], p_joint_perm = p_Dperm, n_resample_col_dropped = n_drop, supports = unname(q2_sup), tested = unname(famA))
res$NC2 <- list(rho_isch = rho_i2, prediction_met = nc2_met)

## 10. Descriptive versions (no alpha)
desc <- list()
add_desc <- function(name, rho, p, n, g) { desc[[name]] <<- data.table(version = name, rho = rho, p_two = unname(p["p_two"]), n_donors = n, n_genes = g,
    null_q025 = unname(p["null_q025"]), null_q975 = unname(p["null_q975"]))
  say(sprintf("descriptive %-34s rho %.4f p %.4g n %d genes %d", name, rho, p["p_two"], n, g)) }
sub_version <- function(name, di, covform, xname, genes = fg, ref = lfc, strat = NULL) {
  md <- meta[di]; Zs <- model.matrix(covform, data = md)
  Zs <- Zs[, c(TRUE, apply(Zs[, -1, drop = FALSE], 2, function(z) length(unique(z)) > 1)), drop = FALSE]
  qz <- qr(Zs); Zs <- Zs[, qz$pivot[seq_len(qz$rank)], drop = FALSE]
  vs <- voomWithQualityWeights(y[, di], design = Zs, plot = FALSE); rownames(vs$weights) <- rownames(vs$E)
  Ps <- fl_prep(vs$E, vs$weights, Zs, md[[xname]], genes)
  IDXs <- perm_matrix(if (is.null(strat)) md$stratum else md[[strat]], B)
  saveRDS(IDXs, file.path(OUT, paste0("perm_index_", name, ".rds")))
  r <- sp(ref, Ps$slope); add_desc(name, r, flp(r, null_rho(null_slopes(Ps, IDXs), ref)), length(di), length(genes))
}
same_version <- function(name, Zn, x, genes = fg, ref = lfc, E = Ef, W = Wf) {
  Ps <- fl_prep(E, W, Zn, x, genes); r <- sp(ref, Ps$slope)
  add_desc(name, r, flp(r, null_rho(null_slopes(Ps, IDX), ref)), ncol(E), length(genes)); invisible(Ps)
}
# Hardy 0 only (strata by center)
meta[, center_stratum := SMCENTER]
sub_version("hardy0_only", which(meta$H == 0 & is.finite(meta$fat_frac)), ~ ageg + factor(SEX) + SMTSISCH + SMRIN + SMCENTER, "fat_rank", strat = "center_stratum")
# family without TNFA/acute-phase genes
exs <- fread(paste0(PL, "exclusion_set_tnfa_acutephase.tsv")); ex_keep <- !(fg %in% exs$ensembl_gene)
r_ex <- sp(lfc[ex_keep], P$slope[ex_keep]); add_desc("family_minus_tnfa_acutephase", r_ex, flp(r_ex, null_rho(NS[, ex_keep], lfc[ex_keep])), nf, sum(ex_keep))
# fibrosis/steatohepatitis notes: (a) covariate (own voom weights on 222), (b) excluded
vfa <- voomWithQualityWeights(y, design = model.matrix(update(covf, ~ . + fib), data = meta), plot = FALSE)
Wfa <- vfa$weights[, fi]; rownames(Wfa) <- rownames(Ef)
same_version("fibrosis_flag_covariate", model.matrix(update(covf, ~ . + fib), data = mf_), mf_$fat_rank, E = vfa$E[, fi], W = Wfa)
sub_version("fibrosis_notes_excluded", which(meta$fib == 0 & is.finite(meta$fat_frac)), covf, "fat_rank")
# composition scores (Guilliams 2022 consensus markers passing filter, not family genes)
mk <- fread(paste0(BASE, "data/external/liver_celltype_markers_combined.tsv"))[source == "Guilliams2022_consensus_ttest"]
immune <- c("B cells","T cells","Macrophages","Mono+mono derived cells","Neutrophils","Resident NK","Circulating NK/NKT","Plasma cells","cDC1s","cDC2s","pDCs","Basophils")
msets <- list(hep = mk[cell_type == "Hepatocytes"]$gene, imm = unique(mk[cell_type %in% immune]$gene), str = unique(mk[cell_type %in% c("Fibroblasts", "Endothelial cells")]$gene))
Zs_all <- t(scale(t(v$E)))
for (nmk in names(msets)) {
  g <- rownames(v$E)[sym[rownames(v$E)] %in% msets[[nmk]] & !(rownames(v$E) %in% fg)]
  meta[, (paste0("comp_", nmk)) := colMeans(Zs_all[g, , drop = FALSE])]
  say("composition score ", nmk, ": marker genes used ", length(g))
}
mf_ <- meta[fi]
covc <- update(covf, ~ . + comp_hep + comp_imm + comp_str)
vc <- voomWithQualityWeights(y, design = model.matrix(covc, data = meta), plot = FALSE)
Wc <- vc$weights[, fi]; rownames(Wc) <- rownames(Ef)
same_version("composition_scores_covariates", model.matrix(covc, data = mf_), mf_$fat_rank, E = vc$E[, fi], W = Wc)
# exonic rate covariate
cove <- update(covf, ~ . + SMEXNCRT)
ve <- voomWithQualityWeights(y, design = model.matrix(cove, data = meta), plot = FALSE)
We <- ve$weights[, fi]; rownames(We) <- rownames(Ef)
same_version("exonic_rate_covariate", model.matrix(cove, data = mf_), mf_$fat_rank, E = ve$E[, fi], W = We)
# exposure variants
meta[, ord_x := ord]
sub_version("note_ordinal_graded_donors", which(!is.na(meta$ord) & is.finite(meta$fat_frac)), covf, "ord_x")
same_version("fat_rank_within_center", Z0, mf_$fat_rank_center)
same_version("log10_fat_fraction", Z0, mf_$logfat)
# gene-set variants
tpm <- fread(paste0(BASE, "data/external/allelic_refs/gtex_v8/gene_tpm_2017-06-05_v8_liver.gct.gz"), skip = 2)
tpm <- tpm[!grepl("_PAR_Y$", Name)]; TM <- as.matrix(tpm[, meta$SAMPID, with = FALSE]); rownames(TM) <- sub("\\..*$", "", tpm$Name); rm(tpm)
medtpm <- apply(TM[fg, ], 1, median); k535 <- medtpm >= 1
r535 <- sp(lfc[k535], P$slope[k535]); add_desc("median_tpm_ge1_within_family", r535, flp(r535, null_rho(NS[, k535], lfc[k535])), nf, sum(k535))
can_g <- can[can$ens %in% rownames(TM)]
univ <- union(rownames(y), can_g$ens)
LT <- log2(TM[univ, fi] + 1)
ft <- lmFit(LT, Xfull)
Pt <- fl_prep(LT, matrix(1, nrow(LT), nf, dimnames = list(rownames(LT), NULL)), Z0, mf_$fat_rank, can_g$ens)
stopifnot(max(abs(Pt$slope - ft$coefficients[can_g$ens, "fat_rank"])) < 1e-6)
rt <- sp(can_g$logFC, Pt$slope); add_desc("all_canonical_limma_trend_log2tpm", rt, flp(rt, null_rho(null_slopes(Pt, IDX), can_g$logFC)), nf, nrow(can_g))
# statistic variant: moderated t
NT <- null_modt(P, IDX, d0, s20); tobs <- eb$t[fg, "fat_rank"]
r_t <- sp(lfc, tobs); add_desc("moderated_t_statistic", r_t, flp(r_t, null_rho(NT, lfc)), nf, length(fg)); rm(NT)
desc_dt <- rbindlist(desc); fwrite(desc_dt, file.path(OUT, "descriptive_versions.tsv"), sep = "\t")
res$desc <- desc_dt

## 11. Per-gene table for all 1,854 calls and Catalog column
crit <- qt(0.975, dft)
se_m <- sqrt(eb$s2.post[match(fg, rownames(fit$coefficients))]) * fit$stdev.unscaled[fg, "fat_rank"]
gt <- data.table(ensembl_id = fg, symbol = fam$symbol, bulk_logFC = lfc, slope = P$slope, ci_low = P$slope - crit * se_m, ci_high = P$slope + crit * se_m,
                 moderated_t = eb$t[fg, "fat_rank"], p_limma = p_lim, p_empirical = p_emp, p_used = p_col, q_BH = q_col)
gt[, state := ifelse(q_BH < 0.05, paste0("associated_", ifelse(slope > 0, "up", "down"), "_", ifelse(sign(slope) == sign(bulk_logFC), "agrees", "disagrees")), "not_detected")]
allc <- data.table(ensembl_id = can$ens, symbol = can$symbol, bulk_logFC = can$logFC)
allc <- merge(allc, gt[, !c("symbol", "bulk_logFC")], by = "ensembl_id", all.x = TRUE, sort = FALSE)
allc[is.na(state), state := ifelse(ensembl_id %in% ens, "not_measurable_below_filter", "not_measurable_no_gencode_v26_match")]
fwrite(allc, file.path(OUT, "per_gene_1854_calls.tsv"), sep = "\t")
res$state_counts <- allc[, .N, by = state]
fwrite(data.table(rho_null = null_p, pc3_null = null_pc3, nc1_null = null_nc1), file.path(OUT, "null_rho_primary_pc3_nc1.tsv.gz"), sep = "\t")
fwrite(data.table(gene = fg, bulk_logFC = lfc, slope = P$slope, nc1_slope = Pnc$slope, q1_slope = Pq1$slope, q2_fat = sf, q2_isch = si),
       file.path(OUT, "family_slopes.tsv"), sep = "\t")
fwrite(data.table(gene = pcg$ens, h8_mu = pcg$h8_mu, slope = Ppc$slope), file.path(OUT, "pc3_slopes.tsv"), sep = "\t")
fwrite(meta[, .(donor, SAMPID, H, SMCENTER, ageg, SEX, SMRIN, SMTSISCH, grade, ord, fib, fat_frac, tissue_px, fat_rank, tissue_rank, isch_rank)],
       file.path(OUT, "donor_exposure_table.tsv"), sep = "\t")

## 12. Verdict by the prespec rule
s6 <- !res$PC3$pass; s7 <- res$NC1$fail
if (s6 || s7) {
  res$verdict <- paste0("not interpretable (", paste(c("S6"[s6], "S7"[s7]), collapse = ", "), ")")
  res$catalog_column_written <- FALSE
} else {
  res$verdict <- if (!famA) "null: alignment not detected" else if (rho_obs > 0) "pass: same direction" else "significant opposite direction"
  res$catalog_column_written <- TRUE
  cat_col <- allc[, .(ensembl_id, gtex_postmortem_fat_slope = slope, gtex_postmortem_fat_ci_low = ci_low, gtex_postmortem_fat_ci_high = ci_high,
                      gtex_postmortem_fat_q = q_BH, gtex_postmortem_fat_state = state)]
  fwrite(cat_col, file.path(OUT, "catalog_column_gtex_postmortem_fat.tsv"), sep = "\t")
}
res$tnfa_wording <- famA && !s6 && !s7 && rho_obs > 0 && desc$family_minus_tnfa_acutephase$p_two >= 0.05
say("VERDICT: ", res$verdict)
save_res()
writeLines(capture.output(sessionInfo()), file.path(RUN, "sessionInfo.txt"))
say("done ", format(Sys.time(), tz = "UTC"))
