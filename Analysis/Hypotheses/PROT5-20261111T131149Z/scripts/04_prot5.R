#!/usr/bin/env Rscript
# PROT5 step 4: plasma hepatocyte-protein panel effect at the 169 serum-ALT leads. Prespecification: ../PRESPEC_PROT5.md
# Unit: loci. Effects aligned to the loci.tsv ALT-raising allele. Protein effects in protein SD per allele.
# Usage: Rscript 04_prot5.R <run_dir> <out_dir>
suppressMessages({ library(data.table); library(ashr) })
set.seed(20261111)
a <- commandArgs(TRUE); RUN <- a[1]; OUT <- a[2]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
LDDIR <- file.path(ROOT, "GWAS/finemapping/data/ld_ref/1kg_eur")
ENZ <- file.path(ROOT, "Analysis/Hypotheses/ENZFAT1-20261110T204401Z/out_S2/locus_effects.tsv")
MVP <- file.path(ROOT, "GWAS/finemapping/data/sumstats/MVP_ALT_EUR_reformatted_hg19.tsv")
X <- file.path(OUT, "extract")
MIN_PROT <- 5L; CIS <- 1e6; EQ_FRAC <- 0.5; SLOPE_EQ <- 0.25
log_cat <- function(...) cat(..., "\n", sep = "")

# ---------------- panels, downloads, stop rule ----------------
pan <- fread(file.path(RUN, "inputs/panels.tsv"), colClasses = c(chr = "character"))
pan <- pan[panel %in% c("hepatocyte", "negative_control") & found == TRUE]
man <- fread(file.path(RUN, "downloads_manifest.tsv"))
pan <- pan[accession %in% man[md5_ok == TRUE, accession]]
pan[, strict := panel == "hepatocyte" & meets_strict_hpa_rule %in% c(TRUE, "True")]
log_cat("panel proteins with verified files: hepatocyte ", pan[panel == "hepatocyte", .N],
        ", control ", pan[panel == "negative_control", .N])
if (pan[panel == "hepatocyte", .N] < 8) {
  writeLines("STOP: fewer than 8 hepatocyte panel proteins available", file.path(OUT, "STOPPED.txt")); quit(status = 0)
}
leads <- fread(file.path(OUT, "leads.tsv"), colClasses = c(locus = "character"))

# ---------------- proxies (1000G EUR, r2 >= 0.8 within 500 kb, signed r to bim A1) ----------------
bim_at <- function(c, pos) {
  b <- fread(cmd = sprintf("awk -F'\\t' 'BEGIN{while((getline l < \"%s\")>0) w[l]=1} ($4 in w)' %s/chr%d_eur.bim",
                           pos, LDDIR, c), header = FALSE, col.names = c("chr", "id", "cm", "pos", "A1", "A2"))
  b
}
prox <- list()
for (c in sort(unique(leads$chr))) {
  L <- leads[chr == c]
  pf <- tempfile(); writeLines(as.character(L$pos), pf)
  b <- bim_at(c, pf)
  L <- merge(L, b[, .(pos, id, A1, A2)], by = "pos")
  L <- L[(A1 == a1 & A2 == a2) | (A1 == a2 & A2 == a1)]
  L <- L[!duplicated(locus)]
  if (!nrow(L)) next
  sf <- file.path(OUT, "ld", sprintf("chr%d_leads.txt", c)); dir.create(dirname(sf), showWarnings = FALSE)
  writeLines(L$id, sf)
  pre <- file.path(OUT, "ld", sprintf("chr%d", c))
  system(sprintf("plink --bfile %s/chr%d_eur --keep-allele-order --r --ld-snp-list %s --ld-window-kb 500 --ld-window 999999 --ld-window-r2 0.8 --out %s > /dev/null",
                 LDDIR, c, sf, pre))
  if (!file.exists(paste0(pre, ".ld"))) next
  ld <- fread(paste0(pre, ".ld"))
  ld <- ld[SNP_A != SNP_B & R^2 >= 0.8]
  pf2 <- tempfile(); writeLines(as.character(unique(ld$BP_B)), pf2)
  bb <- bim_at(c, pf2)[, .(SNP_B = id, pB = pos, B1 = A1, B2 = A2)]
  ld <- merge(ld, bb, by = "SNP_B")
  ld <- merge(ld, L[, .(SNP_A = id, locus, A1, A2, alt_raising_allele)], by = "SNP_A")
  ld[, phase_allele := fifelse(alt_raising_allele == A1, fifelse(R > 0, B1, B2), fifelse(R > 0, B2, B1))]
  prox[[as.character(c)]] <- ld[, .(locus, proxy_pos = pB, B1, B2, r = R, r2 = R^2, dist = abs(pB - BP_A), phase_allele)]
}
prox <- rbindlist(prox)
setorder(prox, locus, -r2, dist)
fwrite(prox, file.path(OUT, "proxies_1000G_r2ge0.8.tsv"), sep = "\t")
need_pos <- unique(rbind(leads[, .(chr, pos)], merge(prox, leads[, .(locus, chr)], by = "locus")[, .(chr, pos = proxy_pos)]))
need_pos[, key := paste(chr, pos, sep = ":")]

# ---------------- read protein extracts ----------------
cn <- c("chr", "pos", "ea", "oa", "beta", "se", "eaf", "p", "variant_id", "rsid", "n")
rd <- function(f) {
  d <- fread(f, colClasses = list(character = c("chromosome", "effect_allele", "other_allele", "variant_id", "rs_id")))
  setnames(d, cn); d[, `:=`(ea = toupper(ea), oa = toupper(oa))]
  d[chr %in% as.character(1:22)][, chr := as.integer(chr)][]
}
pal <- function(x, y) paste0(x, y) %in% c("AT", "TA", "CG", "GC")
W <- list(); NUL <- list(); c0 <- list(); sdy <- list()
nul_ref <- fread(file.path(OUT, "null_snps.tsv"), header = FALSE, col.names = c("chr", "pos", "a1", "a2"))
for (acc in pan$accession) {
  w <- rd(file.path(X, paste0(acc, "_win.tsv.gz")))
  x <- w[chr == 22 & pos == 44324727 & ((ea == "G" & oa == "C") | (ea == "C" & oa == "G"))]
  c0[[acc]] <- data.table(accession = acc, gene = pan[accession == acc, gene], n_rows = nrow(x),
                          rsid = if (nrow(x)) x$rsid[1] else NA_character_, alleles = if (nrow(x)) paste0(x$ea[1], "/", x$oa[1]) else NA_character_)
  W[[acc]] <- w[paste(chr, pos, sep = ":") %in% need_pos$key][, accession := acc]
  nl <- rd(file.path(X, paste0(acc, "_null.tsv.gz")))
  nl <- merge(nl, nul_ref, by = c("chr", "pos"))
  nl <- nl[(ea == a1 & oa == a2) | (ea == a2 & oa == a1)]
  nl[, z := fifelse(ea == a1, beta / se, -beta / se)]
  nl[, f1 := fifelse(ea == a1, eaf, 1 - eaf)]
  sdy[[acc]] <- data.table(accession = acc, sdY = nl[eaf >= 0.01 & eaf <= 0.99, median(se * sqrt(2 * n * eaf * (1 - eaf)), na.rm = TRUE)],
                           n_null = nrow(nl))
  NUL[[acc]] <- nl[, .(key = paste(chr, pos, sep = ":"), z, f1)]
  rm(w, nl); gc(verbose = FALSE)
}
c0 <- rbindlist(c0); sdy <- rbindlist(sdy)
c0[, pass := n_rows == 1 & rsid == "rs738409"]
fwrite(c0, file.path(OUT, "control_C0_alignment.tsv"), sep = "\t"); fwrite(sdy, file.path(OUT, "sdY.tsv"), sep = "\t")
print(sdy)
fwrite(leads[, .(chr, pos)], file.path(OUT, "leads_chrpos.tsv"), sep = "\t", col.names = FALSE)
mvp <- fread(cmd = sprintf("awk -F'\\t' 'BEGIN{while((getline l < \"%s\")>0){split(l,a,\"\\t\"); w[a[1]\":\"a[2]]=1}} NR==1 || (($1\":\"$2) in w)' %s",
                           file.path(OUT, "leads_chrpos.tsv"), MVP), header = TRUE)
setnames(mvp, c("chromosome", "position", "allele1", "allele2", "pval"), c("chr", "pos", "ea", "oa", "p"))
mvp_c0 <- mvp[chr == 22 & pos == 44324727]
c0_mvp <- nrow(mvp_c0) == 1 && ifelse(mvp_c0$ea == "G", mvp_c0$beta, -mvp_c0$beta) > 0
if (!all(c0$pass) || !c0_mvp) {
  print(c0); writeLines("STOP: C0 alignment failed", file.path(OUT, "STOPPED.txt")); quit(status = 0)
}
log_cat("C0 alignment: pass in all ", nrow(c0), " protein files and MVP ALT")

# ---------------- null correlation ----------------
Z <- Reduce(function(x, y) merge(x, y, by = "key"), lapply(names(NUL), function(k) setnames(copy(NUL[[k]])[, .(key, z, f1)], c("key", k, paste0("f_", k)))))
fcol <- grep("^f_", names(Z), value = TRUE)
Zm <- as.matrix(Z[, pan$accession, with = FALSE])
f1 <- Z[[fcol[1]]]
ok <- is.finite(rowSums(Zm)) & f1 >= 0.01 & f1 <= 0.99
R4 <- cor(Zm[ok & apply(abs(Zm), 1, max) < 4, ])
R196 <- cor(Zm[ok & apply(abs(Zm), 1, max) < 1.96, ])
log_cat("null SNPs for R: |z|<4 ", sum(ok & apply(abs(Zm), 1, max) < 4), "; |z|<1.96 ", sum(ok & apply(abs(Zm), 1, max) < 1.96))
dimnames(R4) <- dimnames(R196) <- list(pan$gene, pan$gene)
fwrite(as.data.table(R4, keep.rownames = "gene"), file.path(OUT, "null_correlation_zlt4.tsv"), sep = "\t")
fwrite(as.data.table(R196, keep.rownames = "gene"), file.path(OUT, "null_correlation_zlt1.96.tsv"), sep = "\t")
dimnames(R4) <- dimnames(R196) <- list(pan$accession, pan$accession)
rm(Z, Zm, NUL); gc(verbose = FALSE)

# ---------------- per locus x protein effects ----------------
W <- rbindlist(W)
W <- W[!(duplicated(W[, .(accession, chr, pos, ea, oa)]) | duplicated(W[, .(accession, chr, pos, ea, oa)], fromLast = TRUE))]
gl <- pan[, .(accession, gene, panel, strict, gchr = suppressWarnings(as.integer(chr)), gstart = start, gend = end)]
look <- rbindlist(lapply(seq_len(nrow(leads)), function(i) {
  L <- leads[i]; pr <- prox[locus == L$locus]
  rbindlist(lapply(pan$accession, function(acc) {
    d <- W[accession == acc & chr == L$chr]
    x <- d[pos == L$pos & ((ea == L$a1 & oa == L$a2) | (ea == L$a2 & oa == L$a1))]
    src <- "lead"; ra <- L$alt_raising_allele; r2 <- 1; ppos <- L$pos
    if (nrow(x) != 1) {
      x <- NULL; src <- "missing"
      for (j in seq_len(nrow(pr))) {
        y <- d[pos == pr$proxy_pos[j] & ((ea == pr$B1[j] & oa == pr$B2[j]) | (ea == pr$B2[j] & oa == pr$B1[j]))]
        if (nrow(y) == 1) { x <- y; src <- "proxy"; ra <- pr$phase_allele[j]; r2 <- pr$r2[j]; ppos <- pr$proxy_pos[j]; break }
      }
    }
    if (is.null(x) || !nrow(x)) return(data.table(locus = L$locus, accession = acc, source = "missing"))
    data.table(locus = L$locus, accession = acc, source = src, used_pos = ppos, proxy_r2 = r2,
               raise_allele = ra, ea = x$ea, oa = x$oa, eaf = x$eaf, palindromic = pal(x$ea, x$oa),
               b_raw = fifelse(x$ea == ra, x$beta, -x$beta), se_raw = x$se, p = x$p)
  }), fill = TRUE)
}), fill = TRUE)
look <- merge(look, leads[, .(locus, lchr = chr, lpos = pos)], by = "locus")
look <- merge(look, gl, by = "accession")
look <- merge(look, sdy[, .(accession, sdY)], by = "accession")
look[, `:=`(b = b_raw / sdY, s = se_raw / sdY)]
look[, cis := !is.na(gchr) & gchr == lchr & lpos >= gstart - CIS & lpos <= gend + CIS]
fwrite(look, file.path(OUT, "protein_locus_effects.tsv.gz"), sep = "\t")
log_cat("lookups: ", paste(capture.output(print(look[, .N, by = source])), collapse = " | "))

# ---------------- panel effects ----------------
lin <- function(d, Rm, wsign = NULL) {  # IVW mean, correlation-adjusted variance
  w <- 1 / d$s^2; est <- sum(w * d$b) / sum(w)
  S <- outer(d$s, d$s) * Rm[d$accession, d$accession, drop = FALSE]
  c(est = est, se = sqrt(drop(t(w) %*% S %*% w)) / sum(w), n = nrow(d))
}
diffp <- function(dh, dc, Rm) {
  wh <- 1 / dh$s^2; wc <- 1 / dc$s^2
  cc <- c(wh / sum(wh), -wc / sum(wc)); d <- rbind(dh, dc)
  S <- outer(d$s, d$s) * Rm[d$accession, d$accession, drop = FALSE]
  c(est = sum(cc * d$b), se = sqrt(drop(t(cc) %*% S %*% cc)))
}
use <- look[source != "missing" & !cis & is.finite(b) & is.finite(s)]
pe <- function(Rm, tag) use[, {
  h <- .SD[panel == "hepatocyte"]; k <- .SD[panel == "negative_control"]; st <- .SD[strict == TRUE]
  ph <- if (nrow(h) >= MIN_PROT) lin(h, Rm) else c(est = NA, se = NA, n = nrow(h))
  pc <- if (nrow(k) >= MIN_PROT) lin(k, Rm) else c(est = NA, se = NA, n = nrow(k))
  ps <- if (nrow(st) >= MIN_PROT) lin(st, Rm) else c(est = NA, se = NA, n = nrow(st))
  pd <- if (nrow(h) >= MIN_PROT && nrow(k) >= MIN_PROT) diffp(h, k, Rm) else c(est = NA, se = NA)
  list(P_hep = ph[["est"]], se_hep = ph[["se"]], n_hep = ph[["n"]], P_ctrl = pc[["est"]], se_ctrl = pc[["se"]],
       n_ctrl = pc[["n"]], P_strict = ps[["est"]], se_strict = ps[["se"]], n_strict = ps[["n"]],
       P_diff = pd[["est"]], se_diff = pd[["se"]])
}, by = locus][, Rtag := tag][]
LP <- pe(R4, "zlt4"); LP196 <- pe(R196, "zlt1.96")

# ALT effects
enz <- fread(ENZ, colClasses = c(locus = "character"))
alt <- merge(leads[, .(locus, chr, pos, a1, a2, alt_raising_allele)],
             enz[, .(locus, ref_allele, beta_al_alt_ukb, se_alt_ukb, lfsr_pdff, genes)], by = "locus", all.x = TRUE)
alt[, `:=`(A = fifelse(ref_allele == alt_raising_allele, beta_al_alt_ukb, -beta_al_alt_ukb), A_se = se_alt_ukb)]
m2 <- merge(alt[, .(locus, chr, pos, a1, a2, alt_raising_allele)], mvp[, .(chr, pos, mea = toupper(ea), moa = toupper(oa), mb = beta, ms = se)], by = c("chr", "pos"))
m2 <- m2[(mea == a1 & moa == a2) | (mea == a2 & moa == a1)][!duplicated(locus)]
m2[, A_mvp := fifelse(mea == alt_raising_allele, mb, -mb)]
alt <- merge(alt, m2[, .(locus, A_mvp, A_mvp_se = ms)], by = "locus", all.x = TRUE)
LP <- merge(alt, LP, by = "locus", all.x = TRUE)
LP196 <- merge(alt, LP196, by = "locus", all.x = TRUE)
is_loc <- function(d) !grepl("^ctrl_", d$locus)

# ---------------- slopes ----------------
slope <- function(P, s, A, label, n_min = 3) {
  ok <- is.finite(P) & is.finite(s) & is.finite(A) & s > 0
  P <- P[ok]; s <- s[ok]; A <- A[ok]; n <- length(P)
  if (n < n_min) return(data.table(analysis = label, n_loci = n))
  w <- 1 / s^2; th <- sum(w * A * P) / sum(w * A^2); se0 <- sqrt(1 / sum(w * A^2))
  phi <- sum(w * (P - th * A)^2) / (n - 1); se <- se0 * max(1, sqrt(phi))
  data.table(analysis = label, n_loci = n, slope = th, se = se, lo95 = th - 1.96 * se, hi95 = th + 1.96 * se,
             lo90 = th - 1.645 * se, hi90 = th + 1.645 * se, z = th / se, p = 2 * pnorm(-abs(th / se)), phi = phi,
             min_detectable_slope_80pct = 2.80 * se)
}
L0 <- LP[is_loc(LP)]
res <- list()
res$primary <- slope(L0$P_hep, L0$se_hep, L0$A, "primary: hepatocyte panel on Pan-UKBB ALT")
fat <- L0[lfsr_pdff < 0.05]; nof <- L0[lfsr_pdff >= 0.05]
s1 <- slope(fat$P_hep, fat$se_hep, fat$A, "co-primary arm: liver-fat loci (lfsr_pdff < 0.05)")
s0 <- slope(nof$P_hep, nof$se_hep, nof$A, "co-primary arm: loci without liver-fat effect")
dz <- (s1$slope - s0$slope) / sqrt(s1$se^2 + s0$se^2)
cop <- data.table(analysis = "co-primary: slope(fat) - slope(no fat)", n_loci = s1$n_loci + s0$n_loci,
                  slope = s1$slope - s0$slope, se = sqrt(s1$se^2 + s0$se^2))
cop[, `:=`(lo95 = slope - 1.96 * se, hi95 = slope + 1.96 * se, lo90 = slope - 1.645 * se, hi90 = slope + 1.645 * se,
           z = dz, p = 2 * pnorm(-abs(dz)), min_detectable_slope_80pct = 2.80 * se)]
res$s1 <- s1; res$s0 <- s0; res$cop <- cop
res$ctrl <- slope(L0$P_ctrl, L0$se_ctrl, L0$A, "specificity: negative-control panel on Pan-UKBB ALT")
res$diff <- slope(L0$P_diff, L0$se_diff, L0$A, "specificity: hepatocyte minus control panel on Pan-UKBB ALT")
res$S1 <- slope(L0$P_hep, L0$se_hep, L0$A_mvp, "S1: hepatocyte panel on MVP ALT (no sample overlap)")
res$S1c <- slope(L0$P_ctrl, L0$se_ctrl, L0$A_mvp, "S1: control panel on MVP ALT")
res$S2 <- slope(L0$P_strict, L0$se_strict, L0$A, "S2: strict-HPA hepatocyte panel (15) on Pan-UKBB ALT")
L196 <- LP196[is_loc(LP196)]
res$S3 <- slope(L196$P_hep, L196$se_hep, L196$A, "S3: hepatocyte panel, null R from |z|<1.96")
ok <- is.finite(L0$P_hep) & is.finite(L0$A)
infl <- L0[ok][which.max(A^2 / se_hep^2)]
L4 <- L0[locus != infl$locus]
res$S4 <- slope(L4$P_hep, L4$se_hep, L4$A, sprintf("S4: hepatocyte panel without locus %s (%s; largest A^2/se^2)", infl$locus, infl$genes))
# per-protein slopes
pp <- merge(use, alt[, .(locus, A)], by = "locus")[!grepl("^ctrl_", locus)]
res$prot <- pp[, slope(b, s, A, paste0("per protein: ", gene[1], " (", panel[1], ")")), by = accession][, accession := NULL]
S <- rbindlist(res, fill = TRUE)
fam <- S$analysis %in% c(res$primary$analysis, cop$analysis)
S[, p_holm := NA_real_]; S[fam, p_holm := p.adjust(p, "holm")]
S[, equivalence_0.25 := fifelse(is.finite(lo90), lo90 > -SLOPE_EQ & hi90 < SLOPE_EQ, NA)]
fwrite(S, file.path(OUT, "slopes.tsv"), sep = "\t")
print(S[, .(analysis, n_loci, slope = signif(slope, 3), se = signif(se, 3), lo95 = signif(lo95, 3), hi95 = signif(hi95, 3),
            p = signif(p, 3), p_holm = signif(p_holm, 3), phi = signif(phi, 3))])

# ---------------- per-locus ashr calls ----------------
u <- is_loc(LP) & is.finite(LP$P_hep)
fit <- ash(LP$P_hep[u], LP$se_hep[u], mixcompdist = "halfuniform")
LP[u, `:=`(postmean = get_pm(fit), postsd = get_psd(fit), lfsr = get_lfsr(fit))]
uc <- !is_loc(LP) & is.finite(LP$P_hep)
if (any(uc)) {
  fc <- ash(LP$P_hep[uc], LP$se_hep[uc], g = get_fitted_g(fit), fixg = TRUE, mixcompdist = "halfuniform")
  LP[uc, `:=`(postmean = get_pm(fc), postsd = get_psd(fc), lfsr = get_lfsr(fc))]
}
log_cat("ashr pi0 ", round(get_pi0(fit), 3), "; loci fitted ", sum(u))
LP[, A_z := A / A_se]
LP[, panel_call := fcase(
  !is.finite(P_hep), "untested or underpowered",
  lfsr < 0.05 & postmean > 0, "panel rises",
  lfsr < 0.05 & postmean < 0, "panel falls",
  A_z >= 3 & P_hep - 1.645 * se_hep > -EQ_FRAC * abs(A) & P_hep + 1.645 * se_hep < EQ_FRAC * abs(A), "ALT only",
  default = "untested or underpowered")]
LP[, `:=`(z_hep = P_hep / se_hep, z_ctrl = P_ctrl / se_ctrl, z_diff = P_diff / se_diff,
          powered_for_P_eq_A = se_hep <= abs(A) / 2.80, tost_possible = A_z >= 3 & 1.645 * se_hep < EQ_FRAC * abs(A))]
fwrite(LP, file.path(OUT, "locus_panel_effects.tsv"), sep = "\t")
print(LP[is_loc(LP), .N, by = panel_call])
pw <- LP[is_loc(LP) & is.finite(P_hep), .(loci_with_panel_effect = .N, powered_for_P_eq_A = sum(powered_for_P_eq_A),
                                         tost_possible = sum(tost_possible), median_se_hep = median(se_hep),
                                         median_abs_A = median(abs(A)), median_n_hep = median(n_hep))]
fwrite(pw, file.path(OUT, "power.tsv"), sep = "\t"); print(pw)

# ---------------- controls ----------------
got <- look[locus == "ctrl_GOT1_lead" & gene == "GOT1"]
pn <- LP[locus == "ctrl_PNPLA3_rs738409"]; gt <- LP[locus == "ctrl_GOT1_lead"]; tn <- LP[locus == "34"]
ctl <- rbind(
  data.table(control = "C0 alignment: rs738409 C/G at 22:44324727 in every protein file; MVP ALT G positive",
             value = sprintf("%d/%d protein files; MVP %s", sum(c0$pass), nrow(c0), c0_mvp), pass = all(c0$pass) & c0_mvp),
  data.table(control = "C1 PNPLA3 rs738409 G: hepatocyte panel > 0, p < 0.05",
             value = sprintf("P = %.4f (se %.4f, z %.2f, %d proteins)", pn$P_hep, pn$se_hep, pn$z_hep, pn$n_hep),
             pass = isTRUE(pn$P_hep > 0 & 2 * pnorm(-abs(pn$z_hep)) < 0.05)),
  data.table(control = "C2a GOT1 cis lead 10:101152720 T: GOT1 protein z >= 5.45",
             value = if (nrow(got)) sprintf("b = %.4f SD (se %.4f, z %.2f, %s)", got$b, got$s, got$b / got$s, got$source) else "missing",
             pass = nrow(got) == 1 && isTRUE(got$b / got$s >= 5.45)),
  data.table(control = "C2b GOT1 cis lead: panel without GOT1 |z| < 2",
             value = sprintf("P = %.4f (se %.4f, z %.2f, %d proteins)", gt$P_hep, gt$se_hep, gt$z_hep, gt$n_hep),
             pass = isTRUE(abs(gt$z_hep) < 2)),
  data.table(control = "TNFSF10 locus 34 (report only)",
             value = sprintf("P = %.4f (se %.4f, z %.2f, lfsr %.3g, %s); control panel z %.2f", tn$P_hep, tn$se_hep, tn$z_hep, tn$lfsr, tn$panel_call, tn$z_ctrl),
             pass = NA))
fwrite(ctl, file.path(OUT, "controls.tsv"), sep = "\t"); print(ctl)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
log_cat("PROT5 DONE ", date())
