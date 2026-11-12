#!/usr/bin/env Rscript
# TRAIL4: colocalization at the TNFSF10 and PANX1 serum-ALT loci. Prespecification: ../PRESPEC_TRAIL4.md
# Pairwise coloc.abf (p12 = 1e-5 primary, 5e-6 sensitivity) and SuSiE-coloc with 1000G EUR LD.
# SNP harmonisation follows MQTL2 02_coloc.R (same allele pair, palindromes and duplicated positions dropped).
# Usage: Rscript 01_trail4_coloc.R <out_dir>
suppressMessages({ library(data.table); library(coloc); library(susieR) })
set.seed(20261111)
OUT <- commandArgs(TRUE)[1]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/TRAIL4-20261111T025501Z")
SS <- file.path(ROOT, "GWAS/finemapping/data/sumstats")
LDDIR <- file.path(ROOT, "GWAS/finemapping/data/ld_ref/1kg_eur")
TABIX <- "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/tabix"
PANUKBB <- file.path(ROOT, "GWAS/MR_Data/PanUKBB/biomarkers-30620-ALT-both_sexes-irnt.tsv.bgz")
PQ <- c(UKB = file.path(RUN, "downloads/GCST90470917/GCST90470917.tsv.gz"),
        ZHAO = file.path(RUN, "downloads/GCST90274843/GCST90274843.tsv.gz"))
N <- c(ALT_MVP = 406803, ALB_MVP = 372962, PLT_MVP = 407021, ALT_PANUKBB = 400822,
       TRAIL_UKB = 47745, TRAIL_ZHAO = 14735)
P1 <- 1e-4; P2 <- 1e-4; P12 <- 1e-5; P12_ALT <- 5e-6
MIN_SNPS <- 50L
WIN <- data.table(locus = c("TNFSF10", "PANX1", "PNPLA3"), chr = c(3L, 11L, 22L),
                  lead = c(172294500L, 93870338L, 44324727L))
WIN[, `:=`(s = lead - 500000L, e = lead + 500000L)]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cat("TRAIL4 | coloc", as.character(packageVersion("coloc")), "| susieR", as.character(packageVersion("susieR")),
    "|", date(), "\n")

in_win_awk <- function(ci, pi) paste(sprintf("($%d==%d && $%d>=%d && $%d<=%d)", ci, WIN$chr, pi, WIN$s, pi, WIN$e),
                                     collapse = " || ")
add_locus <- function(d) {
  d[, locus := NA_character_]
  for (i in seq_len(nrow(WIN))) d[chr == WIN$chr[i] & pos >= WIN$s[i] & pos <= WIN$e[i], locus := WIN$locus[i]]
  d[!is.na(locus)]
}

# ---------------- read (GRCh37 throughout; build checked below) ----------------
read_mvp <- function(tr) {
  f <- file.path(SS, sprintf("MVP_%s_EUR_reformatted_hg19.tsv", tr))
  d <- fread(cmd = sprintf("awk -F'\\t' 'NR==1 || %s' %s", in_win_awk(1, 2), f))
  setnames(d, c("chromosome", "position", "allele1", "allele2", "pval"), c("chr", "pos", "ea", "oa", "p"))
  d[, eaf := NA_real_]
  add_locus(d[!is.na(beta) & !is.na(se) & se > 0])
}
read_pq <- function(f) {
  d <- fread(cmd = sprintf("zcat %s | awk -F'\\t' 'NR==1 || %s'", f, in_win_awk(1, 2)))
  setnames(d, c("chromosome", "base_pair_location", "effect_allele", "other_allele", "standard_error",
                "effect_allele_frequency", "p_value"), c("chr", "pos", "ea", "oa", "se", "eaf", "p"))
  if ("rs_id" %in% names(d)) setnames(d, "rs_id", "rsid")
  d[, chr := as.integer(chr)]
  add_locus(d[!is.na(beta) & !is.na(se) & se > 0])
}
read_panukbb <- function() {
  d <- rbindlist(lapply(seq_len(nrow(WIN)), function(i) fread(
    cmd = sprintf("%s %s %d:%d-%d | cut -f1-4,19,25,31,37,43", TABIX, PANUKBB, WIN$chr[i], WIN$s[i], WIN$e[i]),
    header = FALSE, sep = "\t", na.strings = "NA",
    col.names = c("chr", "pos", "oa", "ea", "eaf", "beta", "se", "nlp", "lowconf"))))
  d <- d[lowconf %in% c("false", FALSE) & !is.na(beta) & !is.na(se) & se > 0]
  d[, `:=`(p = 10^(-nlp), lowconf = NULL, nlp = NULL)]
  add_locus(d)
}
D <- list(ALT_MVP = read_mvp("ALT"), ALB_MVP = read_mvp("Albumin"), PLT_MVP = read_mvp("Platelet"),
          ALT_PANUKBB = read_panukbb(), TRAIL_UKB = read_pq(PQ["UKB"]), TRAIL_ZHAO = read_pq(PQ["ZHAO"]))
for (k in names(D)) D[[k]][, `:=`(ea = toupper(ea), oa = toupper(oa))]
print(rbindlist(lapply(names(D), function(k) D[[k]][, .N, by = locus][, trait := k])))

# ---------------- C1 + build checks ----------------
chk <- rbindlist(lapply(names(D), function(k) {
  d <- D[[k]]
  rbindlist(lapply(list(c("rs738409", 22, 44324727, "G", "C"), c("TNFSF10_lead", 3, 172294500, "G", "A"),
                        c("PANX1_lead", 11, 93870338, "C", "T")), function(v) {
    x <- d[chr == as.integer(v[2]) & pos == as.integer(v[3]) &
             ((ea == v[4] & oa == v[5]) | (ea == v[5] & oa == v[4]))]
    data.table(trait = k, variant = v[1], pos = paste(v[2], v[3], sep = ":"), n_found = nrow(x),
               file_alleles = if (nrow(x)) paste0(x$ea[1], "(effect)/", x$oa[1]) else NA_character_,
               rsid_in_file = if (nrow(x) && "rsid" %in% names(x)) as.character(x$rsid[1]) else NA_character_,
               beta_for_named_allele = if (nrow(x)) ifelse(x$ea[1] == v[4], x$beta[1], -x$beta[1]) else NA_real_,
               named_allele = v[4], z = if (nrow(x)) x$beta[1] / x$se[1] * ifelse(x$ea[1] == v[4], 1, -1) else NA_real_)
  }))
}))
fwrite(chk, file.path(OUT, "control_C1_build_allele_checks.tsv"), sep = "\t")
print(chk)
c1 <- chk[variant == "rs738409" & trait %in% c("ALT_MVP", "ALT_PANUKBB")]
stopifnot(nrow(c1) == 2, all(c1$n_found == 1), all(c1$beta_for_named_allele > 0))
stopifnot(chk[trait == "ALT_MVP" & variant == "TNFSF10_lead"]$beta_for_named_allele > 0)
# build: every file must carry rs738409 at its GRCh37 position (pQTL files: also check the rsid column)
build_ok <- chk[variant == "rs738409", .(trait, ok = n_found == 1 & (is.na(rsid_in_file) | rsid_in_file == "rs738409"))]
print(build_ok); stopifnot(all(build_ok$ok))

# ---------------- harmonise and coloc.abf ----------------
pal <- function(a, b) (a %in% c("A", "T") & b %in% c("A", "T")) | (a %in% c("C", "G") & b %in% c("C", "G"))
harm <- function(t1, t2, loc) {
  a <- D[[t1]][locus == loc & nchar(ea) == 1 & nchar(oa) == 1]
  b <- D[[t2]][locus == loc & nchar(ea) == 1 & nchar(oa) == 1]
  a <- a[!(duplicated(pos) | duplicated(pos, fromLast = TRUE))]
  b <- b[!(duplicated(pos) | duplicated(pos, fromLast = TRUE))]
  m <- merge(a[, .(pos, ea1 = ea, oa1 = oa, b1 = beta, se1 = se, p1 = p, f1 = eaf)],
             b[, .(pos, ea2 = ea, oa2 = oa, b2 = beta, se2 = se, p2 = p, f2 = eaf)], by = "pos")
  m <- m[(ea1 == ea2 & oa1 == oa2) | (ea1 == oa2 & oa1 == ea2)]
  m[ea2 != ea1, `:=`(b2 = -b2, f2 = 1 - f2)]
  m[!pal(ea1, oa1)][order(pos)]
}
mk_ds <- function(beta, se, f, n, snp) {
  ds <- list(beta = beta, varbeta = se^2, N = n, type = "quant", snp = snp)
  if (all(is.finite(f))) ds$MAF <- pmin(f, 1 - f) else ds$sdY <- 1
  ds
}
# INT-scaled traits (MVP, Pan-UKBB irnt): sdY = 1. pQTL: sdY estimated by coloc from MAF and N.
is_int <- function(t) !grepl("^TRAIL", t)
run_abf <- function(t1, t2, loc, label) {
  m <- harm(t1, t2, loc)
  if (nrow(m) < MIN_SNPS) return(list(row = data.table(test = label, locus = loc, trait1 = t1, trait2 = t2,
                                                       n_snps = nrow(m), note = "too few SNPs")))
  snp <- paste0(WIN[locus == loc]$chr, ":", m$pos, ":", m$oa1, ":", m$ea1)
  d1 <- if (is_int(t1)) list(beta = m$b1, varbeta = m$se1^2, N = N[[t1]], type = "quant", sdY = 1, snp = snp) else
    mk_ds(m$b1, m$se1, m$f1, N[[t1]], snp)
  d2 <- if (is_int(t2)) list(beta = m$b2, varbeta = m$se2^2, N = N[[t2]], type = "quant", sdY = 1, snp = snp) else
    mk_ds(m$b2, m$se2, m$f2, N[[t2]], snp)
  a <- suppressMessages(coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = P12))
  b <- suppressMessages(coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = P12_ALT))
  ti <- which.max(a$results$SNP.PP.H4)
  i1 <- which.min(m$p1); i2 <- which.min(m$p2)
  pp <- a$summary; pb <- b$summary
  row <- data.table(test = label, locus = loc, trait1 = t1, trait2 = t2, n_snps = nrow(m),
                    min_p1 = m$p1[i1], min_p2 = m$p2[i2],
                    PP.H0 = pp[["PP.H0.abf"]], PP.H1 = pp[["PP.H1.abf"]], PP.H2 = pp[["PP.H2.abf"]],
                    PP.H3 = pp[["PP.H3.abf"]], PP.H4 = pp[["PP.H4.abf"]],
                    s.PP.H0 = pb[["PP.H0.abf"]], s.PP.H1 = pb[["PP.H1.abf"]], s.PP.H2 = pb[["PP.H2.abf"]],
                    s.PP.H3 = pb[["PP.H3.abf"]], s.PP.H4 = pb[["PP.H4.abf"]],
                    top_h4_snp = snp[ti], top_h4_snp_pp = a$results$SNP.PP.H4[ti],
                    top_h4_effect_allele = m$ea1[ti], top_h4_b1 = m$b1[ti], top_h4_se1 = m$se1[ti],
                    top_h4_b2 = m$b2[ti], top_h4_se2 = m$se2[ti],
                    lead1_snp = snp[i1], lead2_snp = snp[i2], note = "")
  list(row = row, m = m, snp = snp, abf = a)
}

# ---------------- LD (1000G EUR, GRCh37; r signed to bim A1 with --keep-allele-order) ----------------
ld_for <- function(loc, tag) {
  w <- WIN[locus == loc]; pre <- file.path(OUT, "ld", paste0(tag))
  dir.create(dirname(pre), showWarnings = FALSE)
  if (!file.exists(paste0(pre, ".ld"))) {
    system(sprintf("plink --bfile %s/chr%d_eur --chr %d --from-bp %d --to-bp %d --snps-only just-acgt --keep-allele-order --make-bed --out %s_w > /dev/null",
                   LDDIR, w$chr, w$chr, w$s, w$e, pre))
    system(sprintf("plink --bfile %s_w --keep-allele-order --r square spaces --out %s > /dev/null", pre, pre))
  }
  bim <- fread(paste0(pre, "_w.bim"), header = FALSE, col.names = c("chr", "id", "cm", "pos", "a1", "a2"))
  list(bim = bim, f = paste0(pre, ".ld"))
}
LDC <- list()
get_ld <- function(loc) { if (is.null(LDC[[loc]])) LDC[[loc]] <<- ld_for(loc, loc); LDC[[loc]] }

run_susie <- function(r, loc) {
  L <- get_ld(loc); bim <- L$bim
  bim[, k := .I]
  bim <- bim[!(duplicated(pos) | duplicated(pos, fromLast = TRUE))]
  m <- r$m
  x <- merge(m[, .(pos, ea1, oa1, z1 = b1 / se1, z2 = b2 / se2)], bim[, .(pos, a1, a2, k)], by = "pos")
  x <- x[(ea1 == a1 & oa1 == a2) | (ea1 == a2 & oa1 == a1)]
  x[ea1 == a2, `:=`(z1 = -z1, z2 = -z2)]
  x <- x[order(k)]
  if (nrow(x) < MIN_SNPS) return(NULL)
  R <- as.matrix(fread(L$f, header = FALSE, select = x$k))[x$k, , drop = FALSE]
  R[!is.finite(R)] <- 0
  dimnames(R) <- list(as.character(x$pos), as.character(x$pos))
  snp <- paste0(WIN[locus == loc]$chr, ":", x$pos, ":", x$a2, ":", x$a1)
  dimnames(R) <- list(snp, snp)
  fit <- function(z, n) {
    lam <- tryCatch(estimate_s_rss(z, R, n = n), error = function(e) NA_real_)
    s <- tryCatch(susie_rss(z = setNames(z, snp), R = R, n = n, L = 10, coverage = 0.95,
                            estimate_residual_variance = FALSE, max_iter = 500),
                  error = function(e) NULL)
    list(s = s, lambda = lam, ncs = if (is.null(s) || is.null(s$sets$cs)) 0L else length(s$sets$cs))
  }
  f1 <- fit(x$z1, N[[r$row$trait1]]); f2 <- fit(x$z2, N[[r$row$trait2]])
  summ <- data.table(test = r$row$test, n_snps_ld = nrow(x), lambda1 = f1$lambda, lambda2 = f2$lambda,
                     ncs1 = f1$ncs, ncs2 = f2$ncs)
  cs_tab <- rbindlist(lapply(list(list("trait1", f1), list("trait2", f2)), function(o) {
    s <- o[[2]]$s; if (is.null(s) || o[[2]]$ncs == 0) return(NULL)
    rbindlist(lapply(seq_along(s$sets$cs), function(j) {
      v <- s$sets$cs[[j]]; top <- v[which.max(s$pip[v])]
      data.table(test = r$row$test, which = o[[1]], cs = names(s$sets$cs)[j], size = length(v),
                 top_snp = snp[top], top_pip = s$pip[top], top_z = c(x$z1, x$z2)[top + if (o[[1]] == "trait2") nrow(x) else 0])
    }))
  }))
  cs <- NULL
  if (f1$ncs >= 1 && f2$ncs >= 1) {
    for (pr in c(P12, P12_ALT)) {
      cc <- tryCatch(suppressMessages(coloc.susie(f1$s, f2$s, p1 = P1, p2 = P2, p12 = pr)), error = function(e) NULL)
      if (!is.null(cc) && !is.null(cc$summary))
        cs <- rbind(cs, data.table(test = r$row$test, p12 = pr, as.data.table(cc$summary)), fill = TRUE)
    }
  }
  # LD between the two traits' lead SNPs (signed to bim A1)
  i1 <- which.max(abs(x$z1)); i2 <- which.max(abs(x$z2))
  summ[, `:=`(lead1 = snp[i1], lead2 = snp[i2], r_lead1_lead2 = R[i1, i2])]
  list(summ = summ, cs_tab = cs_tab, cs = cs, x = x, R = R, snp = snp)
}

call_pp <- function(h0, h1, h2, h3, h4) fifelse(h4 >= 0.8, "colocalized", fifelse(h3 >= 0.8, "distinct",
  fifelse(h0 + h1 >= 0.8, "no signal in trait2", fifelse(h0 + h2 >= 0.8, "no signal in trait1", "unresolved"))))

tests <- list(
  c("C3_ALT_MVP_vs_PanUKBB_PNPLA3", "ALT_MVP", "ALT_PANUKBB", "PNPLA3"),
  c("C2_TRAIL_UKB_vs_ZHAO", "TRAIL_UKB", "TRAIL_ZHAO", "TNFSF10"),
  c("T1_ALT_MVP_vs_TRAIL_UKB", "ALT_MVP", "TRAIL_UKB", "TNFSF10"),
  c("T1s_ALT_PanUKBB_vs_TRAIL_UKB", "ALT_PANUKBB", "TRAIL_UKB", "TNFSF10"),
  c("T1r_ALT_MVP_vs_TRAIL_ZHAO", "ALT_MVP", "TRAIL_ZHAO", "TNFSF10"),
  c("T2_ALT_vs_ALB_TNFSF10", "ALT_MVP", "ALB_MVP", "TNFSF10"),
  c("T2_ALT_vs_PLT_TNFSF10", "ALT_MVP", "PLT_MVP", "TNFSF10"),
  c("T2_ALB_vs_PLT_TNFSF10", "ALB_MVP", "PLT_MVP", "TNFSF10"),
  c("S_ALT_vs_ALB_PANX1", "ALT_MVP", "ALB_MVP", "PANX1"),
  c("S_ALT_vs_PLT_PANX1", "ALT_MVP", "PLT_MVP", "PANX1"),
  c("S_ALB_vs_PLT_PANX1", "ALB_MVP", "PLT_MVP", "PANX1"))
abf <- list(); sus <- list(); cst <- list(); csc <- list()
for (t in tests) {
  cat("\n==", t[1], "\n")
  r <- run_abf(t[2], t[3], t[4], t[1]); abf[[t[1]]] <- r$row; print(r$row)
  if (!is.null(r$m) && r$row$n_snps >= MIN_SNPS) {
    s <- tryCatch(run_susie(r, t[4]), error = function(e) { cat("susie error:", conditionMessage(e), "\n"); NULL })
    if (!is.null(s)) { sus[[t[1]]] <- s$summ; cst[[t[1]]] <- s$cs_tab; csc[[t[1]]] <- s$cs; print(s$summ); print(s$cs) }
    fwrite(r$m, file.path(OUT, paste0("harmonised_", t[1], ".tsv.gz")), sep = "\t")
  }
}
A <- rbindlist(abf, fill = TRUE)
A[, call_primary := call_pp(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)]
A[, call_p12_5e6 := call_pp(s.PP.H0, s.PP.H1, s.PP.H2, s.PP.H3, s.PP.H4)]
fwrite(A, file.path(OUT, "coloc_abf.tsv"), sep = "\t")
fwrite(rbindlist(sus, fill = TRUE), file.path(OUT, "susie_summary.tsv"), sep = "\t")
fwrite(rbindlist(cst, fill = TRUE), file.path(OUT, "susie_credible_sets.tsv"), sep = "\t")
fwrite(rbindlist(csc, fill = TRUE), file.path(OUT, "coloc_susie.tsv"), sep = "\t")

# ---------------- direction at TNFSF10 (aligned to the ALT-raising allele G at 3:172294500) ----------------
dir_rows <- rbindlist(lapply(c("ALT_MVP", "ALT_PANUKBB", "ALB_MVP", "PLT_MVP", "TRAIL_UKB", "TRAIL_ZHAO"), function(k) {
  x <- D[[k]][chr == 3 & pos == 172294500 & ((ea == "G" & oa == "A") | (ea == "A" & oa == "G"))]
  if (!nrow(x)) return(data.table(trait = k, variant = "3:172294500", found = FALSE))
  s <- ifelse(x$ea[1] == "G", 1, -1)
  data.table(trait = k, variant = "3:172294500", found = TRUE, file_effect_allele = x$ea[1], file_other_allele = x$oa[1],
             file_beta = x$beta[1], beta_per_G = s * x$beta[1], se = x$se[1], z_per_G = s * x$beta[1] / x$se[1],
             p = x$p[1], eaf_G = if (is.finite(x$eaf[1])) ifelse(s > 0, x$eaf[1], 1 - x$eaf[1]) else NA_real_)
}), fill = TRUE)
# lead pQTL variant in the window and its LD with the ALT lead
for (k in c("TRAIL_UKB", "TRAIL_ZHAO")) {
  x <- D[[k]][locus == "TNFSF10"][which.min(p)]
  dir_rows <- rbind(dir_rows, data.table(trait = k, variant = paste0("pQTL lead 3:", x$pos), found = TRUE,
                                         file_effect_allele = x$ea, file_other_allele = x$oa, file_beta = x$beta,
                                         se = x$se, p = x$p, eaf_G = NA_real_), fill = TRUE)
}
L <- get_ld("TNFSF10"); bim <- L$bim[, k := .I]
lead_k <- bim[pos == 172294500]$k
if (length(lead_k) == 1) {
  rr <- as.matrix(fread(L$f, header = FALSE, skip = lead_k - 1, nrows = 1))[1, ]
  bim[, r_with_ALT_lead := rr]
  cat("ALT lead in panel:", paste(bim[k == lead_k, .(id, a1, a2)], collapse = " "), "\n")
  pl <- sapply(c("TRAIL_UKB", "TRAIL_ZHAO"), function(k) D[[k]][locus == "TNFSF10"][which.min(p)]$pos)
  fwrite(bim[pos %in% c(172294500, pl)], file.path(OUT, "ld_alt_lead_vs_pqtl_leads.tsv"), sep = "\t")
  print(bim[pos %in% c(172294500, pl)])
} else cat("ALT lead 3:172294500 not in 1000G EUR panel window (found", length(lead_k), ")\n")
fwrite(dir_rows, file.path(OUT, "direction_TNFSF10.tsv"), sep = "\t")
print(dir_rows)
cat("\nPANX1 plasma pQTL: untestable (no open plasma PANX1 pQTL source; see PRESPEC).\n")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("done", date(), "\n")
