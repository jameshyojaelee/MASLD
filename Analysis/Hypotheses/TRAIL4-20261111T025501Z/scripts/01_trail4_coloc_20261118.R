#!/usr/bin/env Rscript
# TRAIL4 (2026-11-18 update): colocalization at the TNFSF10 and PANX1 serum-ALT loci.
# Updates over 01_trail4_coloc.R:
#   1. Parse per-SNP n vector in read_pq() (GCST90470917, GCST90274843) for coloc.abf().
#   2. Compute exact 4-trait single-causal-variant Wakefield ABF joint posterior
#      (ALT_MVP, TRAIL_UKB, PLT_MVP, ALB_MVP) across the harmonized common SNP grid.
#   3. Replace N=379 1000G EUR LD in run_susie() with polyfun_eur (N=337,000 UKBB EUR),
#      PSD-projected, max_iter=5000, set.seed(42).
#   4. Check 3:172294500 (rs79287178) and all r2 >= 0.8 proxies (plus all TNFSF10 CDS
#      variants in LD/sumstats) against TNFSF10 coding exons -> epitope_check.tsv.
# Usage: Rscript 01_trail4_coloc_20261118.R <out_dir>
suppressMessages({ library(data.table); library(coloc); library(susieR) })
set.seed(20261111)
OUT <- commandArgs(TRUE)[1]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RUN <- file.path(ROOT, "Analysis/Hypotheses/TRAIL4-20261111T025501Z")
SS <- file.path(ROOT, "GWAS/finemapping/data/sumstats")
POLYFUN_DIR <- file.path(ROOT, "GWAS/finemapping/data/ld_ref/polyfun_eur")
ONEKG_DIR <- file.path(ROOT, "GWAS/finemapping/data/ld_ref/1kg_eur")
GTF_PATH <- file.path(ROOT, "data/gsmap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf")
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
cat("TRAIL4_20261118 | coloc", as.character(packageVersion("coloc")),
    "| susieR", as.character(packageVersion("susieR")), "|", date(), "\n")

in_win_awk <- function(ci, pi) paste(sprintf("($%d==%d && $%d>=%d && $%d<=%d)", ci, WIN$chr, pi, WIN$s, pi, WIN$e),
                                     collapse = " || ")
add_locus <- function(d) {
  d[, locus := NA_character_]
  for (i in seq_len(nrow(WIN))) d[chr == WIN$chr[i] & pos >= WIN$s[i] & pos <= WIN$e[i], locus := WIN$locus[i]]
  d[!is.na(locus)]
}

# ---------------- read (GRCh37 throughout; build checked below) ----------------
read_mvp <- function(tr, n_scalar) {
  f <- file.path(SS, sprintf("MVP_%s_EUR_reformatted_hg19.tsv", tr))
  d <- fread(cmd = sprintf("awk -F'\\t' 'NR==1 || %s' %s", in_win_awk(1, 2), f))
  setnames(d, c("chromosome", "position", "allele1", "allele2", "pval"), c("chr", "pos", "ea", "oa", "p"))
  d[, `:=`(eaf = NA_real_, n = as.numeric(n_scalar))]
  add_locus(d[!is.na(beta) & !is.na(se) & se > 0])
}
read_pq <- function(f) {
  d <- fread(cmd = sprintf("zcat %s | awk -F'\\t' 'NR==1 || %s'", f, in_win_awk(1, 2)))
  setnames(d, c("chromosome", "base_pair_location", "effect_allele", "other_allele", "standard_error",
                "effect_allele_frequency", "p_value"), c("chr", "pos", "ea", "oa", "se", "eaf", "p"))
  if ("rs_id" %in% names(d)) setnames(d, "rs_id", "rsid")
  stopifnot("n" %in% names(d))
  d[, `:=`(chr = as.integer(chr), n = as.numeric(n))]
  add_locus(d[!is.na(beta) & !is.na(se) & se > 0 & !is.na(n) & n > 0])
}
read_panukbb <- function() {
  d <- rbindlist(lapply(seq_len(nrow(WIN)), function(i) fread(
    cmd = sprintf("%s %s %d:%d-%d | cut -f1-4,19,25,31,37,43", TABIX, PANUKBB, WIN$chr[i], WIN$s[i], WIN$e[i]),
    header = FALSE, sep = "\t", na.strings = "NA",
    col.names = c("chr", "pos", "oa", "ea", "eaf", "beta", "se", "nlp", "lowconf"))))
  d <- d[lowconf %in% c("false", FALSE) & !is.na(beta) & !is.na(se) & se > 0]
  d[, `:=`(p = 10^(-nlp), lowconf = NULL, nlp = NULL, n = as.numeric(N[["ALT_PANUKBB"]]))]
  add_locus(d)
}

D <- list(ALT_MVP = read_mvp("ALT", N[["ALT_MVP"]]),
          ALB_MVP = read_mvp("Albumin", N[["ALB_MVP"]]),
          PLT_MVP = read_mvp("Platelet", N[["PLT_MVP"]]),
          ALT_PANUKBB = read_panukbb(),
          TRAIL_UKB = read_pq(PQ["UKB"]),
          TRAIL_ZHAO = read_pq(PQ["ZHAO"]))
for (k in names(D)) D[[k]][, `:=`(ea = toupper(ea), oa = toupper(oa))]
print(rbindlist(lapply(names(D), function(k) D[[k]][, .(N_rows = .N, min_n = min(n), med_n = median(n), max_n = max(n)), by = locus][, trait := k])))

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
               named_allele = v[4], z = if (nrow(x)) x$beta[1] / x$se[1] * ifelse(x$ea[1] == v[4], 1, -1) else NA_real_,
               n_at_variant = if (nrow(x)) x$n[1] else NA_real_)
  }))
}))
fwrite(chk, file.path(OUT, "control_C1_build_allele_checks.tsv"), sep = "\t")
print(chk)
c1 <- chk[variant == "rs738409" & trait %in% c("ALT_MVP", "ALT_PANUKBB")]
stopifnot(nrow(c1) == 2, all(c1$n_found == 1), all(c1$beta_for_named_allele > 0))
stopifnot(chk[trait == "ALT_MVP" & variant == "TNFSF10_lead"]$beta_for_named_allele > 0)
build_ok <- chk[variant == "rs738409", .(trait, ok = n_found == 1 & (is.na(rsid_in_file) | rsid_in_file == "rs738409"))]
print(build_ok); stopifnot(all(build_ok$ok))

# ---------------- harmonise and pairwise coloc.abf ----------------
pal <- function(a, b) (a %in% c("A", "T") & b %in% c("A", "T")) | (a %in% c("C", "G") & b %in% c("C", "G"))
harm <- function(t1, t2, loc) {
  a <- D[[t1]][locus == loc & nchar(ea) == 1 & nchar(oa) == 1]
  b <- D[[t2]][locus == loc & nchar(ea) == 1 & nchar(oa) == 1]
  a <- a[!(duplicated(pos) | duplicated(pos, fromLast = TRUE))]
  b <- b[!(duplicated(pos) | duplicated(pos, fromLast = TRUE))]
  m <- merge(a[, .(pos, ea1 = ea, oa1 = oa, b1 = beta, se1 = se, p1 = p, f1 = eaf, n1 = n)],
             b[, .(pos, ea2 = ea, oa2 = oa, b2 = beta, se2 = se, p2 = p, f2 = eaf, n2 = n)], by = "pos")
  m <- m[(ea1 == ea2 & oa1 == oa2) | (ea1 == oa2 & oa1 == ea2)]
  m[ea2 != ea1, `:=`(b2 = -b2, f2 = 1 - f2)]
  m <- m[!pal(ea1, oa1)][order(pos)]
  # Ensure MAF is strictly in (0, 1) when finite (required by coloc::check_dataset)
  if (all(is.finite(m$f1))) m <- m[f1 > 0 & f1 < 1]
  if (all(is.finite(m$f2))) m <- m[f2 > 0 & f2 < 1]
  m
}
mk_ds <- function(beta, se, f, n_vec, snp) {
  ds <- list(beta = beta, varbeta = se^2, N = n_vec, type = "quant", snp = snp)
  if (all(is.finite(f))) ds$MAF <- pmin(f, 1 - f) else ds$sdY <- 1
  ds
}
is_int <- function(t) !grepl("^TRAIL", t)

run_abf <- function(t1, t2, loc, label) {
  m <- harm(t1, t2, loc)
  if (nrow(m) < MIN_SNPS) return(list(row = data.table(test = label, locus = loc, trait1 = t1, trait2 = t2,
                                                       n_snps = nrow(m), note = "too few SNPs")))
  snp <- paste0(WIN[locus == loc]$chr, ":", m$pos, ":", m$oa1, ":", m$ea1)
  d1 <- if (is_int(t1)) list(beta = m$b1, varbeta = m$se1^2, N = m$n1, type = "quant", sdY = 1, snp = snp) else
    mk_ds(m$b1, m$se1, m$f1, m$n1, snp)
  d2 <- if (is_int(t2)) list(beta = m$b2, varbeta = m$se2^2, N = m$n2, type = "quant", sdY = 1, snp = snp) else
    mk_ds(m$b2, m$se2, m$f2, m$n2, snp)
  sdy1 <- if ("sdY" %in% names(d1)) d1$sdY else suppressWarnings(coloc:::sdY.est(d1$varbeta, d1$MAF, d1$N))
  sdy2 <- if ("sdY" %in% names(d2)) d2$sdY else suppressWarnings(coloc:::sdY.est(d2$varbeta, d2$MAF, d2$N))
  a <- suppressMessages(suppressWarnings(coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = P12)))
  b <- suppressMessages(suppressWarnings(coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = P12_ALT)))
  ti <- which.max(a$results$SNP.PP.H4)
  i1 <- which.min(m$p1); i2 <- which.min(m$p2)
  pp <- a$summary; pb <- b$summary
  row <- data.table(test = label, locus = loc, trait1 = t1, trait2 = t2, n_snps = nrow(m),
                    n1_min = min(m$n1), n1_med = median(m$n1), n1_max = max(m$n1),
                    n2_min = min(m$n2), n2_med = median(m$n2), n2_max = max(m$n2),
                    sdY1 = sdy1, sdY2 = sdy2,
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

# ---------------- exact 4-trait single-causal-variant Wakefield ABF ----------------
logsum_pair <- function(a, b) {
  if (!is.finite(a)) return(b)
  if (!is.finite(b)) return(a)
  m <- max(a, b)
  m + log1p(exp(-abs(a - b)))
}

# Enumerate all set partitions of {1..m}, returning list of integer cluster-label vectors
# where 0 = inactive (unassigned) and 1..K = cluster IDs.
enum_subset_partitions <- function(m) {
  out <- list()
  # First choose active subset via bitmask 0..(2^m - 1)
  rec_part <- function(idx_vec, cur, next_lbl) {
    if (length(idx_vec) == 0L) { out[[length(out) + 1L]] <<- cur; return(invisible(NULL)) }
    i <- idx_vec[1L]; rest <- idx_vec[-1L]
    for (c_lbl in seq_len(next_lbl - 1L)) {
      nxt <- cur; nxt[i] <- c_lbl; rec_part(rest, nxt, next_lbl)
    }
    nxt <- cur; nxt[i] <- next_lbl; rec_part(rest, nxt, next_lbl + 1L)
  }
  for (mask in 0:(2^m - 1L)) {
    active <- which(bitwAnd(mask, 2^(0:(m - 1L))) > 0L)
    cur <- integer(m)
    rec_part(active, cur, 1L)
  }
  out
}

# Exact log(sum_{i_1..i_K distinct} prod_{r=1}^K exp(L[i_r, r])) via cancellation-free DP
exact_distinct_logsum <- function(L_mat) {
  Q <- nrow(L_mat); K <- ncol(L_mat)
  if (K == 0L) return(0)
  if (K == 1L) return(coloc:::logsum(L_mat[, 1L]))
  n_states <- 2^K
  dp <- rep(-Inf, n_states)
  dp[1L] <- 0  # mask 0 -> index 1
  masks_desc <- rev(seq_len(n_states - 1L))
  # Precompute set bits for each mask 1..(2^K-1)
  bits_list <- lapply(seq_len(n_states - 1L), function(mk) which(bitwAnd(mk, 2^(0:(K - 1L))) > 0L))
  for (i in seq_len(Q)) {
    row_vals <- L_mat[i, ]
    for (mk in masks_desc) {
      val <- dp[mk + 1L]
      for (r in bits_list[[mk]]) {
        prev <- dp[mk - bitwShiftL(1L, r - 1L) + 1L]
        if (is.finite(prev)) val <- logsum_pair(val, prev + row_vals[r])
      }
      dp[mk + 1L] <- val
    }
  }
  dp[n_states]
}

run_multitrait_abf <- function(traits = c("ALT_MVP", "TRAIL_UKB", "PLT_MVP", "ALB_MVP"), loc = "TNFSF10") {
  stopifnot(length(traits) == 4L)
  tabs <- lapply(traits, function(tr) {
    d <- D[[tr]][locus == loc & nchar(ea) == 1 & nchar(oa) == 1]
    d <- d[!(duplicated(pos) | duplicated(pos, fromLast = TRUE))]
    d[, .(pos, ea, oa, beta, se, p, eaf, n)]
  })
  names(tabs) <- traits
  g <- tabs[[1L]][, .(pos, ea1 = ea, oa1 = oa, b_1 = beta, se_1 = se, p_1 = p, f_1 = eaf, n_1 = n)]
  for (k in 2:4) {
    tk <- tabs[[k]][, .(pos, ea_k = ea, oa_k = oa, b_k = beta, se_k = se, p_k = p, f_k = eaf, n_k = n)]
    g <- merge(g, tk, by = "pos")
    g <- g[(ea1 == ea_k & oa1 == oa_k) | (ea1 == oa_k & oa1 == ea_k)]
    flip <- g$ea_k != g$ea1
    g$b_k[flip] <- -g$b_k[flip]
    g$f_k[flip] <- 1 - g$f_k[flip]
    setnames(g, c("b_k", "se_k", "p_k", "f_k", "n_k"),
             paste0(c("b_", "se_", "p_", "f_", "n_"), k))
    g[, c("ea_k", "oa_k") := NULL]
  }
  g <- g[!pal(ea1, oa1)][order(pos)]
  # Check MAF bounds for TRAIL_UKB (trait 2)
  f2 <- g$f_2
  if (any(is.finite(f2))) g <- g[f_2 > 0 & f_2 < 1]
  Q <- nrow(g)
  snp <- paste0(WIN[locus == loc]$chr, ":", g$pos, ":", g$oa1, ":", g$ea1)

  # Compute per-trait lABF on the common grid
  lABF <- matrix(0, nrow = Q, ncol = 4L)
  colnames(lABF) <- traits
  sdY_vec <- numeric(4L); names(sdY_vec) <- traits
  for (k in 1:4) {
    tr <- traits[k]
    bk <- g[[paste0("b_", k)]]; sek <- g[[paste0("se_", k)]]
    fk <- g[[paste0("f_", k)]]; nk <- g[[paste0("n_", k)]]
    if (is_int(tr)) {
      sdy <- 1
    } else {
      sdy <- suppressWarnings(coloc:::sdY.est(sek^2, pmin(fk, 1 - fk), nk))
    }
    sdY_vec[k] <- sdy
    lABF[, k] <- coloc:::approx.bf.estimates(z = bk / sek, V = sek^2, type = "quant", sdY = sdy)$lABF
  }

  # Verify 2-trait DP vs coloc:::combine.abf on (ALT_MVP, TRAIL_UKB)
  p12_check <- suppressMessages(coloc:::combine.abf(lABF[, 1L], lABF[, 2L], P1, P2, P12, quiet = TRUE))
  dp2_h3 <- exact_distinct_logsum(lABF[, 1:2])
  dp2_h4 <- coloc:::logsum(lABF[, 1L] + lABF[, 2L])
  dp2_all <- c(0, log(P1) + coloc:::logsum(lABF[, 1L]), log(P2) + coloc:::logsum(lABF[, 2L]),
               log(P1) + log(P2) + dp2_h3, log(P12) + dp2_h4)
  dp2_pp <- exp(dp2_all - coloc:::logsum(dp2_all))
  max_dp_err <- max(abs(dp2_pp - p12_check))
  cat(sprintf("Multi-trait DP vs coloc:::combine.abf 2-trait sanity check max|diff| = %.3e\n", max_dp_err))
  stopifnot(max_dp_err < 1e-10)

  parts <- enum_subset_partitions(4L)
  stopifnot(length(parts) == 52L)  # Bell(5) = 52
  hyp_tabs <- rbindlist(lapply(seq_along(parts), function(idx) {
    lbl <- parts[[idx]]
    K <- max(lbl)
    if (K == 0L) {
      log_bf <- 0
      cl_sizes <- integer(0)
      label_str <- "NULL (no causal trait)"
    } else {
      L_cl <- sapply(seq_len(K), function(c_id) rowSums(lABF[, which(lbl == c_id), drop = FALSE]))
      if (!is.matrix(L_cl)) L_cl <- matrix(L_cl, nrow = Q, ncol = K)
      log_bf <- exact_distinct_logsum(L_cl)
      cl_sizes <- sapply(seq_len(K), function(c_id) sum(lbl == c_id))
      cl_strs <- sapply(seq_len(K), function(c_id) paste0("{", paste(traits[lbl == c_id], collapse = ","), "}"))
      label_str <- paste(cl_strs, collapse = " | ")
    }
    # Log-prior for a partition with K clusters of sizes cl_sizes:
    # sum_{r=1}^K [log(p1) + (m_r - 1)*log(gamma)] where gamma = conditional sharing probability
    log_prior <- function(p1, gamma) {
      if (K == 0L) return(0)
      sum(log(p1) + (cl_sizes - 1L) * log(gamma))
    }
    data.table(
      hyp_id = idx,
      n_clusters = K,
      n_active_traits = sum(lbl > 0L),
      cluster_pattern = paste(sort(cl_sizes, decreasing = TRUE), collapse = "+"),
      hypothesis = label_str,
      log_BF_sum = log_bf,
      lp_1e5 = log_prior(P1, P12 / P1),
      lp_5e6 = log_prior(P1, P12_ALT / P1),
      lp_hypr_002 = log_prior(P1, 0.02)
    )
  }))
  hyp_tabs[, `:=`(
    PP_p12_1e5 = exp((log_BF_sum + lp_1e5) - coloc:::logsum(log_BF_sum + lp_1e5)),
    PP_p12_5e6 = exp((log_BF_sum + lp_5e6) - coloc:::logsum(log_BF_sum + lp_5e6)),
    PP_hypr_002 = exp((log_BF_sum + lp_hypr_002) - coloc:::logsum(log_BF_sum + lp_hypr_002))
  )]
  hyp_tabs <- hyp_tabs[order(-PP_p12_1e5)]

  # SNP-level posterior on the 4-trait joint colocalization hypothesis {ALT, TRAIL, PLT, ALB}
  sum4 <- rowSums(lABF)
  snp_pp4 <- exp(sum4 - coloc:::logsum(sum4))
  ti4 <- which.max(snp_pp4)

  # Also marginal PP that {ALT_MVP, TRAIL_UKB, PLT_MVP} share the same causal variant
  # (either all 4 share, OR {ALT, TRAIL, PLT} share + ALB distinct/null)
  idx_all4 <- which(hyp_tabs$cluster_pattern == "4")
  idx_3_plt <- which(grepl("\\{ALT_MVP,TRAIL_UKB,PLT_MVP\\}", hyp_tabs$hypothesis))
  summary_dt <- data.table(
    locus = loc,
    traits = paste(traits, collapse = ","),
    n_common_snps = Q,
    sdY_TRAIL_UKB_est = sdY_vec[["TRAIL_UKB"]],
    top_4trait_snp = snp[ti4],
    top_4trait_snp_pp = snp_pp4[ti4],
    PP_4trait_joint_p12_1e5 = hyp_tabs[idx_all4, PP_p12_1e5],
    PP_4trait_joint_p12_5e6 = hyp_tabs[idx_all4, PP_p12_5e6],
    PP_4trait_joint_hypr002 = hyp_tabs[idx_all4, PP_hypr_002],
    PP_ALT_TRAIL_PLT_coloc_marginal_1e5 = sum(hyp_tabs[idx_3_plt, PP_p12_1e5]),
    PP_ALT_TRAIL_PLT_coloc_marginal_5e6 = sum(hyp_tabs[idx_3_plt, PP_p12_5e6]),
    PP_ALT_TRAIL_PLT_coloc_marginal_hypr002 = sum(hyp_tabs[idx_3_plt, PP_hypr_002]),
    top_hypothesis = hyp_tabs$hypothesis[1L],
    top_hypothesis_PP_1e5 = hyp_tabs$PP_p12_1e5[1L]
  )
  list(hypotheses = hyp_tabs, summary = summary_dt, grid = g[, .(pos, ea1, oa1)], snp = snp, snp_pp4 = snp_pp4)
}

# ---------------- PolyFun EUR LD (N=337,000, PSD-projected) ----------------
POLYFUN_BLOCKS <- fread(file.path(POLYFUN_DIR, "approx_LD_blocks.txt"))
LDC <- new.env(parent = emptyenv())

load_polyfun_window <- function(loc) {
  if (!is.null(LDC[[loc]])) return(LDC[[loc]])
  w <- WIN[locus == loc]
  m <- POLYFUN_BLOCKS[chr == w$chr & start < w$e & stop > w$s]
  stopifnot(nrow(m) >= 1L)
  bim_list <- list(); mats <- list()
  for (i in seq_len(nrow(m))) {
    bdir <- file.path(POLYFUN_DIR, paste0("chr", w$chr), paste0(m$start[i], ".", m$stop[i]))
    bf <- file.path(bdir, paste0(m$start[i], ".", m$stop[i], ".bim"))
    lf <- file.path(bdir, paste0(m$start[i], ".", m$stop[i], ".ld"))
    stopifnot(file.exists(bf), file.exists(lf))
    b <- fread(bf, header = FALSE, col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
    b[, `:=`(a1 = toupper(a1), a2 = toupper(a2), block_id = paste0(m$start[i], ".", m$stop[i]))]
    bim_list[[i]] <- b
    mats[[i]] <- as.matrix(fread(lf, header = FALSE))
  }
  bim <- rbindlist(bim_list)
  if (length(mats) == 1L) {
    ld <- mats[[1L]]
  } else {
    tot <- sum(vapply(mats, nrow, integer(1)))
    ld <- matrix(0, tot, tot); off <- 0L
    for (mat in mats) { ix <- off + seq_len(nrow(mat)); ld[ix, ix] <- mat; off <- off + nrow(mat) }
  }
  rm(mats)
  bim[, SNP := paste(chr, pos, a1, a2, sep = ":")]
  keep <- !duplicated(bim$SNP)
  if (!all(keep)) { ld <- ld[keep, keep]; bim <- bim[keep] }
  bim[, k := .I]
  bim[, pair_key := paste(pos, pmin(a1, a2), pmax(a1, a2), sep = ":")]
  res <- list(bim = bim, ld = ld)
  LDC[[loc]] <- res
  res
}

psd_project_corr <- function(R) {
  R <- 0.5 * (R + t(R))
  eig <- eigen(R, symmetric = TRUE)
  pre_min_eig <- min(eig$values)
  n_neg <- sum(eig$values < 0)
  d <- pmax(eig$values, 1e-4)
  R_psd <- eig$vectors %*% (d * t(eig$vectors))
  inv_sd <- 1 / sqrt(diag(R_psd))
  R_out <- (R_psd * inv_sd) * rep(inv_sd, each = nrow(R_psd))
  R_out <- 0.5 * (R_out + t(R_out))
  diag(R_out) <- 1.0
  post_min_eig <- min(eigen(R_out, symmetric = TRUE, only.values = TRUE)$values)
  list(R = R_out, pre_min_eig = pre_min_eig, n_eig_neg = n_neg, post_min_eig = post_min_eig)
}

run_susie <- function(r, loc) {
  L <- load_polyfun_window(loc); bim <- copy(L$bim)
  n_bim_tab <- bim[, .(n_bim = .N), by = pair_key]
  bim <- merge(bim, n_bim_tab, by = "pair_key", sort = FALSE)[n_bim == 1L]
  m <- copy(r$m)
  m[, pair_key := paste(pos, pmin(ea1, oa1), pmax(ea1, oa1), sep = ":")]
  x <- merge(m[, .(pos, pair_key, ea1, oa1, z1 = b1 / se1, z2 = b2 / se2, n1, n2)],
             bim[, .(pair_key, a1, a2, k)], by = "pair_key")
  x <- x[(ea1 == a1 & oa1 == a2) | (ea1 == a2 & oa1 == a1)]
  x[ea1 == a2, `:=`(z1 = -z1, z2 = -z2)]
  x <- x[order(k)]
  if (nrow(x) < MIN_SNPS) return(NULL)

  R <- as.matrix(L$ld[x$k, x$k, drop = FALSE])
  if (sum(is.na(R)) != 0) {
    for (iter in seq_len(100)) {
      na_per_var <- colSums(is.na(R))
      if (max(na_per_var) == 0) break
      bad <- na_per_var > (nrow(R) * 0.5)
      if (sum(bad) == 0) bad <- na_per_var > 0
      R <- R[!bad, !bad, drop = FALSE]; x <- x[!bad]
      if (nrow(R) < MIN_SNPS) break
    }
    if (nrow(R) < MIN_SNPS) return(NULL)
  }

  proj <- psd_project_corr(R)
  R <- proj$R
  snp <- paste0(WIN[locus == loc]$chr, ":", x$pos, ":", x$a2, ":", x$a1)
  dimnames(R) <- list(snp, snp)

  n1_use <- round(median(x$n1))
  n2_use <- round(median(x$n2))
  fit <- function(z, n_val) {
    lam <- tryCatch(estimate_s_rss(z, R, n = n_val), error = function(e) NA_real_)
    set.seed(42)
    s <- tryCatch(susie_rss(z = setNames(z, snp), R = R, n = n_val, L = 10, coverage = 0.95,
                            estimate_residual_variance = FALSE, max_iter = 5000),
                  error = function(e) NULL)
    list(s = s, lambda = lam, converged = if (is.null(s)) FALSE else isTRUE(s$converged),
         ncs = if (is.null(s) || is.null(s$sets$cs)) 0L else length(s$sets$cs))
  }
  f1 <- fit(x$z1, n1_use); f2 <- fit(x$z2, n2_use)
  summ <- data.table(test = r$row$test, ld_panel = "polyfun_eur", n_ld_panel = 337000L,
                     n_snps_ld = nrow(x), pre_min_eig = proj$pre_min_eig, n_eig_neg = proj$n_eig_neg,
                     post_min_eig = proj$post_min_eig, n1_susie = n1_use, n2_susie = n2_use,
                     lambda1 = f1$lambda, lambda2 = f2$lambda,
                     converged1 = f1$converged, converged2 = f2$converged,
                     ncs1 = f1$ncs, ncs2 = f2$ncs)
  cs_tab <- rbindlist(lapply(list(list("trait1", f1), list("trait2", f2)), function(o) {
    s <- o[[2]]$s; if (is.null(s) || o[[2]]$ncs == 0) return(NULL)
    rbindlist(lapply(seq_along(s$sets$cs), function(j) {
      v <- s$sets$cs[[j]]; top <- v[which.max(s$pip[v])]
      data.table(test = r$row$test, which = o[[1]], cs = names(s$sets$cs)[j], size = length(v),
                 top_snp = snp[top], top_pip = s$pip[top],
                 top_z = c(x$z1, x$z2)[top + if (o[[1]] == "trait2") nrow(x) else 0])
    }))
  }))
  cs <- NULL
  if (f1$ncs >= 1 && f2$ncs >= 1) {
    for (pr in c(P12, P12_ALT)) {
      set.seed(42)
      cc <- tryCatch(suppressMessages(coloc.susie(f1$s, f2$s, p1 = P1, p2 = P2, p12 = pr)), error = function(e) NULL)
      if (!is.null(cc) && !is.null(cc$summary))
        cs <- rbind(cs, data.table(test = r$row$test, p12 = pr, as.data.table(cc$summary)), fill = TRUE)
    }
  }
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

# ---------------- exact 4-trait Wakefield ABF joint posterior ----------------
cat("\n== Multi-trait 4-trait Wakefield ABF (ALT_MVP, TRAIL_UKB, PLT_MVP, ALB_MVP) at TNFSF10\n")
mt <- run_multitrait_abf()
fwrite(mt$hypotheses, file.path(OUT, "multitrait_abf_hypotheses.tsv"), sep = "\t")
fwrite(mt$summary, file.path(OUT, "multitrait_abf_summary.tsv"), sep = "\t")
print(mt$summary)
print(mt$hypotheses[1:10, .(hyp_id, n_clusters, cluster_pattern, hypothesis, PP_p12_1e5, PP_p12_5e6, PP_hypr_002)])

# ---------------- direction at TNFSF10 ----------------
dir_rows <- rbindlist(lapply(c("ALT_MVP", "ALT_PANUKBB", "ALB_MVP", "PLT_MVP", "TRAIL_UKB", "TRAIL_ZHAO"), function(k) {
  x <- D[[k]][chr == 3 & pos == 172294500 & ((ea == "G" & oa == "A") | (ea == "A" & oa == "G"))]
  if (!nrow(x)) return(data.table(trait = k, variant = "3:172294500", found = FALSE))
  s <- ifelse(x$ea[1] == "G", 1, -1)
  data.table(trait = k, variant = "3:172294500", found = TRUE, file_effect_allele = x$ea[1], file_other_allele = x$oa[1],
             file_beta = x$beta[1], beta_per_G = s * x$beta[1], se = x$se[1], z_per_G = s * x$beta[1] / x$se[1],
             p = x$p[1], eaf_G = if (is.finite(x$eaf[1])) ifelse(s > 0, x$eaf[1], 1 - x$eaf[1]) else NA_real_,
             n = x$n[1])
}), fill = TRUE)
for (k in c("TRAIL_UKB", "TRAIL_ZHAO")) {
  x <- D[[k]][locus == "TNFSF10"][which.min(p)]
  dir_rows <- rbind(dir_rows, data.table(trait = k, variant = paste0("pQTL lead 3:", x$pos), found = TRUE,
                                         file_effect_allele = x$ea, file_other_allele = x$oa, file_beta = x$beta,
                                         se = x$se, p = x$p, eaf_G = NA_real_, n = x$n), fill = TRUE)
}
fwrite(dir_rows, file.path(OUT, "direction_TNFSF10.tsv"), sep = "\t")
print(dir_rows)

# ---------------- Epitope / coding-exon check at TNFSF10 ----------------
cat("\n== Epitope / TNFSF10 coding exon check for 3:172294500 (rs79287178) and r2 >= 0.8 proxies\n")
gtf_tnf <- fread(cmd = sprintf("awk -F'\\t' '$1==\"chr3\" && $4>=172200000 && $5<=172300000' %s | grep 'gene_name \"TNFSF10\"'", GTF_PATH),
                 header = FALSE, sep = "\t")
setnames(gtf_tnf, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attr"))
gene_box <- gtf_tnf[feature == "gene", .(start = min(start), end = max(end))]
cds_box <- unique(gtf_tnf[feature == "CDS", .(start, end)])
exon_box <- unique(gtf_tnf[feature == "exon", .(start, end)])
cat("TNFSF10 gene (hg19):", gene_box$start, "-", gene_box$end,
    "| CDS intervals:", nrow(cds_box), "| Exon intervals:", nrow(exon_box), "\n")

in_intervals <- function(pos_vec, iv) {
  sapply(pos_vec, function(p) any(p >= iv$start & p <= iv$end))
}
dist_intervals <- function(pos_vec, iv) {
  sapply(pos_vec, function(p) {
    d <- pmax(0L, pmax(iv$start - p, p - iv$end))
    min(d)
  })
}

# PolyFun LD with 3:172294500
L_tnf <- load_polyfun_window("TNFSF10")
pf_bim <- copy(L_tnf$bim)
lead_k <- pf_bim[pos == 172294500]$k
stopifnot(length(lead_k) == 1L)
pf_bim[, r_polyfun := as.numeric(L_tnf$ld[lead_k, k])]
pf_bim[, r2_polyfun := r_polyfun^2]

# 1kg_eur LD with 3:172294500 (rs79287178) via PLINK
onekg_bim <- fread(sprintf("%s/chr3_eur.bim", ONEKG_DIR), header = FALSE,
                   col.names = c("chr", "id", "cm", "pos", "a1", "a2"))
onekg_lead_id <- onekg_bim[pos == 172294500, id]
stopifnot(length(onekg_lead_id) == 1L)
ld_tmp <- file.path(OUT, "tmp_1kg_tnfsf10")
w_tnf <- WIN[locus == "TNFSF10"]
rc <- system(sprintf("plink --bfile %s/chr3_eur --chr 3 --from-bp %d --to-bp %d --r2 --ld-snp %s --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0.0 --out %s > /dev/null 2>&1",
                     ONEKG_DIR, w_tnf$s, w_tnf$e, onekg_lead_id, ld_tmp))
stopifnot(rc == 0L, file.exists(paste0(ld_tmp, ".ld")))
onekg_ld <- fread(paste0(ld_tmp, ".ld"))
unlink(paste0(ld_tmp, c(".ld", ".log", ".nosex")))
onekg_r2_map <- onekg_ld[, .(pos = BP_B, r2_1kg = R2)]

# Combine all variants from PolyFun and 1kg_eur in the TNFSF10 window
ep <- merge(pf_bim[, .(pos, pf_rsid = rsid, pf_a1 = a1, pf_a2 = a2, r_polyfun, r2_polyfun)],
            onekg_r2_map, by = "pos", all = TRUE)
# Annotate TRAIL_UKB and ALT_MVP stats
tr_ukb <- D[["TRAIL_UKB"]][locus == "TNFSF10", .(pos, trail_ukb_p = p, trail_ukb_beta = beta, trail_ukb_ea = ea, trail_ukb_oa = oa, trail_ukb_eaf = eaf)][!duplicated(pos)]
alt_m <- D[["ALT_MVP"]][locus == "TNFSF10", .(pos, alt_mvp_p = p, alt_mvp_beta = beta)][!duplicated(pos)]
ep <- merge(ep, tr_ukb, by = "pos", all.x = TRUE)
ep <- merge(ep, alt_m, by = "pos", all.x = TRUE)
ep[, `:=`(
  chr = 3L,
  is_lead = (pos == 172294500L),
  r2_max = pmax(r2_polyfun, r2_1kg, na.rm = TRUE),
  proxy_r2_ge_0.8 = pmax(r2_polyfun, r2_1kg, na.rm = TRUE) >= 0.8,
  in_tnfsf10_gene = (pos >= gene_box$start & pos <= gene_box$end),
  in_tnfsf10_exon = in_intervals(pos, exon_box),
  in_tnfsf10_cds = in_intervals(pos, cds_box),
  in_tnfsf10_cds_or_splice2bp = in_intervals(pos, data.table(start = cds_box$start - 2L, end = cds_box$end + 2L)),
  dist_to_nearest_cds_bp = dist_intervals(pos, cds_box),
  dist_to_tnfsf10_tss_bp = abs(pos - gene_box$end) # minus strand gene -> TSS is at max coord
)]
ep[, epitope_artifact_flag := proxy_r2_ge_0.8 & in_tnfsf10_cds_or_splice2bp]

# Also check if ANY variant inside TNFSF10 CDS in TRAIL_UKB or ALT_MVP exists outside LD panels
cds_sumstats <- D[["TRAIL_UKB"]][locus == "TNFSF10"][in_intervals(pos, data.table(start = cds_box$start - 2L, end = cds_box$end + 2L))]
cat("TRAIL_UKB variants inside TNFSF10 CDS (+/-2bp):", nrow(cds_sumstats), "\n")
if (nrow(cds_sumstats) > 0) print(cds_sumstats[, .(chr, pos, ea, oa, eaf, beta, se, p, rsid)])

# Keep: lead + all r2>=0.8 + top 15 by r2_max + all variants in TNFSF10 exon/CDS
ep_out <- ep[is_lead | proxy_r2_ge_0.8 | in_tnfsf10_exon | in_tnfsf10_cds | r2_max >= 0.2][order(-r2_max, pos)]
fwrite(ep_out, file.path(OUT, "epitope_check.tsv"), sep = "\t")
cat("Epitope check summary:\n")
cat("  Variants with r2 >= 0.8 (PolyFun or 1KG):", sum(ep$proxy_r2_ge_0.8, na.rm = TRUE), "\n")
cat("  Of those, in TNFSF10 exon/CDS (+/-2bp):", sum(ep$epitope_artifact_flag, na.rm = TRUE), "\n")
cat("  Max r2 with rs79287178 among all TNFSF10 exon variants:",
    if (any(ep$in_tnfsf10_exon, na.rm = TRUE)) max(ep[in_tnfsf10_exon == TRUE]$r2_max, na.rm = TRUE) else NA_real_, "\n")
cat("  Max r2 with rs79287178 among all TNFSF10 CDS variants:",
    if (any(ep$in_tnfsf10_cds, na.rm = TRUE)) max(ep[in_tnfsf10_cds == TRUE]$r2_max, na.rm = TRUE) else NA_real_, "\n")
print(ep_out[1:min(15L, nrow(ep_out)), .(pos, pf_rsid, r2_polyfun = round(r2_polyfun, 4),
                                          r2_1kg = round(r2_1kg, 4), in_tnfsf10_exon, in_tnfsf10_cds,
                                          dist_to_nearest_cds_bp, epitope_artifact_flag, trail_ukb_p, alt_mvp_p)])

cat("\nPANX1 plasma pQTL: untestable (no open plasma PANX1 pQTL source; see PRESPEC).\n")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("done", date(), "\n")
