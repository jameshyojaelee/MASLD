#!/usr/bin/env Rscript
# H53 job 2 (prespec "Compute", job 2): caQTL lbf, coloc.bf_bf for all leads and peaks,
# Currin route, NC1, NC2, three-way, evidence states, T1, descriptive outputs,
# RESULT_H53.md with the verdict computed by the prespec "Reading" rule.
# Usage: Rscript h53_job2.R [--null-debug]   (--null-debug: chr21-22 only, outcome
#        code exercised on circularly shifted caQTL values; no real outcome written)
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H53-20261015T112655Z/scripts/h53_common.R")
args <- commandArgs(trailingOnly = TRUE)
NULL_DEBUG <- "--null-debug" %in% args
NCORE <- as.integer(Sys.getenv("H53_WORKERS", "6"))
DIAG <- file.path(RUN, "diagnostics")
RES  <- if (NULL_DEBUG) file.path(RUN, "debug_null") else file.path(RUN, "results")
dir.create(RES, showWarnings = FALSE)
CHRS <- if (NULL_DEBUG) c("21", "22") else as.character(1:22)
stamp("H53 job 2 start; SLURM_JOB_ID =", Sys.getenv("SLURM_JOB_ID"), "; NULL_DEBUG =", NULL_DEBUG)
cat("Config: WIN", WIN, "MAF_MIN", MAF_MIN, "W", W_PRIOR, "p1", P1, "p2", P2, "p12", P12, "(alt", P12_ALT, ") MIN_SNPS", MIN_SNPS,
    "overlap.min", OVERLAP_MIN, "entry PP.H4 >", PP_ENTRY, "and LD r2 >=", LD_ENTRY, "seed", SEED, "\n")

# ---- 0. checks -------------------------------------------------------------------
verify_inputs(cores = 16)
sel <- readLines(file.path(DIAG, "selection_sha256.txt"))
for (l in sel) {
  h <- sub(" .*", "", l); f <- sub("^[0-9a-f]+  ", "", l)
  if (sha256(file.path(RUN, f), 1) != h) stop("selection table changed since job 1: ", f)
}
gate <- jsonlite::read_json(file.path(DIAG, "job1_gate.json"))
if (!isTRUE(gate$S1_pass) || !isTRUE(gate$S2_pass)) stop("job 1 gate did not pass (S1/S2); job 2 must not run")
cat("CHECKPOINT input-hash, selection-hash and job-1 gate checks PASSED\n")

SL  <- fread(file.path(RUN, "results/study_leads_job1.tsv"), colClasses = list(character = c("chr")))
LOC <- fread(file.path(RUN, "results/loci_job1.tsv"), colClasses = list(character = c("chr")))
MVP <- fread(file.path(RUN, "results/selection_mvp_alt_loci.tsv"), colClasses = list(character = c("chr")))
a <- load_arms()
SL[, lead_key := paste(gwas, hit1, sep = "|")]
LOC[, lead_key := paste(gwas, hit1, sep = "|")]
TL <- SL[reproduced == TRUE & chr %in% CHRS]
reg <- fread(file.path(FM, "config/gwas_registry.tsv"))

read_caqtl_chr <- function(ch) {
  f <- list.files(CQDIR, pattern = paste0("^caQTL_hg38_chr", ch, "_.*txt.gz$"), full.names = TRUE)
  q <- fread(cmd = paste("zcat", shQuote(f)), sep = "\t", showProgress = FALSE)
  setkey(q, peak)
  q
}
peak_centres <- function(q) {
  pc <- q[, .(variant = variant[1], dist = distance_from_peakCenter[1], min_p = min(pvalue), n_rows = .N), by = peak]
  pc[, centre := as.integer(tstrsplit(variant, ":", fixed = TRUE)[[2]]) - dist]
  pc[]
}
load_gwas_chr <- function(st, ch) {   # pilot-coded Currin-route GWAS input
  f <- file.path(FM, reg[study_name == st, sumstats_path])
  g <- fread(cmd = paste0("awk -F'\\t' 'NR==1 || $1==\"", ch, "\"' ", shQuote(f)), showProgress = FALSE)
  g <- g[nchar(allele1) == 1 & nchar(allele2) == 1 & se > 0]
  g[, pos38 := lift19to38(ch, position)]; g <- g[!is.na(pos38)]
  g[, akey := paste(pos38, pmin(allele1, allele2), pmax(allele1, allele2), sep = ":")]
  g <- g[order(pval)][!duplicated(akey)]
  setkey(g, pos38); g[]
}
currin_route <- function(g, qraw, st, hit1, lead_snp, fit, gk) {   # coloc.abf as Currin, descriptive
  q <- copy(qraw); sp <- tstrsplit(q$variant, ":", fixed = TRUE)
  q[, `:=`(pos38 = as.integer(sp[[2]]), ref = sp[[3]], alt = sp[[4]])]
  q <- q[MAF >= 0.01 & nchar(ref) == 1 & nchar(alt) == 1 & is.finite(varbeta) & varbeta > 0]
  if (!nrow(q)) return(list(cr_n = 0L))
  q[, akey := paste(pos38, pmin(ref, alt), pmax(ref, alt), sep = ":")]
  lead <- q[which.min(pvalue)]
  rq <- q[abs(pos38 - lead$pos38) <= WIN]; rg <- g[abs(pos38 - lead$pos38) <= WIN]
  m <- merge(rg[, .(akey, gb = beta, gv = se^2)], rq[, .(akey, qb = beta, qv = varbeta, MAF)], by = "akey")
  m <- m[is.finite(gv) & is.finite(qv) & qv > 0 & gv > 0]
  if (nrow(m) < 50) return(list(cr_n = nrow(m)))
  rg1 <- reg[study_name == st]
  d1 <- if (rg1$trait_type == "binary") list(beta = m$gb, varbeta = m$gv, snp = m$akey, type = "cc", s = rg1$N_cases / rg1$N_tot, N = rg1$N_tot)
        else list(beta = m$gb, varbeta = m$gv, snp = m$akey, type = "quant", sdY = 1, N = rg1$N_tot)
  d2 <- list(beta = m$qb, varbeta = m$qv, snp = m$akey, type = "quant", sdY = 1, N = N_CAQTL, MAF = m$MAF)
  r <- quiet(coloc.abf(d1, d2, p12 = P12_ALT))   # coloc default 1e-5, as PC1 and the pilot
  r6 <- quiet(coloc.abf(d1, d2, p12 = P12))
  k <- paste(lead$pos38, min(lead$ref, lead$alt), max(lead$ref, lead$alt), sep = ":")
  lsnp <- gk$snp_id[match(k, gk$key)]
  rsnp <- if (hit1 %in% rownames(fit$ld_rows)) hit1 else lead_snp
  ld <- if (!is.na(lsnp) && rsnp %in% rownames(fit$ld_rows) && lsnp %in% colnames(fit$ld_rows)) fit$ld_rows[rsnp, lsnp]^2 else NA_real_
  out <- list(cr_n = nrow(m), cr_PP.H3 = r$summary[["PP.H3.abf"]], cr_PP.H4 = r$summary[["PP.H4.abf"]], cr_PP.H4_p12_5e6 = r6$summary[["PP.H4.abf"]],
              cr_caqtl_lead = lead$variant, cr_lead_ld_r2_hit1 = ld)
  if (out$cr_PP.H4 > PP_ENTRY && isTRUE(ld >= LD_ENTRY)) {
    s <- quiet(sensitivity(r, rule = "H4 > 0.8", doplot = FALSE, plot.manhattans = FALSE))
    w <- which(s$pass)
    out$cr_sens_p12_min_pass <- if (length(w)) min(s$p12[w]) else NA_real_
    out$cr_sens_p12_max_pass <- if (length(w)) max(s$p12[w]) else NA_real_
  }
  out
}

# ---- 1. pass 1: primary route + Currin route for every reproduced study-lead ------
t_all <- Sys.time()
PAIRS <- list(); PEAKDIAG <- list()
for (ch in CHRS) {
  leads <- TL[chr == ch]
  if (!nrow(leads)) next
  t0 <- Sys.time()
  q <- read_caqtl_chr(ch); pcen <- peak_centres(q)
  GW <- lapply(setNames(unique(leads$gwas), unique(leads$gwas)), load_gwas_chr, ch = ch)
  stamp("chr", ch, ": caQTL rows", nrow(q), "peaks", nrow(pcen), "| leads", nrow(leads), "| read", format(Sys.time() - t0))
  res <- mclapply(seq_len(nrow(leads)), function(i) {
    r <- leads[i]; fit <- readRDS(r$fit_file); gk <- gwas_keys(fit); grow <- gwas_row(fit, gk, r$idx1)
    alpha1 <- fit$alpha[r$idx1, ]; names(alpha1) <- gk$key
    cand <- pcen[abs(centre - r$hit1_pos38) <= WIN]
    if (is.na(r$hit1_pos38) || !nrow(cand)) return(data.table(lead_key = r$lead_key, peak = NA_character_, n_cand = 0L))
    rbindlist(lapply(seq_len(nrow(cand)), function(j) {
      pk <- cand$peak[j]; qraw <- q[.(pk)]
      qp <- caqtl_prepare(copy(qraw))
      ncent <- uniqueN(qp$pos38 - qp$distance_from_peakCenter)
      o <- if (NULL_DEBUG) {   # outcome code exercised on shifted (null) caQTL values only
        sh <- intersect(gk$key, qp$key); n <- length(sh)
        lb <- qp$lbf; pv <- qp$pvalue
        if (n > 3) { ix <- match(sh, qp$key); k <- floor(n / 2); rot <- (seq_len(n) + k - 1) %% n + 1; lb[ix] <- qp$lbf[ix[rot]]; pv[ix] <- qp$pvalue[ix[rot]] }
        pair_test(grow, gk, qp, fit$ld_rows, r$gwas_lead_snp, qp_lbf = lb, qp_p = pv)
      } else pair_test(grow, gk, qp, fit$ld_rows, r$gwas_lead_snp)
      shared <- setdiff(intersect(gk$key, qp$key), "null")
      o$alpha_mass_shared <- sum(alpha1[shared])
      cr <- currin_route(GW[[r$gwas]], qraw, r$gwas, r$hit1, r$gwas_lead_snp, fit, gk)
      o <- c(o, cr); o <- o[!vapply(o, is.null, logical(1))]
      cbind(data.table(lead_key = r$lead_key, peak = pk, centre = cand$centre[j], n_cand = nrow(cand),
                       n_centres_in_peak = ncent, peak_min_p_all = cand$min_p[j]), as.data.table(o))
    }), fill = TRUE)
  }, mc.cores = NCORE, mc.preschedule = FALSE)
  err <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(err)) { print(res[err]); stop("worker error on chr ", ch) }
  PAIRS[[ch]] <- rbindlist(res, fill = TRUE)
  stamp("chr", ch, "done:", nrow(PAIRS[[ch]]), "(lead, peak) rows in", format(Sys.time() - t0))
  rm(q, GW, res); gc(FALSE)
}
PR <- rbindlist(PAIRS, fill = TRUE)
PR <- merge(PR, SL[, .(lead_key, arm, gwas, hit1, chr, refit_gene, refit_ensembl, idx1, idx2, locus, confirmatory, hit1_is_top_alpha)], by = "lead_key", all.x = TRUE)
fwrite(PR, file.path(RES, "lead_peak_pairs.tsv.gz"), sep = "\t")
stamp("pass 1 done in", format(Sys.time() - t_all))

# ---- 2. diagnostics printed before any outcome table ------------------------------
cat("\n== Diagnostics (before outcome tables) ==\n")
cat("peaks whose centre (position - distance) is not constant:", PR[!is.na(peak) & n_centres_in_peak > 1, uniqueN(peak)], "\n")
pl <- PR[!is.na(peak)]
cat("candidate peaks per lead (median, IQR) by arm:\n"); print(pl[, .(n = n_cand[1]), by = .(lead_key, arm)][, as.list(quantile(n, c(.25, .5, .75))), by = arm])
cat("testable peaks per lead (median, IQR) by arm:\n"); print(pl[, .(n = sum(testable)), by = .(lead_key, arm)][, as.list(quantile(n, c(.25, .5, .75))), by = arm])
cat("shared SNPs per (lead, peak) (median, IQR):\n"); print(quantile(pl$n_shared, c(.25, .5, .75)))
cat("prop_gwas / prop_caqtl among pairs with >=100 shared SNPs (median, IQR) by arm:\n")
print(pl[n_shared >= MIN_SNPS, .(prop_gwas_q = paste(round(quantile(prop_gwas, c(.25, .5, .75)), 3), collapse = "/"),
                                  prop_caqtl_q = paste(round(quantile(prop_caqtl, c(.25, .5, .75)), 3), collapse = "/")), by = arm])
cat("pairs dropped by coloc.bf_bf overlap rule, by arm and side:\n"); print(pl[, .N, by = .(arm, reason)][order(arm, reason)])
cat("idx1 alpha mass on shared SNPs (median, IQR) by arm (testable pairs):\n"); print(pl[testable == TRUE, as.list(quantile(alpha_mass_shared, c(.25, .5, .75))), by = arm])
cat("r2 hit1 differs from top-alpha SNP of idx1 (leads), by arm:\n"); print(SL[reproduced == TRUE, .(n = .N, differs = sum(!hit1_is_top_alpha)), by = arm])
cat("caQTL leads with imputation_R2 < 0.3 (testable pairs):", pl[testable == TRUE, round(mean(caqtl_lead_impR2 < 0.3), 4)], "\n")
cat("caQTL lead missing from PolyFun LD panel (testable pairs):", pl[testable == TRUE, round(mean(is.na(lead_ld_r2)), 4)],
    "(shared SNPs are refit SNPs, all present in the panel by construction)\n")

# ---- 3. evidence states per GWAS signal (study-lead) ------------------------------
pl[, meets_entry := testable & PP.H4 > PP_ENTRY & lead_ld_r2 >= LD_ENTRY]
pl[is.na(meets_entry), meets_entry := FALSE]
pl[, meets_entry_1e5 := testable & PP.H4_p12_1e5 > PP_ENTRY & lead_ld_r2 >= LD_ENTRY]
pl[is.na(meets_entry_1e5), meets_entry_1e5 := FALSE]
if (NULL_DEBUG) {   # force 4 pseudo-entries (highest shifted PP.H4, one per lead) so the entry code path runs; not a result
  fk <- pl[testable == TRUE][order(-PP.H4)][!duplicated(lead_key)][seq_len(min(4, .N)), paste(lead_key, peak)]
  pl[paste(lead_key, peak) %in% fk, `:=`(meets_entry = TRUE, PP.H4 = 0.99, lead_ld_r2 = 0.9)]
  cat("NULL-DEBUG forced pseudo-entries:", fk, "\n")
}
state_of <- function(d) {   # d: rows of one study-lead
  t <- d[testable == TRUE]
  if (!nrow(t)) return("not tested in this source")
  if (any(t$PP.H4 > PP_ENTRY & t$lead_ld_r2 >= LD_ENTRY, na.rm = TRUE)) return("colocalized accessible region")
  if (any(t$PP.H4 > PP_ENTRY)) return("colocalized, lead LD not confirmed")
  if (max(t$PP.H3) > PP_ENTRY) return("different variant (unsupported)")
  if (all(t$PP.H0 + t$PP.H1 > PP_ENTRY)) return("no caQTL detectable at 138 donors (untestable)")
  "indeterminate"
}
lead_sum <- pl[, .(state = state_of(.SD), n_cand = n_cand[1], n_testable = sum(testable),
                   any_entry = any(meets_entry), any_entry_1e5 = any(meets_entry_1e5),
                   any_pp4_gt_0.8_noLD = any(testable & PP.H4 > PP_ENTRY, na.rm = TRUE),
                   any_pp4_gt_0.5 = any(testable & PP.H4 > PP_SECOND, na.rm = TRUE),
                   any_pp4_gt_0.5_LD = any(testable & PP.H4 > PP_SECOND & lead_ld_r2 >= LD_ENTRY, na.rm = TRUE),
                   best_pp4 = suppressWarnings(max(PP.H4[testable], na.rm = TRUE)),
                   mlog10_min_caqtl_p = if (any(testable)) pmin(300, -log10(max(1e-300, min(min_caqtl_p_shared[testable])))) else NA_real_,
                   entry_peaks = paste(peak[meets_entry], collapse = ","),
                   entry_centres_hg38 = if (any(meets_entry)) paste(paste0("chr", chr[1], ":", centre[meets_entry]), collapse = ",") else "",
                   entry_pp4 = paste(signif(PP.H4[meets_entry], 4), collapse = ","),
                   entry_ld_r2 = paste(signif(lead_ld_r2[meets_entry], 3), collapse = ","),
                   entry_caqtl_lead = paste(caqtl_lead_variant[meets_entry], collapse = ","),
                   cr_any_entry = any(cr_PP.H4 > PP_ENTRY & cr_lead_ld_r2_hit1 >= LD_ENTRY, na.rm = TRUE)),
               by = lead_key]
lead_sum[is.infinite(best_pp4), best_pp4 := NA_real_]
lead_sum[state == "not tested in this source", n_testable := 0L]
nt <- SL[!lead_key %in% lead_sum$lead_key, .(lead_key)]
SLo <- merge(SL, lead_sum, by = "lead_key", all.x = TRUE)
SLo[is.na(state), `:=`(state = "not tested in this source", n_testable = 0L, any_entry = FALSE)]
SLo[reproduced == FALSE, `:=`(state = "not tested in this source", n_testable = 0L, any_entry = FALSE)]
SLo[, not_tested_reason := fifelse(!reproduced, "refit not reproduced",
                           fifelse(is.na(hit1_pos38), "hit1 liftover failure",
                           fifelse(state == "not tested in this source", "no testable peak (none within 500 kb, <100 shared SNPs, or overlap rule)", NA_character_)))]

# ---- 4. NC1 (T1 loci) and NC2 -----------------------------------------------------
T1L <- LOC[t1_set %in% c("T1_liver_only", "T1_platelet_only")]
T1L <- merge(T1L[, .(locus, lead_key, t1_set)], SLo, by = c("locus", "lead_key"), all.x = TRUE)
nc1_pairs <- pl[lead_key %in% T1L$lead_key & testable == TRUE, .(lead_key, peak, n_shared)]
setorder(nc1_pairs, lead_key, peak)
set.seed(SEED)
nc1_pairs[, u := runif(.N, 0.25, 0.75)]
nc1_pairs[, offset := pmax(1L, as.integer(floor(u * n_shared)))]
fwrite(nc1_pairs, file.path(RES, "nc1_offsets.tsv"), sep = "\t")

# ---- 5. pass 2: NC1 shifts and three-way (re-reads caQTL rows per chromosome) -----
ent <- pl[meets_entry == TRUE, .(lead_key, peak, chr)]
gene_rows <- rbind(a$L[, .(gwas, hit1, gene, ensembl, idx2, r2_pp4 = PP.H4)], a$P[, .(gwas, hit1, gene, ensembl, idx2, r2_pp4 = PP.H4)])
gene_rows[, lead_key := paste(gwas, hit1, sep = "|")]
eqtl_row_keys <- function(ensembl, ch, idx2) {
  e <- readRDS(file.path(EQTL_SUSIE_DIR, paste0("chr", ch), paste0(ensembl, "_susie.rds")))
  v <- e$lbf_variable[idx2, ]
  ids <- colnames(e$lbf_variable)
  p38 <- lift19to38(ch, as.integer(sub("^[0-9]+:", "", ids)))
  k <- ifelse(is.na(p38), paste0("e|", ids), paste0("p", p38))
  d <- duplicated(k); if (any(d)) k[d] <- paste0("edup|", ids[d])
  names(v) <- k; v
}
NC1 <- list(); TW <- list()
for (ch in CHRS) {
  nl <- nc1_pairs[lead_key %in% SL[chr == ch, lead_key]]
  el <- ent[chr == ch]
  if (!nrow(nl) && !nrow(el)) next
  t0 <- Sys.time(); q <- read_caqtl_chr(ch)
  jobs <- unique(c(nl$lead_key, el$lead_key))
  res <- mclapply(jobs, function(lk) {
    r <- SL[lead_key == lk]; fit <- readRDS(r$fit_file); gk <- gwas_keys(fit); grow <- gwas_row(fit, gk, r$idx1)
    pos19 <- setNames(gk$position, gk$key)
    out_nc <- rbindlist(lapply(which(nl$lead_key == lk), function(j) {
      pk <- nl$peak[j]; qp <- caqtl_prepare(copy(q[.(pk)]))
      sh <- setdiff(intersect(gk$key, qp$key), "null")
      sh <- sh[order(pos19[sh])]                       # hg19 position order
      ix <- match(sh, qp$key); n <- length(ix); k <- nl$offset[j]
      rot <- (seq_len(n) + k - 1L) %% n + 1L
      lb <- qp$lbf; pv <- qp$pvalue
      lb[ix] <- qp$lbf[ix[rot]]; pv[ix] <- qp$pvalue[ix[rot]]
      o <- pair_test(grow, gk, qp, fit$ld_rows, r$gwas_lead_snp, p12s = P12, qp_lbf = lb, qp_p = pv)
      data.table(lead_key = lk, peak = pk, nc1_testable = o$testable, nc1_prop_caqtl = o$prop_caqtl,
                 nc1_PP.H4 = if (is.null(o$PP.H4)) NA_real_ else o$PP.H4,
                 nc1_lead_ld_r2 = if (is.null(o$lead_ld_r2)) NA_real_ else o$lead_ld_r2)
    }), fill = TRUE)
    gl <- gene_rows[lead_key == lk]
    out_tw <- rbindlist(lapply(el[lead_key == lk, peak], function(pk) {
      qp <- caqtl_prepare(copy(q[.(pk)]))
      qp <- qp[order(pvalue)]
      kpos <- paste0("p", qp$pos38); d <- duplicated(kpos); kpos[d] <- paste0("qdup|", qp$variant[d])
      crow <- c(qp$lbf, 0); names(crow) <- c(kpos, "null")
      rbindlist(lapply(seq_len(nrow(gl)), function(g) {
        erow <- eqtl_row_keys(gl$ensembl[g], ch, gl$idx2[g])
        sh <- setdiff(intersect(names(erow), names(crow)), "null")
        b1 <- matrix(erow, 1, dimnames = list(NULL, names(erow))); b2 <- matrix(crow, 1, dimnames = list(NULL, names(crow)))
        pp4 <- NA_real_; dropped <- NA
        if (length(sh) >= MIN_SNPS) {
          rr <- quiet(coloc.bf_bf(b1, b2, p1 = P1, p2 = P2, p12 = P12))$summary
          pp4 <- rr$PP.H4.abf[1]; dropped <- is.na(pp4)
        }
        data.table(lead_key = lk, peak = pk, gene = gl$gene[g], ensembl = gl$ensembl[g], idx2 = gl$idx2[g],
                   r2_eqtl_pp4 = gl$r2_pp4[g], tw_n_shared = length(sh), tw_dropped = dropped, tw_PP.H4 = pp4,
                   tw_position_only_match = TRUE)
      }))
    }), fill = TRUE)
    list(nc = out_nc, tw = out_tw)
  }, mc.cores = NCORE, mc.preschedule = FALSE)
  err <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(err)) { print(res[err]); stop("pass-2 worker error on chr ", ch) }
  NC1[[ch]] <- rbindlist(lapply(res, `[[`, "nc"), fill = TRUE)
  TW[[ch]] <- rbindlist(lapply(res, `[[`, "tw"), fill = TRUE)
  stamp("pass 2 chr", ch, "done in", format(Sys.time() - t0)); rm(q, res); gc(FALSE)
}
NC1 <- rbindlist(NC1, fill = TRUE); TW <- rbindlist(TW, fill = TRUE)
if (!nrow(TW)) TW <- data.table(lead_key = character(), peak = character(), gene = character(), ensembl = character(), idx2 = integer(),
                                r2_eqtl_pp4 = numeric(), tw_n_shared = integer(), tw_dropped = logical(), tw_PP.H4 = numeric(), tw_position_only_match = logical())
fwrite(NC1, file.path(RES, "nc1_shifted_pairs.tsv"), sep = "\t"); fwrite(TW, file.path(RES, "three_way.tsv"), sep = "\t")

# ---- 6. outcomes -------------------------------------------------------------------
wilson <- function(x, n, z = qnorm(0.975)) {
  if (n == 0) return(c(NA, NA, NA)); p <- x / n
  den <- 1 + z^2 / n; cen <- (p + z^2 / (2 * n)) / den
  hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(p, cen - hw, cen + hw)
}
newcombe <- function(x1, n1, x2, n2) {   # hybrid score interval (Newcombe 1998, method 10)
  w1 <- wilson(x1, n1); w2 <- wilson(x2, n2); d <- w1[1] - w2[1]
  c(d, d - sqrt((w1[1] - w1[2])^2 + (w2[3] - w2[1])^2), d + sqrt((w1[3] - w1[1])^2 + (w2[1] - w2[2])^2))
}
fmt <- function(v, d = 3) paste0(formatC(v[1], format = "f", digits = d), " [", formatC(v[2], format = "f", digits = d), ", ", formatC(v[3], format = "f", digits = d), "]")

# NC1 outcomes per locus
nc1l <- NC1[, .(nc1_any_noLD = any(nc1_testable & nc1_PP.H4 > PP_ENTRY, na.rm = TRUE),
                nc1_any_full = any(nc1_testable & nc1_PP.H4 > PP_ENTRY & nc1_lead_ld_r2 >= LD_ENTRY, na.rm = TRUE),
                nc1_prop_unchanged = all(abs(nc1_prop_caqtl - pl$prop_caqtl[match(paste(lead_key, peak), paste(pl$lead_key, pl$peak))]) < 1e-9)),
            by = lead_key]
T1L <- merge(T1L, nc1l, by = "lead_key", all.x = TRUE)
T1L[, Y := as.integer(any_entry)]
T1L[, arm01 := as.integer(t1_set == "T1_liver_only")]
T1in <- T1L[reproduced == TRUE & n_testable >= 1]
n_test <- T1in[, .N, by = t1_set]
S3_PASS <- isTRUE(n_test[t1_set == "T1_liver_only", N] >= 120) && isTRUE(n_test[t1_set == "T1_platelet_only", N] >= 60)
nc1_rates <- T1in[, .(n = .N, rate_noLD = mean(nc1_any_noLD %in% TRUE), rate_full = mean(nc1_any_full %in% TRUE),
                      prop_caqtl_unchanged = all(nc1_prop_unchanged %in% TRUE)), by = t1_set]
S4_PASS <- all(nc1_rates$rate_noLD <= 0.10)
nc2 <- pl[testable == TRUE & min_caqtl_p_shared > 1e-3, .(n_pairs = .N, n_pp4_gt_0.8 = sum(PP.H4 > PP_ENTRY), frac = mean(PP.H4 > PP_ENTRY))]

# Primary statistic: 238 confirmatory liver loci
CL <- merge(LOC[arm == "liver" & confirmatory == TRUE, .(locus, lead_key, t1_set)], SLo, by = c("locus", "lead_key"))
stopifnot(NULL_DEBUG || nrow(CL) == 238)
x_prim <- sum(CL$any_entry %in% TRUE)
prim_all <- wilson(x_prim, nrow(CL))
prim_test <- wilson(x_prim, CL[n_testable >= 1, .N])
prim_repr <- wilson(x_prim, CL[reproduced == TRUE, .N])

# T1
t1_res <- list(status = NA_character_)
fit_t1 <- NULL
if (!S3_PASS) {
  t1_res$status <- "untestable (S3: too few loci with a testable peak)"
} else if (!S4_PASS) {
  t1_res$status <- "untestable (S4: NC1 rate above 0.10)"
} else {
  fit_t1 <- glm(Y ~ arm01 + mlog10_lead_p + n_gws_window + n_testable + mlog10_min_caqtl_p, data = T1in, family = binomial())
  co <- coef(fit_t1)
  if (!fit_t1$converged || any(abs(co) > 10)) {
    t1_res$status <- "not estimable (S5: non-convergence or separation)"
  } else {
    sm <- summary(fit_t1)$coefficients
    fit0 <- update(fit_t1, . ~ . - arm01)
    lr <- anova(fit0, fit_t1, test = "LRT")
    ci <- suppressMessages(confint(fit_t1, "arm01"))
    t1_res <- list(status = "estimated", beta = co[["arm01"]], OR = exp(co[["arm01"]]), z = sm["arm01", "z value"],
                   p_wald = sm["arm01", "Pr(>|z|)"], p_LR = lr$`Pr(>Chi)`[2], OR_profile_lo = exp(ci[[1]]), OR_profile_hi = exp(ci[[2]]),
                   n_liver = T1in[arm01 == 1, .N], n_platelet = T1in[arm01 == 0, .N])
  }
}
if (NULL_DEBUG) {   # exercise the model code on the shifted data regardless of gates
  dbg <- tryCatch(glm(Y ~ arm01 + mlog10_lead_p + n_gws_window + n_testable + mlog10_min_caqtl_p, data = T1in, family = binomial()), error = function(e) e)
  cat("NULL-DEBUG glm:\n"); print(if (inherits(dbg, "error")) dbg else summary(dbg)$coefficients)
}
x1 <- T1in[arm01 == 1, sum(Y)]; n1 <- T1in[arm01 == 1, .N]; x2 <- T1in[arm01 == 0, sum(Y)]; n2 <- T1in[arm01 == 0, .N]
unadj <- newcombe(x1, n1, x2, n2)
fisher_p <- if (n1 > 0 && n2 > 0) fisher.test(matrix(c(x1, n1 - x1, x2, n2 - x2), 2))$p.value else NA_real_
bal <- T1in[, lapply(.SD, function(v) c(mean(v), sd(v))), by = arm01, .SDcols = c("mlog10_lead_p", "n_gws_window", "n_testable", "mlog10_min_caqtl_p")]
smd <- sapply(c("mlog10_lead_p", "n_gws_window", "n_testable", "mlog10_min_caqtl_p"), function(v) {
  a1 <- T1in[arm01 == 1, get(v)]; a0 <- T1in[arm01 == 0, get(v)]; (mean(a1) - mean(a0)) / sqrt((var(a1) + var(a0)) / 2) })

# lead-level (max over studies and peaks in the locus) with locus-cluster bootstrap
SLo[, locus_any := any(any_entry %in% TRUE), by = locus]
LL <- unique(SLo[locus %in% T1L$locus, .(locus, locus_any)])
LL <- merge(LL, unique(T1L[, .(locus, arm01, t1_set)]), by = "locus")
LL <- LL[locus %in% SLo[n_testable >= 1, locus]]
set.seed(SEED)
bs <- replicate(1000, { b1 <- sample(LL[arm01 == 1, locus_any], replace = TRUE); b0 <- sample(LL[arm01 == 0, locus_any], replace = TRUE); mean(b1) - mean(b0) })
lead_level <- c(LL[arm01 == 1, mean(locus_any)] - LL[arm01 == 0, mean(locus_any)], quantile(bs, c(.025, .975)))

# MVP_ALT_EUR same-biobank sensitivity
MV <- merge(MVP[mvp_set == "MVP_ALT_liver_only_confirmatory", .(gwas, hit1)], SLo, by = c("gwas", "hit1"))
MV[, `:=`(arm01 = 1L, Y = as.integer(any_entry %in% TRUE))]
MVd <- rbind(MV[reproduced == TRUE & n_testable >= 1, .(Y, arm01, mlog10_lead_p, n_gws_window, n_testable, mlog10_min_caqtl_p)],
             T1in[arm01 == 0, .(Y, arm01, mlog10_lead_p, n_gws_window, n_testable, mlog10_min_caqtl_p)])
mv_fit <- tryCatch(glm(Y ~ arm01 + mlog10_lead_p + n_gws_window + n_testable + mlog10_min_caqtl_p, data = MVd, family = binomial()), error = function(e) NULL)
mv_or <- if (!is.null(mv_fit) && mv_fit$converged) exp(coef(mv_fit)[["arm01"]]) else NA_real_

# 53 shared loci + refit 238 vs 101
SH <- CL[t1_set == "shared_liver_descriptive"]
sh_rate <- wilson(SH[reproduced == TRUE & n_testable >= 1, sum(any_entry)], SH[reproduced == TRUE & n_testable >= 1, .N])
T238 <- rbind(CL[reproduced == TRUE & n_testable >= 1, .(Y = as.integer(any_entry), arm01 = 1L, mlog10_lead_p, n_gws_window, n_testable, mlog10_min_caqtl_p)],
              T1in[arm01 == 0, .(Y, arm01, mlog10_lead_p, n_gws_window, n_testable, mlog10_min_caqtl_p)])
f238 <- tryCatch(glm(Y ~ arm01 + mlog10_lead_p + n_gws_window + n_testable + mlog10_min_caqtl_p, data = T238, family = binomial()), error = function(e) NULL)
or238 <- if (!is.null(f238) && f238$converged) exp(coef(f238)[["arm01"]]) else NA_real_

# post-pilot chr1/chr19 loci and direct MASLD traits
PP19 <- merge(LOC[t1_set == "post_pilot_chr1_chr19" & arm == "liver", .(locus, lead_key)], SLo, by = c("locus", "lead_key"))
tier1 <- fread(file.path(FM, "config/gwas_trait_tier.tsv"))[tier == 1, study_name]
DM <- SLo[arm == "liver" & gwas %in% tier1]

# three-way and gene-level Catalog annotation (liver arm genes)
TW[, three_way := !is.na(tw_PP.H4) & tw_PP.H4 > PP_ENTRY]
liv <- SLo[arm == "liver"]
gl_liv <- merge(a$L[, .(gwas, hit1, gene, ensembl, idx1, idx2, r2_pp4 = PP.H4)][, lead_key := paste(gwas, hit1, sep = "|")],
                liv[, .(lead_key, state, any_entry, entry_peaks, entry_centres_hg38, entry_pp4, entry_ld_r2)], by = "lead_key")
twg <- TW[, .(three_way = any(three_way), tw_best_pp4 = suppressWarnings(max(tw_PP.H4, na.rm = TRUE))), by = .(lead_key, ensembl)]
gl_liv <- merge(gl_liv, twg, by = c("lead_key", "ensembl"), all.x = TRUE)
gl_liv[is.na(three_way), three_way := FALSE]

# ABC (liver / hepatocyte / HepG2 elements, hg19 -> hg38), centre +-250 bp
abc <- fread(file.path(PROJ, "data/external/abc_liver/abc_liver_hepg2_enhancers.tsv.gz"))
abc[, chrn := sub("^chr", "", chr)]
abc <- abc[chrn %in% as.character(1:22)]
abc[, `:=`(s38 = lift19to38(chrn[1], start), e38 = lift19to38(chrn[1], end)), by = chrn]
abc <- abc[!is.na(s38) & !is.na(e38) & e38 > s38]
ep <- unique(pl[meets_entry == TRUE, .(lead_key, peak, chr, centre)])
abc_hits <- rbindlist(lapply(seq_len(nrow(ep)), function(i) {
  h <- abc[chrn == ep$chr[i] & s38 <= ep$centre[i] + 250 & e38 >= ep$centre[i] - 250]
  if (!nrow(h)) return(NULL)
  data.table(lead_key = ep$lead_key[i], peak = ep$peak[i], abc_gene = unique(h$TargetGene))
}))
if (is.null(abc_hits) || !nrow(abc_hits)) abc_hits <- data.table(lead_key = character(), peak = character(), abc_gene = character())
gl_liv[, abc_predicted_link := paste(lead_key, gene) %in% paste(abc_hits$lead_key, abc_hits$abc_gene)]

# Currin / Broadaway source flags and driver peaks
s23 <- fread(file.path(RUN, "downloads_untrusted/currin_S23A.tsv"))
s9  <- fread(file.path(RUN, "downloads_untrusted/currin_S9A.tsv"))
trait_map <- c(UKBB_GGT = "GGT", UKBB_ALT = "ALT", MVP_ALT_EUR = "ALT")
pe <- pl[meets_entry == TRUE, .(lead_key, gwas, peak)]
pe[, currin_trait := trait_map[gwas]]
pe[, reported_currin := !is.na(currin_trait) & paste(currin_trait, peak) %in% paste(s23$GWAS_trait, s23$caPeak)]
s9l <- s9[, .(peak = unlist(strsplit(caPeaks, ","))), by = .(CoordinatedSet_ID, driver_peak)]
pe <- merge(pe, s9l, by = "peak", all.x = TRUE)
pe[, driver_status := fifelse(is.na(CoordinatedSet_ID), "not in a Currin coordinated set",
                      fifelse(is.na(driver_peak) | driver_peak %in% c("na", ""), "coordinated set, no driver called",
                      fifelse(driver_peak == peak, "Currin driver peak", "coordinated set, not the driver")))]
src <- pe[, .(source_flag = if (any(reported_currin)) "reported by Currin 2025" else
                  if (all(is.na(currin_trait))) "not in source tables (trait not in Currin 2025)" else "not in source tables",
              broadaway_flag = "not checked (supplement not open access)",
              driver_status = paste(unique(driver_status), collapse = ";")), by = lead_key]
liv <- merge(liv, src, by = "lead_key", all.x = TRUE)

# gene-level annotation
state_rank <- c("colocalized accessible region", "colocalized, lead LD not confirmed", "different variant (unsupported)",
                "no caQTL detectable at 138 donors (untestable)", "indeterminate", "not tested in this source")
gl_liv <- merge(gl_liv, liv[, .(lead_key, source_flag)], by = "lead_key", all.x = TRUE)
G <- gl_liv[, .(caqtl_region_state = state_rank[min(match(state, state_rank))],
                caqtl_region_hg38 = paste(unique(entry_centres_hg38[entry_centres_hg38 != ""]), collapse = ";"),
                caqtl_pp_h4 = paste(unique(entry_pp4[entry_pp4 != ""]), collapse = ";"),
                caqtl_lead_ld_r2 = paste(unique(entry_ld_r2[entry_ld_r2 != ""]), collapse = ";"),
                caqtl_three_way_eqtl = any(three_way),
                caqtl_abc_predicted_link = any(abc_predicted_link & any_entry),
                caqtl_source_flag = paste(unique(na.omit(source_flag)), collapse = ";"),
                n_gwas_signals = uniqueN(lead_key)), by = .(gene, ensembl)]
G[caqtl_region_state != "colocalized accessible region", caqtl_three_way_eqtl := FALSE]
missing_primary <- setdiff(a$pos_genes, G$ensembl)
G <- rbind(G, data.table(gene = NA_character_, ensembl = missing_primary, caqtl_region_state = "not tested in this source (no tier-1/2 EUR pair)"), fill = TRUE)
M1 <- G[caqtl_three_way_eqtl == TRUE, uniqueN(ensembl)]
M2 <- G[caqtl_region_state == "colocalized accessible region" & caqtl_three_way_eqtl == FALSE, uniqueN(ensembl)]
M_abc <- G[caqtl_region_state == "colocalized accessible region" & caqtl_three_way_eqtl == FALSE & caqtl_abc_predicted_link == TRUE, uniqueN(ensembl)]

N_entries <- liv[state == "colocalized accessible region", .N]
K_currin <- liv[state == "colocalized accessible region" & source_flag == "reported by Currin 2025", .N]
if (!S4_PASS) { liv[, entry_flag := "rule not specific at this locus set"]; G[, entry_flag := "rule not specific at this locus set"] }

# ---- 7. write tables ----------------------------------------------------------------
fwrite(liv, file.path(RES, "entry_table_liver_study_leads.tsv"), sep = "\t")
fwrite(SLo[arm == "platelet"], file.path(RES, "entry_table_platelet_leads.tsv"), sep = "\t")
fwrite(G, file.path(RES, "catalog_gene_caqtl_annotation.tsv"), sep = "\t")
fwrite(T1in, file.path(RES, "t1_model_data.tsv"), sep = "\t")
fwrite(abc_hits, file.path(RES, "abc_overlaps.tsv"), sep = "\t")
cr <- pl[, .(cr_entry = any(cr_PP.H4 > PP_ENTRY & cr_lead_ld_r2_hit1 >= LD_ENTRY, na.rm = TRUE)), by = .(lead_key, arm)]
cr_sens <- pl[!is.na(cr_sens_p12_min_pass), .(lead_key, peak, cr_PP.H4, cr_sens_p12_min_pass, cr_sens_p12_max_pass)]
fwrite(cr_sens, file.path(RES, "currin_route_sensitivity.tsv"), sep = "\t")

# ---- 8. figures (supplementary panels; observations on real axes) -------------------
suppressPackageStartupMessages(library(ggplot2))
th <- theme_classic(base_size = 6) + theme(text = element_text(size = 6, colour = "black", face = "plain"),
        axis.text = element_text(size = 6, colour = "black"), axis.title = element_text(size = 6, colour = "black"),
        legend.position = "none", plot.title = element_blank())
fr <- rbind(data.table(set = "Liver-only loci", x = x1, n = n1, t(wilson(x1, n1))),
            data.table(set = "Platelet-only loci", x = x2, n = n2, t(wilson(x2, n2))))
setnames(fr, c("V1", "V2", "V3"), c("p", "lo", "hi"))
fr[, set := factor(set, levels = c("Liver-only loci", "Platelet-only loci"))]
p1 <- ggplot(fr, aes(set, p)) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.3, colour = "black") +
  geom_point(aes(fill = set), shape = 21, size = 1.8, stroke = 0.3, colour = "black") +
  scale_fill_manual(values = c("Liver-only loci" = "#C2185B", "Platelet-only loci" = "#9E9E9E")) +
  scale_y_continuous(limits = c(0, 1), name = "Fraction of loci with a colocalized\naccessible region (Wilson 95% CI)") + xlab(NULL) + th
pdf(file.path(RES, "supp_panel_locus_fraction_by_arm.pdf"), width = 1.8, height = 2.0, useDingbats = FALSE); print(p1); dev.off()
if (identical(t1_res$status, "estimated")) {
  od <- data.table(model = "T1 (adjusted)", OR = t1_res$OR, lo = t1_res$OR_profile_lo, hi = t1_res$OR_profile_hi)
  p2 <- ggplot(od, aes(OR, model)) + geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "#9E9E9E") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.3, colour = "black") +
    geom_point(size = 1.5, colour = "black") + scale_x_log10(name = "Odds ratio, liver-only vs platelet-only\n(profile 95% CI)") + ylab(NULL) + th
  pdf(file.path(RES, "supp_panel_t1_odds_ratio.pdf"), width = 2.2, height = 1.0, useDingbats = FALSE); print(p2); dev.off()
}

# ---- 9. verdict (prespec "Reading") and RESULT_H53.md ------------------------------
t1_sig <- identical(t1_res$status, "estimated") && t1_res$p_wald < 0.05
t1_word <- if (!identical(t1_res$status, "estimated")) paste("T1", t1_res$status) else
           if (t1_sig && t1_res$OR > 1) "T1 positive (liver-only loci share a liver caQTL more often at equal GWAS and caQTL strength)" else
           if (t1_sig) "T1 significant in the opposite direction (platelet-only loci share more often)" else
           "T1 not significant (a difference smaller than about 0.17-0.20 cannot be excluded)"
verdict <- if (N_entries == 0) "NULL (a): zero entries with controls passing; the abstract must drop 'regulatory DNA elements'" else
           if (!S4_PASS) paste0("ENTRIES WRITTEN WITH FLAG 'rule not specific at this locus set' (S4); T1 untestable; ", N_entries, " GWAS signals") else
           paste0("PASS: ", N_entries, " liver GWAS signals have a colocalized accessible region; ", M1,
                       " Catalog genes three-way supported (M1), ", M2, " reached only through the r2 pair (M2); ", t1_word)
if (NULL_DEBUG) verdict <- paste("NULL-DEBUG RUN (shifted caQTL values, chr21-22 only; not a result):", verdict)

L <- c("# RESULT H53", "",
  paste0("Run: ", RUN, " | job ", Sys.getenv("SLURM_JOB_ID"), " | written ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("Prespec sha256 ", PRESPEC_SHA, " (verified). Selection sha256: ", paste(sub(" .*", "", sel), collapse = ", ")), "",
  paste0("**Verdict (computed by the prespec Reading rule):** ", verdict), "",
  "## Stop rules", "",
  paste0("- S1 reproduction: PASS (job 1). T1 arms: ", paste(sprintf("%s %d/%d", unlist(lapply(gate$S1, `[[`, "t1_set")), unlist(lapply(gate$S1, `[[`, "reproduced")), unlist(lapply(gate$S1, `[[`, "n"))), collapse = "; ")),
  paste0("- S2 positive controls: PASS (job 1; see results/controls_pc1_pc3.tsv, controls_pc2_pc2b.tsv)."),
  paste0("- S3 testable loci: ", if (S3_PASS) "PASS" else "TRIGGERED", " (liver-only ", n1, ", platelet-only ", n2, "; floor 120 / 60)."),
  paste0("- S4 NC1 locus rate at PP.H4 > 0.8, no LD rule: ", paste(sprintf("%s %.4f (n %d)", nc1_rates$t1_set, nc1_rates$rate_noLD, nc1_rates$n), collapse = "; "),
         "; full rule: ", paste(sprintf("%s %.4f", nc1_rates$t1_set, nc1_rates$rate_full), collapse = "; "), " -> ", if (S4_PASS) "PASS (<= 0.10)" else "TRIGGERED (> 0.10)",
         ". prop_caqtl unchanged by the shift: ", paste(nc1_rates$prop_caqtl_unchanged, collapse = "/"), ". Prediction was <= 0.03 per arm."),
  paste0("- S5 model: ", t1_res$status, "."), "",
  "## Primary statistic (estimation, no test)", "",
  paste0("Confirmatory liver loci (chr1, chr19 excluded) whose selected lead has >= 1 peak meeting the entry rule (bf_bf PP.H4 > 0.8 and lead LD r2 >= 0.5): ",
         x_prim, "."),
  paste0("- of all ", nrow(CL), " loci: ", fmt(prim_all)),
  paste0("- of ", CL[n_testable >= 1, .N], " loci with a testable peak: ", fmt(prim_test)),
  paste0("- of ", CL[reproduced == TRUE, .N], " loci with a reproduced lead: ", fmt(prim_repr)),
  "Conditional on prior liver-eQTL colocalization (SuSiE PP.H4 > 0.5); not the share of liver GWAS loci with a caQTL.", "",
  "## T1 (the only test; family = {T1}, no correction)", "",
  if (identical(t1_res$status, "estimated")) paste0("Logistic regression, ", t1_res$n_liver, " liver-only vs ", t1_res$n_platelet,
     " platelet-only confirmatory loci with >= 1 testable peak; covariates -log10 lead p, n GWS variants in 500 kb, n testable peaks, -log10 min caQTL p. ",
     sprintf("arm OR %.3f (profile 95%% CI %.3f to %.3f), Wald z %.3f, p %.4g; LR p %.4g.", t1_res$OR, t1_res$OR_profile_lo, t1_res$OR_profile_hi, t1_res$z, t1_res$p_wald, t1_res$p_LR))
  else paste0("T1 ", t1_res$status, "."),
  paste0("Unadjusted: liver-only ", x1, "/", n1, ", platelet-only ", x2, "/", n2, "; difference ", fmt(unadj), " (Newcombe hybrid score); Fisher exact p ", signif(fisher_p, 4), " (descriptive)."),
  paste0("Covariate standardized mean differences (liver - platelet): ", paste(sprintf("%s %.2f", names(smd), smd), collapse = "; "), "."), "",
  "## Controls", "",
  paste0("- NC2 (testable pairs with min caQTL p on shared SNPs > 1e-3): ", nc2$n_pp4_gt_0.8, "/", nc2$n_pairs, " = ", signif(nc2$frac, 3), " with PP.H4 > 0.8 (prediction <= 0.005)."),
  "- PC1/PC2/PC2b/PC3: see job 1 log and results/controls_*.tsv.", "",
  "## Catalog deliverable (descriptive)", "",
  paste0("- Liver GWAS signals (study-leads, all chromosomes) by evidence state: ", paste(liv[, .N, by = state][, sprintf("%s %d", state, N)], collapse = "; "), "."),
  paste0("- Signals with a colocalized accessible region: ", N_entries, " (", K_currin, " already reported by Currin 2025 for the same trait and caPeak; Broadaway 2024 not checked)."),
  paste0("- Same at p12 = 1e-5: ", liv[any_entry_1e5 %in% TRUE, .N], ". PP.H4 > 0.5 with LD rule: ", liv[any_pp4_gt_0.5_LD %in% TRUE, .N], "."),
  paste0("- Genes: M1 three-way supported ", M1, "; M2 reached only through the r2 GWAS-eQTL pair ", M2, " (of which ABC-only predicted link ", M_abc, ")."),
  paste0("- Currin route (coloc.abf, +-500 kb, lead LD r2 >= 0.5 to hit1) liver signals with an entry: ", cr[arm == "liver", sum(cr_entry)], "; platelet: ", cr[arm == "platelet", sum(cr_entry)], "."), "",
  "## Descriptive comparisons (no p-value calls)", "",
  paste0("- Lead-level (max over studies and peaks) difference, locus-cluster bootstrap (1,000; seed ", SEED, "): ", fmt(lead_level), "."),
  paste0("- MVP_ALT_EUR liver-only (", MV[reproduced == TRUE & n_testable >= 1, .N], " loci) vs platelet-only: adjusted OR ", signif(mv_or, 4), "."),
  paste0("- 53 shared confirmatory liver loci: rate ", fmt(sh_rate), "; T1 model with them kept (238 vs 101): OR ", signif(or238, 4), "."),
  paste0("- Post-pilot chr1/chr19 liver loci: ", PP19[any_entry %in% TRUE, .N], "/", nrow(PP19), " with an entry."),
  paste0("- Direct MASLD-trait signals (tier 1): ", paste(DM[, .N, by = state][, sprintf("%s %d", state, N)], collapse = "; "), "."), "",
  "Tables: results/entry_table_liver_study_leads.tsv, catalog_gene_caqtl_annotation.tsv, lead_peak_pairs.tsv.gz, three_way.tsv, nc1_*.tsv, t1_model_data.tsv.",
  "Figures: results/supp_panel_locus_fraction_by_arm.pdf, supp_panel_t1_odds_ratio.pdf.")
writeLines(L, file.path(RES, "RESULT_H53.md"))
if (!NULL_DEBUG) file.copy(file.path(RES, "RESULT_H53.md"), file.path(RUN, "RESULT_H53.md"), overwrite = FALSE)
sink(file.path(RES, "t1_glm_summary.txt")); if (!is.null(fit_t1)) print(summary(fit_t1)); print(bal); sink()
writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo_job2.txt"))
cat(L, sep = "\n")
stamp("H53 job 2 done")
