#!/usr/bin/env Rscript
# format_mvp_for_coloc.R
# Reformat one MVP R4 (dbGaP phs002453) GWAS sumstat into the finemapping COLOC
# schema. MVP R4 GIA coordinates are GRCh38/hg38 (verified: rs738409 @ chr22:43,928,847),
# but the Broadaway eQTL + LD panels + the canonical registry are hg19, and 06_susie_coloc.R
# merges on chr:pos -> we MUST liftOver hg38 -> hg19.
#
# MVP raw columns: SNP_ID chrom pos ref alt ea af num_samples beta sebeta pval r2 q_pval i2 direction
#   (binary phecodes instead carry: ... or ci pval ... -> beta=log(OR), se from 95% CI)
#   - `ea` = effect allele; `beta` is wrt ea.
# Output schema (tab, header): chromosome position[hg19] allele1 allele2 beta se pval
#   - allele1 = effect allele (06_susie_coloc.R aligns eQTL EA to gwas allele1) = MVP ea
#   - allele2 = the non-effect allele (ref or alt, whichever != ea)
#
# Usage: Rscript format_mvp_for_coloc.R <input.txt.gz> <output_reformatted_hg19.tsv>

suppressPackageStartupMessages({
  library(data.table); library(rtracklayer); library(GenomicRanges)
})
CHAIN <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/broadaway_eqtl/hg38ToHg19.over.chain"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: format_mvp_for_coloc.R <input.gz> <output.tsv>")
in_gz  <- args[1]
out_ts <- args[2]

cat("[format_mvp] in :", in_gz, "\n")
cat("[format_mvp] out:", out_ts, "\n")

dt <- fread(in_gz, showProgress = FALSE)
n0 <- nrow(dt)
cat("[format_mvp] raw rows:", n0, " cols:", paste(names(dt), collapse=","), "\n")

# --- column resolution (robust to minor header variation) ---
pick <- function(cands) { hit <- intersect(cands, names(dt)); if (length(hit)) hit[1] else NA_character_ }
c_chr  <- pick(c("chrom","CHR","chr","chromosome","#chrom"))
c_pos  <- pick(c("pos","POS","position","base_pair_location","BP"))
c_ref  <- pick(c("ref","REF","other_allele","A2","ALLELE0"))
c_alt  <- pick(c("alt","ALT","A1","ALLELE1"))
c_ea   <- pick(c("ea","EA","effect_allele","A1"))
c_beta <- pick(c("beta","BETA","b","effect"))           # quantitative traits
c_se   <- pick(c("sebeta","se","SE","standard_error","stderr"))
c_or   <- pick(c("or","OR","odds_ratio"))               # binary traits (SAIGE -> OR + CI)
c_ci   <- pick(c("ci","CI","or_ci","confidence_interval"))
c_p    <- pick(c("pval","PVAL","P","p_value","pvalue"))

is_binary <- is.na(c_beta) && !is.na(c_or)
need <- c(chrom=c_chr, pos=c_pos, ref=c_ref, alt=c_alt, ea=c_ea, pval=c_p)
if (is_binary) need <- c(need, or=c_or, ci=c_ci) else need <- c(need, beta=c_beta, se=c_se)
if (any(is.na(need))) stop(paste("Missing required column(s):",
                                 paste(names(need)[is.na(need)], collapse=", ")))
cat("[format_mvp] trait type:", if (is_binary) "BINARY (OR+CI -> log-OR, se from 95% CI)" else "QUANTITATIVE (beta/se)", "\n")

if (is_binary) {
  # beta = log(OR); se from 95% CI "lower,upper": (log U - log L)/(2*1.959964)
  ci_str <- as.character(dt[[c_ci]])
  parts  <- tstrsplit(ci_str, ",", fixed = TRUE)
  ci_lo  <- suppressWarnings(as.numeric(parts[[1]]))
  ci_hi  <- suppressWarnings(as.numeric(parts[[2]]))
  orv    <- suppressWarnings(as.numeric(dt[[c_or]]))
  d <- data.table(
    chrom = as.character(dt[[c_chr]]),
    position = suppressWarnings(as.integer(dt[[c_pos]])),
    ref = toupper(as.character(dt[[c_ref]])),
    alt = toupper(as.character(dt[[c_alt]])),
    ea  = toupper(as.character(dt[[c_ea]])),
    beta = log(orv),
    se   = (log(ci_hi) - log(ci_lo)) / (2 * 1.959963985),
    pval = suppressWarnings(as.numeric(dt[[c_p]])))
} else {
  d <- dt[, .(chrom = as.character(get(c_chr)),
              position = suppressWarnings(as.integer(get(c_pos))),
              ref = toupper(as.character(get(c_ref))),
              alt = toupper(as.character(get(c_alt))),
              ea  = toupper(as.character(get(c_ea))),
              beta = suppressWarnings(as.numeric(get(c_beta))),
              se   = suppressWarnings(as.numeric(get(c_se))),
              pval = suppressWarnings(as.numeric(get(c_p))))]
}

# chromosome -> integer 1..22 (drop X/Y/MT and 'chr' prefix)
d[, chrom := sub("^chr", "", chrom, ignore.case = TRUE)]
d[, chromosome := suppressWarnings(as.integer(chrom))]
d <- d[!is.na(chromosome) & chromosome >= 1 & chromosome <= 22]
n_chr <- nrow(d)                                    # C2 attrition: post autosome 1..22 filter

# effect allele = allele1 (= ea); other allele = ref/alt that isn't ea
d[, allele1 := ea]
d[, allele2 := fifelse(ea == alt, ref, fifelse(ea == ref, alt, NA_character_))]

# QC: finite stats, valid alleles, ea must be one of {ref,alt}
d <- d[!is.na(position) & is.finite(beta) & is.finite(se) & se > 0 &
       is.finite(pval) & pval >= 0 & pval <= 1 &
       !is.na(allele2) & nchar(allele1) >= 1 & nchar(allele2) >= 1 &
       allele1 != allele2]
n_QC <- nrow(d)                                     # C2 attrition: post finite-stat / valid-allele QC

# --- liftOver hg38 -> hg19 (match Broadaway eQTL + LD panels) ---
chain  <- import.chain(CHAIN)
gr     <- GRanges(seqnames = paste0("chr", d$chromosome),
                  ranges = IRanges(start = d$position, width = 1))
lifted <- liftOver(gr, chain)
nmap   <- lengths(lifted)
pos19  <- rep(NA_integer_,   nrow(d))
chr19  <- rep(NA_character_, nrow(d))
u      <- which(nmap == 1L)                 # unique 1:1 maps only
if (length(u)) { lu <- unlist(lifted[u]); pos19[u] <- start(lu); chr19[u] <- as.character(seqnames(lu)) }
d[, `:=`(pos_hg19 = pos19, chr_hg19 = chr19)]
n_pre <- nrow(d)
d <- d[!is.na(pos_hg19) & chr_hg19 == paste0("chr", chromosome)]   # keep same-chr 1:1 lifts
n_lifted <- nrow(d)                                 # C2 attrition: post same-chr 1:1 hg38->hg19 liftOver
d[, position := as.integer(pos_hg19)]
d[, c("pos_hg19","chr_hg19") := NULL]
cat(sprintf("[format_mvp] liftOver hg38->hg19: %d/%d mapped (%.1f%%)\n",
            nrow(d), n_pre, 100*nrow(d)/max(n_pre,1)))

# exact-variant dedup (chr:pos:a1:a2 on hg19), keep most significant
setorder(d, chromosome, position, pval)
d <- d[!duplicated(paste(chromosome, position, allele1, allele2))]

out <- d[, .(chromosome, position, allele1, allele2, beta, se, pval)]
fwrite(out, out_ts, sep = "\t")
cat("[format_mvp] kept rows:", nrow(out), " (", round(100*nrow(out)/n0,1), "% )\n")
cat("[format_mvp] chr range:", min(out$chromosome), "-", max(out$chromosome),
    " | n genome-wide sig (p<5e-8):", sum(out$pval < 5e-8), "\n")

# --- C1 standing anchor guard: rs738409 / PNPLA3 (hg19 chr22:44324727) ------------------
# The PNPLA3 I148M missense variant is the canonical MASLD effect anchor. Its MVP
# effect-allele coding must be INVARIANT across every reformatted stratum (allele1/EA = C,
# allele2/OA = G). This guard hard-fails the reformat if the anchor allele orientation drifts,
# catching an upstream ea/ref/alt column mix-up BEFORE it silently corrupts downstream
# allele harmonization (06_susie_coloc.R flips eqtl_beta on the EA/NEA orientation).
# NOTE: rs738409 is itself a C/G palindrome, so it is later dropped from the in-coloc merge
# by 06:197-199 (strand-ambiguous filter). This guard therefore validates MVP-level
# formatting/harmonization, not the in-coloc allele set. The beta SIGN is trait-dependent
# (liver-injury/disease traits negative for the C allele; Platelet positive = genuine
# biology), so the guard asserts allele identity only and merely REPORTS the sign.
ANCHOR_CHR <- 22L; ANCHOR_POS <- 44324727L
ANCHOR_EA  <- "C"; ANCHOR_OA  <- "G"
anc <- out[chromosome == ANCHOR_CHR & position == ANCHOR_POS]
if (nrow(anc) >= 1L) {
  a  <- anc[1]
  ok <- (a$allele1 == ANCHOR_EA && a$allele2 == ANCHOR_OA)
  cat(sprintf("[format_mvp] anchor rs738409 chr22:44324727  EA=%s OA=%s beta=%.4g sign=%s  [%s]\n",
              a$allele1, a$allele2, a$beta, ifelse(a$beta < 0, "-", "+"),
              ifelse(ok, "OK", "DRIFT")))
  if (!ok) stop(sprintf(
    "[format_mvp] ANCHOR DRIFT: rs738409 EA/OA = %s/%s but expected %s/%s -- allele harmonization would be corrupted. Aborting reformat of %s",
    a$allele1, a$allele2, ANCHOR_EA, ANCHOR_OA, out_ts))
} else {
  cat("[format_mvp] anchor rs738409 chr22:44324727 absent from output ",
      "(partial / non-genome-wide input?) -- anchor guard skipped\n")
}

# --- C2 coordinate/build liftOver attrition record + loss guard (2026-07-05) --
# Persist a machine-readable per-stratum attrition record so the standalone C2
# audit (qa_campaign/audit/C2_liftover_attrition_audit.R) can assemble
# GWAS/MR_Data/MVP/liftover_attrition.tsv WITHOUT re-reading raw gz or re-running
# liftOver. The five stage counts are captured inline where each filter is applied:
#   n_raw   = n0        raw fread rows
#   n_chr   = n_chr     after autosome 1..22 filter (line ~84)
#   n_QC    = n_QC      after finite-stat / valid-allele QC (line ~95)
#   n_lifted= n_lifted  after same-chr 1:1 hg38->hg19 liftOver (line ~109)
#   n_dedup = nrow(out) after chr:pos:a1:a2 exact-variant dedup (line ~117)
# Race-safe: each stratum writes its OWN part file (parallel reformat jobs never
# share a handle); the audit consolidates all parts + logs + reformatted files.
# GUARD: warn when total coordinate/build loss > 1%, hard-fail > 5% (a >5% loss
# almost always means a wrong chain / build mismatch, not benign attrition).
attr <- data.table(
  stratum        = sub("_reformatted_hg19\\.tsv$", "", basename(out_ts)),
  n_raw          = n0,
  n_chr          = n_chr,
  n_QC           = n_QC,
  n_lifted       = n_lifted,
  n_dedup        = nrow(out),
  anchor_present = nrow(anc) >= 1L,
  anchor_EA      = if (nrow(anc) >= 1L) anc[1]$allele1 else NA_character_,
  anchor_OA      = if (nrow(anc) >= 1L) anc[1]$allele2 else NA_character_)
attr[, `:=`(
  loss_chr_pct   = 100 * (n_raw    - n_chr)    / pmax(n_raw,    1L),
  loss_qc_pct    = 100 * (n_chr    - n_QC)     / pmax(n_chr,    1L),
  loss_lift_pct  = 100 * (n_QC     - n_lifted) / pmax(n_QC,     1L),
  loss_dedup_pct = 100 * (n_lifted - n_dedup)  / pmax(n_lifted, 1L),
  loss_total_pct = 100 * (n_raw    - n_dedup)  / pmax(n_raw,    1L))]
LOSS_WARN_PCT <- 1.0; LOSS_FAIL_PCT <- 5.0
attr[, flag := fifelse(loss_total_pct > LOSS_FAIL_PCT | !anchor_present, "FAIL",
                fifelse(loss_total_pct > LOSS_WARN_PCT, "WARN", "PASS"))]
parts_dir <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/attrition_parts"
try({
  dir.create(parts_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(attr, file.path(parts_dir, paste0(attr$stratum, ".tsv")), sep = "\t")
}, silent = TRUE)
cat(sprintf("[format_mvp] ATTRITION\t%s\traw=%d chr=%d QC=%d lifted=%d dedup=%d loss_total=%.4f%% anchor=%s flag=%s\n",
            attr$stratum, attr$n_raw, attr$n_chr, attr$n_QC, attr$n_lifted, attr$n_dedup,
            attr$loss_total_pct, ifelse(attr$anchor_present, "present", "ABSENT"), attr$flag))
if (attr$loss_total_pct > LOSS_WARN_PCT)
  cat(sprintf("[format_mvp] WARN: %s coordinate/build attrition %.4f%% exceeds warn threshold %.1f%%\n",
              attr$stratum, attr$loss_total_pct, LOSS_WARN_PCT))
if (attr$loss_total_pct > LOSS_FAIL_PCT)
  stop(sprintf("[format_mvp] ATTRITION FAIL: %s total loss %.4f%% exceeds fail threshold %.1f%% -- likely wrong liftOver chain or build mismatch. Aborting reformat of %s",
               attr$stratum, attr$loss_total_pct, LOSS_FAIL_PCT, out_ts))
