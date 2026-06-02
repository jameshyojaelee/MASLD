#!/usr/bin/env Rscript
# 00d_meta_sveinbjornsson.R
# ---------------------------------------------------------------------------
# Fixed-effects inverse-variance-weighted (IVW) meta-analysis of the Sveinbjornsson-2022
# per-cohort case-control arms into one well-powered GWAS for COLOC. The per-cohort arms are
# individually underpowered (e.g. cirrhosis: deCODE 691 / UKB 1,425 / Intermountain 2,301 /
# FinnGen 392 cases); the meta recovers the paper's actual power (cirrhosis 4,809 cases; HCC 861).
#
# Effect handling: deCODE-format files give Effect=OR (no SE) -> beta=log(OR), SE back-computed from
# beta+pval; the FinnGen-format file gives Beta(log-OR)+Standard_Error directly. QC per cohort
# (sentinel-OR drop, MAF>=1%, INFO>=0.3, |log-OR|<=5, strand-ambiguous drop), allele-aligned to the
# alphabetically-first allele, then IVW-combined per variant. hg38 -> hg19 liftover at the end.
#
# Usage: META_PHENO=CIRRHOSIS micromamba run -n rnaseq Rscript src/00d_meta_sveinbjornsson.R   (or HCC)
# Output: data/sumstats/Sveinbjornsson2022_<Pheno>_meta_EUR_reformatted_hg19.tsv
# ---------------------------------------------------------------------------
suppressMessages({ library(data.table); library(rtracklayer); library(GenomicRanges) })
setDTthreads(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SV   <- file.path(BASE, "GWAS/MR_Data/Sveinbjornsson2022")
OUT  <- file.path(BASE, "GWAS/finemapping/data/sumstats")
chain <- import.chain(file.path(BASE, "data/broadaway_eqtl/hg38ToHg19.over.chain"))
PHENO <- toupper(Sys.getenv("META_PHENO", "CIRRHOSIS"))

# de-rotated CMDKP N (see docs/review/2026-05-31-parallel-findings.md; verify vs paper supplement)
COHORTS <- list(
  CIRRHOSIS = list(
    list(f = "Cirrhosis_deCODE_sumstat.txt",        fmt = "decode",  ca = 691,  co = 396414),
    list(f = "Cirrhosis_INTERMOUNTAIN_sumstat.txt", fmt = "decode",  ca = 2301, co = 24715),
    list(f = "Cirrhosis_UKBB_sumstat.txt",          fmt = "decode",  ca = 1425, co = 366616),
    list(f = "Cirrhosis_FINNGEN_sumstats.txt",      fmt = "finngen", ca = 392,  co = 176196)),
  HCC = list(
    list(f = "HCC_deCODE_sumstat.txt",        fmt = "decode", ca = 406, co = 308683),
    list(f = "HCC_INTERMOUNTAIN_sumstat.txt", fmt = "decode", ca = 81,  co = 26333),
    list(f = "HCC_UKBB_sumstat.txt",          fmt = "decode", ca = 374, co = 406264))
)[[PHENO]]
stopifnot("unknown META_PHENO" = !is.null(COHORTS))
n_cases <- sum(sapply(COHORTS, `[[`, "ca")); n_tot <- n_cases + sum(sapply(COHORTS, `[[`, "co"))
cat(sprintf("=== Sveinbjornsson-2022 %s meta: %d cohorts, %d cases / %d total ===\n",
            PHENO, length(COHORTS), n_cases, n_tot))

liftover_hg19 <- function(dt) {
  gr <- GRanges(paste0("chr", dt$chromosome), IRanges(dt$position, width = 1))
  lifted <- liftOver(gr, chain); keep <- lengths(lifted) == 1
  cat(sprintf("  liftover: %d / %d mapped 1:1\n", sum(keep), length(keep)))
  dt <- dt[keep]; dt$position <- start(unlist(lifted[keep])); dt
}

parse_cohort <- function(path, fmt, label) {
  if (!file.exists(path)) { cat("  MISSING:", path, "\n"); return(NULL) }
  dt <- fread(path)
  if (fmt == "decode") {
    setnames(dt, c("Pval","Effect","Chrom","Pos","Amin","Amaj","MAF_PC","Info"),
                 c("pval","eff","chr","pos","ea","oa","maf","info"), skip_absent = TRUE)
    dt <- dt[, .(chr, pos, ea = toupper(as.character(ea)), oa = toupper(as.character(oa)),
                 pval = as.numeric(pval), beta = log(as.numeric(eff)),
                 maf = as.numeric(maf) / 100, info = as.numeric(info))]
    z <- qnorm(1 - pmin(pmax(dt$pval, 1e-300), 1) / 2); dt[, se := abs(beta) / z]
  } else {  # finngen: Beta(log-OR) + Standard_Error; effect allele = effectAll; derive other
    setnames(dt, c("Pval","Standard_Error","Beta","Chrom","Pos","effectAll","canonRef","canonAlt","Freq_PC","Info"),
                 c("pval","se","beta","chr","pos","ea","cref","calt","maf","info"), skip_absent = TRUE)
    dt <- dt[, .(chr, pos, ea = toupper(as.character(ea)), cref = toupper(as.character(cref)),
                 calt = toupper(as.character(calt)), pval = as.numeric(pval), beta = as.numeric(beta),
                 se = as.numeric(se), maf = as.numeric(maf) / 100, info = as.numeric(info))]
    dt[, oa := fifelse(ea == calt, cref, calt)][, `:=`(cref = NULL, calt = NULL)]
    dt[, maf := pmin(maf, 1 - maf)]
  }
  dt[, chr := as.integer(sub("^chr", "", as.character(chr), ignore.case = TRUE))]
  dt[, pos := as.integer(pos)]
  dt <- dt[!is.na(chr) & chr >= 1 & chr <= 22 & !is.na(pos) &
             is.finite(beta) & is.finite(se) & se > 0 & beta != 0]
  n0 <- nrow(dt)
  top <- dt[, .N, by = beta][order(-N)][1]                    # placeholder-OR sentinel
  if (!is.na(top$N) && top$N / nrow(dt) > 0.01) dt <- dt[beta != top$beta]
  dt <- dt[(is.na(maf) | maf >= 0.01) & (is.na(info) | info >= 0.3) & abs(beta) <= 5]
  dt <- dt[!((ea %in% c("A","T") & oa %in% c("A","T")) | (ea %in% c("C","G") & oa %in% c("C","G")))]
  dt <- dt[nchar(ea) == 1 & nchar(oa) == 1 & ea != oa]
  dt[, A1 := pmin(ea, oa)][, A2 := pmax(ea, oa)]              # canonical orientation
  dt[, beta := fifelse(ea == A1, beta, -beta)]
  dt[, key := paste(chr, pos, A1, A2, sep = ":")]
  dt <- dt[!duplicated(key)]
  cat(sprintf("  [%s] %d -> %d variants (QC)\n", label, n0, nrow(dt)))
  dt[, .(key, chr, pos, A1, A2, beta, se)]
}

allc <- rbindlist(lapply(COHORTS, function(c) parse_cohort(file.path(SV, c$f), c$fmt, c$f)))
if (is.null(allc) || !nrow(allc)) stop("no variants parsed")
allc[, w := 1 / se^2]
meta <- allc[, .(beta = sum(beta * w) / sum(w), se = sqrt(1 / sum(w)), n_coh = .N,
                 chromosome = chr[1], position = pos[1], allele1 = A1[1], allele2 = A2[1]), by = key]
meta[, z := beta / se][, pval := 2 * pnorm(-abs(z))]
lambda <- median(meta$z^2, na.rm = TRUE) / qchisq(0.5, 1)
cat(sprintf("\nMeta variants: %d (n_cohorts: 1=%d 2=%d 3=%d 4=%d) | lambda_GC=%.3f | min p=%.2e\n",
            nrow(meta), sum(meta$n_coh==1), sum(meta$n_coh==2), sum(meta$n_coh==3), sum(meta$n_coh==4),
            lambda, min(meta$pval, na.rm = TRUE)))
out <- liftover_hg19(meta[, .(chromosome, position, allele1, allele2, beta, se, pval)])
study <- sprintf("Sveinbjornsson2022_%s_meta_EUR", tools::toTitleCase(tolower(PHENO)))
of <- file.path(OUT, paste0(study, "_reformatted_hg19.tsv"))
fwrite(out, of, sep = "\t")
cat(sprintf("\nWritten: %s (%d variants)\n", of, nrow(out)))
cat(sprintf("Registry row to add (tab-separated):\n%s\tdata/sumstats/%s_reformatted_hg19.tsv\tdata/lead_snps/%s_leadSNPs.tsv\tEUR\tbinary\t%d\t%d\tukbb_eur\t0.5\n",
            study, study, study, n_tot, n_cases))
