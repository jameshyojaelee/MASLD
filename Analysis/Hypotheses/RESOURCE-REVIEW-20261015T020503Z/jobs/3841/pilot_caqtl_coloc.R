#!/usr/bin/env Rscript
# Gate L10-04 pilot (feasibility only, not a result).
# 1. Positive control: Currin 2025 Genome Res (PMID 40393811) reports two GGT
#    GWAS signals colocalized with liver caQTL peaks caPeak21014 (chr1) and
#    caPeak156441 (chr19). Recompute coloc.abf(UKBB_GGT, Currin caQTL) here.
# 2. Pilot estimand on chr1 + chr19 only: for UKBB_GGT EUR tier-1/2 colocalized
#    leads and MVP_Platelet_EUR colocalized leads, max PP.H4 over caQTL peaks
#    whose centre lies within +-500 kb of the GWAS lead.
suppressPackageStartupMessages({
  library(data.table); library(coloc); library(rtracklayer); library(GenomicRanges)
})
set.seed(20261015)
cat("coloc", as.character(packageVersion("coloc")), "\n")
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CQ <- file.path(BASE, "data/external/currin_2025_caqtl")
SS <- file.path(BASE, "GWAS/finemapping/data/sumstats")
SCR <- file.path(BASE, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L10-04")
JOB <- getwd()
WIN <- 5e5
N_CAQTL <- 138
N_GWAS <- c(UKBB_GGT = 343850, MVP_Platelet_EUR = 407021)  # gwas_registry.tsv N_tot

chain <- import.chain(file.path(JOB, "hg19ToHg38.over.chain"))

lift <- function(chr, pos) {
  gr <- GRanges(paste0("chr", chr), IRanges(pos, width = 1), idx = seq_along(pos))
  l <- unlist(liftOver(gr, chain))
  l <- l[!duplicated(l$idx)]
  out <- rep(NA_integer_, length(pos)); out[l$idx] <- start(l)
  ok <- rep(FALSE, length(pos)); ok[l$idx] <- as.character(seqnames(l)) == paste0("chr", chr)
  out[!ok] <- NA_integer_
  out
}

read_caqtl <- function(ch) {
  f <- list.files(CQ, pattern = paste0("^caQTL_hg38_chr", ch, "_.*txt.gz$"), full.names = TRUE)
  d <- fread(cmd = paste("zcat", shQuote(f)), sep = "\t",
             select = c("variant", "peak", "distance_from_peakCenter", "pvalue", "beta", "varbeta", "MAF"))
  sp <- tstrsplit(d$variant, ":", fixed = TRUE)
  d[, `:=`(pos38 = as.integer(sp[[2]]), ref = sp[[3]], alt = sp[[4]])]
  d[, centre := pos38 - distance_from_peakCenter]
  d <- d[MAF >= 0.01 & nchar(ref) == 1 & nchar(alt) == 1]
  d[, akey := paste(pos38, pmin(ref, alt), pmax(ref, alt), sep = ":")]
  d
}

read_gwas <- function(study, ch) {
  f <- file.path(JOB, paste0(study, "_chr", ch, ".tsv"))
  g <- fread(f, header = FALSE, col.names = c("chromosome", "position", "allele1", "allele2", "beta", "se", "pval"))
  g <- g[nchar(allele1) == 1 & nchar(allele2) == 1 & se > 0]
  g[, pos38 := lift(ch, position)]
  g <- g[!is.na(pos38)]
  g[, akey := paste(pos38, pmin(allele1, allele2), pmax(allele1, allele2), sep = ":")]
  g <- g[!duplicated(akey)]
  g
}

run_coloc <- function(g, q, study) {
  m <- merge(g[, .(akey, gb = beta, gv = se^2, gp = pval)], q[, .(akey, qb = beta, qv = varbeta, MAF)], by = "akey")
  m <- m[is.finite(gv) & is.finite(qv) & qv > 0 & gv > 0]
  if (nrow(m) < 50) return(list(pp4 = NA_real_, pp3 = NA_real_, n = nrow(m), gmin = NA_real_))
  r <- suppressMessages(coloc.abf(
    dataset1 = list(beta = m$gb, varbeta = m$gv, snp = m$akey, type = "quant", sdY = 1, N = N_GWAS[[study]]),
    dataset2 = list(beta = m$qb, varbeta = m$qv, snp = m$akey, type = "quant", sdY = 1, N = N_CAQTL, MAF = m$MAF)))
  list(pp4 = unname(r$summary["PP.H4.abf"]), pp3 = unname(r$summary["PP.H3.abf"]), n = nrow(m), gmin = min(m$gp))
}

res_pc <- list(); res_loci <- list()
eur <- fread(file.path(SCR, "eur_tier12_pairs.tsv"))
plt <- fread(file.path(SCR, "platelet_eur_pairs.tsv"))

for (ch in c("19", "1")) {
  t0 <- Sys.time()
  q <- read_caqtl(ch)
  cat("chr", ch, "caQTL rows", nrow(q), "peaks", uniqueN(q$peak), "read", format(Sys.time() - t0), "\n")
  setkey(q, peak)
  peaks <- q[, .(centre = centre[1], minp = min(pvalue)), by = peak]
  for (study in c("UKBB_GGT", "MVP_Platelet_EUR")) {
    g <- read_gwas(study, ch)
    cat("  ", study, "chr", ch, "GWAS rows lifted", nrow(g), "\n")
    # positive control (GGT only)
    if (study == "UKBB_GGT") {
      for (pk in c("peak21014", "peak156441")) {
        qq <- q[.(pk), nomatch = 0L]
        if (nrow(qq) == 0) next
        lead <- qq[which.min(pvalue)]
        reg_q <- qq[abs(pos38 - lead$pos38) <= WIN]
        reg_g <- g[abs(pos38 - lead$pos38) <= WIN]
        r <- run_coloc(reg_g, reg_q, study)
        res_pc[[length(res_pc) + 1]] <- data.table(study, peak = pk, chr = ch, caqtl_lead = lead$variant,
          caqtl_lead_p = lead$pvalue, n_snps = r$n, gwas_min_p = r$gmin, PP.H3 = r$pp3, PP.H4 = r$pp4)
      }
    }
    src <- if (study == "UKBB_GGT") eur[gwas == "UKBB_GGT"] else plt
    hits <- unique(src$hit1)
    hits <- hits[sub(":.*", "", hits) == ch]
    for (h in hits) {
      p19 <- as.integer(sub(".*:", "", h)); p38 <- lift(ch, p19)
      if (is.na(p38)) { res_loci[[length(res_loci) + 1]] <- data.table(study, hit1 = h, lead38 = NA, n_peaks = 0L, best_peak = NA, best_pp4 = NA_real_, n_snps = NA_integer_, lifted = FALSE); next }
      cand <- peaks[abs(centre - p38) <= WIN]
      best <- NA_real_; bp <- NA_character_; bn <- NA_integer_
      gr <- g[abs(pos38 - p38) <= WIN]
      for (pk in cand$peak) {
        qq <- q[.(pk), nomatch = 0L]
        r <- run_coloc(gr, qq[abs(pos38 - p38) <= WIN], study)
        if (!is.na(r$pp4) && (is.na(best) || r$pp4 > best)) { best <- r$pp4; bp <- pk; bn <- r$n }
      }
      res_loci[[length(res_loci) + 1]] <- data.table(study, hit1 = h, lead38 = p38, n_peaks = nrow(cand),
        best_peak = bp, best_pp4 = best, n_snps = bn, lifted = TRUE)
    }
    rm(g); gc()
  }
  rm(q); gc()
  cat("chr", ch, "done", format(Sys.time() - t0), "\n")
}
pc <- rbindlist(res_pc); loci <- rbindlist(res_loci)
fwrite(pc, file.path(JOB, "positive_control_coloc.tsv"), sep = "\t")
fwrite(loci, file.path(JOB, "pilot_locus_max_pp4.tsv"), sep = "\t")
cat("\nPOSITIVE CONTROL\n"); print(pc)
cat("\nPILOT (chr1+chr19)\n")
print(loci[, .(n_leads = .N, n_lifted = sum(lifted), n_with_peak = sum(n_peaks > 0),
               frac_pp4_gt_0.5 = sum(best_pp4 > 0.5, na.rm = TRUE) / .N,
               n_pp4_gt_0.5 = sum(best_pp4 > 0.5, na.rm = TRUE),
               n_pp4_gt_0.8 = sum(best_pp4 > 0.8, na.rm = TRUE),
               median_peaks = median(n_peaks)), by = study])
writeLines(capture.output(sessionInfo()), file.path(JOB, "sessionInfo.txt"))
