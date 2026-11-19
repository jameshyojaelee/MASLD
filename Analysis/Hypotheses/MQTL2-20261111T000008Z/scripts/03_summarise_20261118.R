#!/usr/bin/env Rscript
# MQTL2 step 03 (2026-11-18 remediation, N15):
# - Dynamic GCKR & TM6SF2 control lookups.
# - Reads ENZFAT1 out_S2_20261118/locus_effects.tsv (locus, lfsr_pdff, postmean_pdff).
# - Reports sdY_est distribution + all original and sensitivity Fisher tests:
#   primary_PPH4_0.8, sens1_PPH4_0.5, sens2_p12_5e-6, sens3_drop_GCKR_TM6SF2,
#   alt_gws_PPH4_0.8, alt_gws_PPH4_0.5,
#   steatosis_concordant_PPH4_0.8, steatosis_concordant_PPH4_0.5,
#   alt_gws_steatosis_concordant_PPH4_0.8, alt_gws_steatosis_concordant_PPH4_0.5.
# Usage: Rscript 03_summarise_20261118.R <out_dir> <loci_file>
suppressMessages(library(data.table))
set.seed(20261110)
args <- commandArgs(TRUE); OUT <- args[1]; LOCI_FILE <- args[2]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
TH <- 0.8; TH_S <- 0.5; ALT_GWS_NLP <- 7.301
sink(file.path(OUT, "summary.txt"), split = TRUE)

loci <- fread(LOCI_FILE, colClasses = list(character = "locus_id"))
res <- fread(file.path(OUT, "coloc_all_tests.tsv.gz"), colClasses = list(character = "locus_id"))
chk <- fread(file.path(OUT, "file_checks.tsv"))
alt_sum <- fread(file.path(OUT, "alt_window_summary.tsv"), colClasses = list(character = "locus_id"))

gckr_id <- loci[grepl("(^|,)GCKR(,|$)", genes), locus_id]
tm6_id <- loci[chromosome == 19L & cs_lead_pos == 19379549L, locus_id]
if (length(tm6_id) == 0L) tm6_id <- loci[grepl("(^|,)(TM6SF2|MAU2)(,|$)", genes), locus_id]
stopifnot(length(gckr_id) == 1L, length(tm6_id) == 1L)
cat("Dynamic control locus IDs: GCKR =", gckr_id, "| TM6SF2 =", tm6_id, "\n")

# ---------------- trait classes (by name, fixed in prespec) ----------------
cls_of <- function(src, t) {
  tl <- tolower(t)
  fifelse(grepl(" to .* ratio", tl), "ratio",
  fifelse(src == "ottensmann", "lipid",
  fifelse(grepl("^(ceramide|glycosyl ceramide) ", tl), "lipid",
  fifelse(tl %in% c("alanine levels", "glutamate levels", "glutamine levels", "tyrosine levels",
                    "leucine levels", "isoleucine levels", "valine levels"), "amino_acid",
  fifelse(tl == "pyruvate levels", "pyruvate",
  fifelse(grepl("^(n-acetyl|n-lactoyl|gamma-glutamyl)", tl), "aa_derivative",
  fifelse(grepl("carnitine", tl), "acylcarnitine",
  fifelse(grepl("cholate|cholic|muricholate|chenodeoxychol|lithochol|ursodeoxychol|hyochol|cholenate", tl),
          "bile_acid", "other"))))))))
}
res[, class := cls_of(source, trait)]
res[, triglyceride := source == "ottensmann" & grepl("^(triacylglycerol|diacylglycerol)", tolower(trait))]
cat("Trait classes (traits):\n"); print(unique(res[, .(accession, class)])[, .N, by = class])

# ---------------- checks: build and alleles ----------------
cat("\n== File checks ==\n"); print(chk[, .N, by = status])
g <- chk[snp == "rs1260326"]
bad_allele <- g[n_other_alleles_at_pos > 0 & n_match == 0]
cat("GCKR rs1260326 present with matching {C,T} alleles in", g[n_match == 1, .N], "files; files with other alleles only:",
    nrow(bad_allele), "\n")
if (nrow(bad_allele) > 0) { print(bad_allele); stop("Allele alignment failed at GCKR; stopping before any test") }
cat("Ottensmann rsid positions (GRCh38) seen:\n"); print(chk[!is.na(rsid_pos), .N, by = .(snp, rsid_pos)])
cat("Chen rows at GRCh38 control positions:\n"); print(chk[, .N, by = .(snp, alleles_at_pos38)])
chk[, eaf_diff := abs(eaf_of_panukbb_alt - panukbb_af_alt)]
bad_eaf <- unique(chk[!is.na(eaf_diff) & eaf_diff > 0.1, accession])
cat("Files failing effect-allele frequency check (|diff| > 0.1):", length(bad_eaf), "\n")
if (length(bad_eaf)) print(chk[accession %in% bad_eaf, .(accession, snp, met_ea, eaf_of_panukbb_alt, panukbb_af_alt)])
cat("EAF difference summary by SNP:\n"); print(chk[, .(n = sum(!is.na(eaf_diff)), median = median(eaf_diff, na.rm = TRUE),
                                                       max = max(eaf_diff, na.rm = TRUE)), by = snp])
res <- res[!accession %in% bad_eaf]
fwrite(chk, file.path(OUT, "file_checks_annotated.tsv"), sep = "\t")

pu <- res[gwas == "PanUKBB_ALT_EUR"]
cat("\nPan-UKBB tests:", nrow(pu), "| traits:", uniqueN(pu$accession), "| loci:", uniqueN(pu$locus_id), "\n")

# ---------------- sdY_est distribution ----------------
cat("\n== Empirical sdY_est distribution (Pan-UKBB pairs) ==\n")
sdy_summ <- rbind(
  pu[, .(group = "ALL", n = .N, n_used_sdY_est = sum(used_sdY_est, na.rm = TRUE),
         min = min(sdY_est, na.rm = TRUE), q25 = quantile(sdY_est, 0.25, na.rm = TRUE),
         median = median(sdY_est, na.rm = TRUE), q75 = quantile(sdY_est, 0.75, na.rm = TRUE),
         max = max(sdY_est, na.rm = TRUE))],
  pu[, .(n = .N, n_used_sdY_est = sum(used_sdY_est, na.rm = TRUE),
         min = min(sdY_est, na.rm = TRUE), q25 = quantile(sdY_est, 0.25, na.rm = TRUE),
         median = median(sdY_est, na.rm = TRUE), q75 = quantile(sdY_est, 0.75, na.rm = TRUE),
         max = max(sdY_est, na.rm = TRUE)), by = .(group = paste0("source:", source))],
  pu[, .(n = .N, n_used_sdY_est = sum(used_sdY_est, na.rm = TRUE),
         min = min(sdY_est, na.rm = TRUE), q25 = quantile(sdY_est, 0.25, na.rm = TRUE),
         median = median(sdY_est, na.rm = TRUE), q75 = quantile(sdY_est, 0.75, na.rm = TRUE),
         max = max(sdY_est, na.rm = TRUE)), by = .(group = paste0("maf_source:", maf_source))]
)
print(sdy_summ)
fwrite(sdy_summ, file.path(OUT, "sdY_est_summary.tsv"), sep = "\t")

# ---------------- controls ----------------
cat("\n== Controls ==\n")
self <- fread(file.path(OUT, "control_self_coloc_TM6SF2.tsv"))
cat("Self-colocalization ALT vs ALT, TM6SF2 locus", tm6_id, ": PP.H4 =", signif(self$PP.H4.abf, 4), "\n")
best <- function(d, by_cls) d[order(-PP.H4)][, .SD[1], by = by_cls]
ctrl_tab <- rbind(
  best(pu[locus_id == gckr_id & triglyceride == TRUE], "gwas")[, control := sprintf("GCKR (%s) x triglyceride, Pan-UKBB ALT", gckr_id)],
  best(res[gwas == "MVP_ALT_EUR_control" & triglyceride == TRUE], "gwas")[, control := sprintf("GCKR (%s) x triglyceride, MVP ALT", gckr_id)],
  best(pu[locus_id == tm6_id & class == "lipid"], "gwas")[, control := sprintf("TM6SF2 (%s) x lipid", tm6_id)],
  best(pu[locus_id == "ctrl_PNPLA3_rs738409" & class == "lipid"], "gwas")[, control := "PNPLA3 window x lipid"],
  fill = TRUE)
ctrl_tab[, pass := PP.H4 >= TH]
print(ctrl_tab[, .(control, trait, n_snps, sdY_est = round(sdY_est, 3), PP.H3 = round(PP.H3, 4), PP.H4 = round(PP.H4, 4), top_pos, pass)])
fwrite(ctrl_tab, file.path(OUT, "controls.tsv"), sep = "\t")
cat("Number of loci-level lipid colocalizations by control locus (PP.H4 >= 0.8):\n")
print(pu[locus_id %in% c(gckr_id, tm6_id, "ctrl_PNPLA3_rs738409") & PP.H4 >= TH, .N, by = .(locus_id, class)])

# ---------------- exposure (ENZFAT1 out_S2_20261118) ----------------
ex <- fread(file.path(ROOT, "Analysis/Hypotheses/ENZFAT1-20261110T204401Z/out_S2_20261118/locus_effects.tsv"),
            select = c("locus", "genes", "lfsr_pdff", "postmean_pdff"), colClasses = list(character = "locus"))
ex <- ex[!grepl("^ctrl_", locus) & is.finite(lfsr_pdff)]
ex <- merge(ex, alt_sum[, .(locus = locus_id, alt_max_nlp)], by = "locus", all.x = TRUE)
ex[, fat := lfsr_pdff < 0.05]
ex[, fat_concordant := lfsr_pdff < 0.05 & postmean_pdff > 0]
ex[, alt_gws := is.finite(alt_max_nlp) & alt_max_nlp >= ALT_GWS_NLP]
cat("\nExposure: loci with lfsr_pdff:", nrow(ex),
    "| fat (lfsr<0.05):", sum(ex$fat),
    "| fat_concordant (lfsr<0.05 & pm>0):", sum(ex$fat_concordant),
    "| reversed fat (lfsr<0.05 & pm<0):", sum(ex$fat & !ex$fat_concordant),
    "| alt_gws (max_nlp >= 7.301):", sum(ex$alt_gws), "\n")
print(ex[fat == TRUE, .(locus, genes, lfsr_pdff = signif(lfsr_pdff, 3), postmean_pdff = signif(postmean_pdff, 3), alt_max_nlp = round(alt_max_nlp, 2), fat_concordant, alt_gws)])

# ---------------- locus-level outcome & Fisher tests ----------------
locus_outcome <- function(d, col, th, cls = "lipid") {
  t <- d[class == cls, .(tested = .N > 0, best = max(get(col), na.rm = TRUE)), by = locus_id]
  t[, coloc := best >= th]
  t
}
primary_test <- function(lo, ex_sub, label, exp_col = "fat", drop = character()) {
  m <- merge(ex_sub[!locus %in% drop], lo, by.x = "locus", by.y = "locus_id", all.x = TRUE)
  n_untest <- sum(is.na(m$tested))
  m <- m[!is.na(tested)]
  ev <- m[[exp_col]]
  tab <- table(factor(ev, c(TRUE, FALSE), c("fat_effect", "no_fat_effect")),
               factor(m$coloc, c(TRUE, FALSE), c("lipid_coloc", "no_lipid_coloc")))
  ft <- fisher.test(tab)
  mtot <- sum(m$coloc); k <- sum(ev); n <- nrow(m)
  pmin1 <- if (mtot >= k) choose(mtot, k) / choose(n, k) else NA_real_
  xmax <- min(k, mtot)
  tmax <- matrix(c(xmax, k - xmax, mtot - xmax, n - k - mtot + xmax), 2)
  pmin2 <- fisher.test(tmax)$p.value
  cat("\n--", label, "--\n"); print(tab)
  cat(sprintf("n loci = %d (untestable/excluded %d); OR = %.2f (95%% CI %.2f to %.2f); two-sided p = %.3g; smallest attainable two-sided p = %.3g\n",
              n, n_untest, ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value, pmin2))
  cat("Positive loci and colocalization:\n")
  print(m[ev == TRUE, .(locus, genes, lfsr_pdff = signif(lfsr_pdff, 3), postmean_pdff = signif(postmean_pdff, 3),
                        alt_max_nlp = round(alt_max_nlp, 2), best_lipid_PP.H4 = round(best, 3), coloc)])
  data.table(analysis = label, n = n, n_untestable = n_untest, a = tab[1, 1], b = tab[1, 2], c = tab[2, 1], d = tab[2, 2],
             OR = unname(ft$estimate), CI_lo = ft$conf.int[1], CI_hi = ft$conf.int[2], p_two_sided = ft$p.value,
             min_attainable_p_two_sided = pmin2, min_attainable_p_one_sided = pmin1)
}

cat("\n== Primary and Sensitivity Fisher Tests ==\n")
lo_08 <- locus_outcome(pu, "PP.H4", TH)
lo_05 <- locus_outcome(pu, "PP.H4", TH_S)
lo_p12 <- locus_outcome(pu, "PP.H4_p12_5e6", TH)
ex_gws <- ex[alt_gws == TRUE]

P <- rbind(
  primary_test(lo_08, ex, "primary_PPH4_0.8", "fat"),
  primary_test(lo_05, ex, "sens1_PPH4_0.5", "fat"),
  primary_test(lo_p12, ex, "sens2_p12_5e-6", "fat"),
  primary_test(lo_08, ex, "sens3_drop_GCKR_TM6SF2", "fat", drop = c(gckr_id, tm6_id)),
  primary_test(lo_08, ex_gws, "alt_gws_PPH4_0.8", "fat"),
  primary_test(lo_05, ex_gws, "alt_gws_PPH4_0.5", "fat"),
  primary_test(lo_08, ex, "steatosis_concordant_PPH4_0.8", "fat_concordant"),
  primary_test(lo_05, ex, "steatosis_concordant_PPH4_0.5", "fat_concordant"),
  primary_test(lo_08, ex_gws, "alt_gws_steatosis_concordant_PPH4_0.8", "fat_concordant"),
  primary_test(lo_05, ex_gws, "alt_gws_steatosis_concordant_PPH4_0.5", "fat_concordant")
)
fwrite(P, file.path(OUT, "primary_and_sensitivity.tsv"), sep = "\t")
cat("\n== Summary Fisher table ==\n")
print(P[, .(analysis, n, a, b, c, d, OR = round(OR, 2), CI_lo = round(CI_lo, 2), CI_hi = round(CI_hi, 2), p_two_sided = signif(p_two_sided, 3))])

# ---------------- descriptive per locus ----------------
hits <- pu[PP.H4 >= TH]
per_locus <- hits[, .(n_traits = .N, classes = paste(sort(unique(class)), collapse = ";"),
                      best_PP.H4 = max(PP.H4),
                      traits = paste(head(trait[order(-PP.H4)], 8), collapse = " | ")), by = locus_id]
per_class <- dcast(pu[, .(best = max(PP.H4)), by = .(locus_id, class)], locus_id ~ class, value.var = "best")
setnames(per_class, setdiff(names(per_class), "locus_id"), paste0("bestPPH4_", setdiff(names(per_class), "locus_id")))
D <- merge(loci[, .(locus_id, chromosome, cs_lead_pos, genes, cls_DIG1, defining_enzyme)], alt_sum, by = "locus_id", all.x = TRUE)
D <- merge(D, ex[, .(locus_id = locus, lfsr_pdff, postmean_pdff, fat, fat_concordant, alt_gws)], by = "locus_id", all.x = TRUE)
D <- merge(D, per_locus, by = "locus_id", all.x = TRUE)
D <- merge(D, per_class, by = "locus_id", all.x = TRUE)
D[is.na(n_traits), n_traits := 0L]
setorder(D, -n_traits)
fwrite(D, file.path(OUT, "per_locus.tsv"), sep = "\t")
cat("\n== Descriptive ==\n")
cat("Loci tested (Pan-UKBB):", uniqueN(pu$locus_id), "| loci with >= 1 colocalized trait:", D[n_traits > 0, .N], "\n")
cat("Loci by colocalized class (PP.H4 >= 0.8):\n"); print(hits[, .(n_loci = uniqueN(locus_id)), by = class][order(-n_loci)])
print(D[n_traits > 0, .(locus_id, genes, lfsr_pdff = signif(lfsr_pdff, 2), n_traits, classes, best_PP.H4 = round(best_PP.H4, 3))])

aa <- pu[class %in% c("amino_acid", "pyruvate") & PP.H4 >= TH]
aa <- merge(aa, ex[, .(locus_id = locus, lfsr_pdff, fat)], by = "locus_id", all.x = TRUE)
aa <- merge(aa, loci[, .(locus_id, genes)], by = "locus_id")
aa[, met_beta_alt_raising := top_met_beta * sign(top_alt_beta)]
cat("\nAmino-acid / pyruvate colocalizations (all loci; fat = liver-fat effect):\n")
print(aa[, .(locus_id, genes, fat, trait, PP.H4 = round(PP.H4, 3), top_pos, met_beta_alt_raising = signif(met_beta_alt_raising, 3))])
fwrite(aa, file.path(OUT, "amino_acid_colocalizations.tsv"), sep = "\t")
cat("Amino-acid traits tested per locus (median):", pu[class == "amino_acid", .N, by = locus_id][, median(N)], "\n")
cat("\nPairs with PP.H4 >= 0.8:", nrow(hits), "of", nrow(pu), "| PP.H3 >= 0.8 (distinct variants):", pu[PP.H3 >= TH, .N], "\n")
sink()
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
