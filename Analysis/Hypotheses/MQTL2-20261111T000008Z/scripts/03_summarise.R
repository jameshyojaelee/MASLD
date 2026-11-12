#!/usr/bin/env Rscript
# MQTL2 step 03: controls, primary Fisher test, sensitivity, descriptive per-locus tables. Rules fixed in
# ../PRESPEC_MQTL2.md. Unit = locus. Usage: Rscript 03_summarise.R <out_dir> <loci_file>
suppressMessages(library(data.table))
set.seed(20261110)
args <- commandArgs(TRUE); OUT <- args[1]; LOCI_FILE <- args[2]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
TH <- 0.8; TH_S <- 0.5
sink(file.path(OUT, "summary.txt"), split = TRUE)

loci <- fread(LOCI_FILE, colClasses = list(character = "locus_id"))
res <- fread(file.path(OUT, "coloc_all_tests.tsv.gz"), colClasses = list(character = "locus_id"))
chk <- fread(file.path(OUT, "file_checks.tsv"))
alt_sum <- fread(file.path(OUT, "alt_window_summary.tsv"), colClasses = list(character = "locus_id"))

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

# ---------------- controls ----------------
cat("\n== Controls ==\n")
self <- fread(file.path(OUT, "control_self_coloc_locus148.tsv"))
cat("Self-colocalization ALT vs ALT, locus 148: PP.H4 =", signif(self$PP.H4.abf, 4), "\n")
best <- function(d, by_cls) d[order(-PP.H4)][, .SD[1], by = by_cls]
ctrl_tab <- rbind(
  best(pu[locus_id == "17" & triglyceride == TRUE], "gwas")[, control := "GCKR (17) x triglyceride, Pan-UKBB ALT"],
  best(res[gwas == "MVP_ALT_EUR_control" & triglyceride == TRUE], "gwas")[, control := "GCKR (17) x triglyceride, MVP ALT"],
  best(pu[locus_id == "148" & class == "lipid"], "gwas")[, control := "TM6SF2 (148) x lipid"],
  best(pu[locus_id == "ctrl_PNPLA3_rs738409" & class == "lipid"], "gwas")[, control := "PNPLA3 window x lipid"],
  fill = TRUE)
ctrl_tab[, pass := PP.H4 >= TH]
print(ctrl_tab[, .(control, trait, n_snps, PP.H3, PP.H4, top_pos, pass)])
fwrite(ctrl_tab, file.path(OUT, "controls.tsv"), sep = "\t")
cat("Number of loci-level lipid colocalizations by control locus (PP.H4 >= 0.8):\n")
print(pu[locus_id %in% c("17", "148", "ctrl_PNPLA3_rs738409") & PP.H4 >= TH, .N, by = .(locus_id, class)])

# ---------------- exposure (2 columns only) ----------------
ex <- fread(file.path(ROOT, "Analysis/Hypotheses/ENZFAT1-20261110T204401Z/out_S2/locus_effects.tsv"),
            select = c("locus", "lfsr_pdff"), colClasses = list(character = "locus"))
ex <- ex[!grepl("^ctrl_", locus) & !is.na(lfsr_pdff)]
ex[, fat := lfsr_pdff < 0.05]
cat("\nExposure: loci with lfsr_pdff:", nrow(ex), "| positive:", sum(ex$fat), "(", paste(ex[fat == TRUE, locus], collapse = ","), ")\n")

# ---------------- locus-level outcome ----------------
locus_outcome <- function(d, col, th, cls = "lipid") {
  t <- d[class == cls, .(tested = .N > 0, best = max(get(col), na.rm = TRUE)), by = locus_id]
  t[, coloc := best >= th]
  t
}
primary_test <- function(lo, ex, label, drop = character()) {
  m <- merge(ex[!locus %in% drop], lo, by.x = "locus", by.y = "locus_id", all.x = TRUE)
  n_untest <- sum(is.na(m$tested))
  m <- m[!is.na(tested)]
  tab <- table(factor(m$fat, c(TRUE, FALSE), c("fat_effect", "no_fat_effect")),
               factor(m$coloc, c(TRUE, FALSE), c("lipid_coloc", "no_lipid_coloc")))
  ft <- fisher.test(tab)
  mtot <- sum(m$coloc); k <- sum(m$fat); n <- nrow(m)
  pmin1 <- if (mtot >= k) choose(mtot, k) / choose(n, k) else NA_real_
  # smallest attainable two-sided p given margins: most extreme table in the positive direction
  xmax <- min(k, mtot)
  tmax <- matrix(c(xmax, k - xmax, mtot - xmax, n - k - mtot + xmax), 2)
  pmin2 <- fisher.test(tmax)$p.value
  cat("\n--", label, "--\n"); print(tab)
  cat(sprintf("n loci = %d (untestable/excluded %d); OR = %.2f (95%% CI %.2f to %.2f); two-sided p = %.3g; smallest attainable two-sided p = %.3g\n",
              n, n_untest, ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value, pmin2))
  cat("Positive loci and colocalization:\n"); print(m[fat == TRUE, .(locus, lfsr_pdff, best_lipid_PP.H4 = round(best, 3), coloc)])
  data.table(analysis = label, n = n, n_untestable = n_untest, a = tab[1, 1], b = tab[1, 2], c = tab[2, 1], d = tab[2, 2],
             OR = unname(ft$estimate), CI_lo = ft$conf.int[1], CI_hi = ft$conf.int[2], p_two_sided = ft$p.value,
             min_attainable_p_two_sided = pmin2, min_attainable_p_one_sided = pmin1)
}
cat("\n== Primary: liver-fat effect vs lipid colocalization (PP.H4 >= 0.8, default priors) ==\n")
P <- rbind(
  primary_test(locus_outcome(pu, "PP.H4", TH), ex, "primary_PPH4_0.8"),
  primary_test(locus_outcome(pu, "PP.H4", TH_S), ex, "sens1_PPH4_0.5"),
  primary_test(locus_outcome(pu, "PP.H4_p12_5e6", TH), ex, "sens2_p12_5e-6"),
  primary_test(locus_outcome(pu, "PP.H4", TH), ex, "sens3_drop_17_148", drop = c("17", "148")))
fwrite(P, file.path(OUT, "primary_and_sensitivity.tsv"), sep = "\t")

# ---------------- descriptive per locus ----------------
hits <- pu[PP.H4 >= TH]
per_locus <- hits[, .(n_traits = .N, classes = paste(sort(unique(class)), collapse = ";"),
                      best_PP.H4 = max(PP.H4),
                      traits = paste(head(trait[order(-PP.H4)], 8), collapse = " | ")), by = locus_id]
per_class <- dcast(pu[, .(best = max(PP.H4)), by = .(locus_id, class)], locus_id ~ class, value.var = "best")
setnames(per_class, setdiff(names(per_class), "locus_id"), paste0("bestPPH4_", setdiff(names(per_class), "locus_id")))
D <- merge(loci[, .(locus_id, chromosome, cs_lead_pos, genes, cls_DIG1)], alt_sum, by = "locus_id", all.x = TRUE)
D <- merge(D, ex[, .(locus_id = locus, lfsr_pdff, fat)], by = "locus_id", all.x = TRUE)
D <- merge(D, per_locus, by = "locus_id", all.x = TRUE)
D <- merge(D, per_class, by = "locus_id", all.x = TRUE)
D[is.na(n_traits), n_traits := 0L]
setorder(D, -n_traits)
fwrite(D, file.path(OUT, "per_locus.tsv"), sep = "\t")
cat("\n== Descriptive ==\n")
cat("Loci tested (Pan-UKBB):", uniqueN(pu$locus_id), "| loci with >= 1 colocalized trait:", D[n_traits > 0, .N], "\n")
cat("Loci by colocalized class (PP.H4 >= 0.8):\n"); print(hits[, .(n_loci = uniqueN(locus_id)), by = class][order(-n_loci)])
print(D[n_traits > 0, .(locus_id, genes, lfsr_pdff = signif(lfsr_pdff, 2), n_traits, classes, best_PP.H4 = round(best_PP.H4, 3))])

# amino acid and pyruvate at loci without liver-fat effect; metabolite effect for ALT-raising allele at top SNP
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
