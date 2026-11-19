#!/usr/bin/env Rscript
# MQTL2 step 01 (2026-11-18 remediation, N15): convert DIG1 out_release_vs_disease_loci_20261118/loci.tsv
# to the locus-table format read by 02_coloc_20261118.R, plus one control-only PNPLA3 window.
# Usage: Rscript 01_make_loci_20261118.R <out_file>
suppressMessages(library(data.table))
args <- commandArgs(TRUE)
out_file <- args[1]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WIN <- 5e5

L <- fread(file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z/out_release_vs_disease_loci_20261118/loci.tsv"),
           select = c("locus", "chr", "pos", "a1", "a2", "defining_enzyme", "enzyme_raising_allele", "alt_raising_allele", "genes", "cls"))
n_dig1 <- nrow(L)
stopifnot(n_dig1 > 100L, !anyDuplicated(L$locus))
L[, locus_id := as.character(locus)]
ctrl <- data.table(
  locus = NA_integer_, chr = 22L, pos = 44324727L, a1 = "G", a2 = "C",
  defining_enzyme = "MVP_ALT_EUR", enzyme_raising_allele = "G", alt_raising_allele = "G",
  genes = "PNPLA3", cls = "control_only", locus_id = "ctrl_PNPLA3_rs738409"
)
L <- rbind(L, ctrl)

loci <- L[, .(
  study = "PanUKBB_ALT_EUR", chromosome = chr, locus = locus_id,
  cs_lead_pos = pos, cs_start = pos, cs_end = pos, n_cs_variants = 1L, max_pip = NA_real_,
  win_start = pmax(1L, pos - WIN), win_end = pos + WIN, locus_id,
  sumstats_path = "GWAS/MR_Data/PanUKBB/biomarkers-30620-ALT-both_sexes-irnt.tsv.bgz",
  trait_type = "quantitative", N_tot = 400822L, N_cases = 0L, ancestry = "EUR", coloc_type = "quant",
  genes, cls_DIG1 = cls, defining_enzyme, lead_a1 = a1, lead_a2 = a2,
  enzyme_raising_allele, alt_raising_allele
)]
fwrite(loci, out_file, sep = "\t")
cat("wrote", nrow(loci), "loci (", n_dig1, "+ 1 control) to", out_file, "\n")
