#!/usr/bin/env Rscript
# 00c_prep_coloc_targeted_leads.R
# Build per-GWAS lead-SNP files for the targeted COLOC fine-mapping re-run.
#
# For each of the 22 target GWAS:
#   1. Take all (gene, gwas) COLOC pairs with PP.H4 >= 0.5
#      (using pp4 = max(PP.H4.susie, PP.H4.abf)).
#   2. Drop pairs whose top_snp is already in combined_finemapping.csv.
#   3. Look up the GWAS p-value at each top_snp.
#   4. Sanity gate: require p <= 1e-5 (loci weaker than this are too
#      underpowered to fine-map informatively even with COLOC support).
#   5. Cluster within 500 kb per chromosome; keep the lowest-p variant per
#      cluster as the locus lead.
#   6. Write data/lead_snps_coloc_targeted/<gwas>_leadSNPs.tsv (CHR, BP, locus).

suppressPackageStartupMessages(library(data.table))

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE, "GWAS/finemapping")
OUT_DIR <- file.path(FM_DIR, "data/lead_snps_coloc_targeted")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

EUR_17 <- c(
  "2019_31311600_NAFLD_EUR", "2020_32298765_NAFLD_EUR",
  "2021_34128465_PDFF_EUR",  "2021_34841290_NAFLD_EUR",
  "2021_34957434_PDFF_EUR",  "2022_36402844_PDFF_EUR",
  "2023_36280732_NAFLD_deCode_EUR",
  "2023_36280732_NAFLD_Intermountain_EUR",
  "2023_36280732_NAFLD_UKBB_EUR",
  "FinnGen_HCC", "FinnGen_NAFLD", "FinnGen_NASH",
  "Ghouse_Cirrhosis", "Ghouse_HCC",
  "UKBB_ALT", "UKBB_AST", "UKBB_GGT")
BBJ_5  <- c("2020_32514122_Cirrhosis_EAS", "2020_32514122_HCC_EAS",
            "BBJ_ALT", "BBJ_AST", "BBJ_GGT")
TARGET <- c(EUR_17, BBJ_5)
P_GATE <- 1e-5

gwas_registry <- fread(file.path(FM_DIR, "config/gwas_registry.tsv"))
sumstats_map  <- gwas_registry[study_name %in% TARGET,
                                .(study_name, sumstats_path)]

cat("[prep] Loading COLOC results ...\n")
sc <- fread(file.path(FM_DIR,
            "results/susie_coloc/susie_coloc_all_gwas.csv"),
            select = c("gene", "gwas_name", "top_snp",
                       "PP.H4.susie", "PP.H4.abf"))
sc <- sc[gwas_name %in% TARGET & !is.na(top_snp)]
sc[, pp4 := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
sc[, snp_chr := as.integer(tstrsplit(top_snp, ":", fixed = TRUE)[[1]])]
sc[, snp_pos := as.integer(tstrsplit(top_snp, ":", fixed = TRUE)[[2]])]
sc <- sc[pp4 >= 0.5 & !is.na(snp_chr) & !is.na(snp_pos)]
cat(sprintf("  %d COLOC-positive (pp4 >= 0.5) rows across %d GWAS\n",
            nrow(sc), uniqueN(sc$gwas_name)))

cat("[prep] Loading existing fine-mapping coverage ...\n")
fmap <- fread(file.path(FM_DIR,
              "results/combined_finemapping.csv"),
              select = c("study", "chromosome", "position"))
fmap <- fmap[study %in% TARGET & chromosome %in% 1:22]
fmap_keys <- unique(paste(fmap$study, fmap$chromosome, fmap$position, sep = "_"))
sc[, key := paste(gwas_name, snp_chr, snp_pos, sep = "_")]
sc <- sc[!(key %in% fmap_keys)]
cat(sprintf("  After dropping already-finemapped variants: %d rows\n", nrow(sc)))

# --- per-GWAS lead-SNP file generation ---
total_leads        <- 0L
total_dropped_p    <- 0L
total_dropped_miss <- 0L
summary_rows       <- list()

for (g in sort(unique(sc$gwas_name))) {
  spath <- sumstats_map[study_name == g, sumstats_path][1]
  if (is.na(spath) || !file.exists(file.path(FM_DIR, spath))) {
    cat(sprintf("\n[%s] sumstats not found -- skipping\n", g))
    next
  }
  cat(sprintf("\n[%s] reading %s ...\n", g, basename(spath)))
  ss <- fread(file.path(FM_DIR, spath),
              select = c("chromosome", "position", "pval"))
  # Coerce types (sumstats files vary: some have chr as character "X","Y";
  # some have pval read as character because of non-numeric outlier rows).
  ss[, chromosome := suppressWarnings(as.integer(as.character(chromosome)))]
  ss[, position   := as.integer(position)]
  ss[, pval       := suppressWarnings(as.numeric(pval))]
  ss <- ss[!is.na(chromosome) & chromosome %in% 1:22 & !is.na(pval)]
  setkey(ss, chromosome, position)

  cand <- sc[gwas_name == g, .(gene, snp_chr, snp_pos, pp4)]
  cand <- merge(cand, ss,
                by.x = c("snp_chr", "snp_pos"),
                by.y = c("chromosome", "position"),
                all.x = TRUE)
  n_total       <- nrow(cand)
  n_miss        <- sum(is.na(cand$pval))
  cand          <- cand[!is.na(pval)]
  n_after_miss  <- nrow(cand)
  cand          <- cand[pval <= P_GATE]
  n_pass        <- nrow(cand)
  n_drop_p      <- n_after_miss - n_pass

  total_dropped_miss <- total_dropped_miss + n_miss
  total_dropped_p    <- total_dropped_p    + n_drop_p

  cat(sprintf("  candidates=%d  missing-from-sumstats=%d  p>1e-5=%d  passing=%d\n",
              n_total, n_miss, n_drop_p, n_pass))

  if (n_pass == 0) {
    cat(sprintf("  [%s] 0 leads after gate -- no file written\n", g))
    summary_rows[[g]] <- data.table(study = g, n_candidates = n_total,
                                     n_missing = n_miss, n_drop_p = n_drop_p,
                                     n_leads_after_cluster = 0L)
    next
  }

  # Cluster within 500 kb per chromosome; keep lowest-p variant per cluster
  cand <- cand[order(snp_chr, snp_pos)]
  cand[, gap        := snp_pos - shift(snp_pos, 1L), by = snp_chr]
  cand[, new_locus  := is.na(gap) | gap > 500000]
  cand[, locus_id   := cumsum(new_locus), by = snp_chr]
  leads <- cand[, .SD[which.min(pval)], by = .(snp_chr, locus_id)]
  leads_out <- leads[, .(CHR = snp_chr,
                          BP  = snp_pos,
                          locus = paste(snp_chr, snp_pos, sep = "."))]

  out_path <- file.path(OUT_DIR, paste0(g, "_leadSNPs.tsv"))
  fwrite(leads_out, out_path, sep = "\t")
  total_leads <- total_leads + nrow(leads_out)
  cat(sprintf("  [%s] wrote %d leads -> %s\n",
              g, nrow(leads_out), basename(out_path)))

  summary_rows[[g]] <- data.table(study = g, n_candidates = n_total,
                                   n_missing = n_miss, n_drop_p = n_drop_p,
                                   n_leads_after_cluster = nrow(leads_out))
}

cat("\n========= SUMMARY =========\n")
summary_dt <- rbindlist(summary_rows)
print(summary_dt[order(-n_leads_after_cluster)])
cat(sprintf("\nTotal: %d unique loci across %d GWAS\n",
            total_leads,
            sum(summary_dt$n_leads_after_cluster > 0)))
cat(sprintf("Dropped: %d candidates (missing in sumstats: %d, p > %.0e: %d)\n",
            total_dropped_miss + total_dropped_p,
            total_dropped_miss, P_GATE, total_dropped_p))

fwrite(summary_dt, file.path(OUT_DIR, "_summary.csv"))
cat("Wrote summary to", file.path(OUT_DIR, "_summary.csv"), "\n")
