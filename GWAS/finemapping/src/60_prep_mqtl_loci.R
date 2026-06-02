#!/usr/bin/env Rscript
# 60_prep_mqtl_loci.R
# ---------------------------------------------------------------------------
# Build the bounded MASLD credible-set locus table for the metabolite/lipid-QTL
# COLOC layer, plus an hg19 ENSG -> cis-window map (from Broadaway eQTLs) for
# many-to-many locus->gene assignment.
#
# COMPLEMENTARY (supporting) layer — NOT a primary genetic-causal tier.
#
# Gates (MASLD relevance):
#   (a) MASLD-specific GWAS only (FinnGen NAFLD/NASH/HCC, 2021 NAFLD EUR,
#       2023 deCODE NAFLD EUR, 2020 EAS Cirrhosis/HCC, Ghouse Cirrhosis/HCC).
#   (b) Credible-set loci only (either_in_cs==TRUE | recommended_pip>0.1).
#
# Output:
#   data/external/chen2023_mqtl/masld_credset_loci.tsv  (one row per MASLD locus)
#   data/external/chen2023_mqtl/broadaway_gene_hg19_window.tsv  (ENSG cis windows)
# ---------------------------------------------------------------------------
suppressMessages({library(data.table)})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR  <- file.path(BASE_DIR, "GWAS/finemapping")
OUT_DIR <- file.path(BASE_DIR, "data/external/chen2023_mqtl")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- MASLD-specific GWAS (per task spec / registry) ---
MASLD_GWAS <- c(
  "FinnGen_NAFLD", "FinnGen_NASH", "FinnGen_HCC",
  "2021_34841290_NAFLD_EUR", "2023_36280732_NAFLD_deCode_EUR",
  "2020_32514122_Cirrhosis_EAS", "2020_32514122_HCC_EAS",
  "Ghouse_Cirrhosis", "Ghouse_HCC"
)
LOCUS_WINDOW <- 5e5  # +/- 500kb around credible-set centre (matches registry window_mb=0.5)

cat("=== 60_prep_mqtl_loci ===\n")
cf <- fread(file.path(FM_DIR, "results/combined_finemapping.csv"))
g <- cf[study %in% MASLD_GWAS &
        (either_in_cs == TRUE | (!is.na(recommended_pip) & recommended_pip > 0.1)) &
        !is.na(chromosome) & !is.na(position)]
cat("Gated MASLD credible-set variants:", nrow(g), "\n")

# Collapse credible-set variants into loci: per (study, chr, locus) take min/max pos.
# Some studies share the same 'locus' label only within-study, so key on study+chr+locus.
g[, locus := as.character(locus)]
loci <- g[, .(
  cs_lead_pos   = position[which.max(ifelse(is.na(recommended_pip), -1, recommended_pip))],
  cs_start      = min(position),
  cs_end        = max(position),
  n_cs_variants = .N,
  max_pip       = max(recommended_pip, na.rm = TRUE)
), by = .(study, chromosome, locus)]
loci[, win_start := pmax(1, cs_start - LOCUS_WINDOW)]
loci[, win_end   := cs_end + LOCUS_WINDOW]
loci[, locus_id  := paste(study, chromosome, locus, sep = "__")]

# Merge GWAS metadata from registry (N, type, ancestry)
reg <- fread(file.path(FM_DIR, "config/gwas_registry.tsv"))
loci <- merge(loci, reg[, .(study = study_name, sumstats_path, trait_type, N_tot, N_cases, ancestry)],
              by = "study", all.x = TRUE)
loci[, coloc_type := ifelse(trait_type == "binary", "cc", "quant")]

cat("Distinct MASLD loci (study x chr x locus):", nrow(loci), "\n")
cat("Distinct genomic windows (chr + ~1Mb bin):",
    loci[, uniqueN(paste(chromosome, round(cs_lead_pos/1e6)))], "\n")
cat("Per-study loci:\n"); print(loci[, .N, by = study][order(-N)])

fwrite(loci, file.path(OUT_DIR, "masld_credset_loci.tsv"), sep = "\t")
cat("Wrote", file.path(OUT_DIR, "masld_credset_loci.tsv"), "\n\n")

# --- hg19 ENSG cis-window map from Broadaway eQTLs (for locus->gene) ---
# Only needed by the COMBINER (61b), not by the COLOC array. Skip with
# SKIP_GENE_WINDOW=1 if it already exists. This loop reads ~12GB of eQTL
# files (chr1 alone is 1.2GB) so it is the slow part of this prep.
GW_FILE <- file.path(OUT_DIR, "broadaway_gene_hg19_window.tsv")
if (Sys.getenv("SKIP_GENE_WINDOW", "0") == "1" && file.exists(GW_FILE)) {
  cat("SKIP_GENE_WINDOW=1 and", GW_FILE, "exists — skipping gene-window rebuild.\n")
} else {
  cat("Building hg19 ENSG cis-window map from Broadaway eQTLs (reads ~12GB; slow)...\n")
  EQTL_DIR <- file.path(BASE_DIR, "data/broadaway_eqtl")
  gene_pos <- vector("list", 22L)
  for (chr in 1:22) {
    f <- file.path(EQTL_DIR, paste0("chr", chr, "_marginal_summary_results.tsv"))
    if (!file.exists(f)) next
    e <- fread(f, select = c("ENSG", "GeneSymbol", "CHR", "POS"))
    gp <- e[, .(chr = CHR[1], cis_start = min(POS), cis_end = max(POS),
                cis_mid = as.integer(median(POS)), gene = GeneSymbol[1]), by = ENSG]
    gene_pos[[chr]] <- gp
    cat("  chr", chr, "done (", nrow(gp), "genes)\n")
  }
  gene_pos <- rbindlist(gene_pos)
  fwrite(gene_pos, GW_FILE, sep = "\t")
  cat("Wrote", nrow(gene_pos), "ENSG cis-windows (hg19) to broadaway_gene_hg19_window.tsv\n")
}
cat("=== done ===\n")
