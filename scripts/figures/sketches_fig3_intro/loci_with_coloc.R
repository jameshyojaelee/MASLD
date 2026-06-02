#!/usr/bin/env Rscript
# Count # GWAS loci (lead-SNP windows) with at least one COLOC-positive gene
# at PP4>=0.5 — separately for SuSiE COLOC and ABF (regular) COLOC. Canonical
# panel choice (PolyFun for UKBB EUR, 1KG otherwise).
suppressPackageStartupMessages({ library(data.table) })

WIN <- 0.5e6
ukbb_studies <- c("UKBB_ALT","UKBB_AST","UKBB_GGT",
                  "2023_36280732_NAFLD_UKBB_EUR",
                  "2021_34128465_PDFF_EUR",
                  "2021_34957434_PDFF_EUR",
                  "2022_36402844_PDFF_EUR")
reg <- fread("GWAS/finemapping/config/gwas_registry.tsv",
             select = c("study_name","ancestry","leadsnps_path"))
setnames(reg, "study_name", "study")

# ---- Build per-study lead SNP table ----
lead_list <- list()
for (i in seq_len(nrow(reg))) {
  s <- reg$study[i]
  p <- file.path("GWAS/finemapping", reg$leadsnps_path[i])
  if (!file.exists(p)) next
  ld <- fread(p)
  chr_col <- intersect(c("chr","chromosome","CHR"), names(ld))[1]
  pos_col <- intersect(c("pos","position","BP","POS","bp"), names(ld))[1]
  if (is.na(chr_col)||is.na(pos_col)) next
  ld[, lead_chr := as.integer(sub("^chr","", get(chr_col)))]
  ld[, lead_pos := as.integer(get(pos_col))]
  ld <- ld[!is.na(lead_chr)&!is.na(lead_pos),
           .(study=s, lead_chr, lead_pos,
             lead_id = paste0(s, "_", lead_chr, ":", lead_pos))]
  lead_list[[s]] <- ld
}
leads <- rbindlist(lead_list)

# ---- Read per-study COLOC results ----
read_per_study <- function(root) {
  out <- list()
  for (sd in list.dirs(root, recursive = FALSE, full.names = FALSE)) {
    files <- list.files(file.path(root, sd), pattern = "csv$", full.names = TRUE)
    if (!length(files)) next
    d <- rbindlist(lapply(files, fread), fill = TRUE, use.names = TRUE)
    keep <- intersect(c("gene","gwas_name","chr","top_snp",
                        "PP.H4.susie","PP.H4.abf","method"), names(d))
    out[[sd]] <- d[, ..keep][, study := sd]
  }
  rbindlist(out, fill = TRUE)
}

poly  <- read_per_study("GWAS/finemapping/results/susie_coloc_polyfun")
poly  <- poly[study %in% ukbb_studies]
non_ukbb_eur <- reg[ancestry == "EUR" & !study %in% ukbb_studies, study]
kg_eur <- read_per_study("GWAS/finemapping/results/susie_coloc_1kg")
kg_eur <- kg_eur[study %in% non_ukbb_eur]
kg_other <- read_per_study("GWAS/finemapping/results/susie_coloc")
kg_other <- kg_other[study %in% reg[ancestry %in% c("EAS","AFR","SAS"), study]]

all_dt <- rbindlist(list(poly, kg_eur, kg_other), fill = TRUE)

# Parse top_snp into chr:pos
all_dt[, top_pos := as.integer(sub("^[0-9XY]+:", "", top_snp))]
all_dt[, top_chr := as.integer(sub(":.*", "", top_snp))]

# ---- Map each gene's top_snp to its enclosing lead-SNP window ----
attach_lead <- function(dt, leads_dt) {
  out <- list()
  for (s in unique(dt$study)) {
    a <- dt[study == s & !is.na(top_chr) & !is.na(top_pos)]
    b <- leads_dt[study == s]
    if (!nrow(a) || !nrow(b)) next
    matched <- character(nrow(a))
    for (i in seq_len(nrow(a))) {
      hit <- which(b$lead_chr == a$top_chr[i] &
                   abs(b$lead_pos - a$top_pos[i]) <= WIN)
      if (length(hit)) {
        nearest <- hit[which.min(abs(b$lead_pos[hit] - a$top_pos[i]))]
        matched[i] <- b$lead_id[nearest]
      } else matched[i] <- NA_character_
    }
    a[, lead_id := matched]
    out[[s]] <- a
  }
  rbindlist(out, fill = TRUE)
}

mapped <- attach_lead(all_dt, leads)

# ---- Count loci with >=1 COLOC-positive gene per method ----
n_loci_susie <- mapped[!is.na(lead_id) & !is.na(PP.H4.susie) & PP.H4.susie >= 0.5,
                        uniqueN(lead_id)]
n_loci_abf   <- mapped[!is.na(lead_id) & !is.na(PP.H4.abf)   & PP.H4.abf   >= 0.5,
                        uniqueN(lead_id)]
n_loci_input <- mapped[!is.na(lead_id), uniqueN(lead_id)]   # any locus tested

cat("Total loci tested by COLOC (any gene):", n_loci_input, "\n")
cat("Loci with >=1 SuSiE-COLOC gene PP4>=0.5:", n_loci_susie, "\n")
cat("Loci with >=1 ABF-COLOC gene PP4>=0.5:  ", n_loci_abf, "\n")

# Per-ancestry rollup
mapped <- merge(mapped, reg[, .(study, ancestry)], by = "study", all.x = TRUE)
agg <- mapped[!is.na(lead_id),
  .(loci_input  = uniqueN(lead_id),
    loci_abf    = uniqueN(lead_id[!is.na(PP.H4.abf)   & PP.H4.abf   >= 0.5]),
    loci_susie  = uniqueN(lead_id[!is.na(PP.H4.susie) & PP.H4.susie >= 0.5])),
  by = ancestry]
cat("\nPer-ancestry breakdown:\n"); print(agg)
fwrite(agg, "scripts/figures/sketches_fig3_intro/data/loci_with_coloc_per_ancestry.csv")

# Loci with any COLOC test run (no PP4 filter) — denominators
n_loci_abf_any   <- mapped[!is.na(lead_id) & !is.na(PP.H4.abf),   uniqueN(lead_id)]
n_loci_susie_any <- mapped[!is.na(lead_id) & !is.na(PP.H4.susie), uniqueN(lead_id)]

# Cascade with NO PP4 threshold on the COLOC tier — every locus where the
# method produced a COLOC test (any cis gene, any PP4).
cascade <- data.table(
  step = c(
    "Lead loci screened",
    "SuSiE converged loci (canonical)",
    "SuSiE-X high-PIP physical loci (cross-ancestry)",
    "meSuSiE shared CS physical loci (cross-ancestry)",
    "Loci — regular (ABF) COLOC tested",
    "Loci — SuSiE COLOC tested"),
  count = c(308L, 272L, 91L, 110L, n_loci_abf_any, n_loci_susie_any))
fwrite(cascade, "scripts/figures/sketches_fig3_intro/data/method_cascade_loci.csv")

cat("\nLoci cascade:\n"); print(cascade)
