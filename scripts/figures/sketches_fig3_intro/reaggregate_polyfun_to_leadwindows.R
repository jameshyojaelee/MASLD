#!/usr/bin/env Rscript
# Re-aggregate PolyFun per-gene SuSiE results to lead-SNP windows so the EUR
# locus count is comparable to UKBB sghatan / 1KG / TOP-LD.
#
# PolyFun pipeline runs SuSiE per gene's cis-window (gene-centric); other
# panels run SuSiE per GWAS lead-SNP region (locus-centric, ~1Mb window).
# To compare apples-to-apples, map each PolyFun-converged gene's `top_snp` to
# its enclosing lead-SNP window and count unique lead SNPs with >=1 converged
# gene per study.
suppressPackageStartupMessages({ library(data.table) })

BASE <- "GWAS/finemapping"
WIN  <- 0.5e6  # +/- 0.5 Mb (matches registry window_mb=0.5)

reg <- fread(file.path(BASE, "config/gwas_registry.tsv"),
             select = c("study_name","ancestry","leadsnps_path"))
setnames(reg, "study_name", "study")

# ---- Build lead-SNP table per study ----
lead_list <- list()
for (i in seq_len(nrow(reg))) {
  s    <- reg$study[i]
  path <- file.path(BASE, reg$leadsnps_path[i])
  if (!file.exists(path)) next
  ld <- fread(path)
  chr_col <- intersect(c("chr","chromosome","CHR","chrom"), names(ld))[1]
  pos_col <- intersect(c("pos","position","BP","POS","bp","start"), names(ld))[1]
  if (is.na(chr_col) || is.na(pos_col)) next
  ld[, chr := as.integer(sub("^chr","", get(chr_col)))]
  ld[, pos := as.integer(get(pos_col))]
  ld <- ld[!is.na(chr) & !is.na(pos), .(study = s, chr, pos)]
  ld[, lead_id := paste0(chr, ":", pos)]
  lead_list[[s]] <- ld
}
leads <- rbindlist(lead_list)
cat("Lead-SNP table:", nrow(leads), "rows across", uniqueN(leads$study), "studies\n")

# ---- Load PolyFun per-gene COLOC ----
poly <- fread(file.path(BASE,
  "results/susie_coloc_polyfun_bdiag_2026-04-28/susie_coloc_all_gwas_polyfun.csv"),
  select = c("gwas_name","chr","top_snp","method"))
setnames(poly, "gwas_name", "study")
poly[, top_pos := as.integer(sub("^[0-9XY]+:", "", top_snp))]
poly[, top_chr := as.integer(sub(":.*", "", top_snp))]
poly_susie <- poly[method == "susie" & !is.na(top_pos) & !is.na(top_chr)]
cat("\nPolyFun susie-converged rows:", nrow(poly_susie),
    "across", uniqueN(poly_susie$study), "studies\n")

# Try the older PolyFun dir too (May 1-3 chr files) by collapsing all per-study
poly_dir <- file.path(BASE, "results/susie_coloc_polyfun")
poly_old_list <- list()
for (sd in list.dirs(poly_dir, recursive = FALSE)) {
  files <- list.files(sd, pattern = "csv$", full.names = TRUE)
  if (!length(files)) next
  d <- rbindlist(lapply(files, fread), fill = TRUE,
                 use.names = TRUE)[, .(gwas_name, chr, top_snp, method)]
  poly_old_list[[basename(sd)]] <- d
}
poly_old <- rbindlist(poly_old_list, fill = TRUE)
setnames(poly_old, "gwas_name", "study")
poly_old[, top_pos := as.integer(sub("^[0-9XY]+:", "", top_snp))]
poly_old[, top_chr := as.integer(sub(":.*", "", top_snp))]
poly_old_susie <- poly_old[method == "susie" & !is.na(top_pos) & !is.na(top_chr)]
cat("PolyFun (newer per-study dirs) susie-converged rows:", nrow(poly_old_susie),
    "across", uniqueN(poly_old_susie$study), "studies\n")

# Use the larger / more complete set
poly_use <- if (nrow(poly_old_susie) > nrow(poly_susie)) poly_old_susie else poly_susie

# ---- For each susie row, attach the closest lead SNP within 0.5 Mb ----
attach_lead <- function(susie_dt, leads_dt, win = WIN) {
  out <- list()
  for (s in unique(susie_dt$study)) {
    a <- susie_dt[study == s]
    b <- leads_dt[study == s]
    if (!nrow(a) || !nrow(b)) next
    # For each row in a, scan b for matching chr and within window
    a_chr <- a$top_chr; a_pos <- a$top_pos
    b_chr <- b$chr;     b_pos <- b$pos; b_id <- b$lead_id
    matched <- character(nrow(a))
    for (i in seq_len(nrow(a))) {
      hit_idx <- which(b_chr == a_chr[i] & abs(b_pos - a_pos[i]) <= win)
      if (length(hit_idx)) {
        # If multiple, pick closest
        nearest <- hit_idx[which.min(abs(b_pos[hit_idx] - a_pos[i]))]
        matched[i] <- b_id[nearest]
      } else {
        matched[i] <- NA_character_
      }
    }
    a[, lead_id := matched]
    out[[s]] <- a
  }
  rbindlist(out, fill = TRUE)
}

mapped <- attach_lead(poly_use, leads)

# ---- Count unique lead-SNP windows with >=1 converged gene per study ----
poly_loci <- mapped[!is.na(lead_id),
  .(n_loci_polyfun = uniqueN(lead_id)), by = study]
poly_loci <- merge(poly_loci, reg[, .(study, ancestry)], by = "study", all.x = TRUE)

cat("\n=== PolyFun loci re-aggregated to lead-SNP windows ===\n")
print(poly_loci)
cat("\nPer-ancestry totals:\n")
print(poly_loci[, .(n_loci = sum(n_loci_polyfun, na.rm = TRUE),
                     n_studies = uniqueN(study)),
                 by = ancestry])

fwrite(poly_loci, "scripts/figures/sketches_fig3_intro/data/polyfun_loci_per_study.csv")

# ---- Also report unmapped (genes whose top_snp didn't fall in any lead window) ----
unmapped <- mapped[is.na(lead_id), .N, by = study]
cat("\nUnmapped (top_snp outside all lead windows) by study:\n")
print(unmapped)

cat("\nDone.\n")
