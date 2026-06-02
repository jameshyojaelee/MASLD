#!/usr/bin/env Rscript
# Count converged finemapped loci per study, per ancestry, per LD panel.
# Sources:
#   UKBB sghatan (Apr 11)  : combined_finemapping.csv          (default)
#   1KG (Apr 27)           : combined_finemapping_1kg.csv
#   TOP-LD (Apr 27)        : combined_finemapping_topld.csv
#   PolyFun (Apr 28-May 3) : susie_coloc_polyfun/{study}/      (per-chr files;
#                                no aggregate combined file)
suppressPackageStartupMessages({ library(data.table) })

BASE <- "GWAS/finemapping/results"
reg  <- fread("GWAS/finemapping/config/gwas_registry.tsv",
              select = c("study_name","ancestry"))
setnames(reg, "study_name", "study")

count_combined <- function(path, label) {
  if (!file.exists(path)) return(NULL)
  dt <- fread(path, select = c("study","locus","susie_converged"))
  out <- dt[susie_converged == TRUE, .(n_loci = uniqueN(locus)), by = study]
  out[, panel := label]; out
}

# UKBB sghatan, 1KG, TOP-LD
cu <- count_combined(file.path(BASE, "combined_finemapping.csv"),
                     "UKBB sghatan (Apr 11)")
ck <- count_combined(file.path(BASE, "combined_finemapping_1kg.csv"),
                     "1KG (Apr 27)")
ct <- count_combined(file.path(BASE, "combined_finemapping_topld.csv"),
                     "TOP-LD (Apr 27)")

# PolyFun — derive from per-study chr files in the most-recent dir
poly_root <- file.path(BASE, "susie_coloc_polyfun")
study_dirs <- list.dirs(poly_root, recursive = FALSE)
study_dirs <- study_dirs[basename(study_dirs) %in% reg$study]
poly_list <- list()
for (sd in study_dirs) {
  files <- list.files(sd, pattern = "csv$", full.names = TRUE)
  if (!length(files)) next
  d <- rbindlist(lapply(files, fread), fill = TRUE)
  if (!"method" %in% names(d)) next
  # PolyFun output is gene-level; "loci" approximated by unique top_snp where
  # method=="susie" (means SuSiE converged for this gene's locus window).
  # Multiple genes in the same locus would share a top_snp — uniqueN is correct.
  poly_list[[basename(sd)]] <- data.table(
    n_loci = d[method == "susie", uniqueN(top_snp)])
}
cp <- rbindlist(poly_list, idcol = "study")
cp[, panel := "PolyFun (Apr 28 / May 3)"]

all <- rbindlist(list(cu, ck, ct, cp), use.names = TRUE)
all <- merge(all, reg, by = "study", all.x = TRUE)

# Per-study wide
wide <- dcast(all, study + ancestry ~ panel, value.var = "n_loci")
fwrite(wide, "scripts/figures/sketches_fig3_intro/data/loci_per_study_per_panel.csv")
cat("=== Per-study converged loci across LD panels ===\n")
print(wide)

# Per-ancestry rollup
agg <- all[, .(n_loci = sum(n_loci, na.rm = TRUE),
               n_studies = uniqueN(study[!is.na(n_loci) & n_loci > 0])),
           by = .(panel, ancestry)]
wide2 <- dcast(agg, ancestry ~ panel, value.var = "n_loci")
wide2_studies <- dcast(agg, ancestry ~ panel, value.var = "n_studies")
fwrite(wide2, "scripts/figures/sketches_fig3_intro/data/loci_per_ancestry_per_panel.csv")
cat("\n=== Loci by ancestry × LD panel ===\n")
print(wide2)
cat("\n=== Studies covered by ancestry × LD panel ===\n")
print(wide2_studies)

cat("\nFiles written to scripts/figures/sketches_fig3_intro/data/\n")
