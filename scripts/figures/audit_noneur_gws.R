#!/usr/bin/env Rscript
# Audit: of the non-EUR-unique colocalizing genes (Fig2E), how many reach genome-wide
# significance (p < 5e-8) at their colocalizing lead variant in a non-European GWAS?
# Matches coloc top_snp (hg19 chr:pos) -> the SAME reformatted-hg19 sumstats the
# finemapping/COLOC pipeline ran on -> p_value (col 7).
#
# COLOC SET (2026-07-06): SuSiE-PRIMARY. Colocalization is defined as PP.H4.susie > 0.5
# ONLY (the 473-gene SuSiE-COLOC set), NOT the former SuSiE-OR-ABF union. This matches the
# gated Fig2E multi-signal COLOC panel (fig2_ancestry_unique_coloc_gated.py); the non-EUR-unique
# gene set is therefore SMALLER than the retired union audit.
#
# ANCESTRY (2026-07-05, MVP-fix): registry-driven via gwas_ancestry() from
# load_figure_data.R, covering the 50-GWAS portfolio (23 legacy + 27 MVP strata) incl.
# AMR. RETIRES (a) the hardcoded heuristic grepl(BBJ->EAS / AFR->AFR / CSA->SAS /
# else->EUR), which had NO AMR bin and misrouted every MVP_* stratum into EUR, and
# (b) the old lookup that read only GWAS/MR_Data/{BBJ,PanUKBB}/*_harmonised_hg38.tsv.gz
# (structurally blind to MVP p-values). Sumstats paths now come straight from
# gwas_registry.tsv `sumstats_path` for all 50 strata.
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))  # gwas_ancestry() + load_gwas_registry()
reg <- load_gwas_registry()
# MAIN (Tier-1/2, liver-specific) restriction (2026-07-06): audit only the placement=="main"
# strata so the non-EUR-unique gene set matches the main-scoped Fig2E ancestry panel.
MAIN_STUDIES <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[
  placement == "main", study_name]
COLOC_INPUT <- Sys.getenv(
  "FIG2_COLOC_INPUT",
  file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
co <- fread(COLOC_INPUT)
co <- co[gwas_name %in% MAIN_STUDIES]                            # keep MAIN (Tier-1/2) strata only
co[, ancestry := as.character(gwas_ancestry(gwas_name))]         # registry-driven EUR/AFR/AMR/EAS/SAS
co[, pp4 := PP.H4.susie]                                         # SuSiE-PRIMARY (2026-07-06): SuSiE only, no ABF union
sig <- co[is.finite(pp4) & pp4 > 0.5]

# non-EUR-unique genes + their specificity bin
g <- sig[, .(has_eur = any(ancestry=="EUR"),
             ancs = paste(sort(unique(ancestry[ancestry!="EUR"])), collapse="+"),
             n_non = uniqueN(ancestry[ancestry!="EUR"])), by = gene]
nu <- g[has_eur == FALSE]
nu[, bin := fcase(n_non>=2, "multiple non-EUR", ancs=="EAS","EAS only",
                  ancs=="AFR","AFR only", ancs=="AMR","AMR only",
                  ancs=="SAS","SAS only", default="other")]
cat("non-EUR-unique genes:", nrow(nu), "| by bin:\n"); print(nu[, .N, by=bin][order(-N)])

# coloc lead variants for these genes in their non-EUR GWAS
cc <- sig[gene %in% nu$gene & ancestry!="EUR", .(gene, gwas_name, ancestry, top_snp, pp4,
                                                 method = fifelse(!is.na(PP.H4.susie),"susie","abf"))]
cc[, c("tchr","tpos") := tstrsplit(top_snp, ":", keep = 1:2)]
cc[, key := paste0(tchr, ":", tpos)]

# look up GWAS p_value per non-EUR GWAS in its reformatted-hg19 sumstats (registry path).
# Uniform schema (all 50 strata): chromosome pos allele1 allele2 beta se pval
#   -> match key = $1":"$2 (chr:pos hg19), read pval = $7. Plain .tsv (no gz).
fpath <- function(gw) file.path(BASE, "GWAS/finemapping",
                                reg$sumstats_path[match(gw, reg$study_name)])
res <- rbindlist(lapply(unique(cc$gwas_name), function(gw) {
  fp <- fpath(gw)
  if (is.na(fp) || !file.exists(fp)) { message("MISSING sumstats for ", gw, " -> ", fp); return(NULL) }
  keys <- unique(cc[gwas_name==gw]$key); wf <- tempfile(); writeLines(keys, wf)
  cmd <- sprintf("awk -F'\\t' 'NR==FNR{w[$1]=1;next} FNR>1 && (($1\":\"$2) in w){print $1\":\"$2\"\\t\"$7}' %s %s",
                 shQuote(wf), shQuote(fp))
  r <- tryCatch(fread(cmd=cmd, header=FALSE, col.names=c("key","pval")), error=function(e) data.table(key=character(),pval=numeric()))
  if (nrow(r)) r[, gwas_name := gw]; r
}), fill = TRUE)
res <- res[, .(pval = min(pval, na.rm=TRUE)), by=.(gwas_name, key)]   # dedup
cc <- merge(cc, res, by=c("gwas_name","key"), all.x=TRUE)

# per gene: best (min) GWAS p across its non-EUR coloc leads, + max pp4
pg <- cc[, .(min_p = suppressWarnings(min(pval, na.rm=TRUE)),
             max_pp4 = max(pp4),
             best_gwas = gwas_name[which.min(pval)],
             best_snp  = top_snp[which.min(pval)],
             any_susie = any(method=="susie")), by=gene]
pg[is.infinite(min_p), min_p := NA_real_]
pg <- merge(pg, nu[, .(gene, bin)], by="gene")
pg[, tier := fcase(is.na(min_p), "no_p_found",
                   min_p < 5e-8, "GWS (<5e-8)",
                   min_p < 1e-6, "suggestive (5e-8..1e-6)",
                   default = "weak (>=1e-6)")]

cat("\n================ AUDIT RESULT ================\n")
cat("Total non-EUR-unique genes audited:", nrow(pg), "\n\n")
cat("By significance tier of the colocalizing lead variant (min p across its non-EUR GWAS):\n")
print(pg[, .N, by=tier][order(-N)])
cat("\nGenome-wide-significant (p<5e-8) by bin:\n")
print(dcast(pg, bin ~ tier, value.var="gene", fun.aggregate=length))
cat("\n--- the GWS-passing genes (real non-EUR discoveries) ---\n")
print(pg[tier=="GWS (<5e-8)"][order(min_p), .(gene, bin, min_p, max_pp4=round(max_pp4,2), best_gwas, any_susie)])
cat(sprintf("\n--- all %d non-EUR-unique, sorted by p (top 25 shown) ---\n", nrow(pg)))
print(head(pg[order(min_p), .(gene, bin, min_p, max_pp4=round(max_pp4,2), best_gwas, tier)], 25))

# GWS / suggestive counts by ancestry bin (for the report + Fig2E caption)
cat("\nGWS (<5e-8) count by bin:\n")
print(pg[tier=="GWS (<5e-8)", .N, by=bin][order(-N)])
cat("\nsuggestive (5e-8..1e-6) count by bin:\n")
print(pg[tier=="suggestive (5e-8..1e-6)", .N, by=bin][order(-N)])

# Canonical output consumed by scripts/figures/fig2_ancestry_unique_coloc_gated.py
outdir <- Sys.getenv(
  "FIG2_VARIANT_CLASS_DIR",
  file.path(BASE, "RNA-seq/results/coloc_variant_classes"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, "noneur_gws_audit.csv")
fwrite(pg[order(min_p)], outfile)
cat("\nwrote:", outfile, "|", nrow(pg), "non-EUR-unique genes\n")
