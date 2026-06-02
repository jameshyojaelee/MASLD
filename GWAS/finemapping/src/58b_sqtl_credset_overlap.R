#!/usr/bin/env Rscript
# ============================================================================
# 58b_sqtl_credset_overlap.R
# Tier-1 atlas layer: GTEx v8 Liver splicing-QTL credible-set overlap.
#
# Orthogonal genetic axis to the existing eQTL COLOC: does a MASLD GWAS
# credible-set variant ALSO act as a *significant liver sQTL* for a gene
# (i.e., alter its splicing)? GTEx signifpairs contains only significant
# variant-intron pairs, so any match is a significant liver sQTL hit. The
# binary credible-set overlap is the honest signal at GTEx liver N=208;
# a full coloc.abf arm (58c) is held as an optional sensitivity upgrade.
#
# Reuses the verified 55_gwas_atac_variant_overlap.R frontend:
#   combined_finemapping.csv (hg19) -> CS/PIP>0.1 filter -> dedup -> liftOver hg38.
# Match to GTEx Liver sQTL (hg38) by chr:pos + allele-set concordance.
# Gene from sQTL phenotype_id trailing Ensembl; mapped to symbol via GENCODE v49.
#
# Output: RNA-seq/results/multi_evidence/sqtl_atlas_columns.tsv (gene-level)
# Merged by 27a's post-assembly block (drop-then-merge on human_symbol).
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr)
  library(GenomicRanges); library(rtracklayer)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE, "GWAS/finemapping")
SQTL   <- file.path(BASE, "data/external/gtex_v8_liver_sqtl/GTEx_Analysis_v8_sQTL/Liver.v8.sqtl_signifpairs.txt.gz")
GENCODE <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
OUT    <- file.path(BASE, "RNA-seq/results/multi_evidence/sqtl_atlas_columns.tsv")

cat("== 58b_sqtl_credset_overlap ==\n")

# ---- 1-4. 55 frontend: load -> filter CS -> dedup -> liftOver hg19->hg38 ----
fm <- fread(file.path(FM_DIR, "results/combined_finemapping.csv"))
fm_cs <- fm %>% filter(
  (either_in_cs == TRUE) |
  (!is.na(recommended_pip) & recommended_pip > 0.1) |
  (is.na(recommended_pip) & !is.na(max_pip) & max_pip > 0.1))
variants <- as.data.table(fm_cs %>%
  group_by(chromosome, position, allele1, allele2) %>%
  summarise(max_pip = max(max_pip, na.rm = TRUE),
            max_rec_pip = if (any(!is.na(recommended_pip))) max(recommended_pip, na.rm = TRUE) else max(max_pip, na.rm = TRUE),
            n_studies = n_distinct(study),
            study_list = paste(sort(unique(study)), collapse = ";"),
            .groups = "drop"))
cat(sprintf("  credible-set variants (deduped): %d\n", nrow(variants)))

chain <- import.chain(file.path(BASE, "data/broadaway_eqtl/hg19ToHg38.over.chain"))
gr19 <- GRanges(paste0("chr", variants$chromosome),
                IRanges(variants$position, width = 1), strand = "*")
mcols(gr19) <- variants
lifted <- liftOver(gr19, chain)
gr38 <- unlist(lifted[lengths(lifted) == 1])
vh <- as.data.table(mcols(gr38))
vh[, `:=`(chr_hg38 = as.character(seqnames(gr38)), pos_hg38 = start(gr38))]
cat(sprintf("  lifted to hg38: %d\n", nrow(vh)))

# ---- 5. GTEx Liver sQTL: parse variant_id (hg38) + phenotype_id (gene) ----
sq <- fread(SQTL, select = c("variant_id", "phenotype_id", "pval_nominal"))
sq[, c("chr_hg38", "pos_hg38", "ref", "alt") :=
     tstrsplit(variant_id, "_", keep = 1:4)]
sq[, pos_hg38 := as.integer(pos_hg38)]
sq[, ensembl_base := sub("\\..*$", "", sub("^.*:", "", phenotype_id))]
cat(sprintf("  GTEx liver sQTL sig pairs: %d (introns in %d genes)\n",
            nrow(sq), uniqueN(sq$ensembl_base)))

# ---- 6. Overlap by chr:pos + allele-set concordance ----
mrg <- merge(vh[, .(chr_hg38, pos_hg38, allele1, allele2, max_pip, max_rec_pip,
                    n_studies, study_list)],
             sq[, .(chr_hg38, pos_hg38, ref, alt, ensembl_base, pval_nominal)],
             by = c("chr_hg38", "pos_hg38"), allow.cartesian = TRUE)
allele_ok <- function(a1, a2, r, a)
  mapply(function(a1, a2, r, a)
    setequal(toupper(c(a1, a2)), toupper(c(r, a))), a1, a2, r, a)
mrg <- mrg[allele_ok(allele1, allele2, ref, alt)]
cat(sprintf("  credset x sQTL matches (allele-concordant): %d (%d genes)\n",
            nrow(mrg), uniqueN(mrg$ensembl_base)))

# ---- 7. Aggregate per gene ----
g2s <- unique(fread(GENCODE, select = c("ensembl_base", "gene_name")))
agg <- mrg[, .(sqtl_n_introns = .N,
               sqtl_min_pval = min(pval_nominal, na.rm = TRUE),
               sqtl_max_credset_pip = max(max_rec_pip, na.rm = TRUE),
               sqtl_studies = paste(sort(unique(unlist(strsplit(study_list, ";")))), collapse = ";")),
           by = ensembl_base]
agg <- merge(agg, g2s, by = "ensembl_base", all.x = TRUE)
agg <- agg[!is.na(gene_name) & gene_name != ""]
out <- agg[, .(human_symbol = gene_name, sqtl_credset_hit = TRUE,
               sqtl_n_introns, sqtl_min_pval, sqtl_max_credset_pip, sqtl_studies)]
out <- out[order(sqtl_min_pval)][!duplicated(human_symbol)]
# v2 relevance gate: flag hits driven by a MASLD-DISEASE GWAS (NAFLD/NASH/cirrhosis/HCC/
# steatosis endpoints; Ghodsian; deCODE NAFL), vs hits supported only by liver-enzyme
# proxies (ALT/AST/GGT/PDFF). MASLD-driven hits are the defensible T1 evidence.
out[, sqtl_masld_gwas_driven :=
      grepl("NAFLD|NASH|cirrh|HCC|steato|hepatocell|ghodsian|nafl", sqtl_studies, ignore.case = TRUE)]

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(out, OUT, sep = "\t")
cat(sprintf("WROTE %s : %d genes with a liver-sQTL credible-set hit\n", OUT, nrow(out)))
print(head(out[order(sqtl_min_pval)], 10))
