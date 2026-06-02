#!/usr/bin/env Rscript
# Team A4 — External Sex-Stratified eQTL Integrator
# -------------------------------------------------------------------------
# Harmonize Oliva 2020 (Science) GTEx v8 Liver sex-biased cis-eQTL data and
# cross-tabulate against:
#   - canonical SuSiE-COLOC hits (gene_level_coloc.csv, polyfun, 1kg, topld)
#   - sex_causal_scores.csv (sex_class table from RNA-seq/results/stratified_causal/)
#
# Inputs:
#   data/gtex_sb_eqtl/GTEx_Analysis_v8_sbeQTLs/GTEx_Analysis_v8_sbeQTLs.txt
#   data/gencode_v49_gene_metadata.tsv.gz
#   GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv      (sghatan v1)
#   GWAS/finemapping/results/susie_coloc_polyfun/gene_level_coloc_polyfun.csv  (production)
#   RNA-seq/results/stratified_causal/sex_causal_scores.csv
#
# Outputs:
#   data/gtex_sb_eqtl/oliva2020_liver_sbeqtl_harmonized.tsv.gz
#   RNA-seq/results/stratified_causal/oliva_sbeqtl_masld_overlap.csv
#
# Schema columns (sb-eQTL file, see README in tarball):
#   1 ensembl_gene_id    2 hugo_gene_id    3 gene_type
#   4 variant_id         5 rs_id            6 Tissue
#   7 maf
#   8 pval_nominal_sb    9 slope_sb        10 slope_se_sb
#  11 numtested         12 pvals.corrected  13 qval   (Storey, gene-level FDR)
#  14 pval_nominal_f   15 slope_f          16 slope_se_f
#  17 pval_nominal_m   18 slope_m          19 slope_se_m
#  20 pval_nominal    21 slope            22 slope_se
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(R.utils)   # for fwrite gz support
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

SB_FILE <- file.path(PROJ,
  "data/gtex_sb_eqtl/GTEx_Analysis_v8_sbeQTLs/GTEx_Analysis_v8_sbeQTLs.txt")
GENE_META <- file.path(PROJ, "data/gencode_v49_gene_metadata.tsv.gz")
COLOC_CANONICAL <- file.path(PROJ,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
COLOC_POLYFUN <- file.path(PROJ,
  "GWAS/finemapping/results/susie_coloc_polyfun/gene_level_coloc_polyfun.csv")
SEX_TABLE <- file.path(PROJ,
  "RNA-seq/results/stratified_causal/sex_causal_scores.csv")

OUT_HARM <- file.path(PROJ,
  "data/gtex_sb_eqtl/oliva2020_liver_sbeqtl_harmonized.tsv.gz")
OUT_OVER <- file.path(PROJ,
  "RNA-seq/results/stratified_causal/oliva_sbeqtl_masld_overlap.csv")

dir.create(dirname(OUT_HARM), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_OVER), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# 1. Load full sb-eQTL file; filter to Liver
# ------------------------------------------------------------------
cat("[1/5] Loading Oliva 2020 sb-eQTL master file ...\n")
sb <- fread(SB_FILE, sep = "\t", header = TRUE,
            colClasses = list(character = c("ensembl_gene_id", "hugo_gene_id",
                                            "gene_type", "variant_id", "rs_id",
                                            "Tissue")))
cat(sprintf("  Total sb-eQTL rows: %d across %d tissues\n",
            nrow(sb), uniqueN(sb$Tissue)))

liv <- sb[Tissue == "Liver"]
cat(sprintf("  Liver rows: %d\n", nrow(liv)))

# ------------------------------------------------------------------
# 2. Harmonize variant_id to chr:pos:REF:ALT (hg38)
#    GTEx variant_id format: chr10_69222380_T_TCA_b38  ->  chr10:69222380:T:TCA
# ------------------------------------------------------------------
cat("[2/5] Harmonizing variant_id to chr:pos:REF:ALT (hg38) ...\n")
liv[, c("chr", "pos", "ref", "alt", "build") :=
       tstrsplit(variant_id, "_", fixed = TRUE)]
liv[, variant_hg38 := paste(chr, pos, ref, alt, sep = ":")]
liv[, pos := as.integer(pos)]

# ------------------------------------------------------------------
# 3. Annotate Ensembl base ID + symbol from gencode_v49 metadata
# ------------------------------------------------------------------
cat("[3/5] Annotating Ensembl base + symbol from GENCODE v49 ...\n")
liv[, ensembl_base := sub("\\..*$", "", ensembl_gene_id)]

gm <- fread(GENE_META)
# columns: gene_id  gene_name  chromosome  gene_biotype  ensembl_base
gm_sub <- unique(gm[, .(ensembl_base, symbol_gencode = gene_name,
                        biotype_gencode = gene_biotype)])
liv <- merge(liv, gm_sub, by = "ensembl_base", all.x = TRUE, sort = FALSE)

# Prefer GENCODE symbol; fall back to GTEx Hugo symbol
liv[, symbol := ifelse(!is.na(symbol_gencode) & symbol_gencode != "",
                       symbol_gencode, hugo_gene_id)]

# Significance flags (Oliva 2020 uses Storey q for sb-eQTL discovery)
liv[, sb_sig_q25 := !is.na(qval) & qval <= 0.25]
liv[, sb_sig_q10 := !is.na(qval) & qval <= 0.10]
liv[, sb_sig_q05 := !is.na(qval) & qval <= 0.05]
# Sex bias direction: positive slope_sb => male-biased eQTL effect (sb slope is
# the male - female interaction slope per Oliva 2020 README convention; flagged
# as such because female stats are reported separately with their own slope).
liv[, sb_direction := ifelse(is.na(slope_sb), NA_character_,
                             ifelse(slope_sb > 0, "male_biased_effect",
                                                  "female_biased_effect"))]

# Reorder + select columns for harmonized output
harm <- liv[, .(
  ensembl_gene_id, ensembl_base, symbol, biotype_gencode, gene_type,
  variant_id, variant_hg38, rs_id, chr, pos, ref, alt, build, maf = maf,
  pval_sb = pval_nominal_sb, slope_sb, slope_se_sb,
  pvals_corrected = pvals.corrected, qval, numtested,
  pval_f = pval_nominal_f, slope_f, slope_se_f,
  pval_m = pval_nominal_m, slope_m, slope_se_m,
  pval_combined = pval_nominal, slope_combined = slope, slope_se_combined = slope_se,
  sb_direction, sb_sig_q25, sb_sig_q10, sb_sig_q05
)]

fwrite(harm, OUT_HARM, sep = "\t")
cat(sprintf("  Harmonized table -> %s  (rows=%d, cols=%d)\n",
            OUT_HARM, nrow(harm), ncol(harm)))
cat(sprintf("  Significant sb-eQTLs in Liver: q<=0.25 -> %d, q<=0.10 -> %d, q<=0.05 -> %d\n",
            sum(harm$sb_sig_q25, na.rm = TRUE),
            sum(harm$sb_sig_q10, na.rm = TRUE),
            sum(harm$sb_sig_q05, na.rm = TRUE)))

# ------------------------------------------------------------------
# 4. Integrate with canonical COLOC + sex_class table
# ------------------------------------------------------------------
cat("[4/5] Cross-tabulating with COLOC + sex_class ...\n")

read_coloc <- function(path, tag) {
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)
  # Standardize column to ensembl_base for join
  if ("ensembl" %in% names(dt)) dt[, ensembl_base := sub("\\..*$", "", ensembl)]
  setnames(dt, "coloc_best_pp4",       paste0("coloc_best_pp4_",       tag), skip_absent = TRUE)
  setnames(dt, "coloc_best_gwas",      paste0("coloc_best_gwas_",      tag), skip_absent = TRUE)
  setnames(dt, "coloc_best_susie_pp4", paste0("coloc_best_susie_pp4_", tag), skip_absent = TRUE)
  setnames(dt, "coloc_best_susie_gwas",paste0("coloc_best_susie_gwas_",tag), skip_absent = TRUE)
  dt[, c("ensembl_base", "gene",
         paste0("coloc_best_pp4_",        tag),
         paste0("coloc_best_gwas_",       tag),
         paste0("coloc_best_susie_pp4_",  tag),
         paste0("coloc_best_susie_gwas_", tag)), with = FALSE]
}

coloc_main <- read_coloc(COLOC_CANONICAL, "canonical")
coloc_pf   <- read_coloc(COLOC_POLYFUN,   "polyfun")

# Aggregate to one row per gene (max PP4) just in case
agg <- function(dt, tag) {
  if (is.null(dt)) return(NULL)
  dt[, .(
    gene_canonical = first(gene[!is.na(gene) & gene != ""]),
    pp4_abf   = suppressWarnings(max(get(paste0("coloc_best_pp4_",        tag)), na.rm = TRUE)),
    pp4_susie = suppressWarnings(max(get(paste0("coloc_best_susie_pp4_",  tag)), na.rm = TRUE))
  ), by = ensembl_base][, c("pp4_abf","pp4_susie") := .(
    ifelse(is.finite(pp4_abf), pp4_abf, NA_real_),
    ifelse(is.finite(pp4_susie), pp4_susie, NA_real_)
  )][]
}

cm <- agg(coloc_main, "canonical")
cp <- agg(coloc_pf,   "polyfun")
if (!is.null(cm)) setnames(cm, c("gene_canonical","pp4_abf","pp4_susie"),
                                c("gene_canonical_v1","pp4_abf_v1","pp4_susie_v1"))
if (!is.null(cp)) setnames(cp, c("gene_canonical","pp4_abf","pp4_susie"),
                                c("gene_canonical_pf","pp4_abf_pf","pp4_susie_pf"))

# sex_class table
sx <- fread(SEX_TABLE)
sx[, ensembl_base := sub("\\..*$", "", ensembl_id)]
sx_sub <- unique(sx[, .(ensembl_base, sex_class,
                        sex_logFC_F = logFC_F, sex_logFC_M = logFC_M,
                        sex_padj_F  = padj_F,  sex_padj_M  = padj_M,
                        sex_causal_score, sex_causal_signed)])

# Reduce harmonized to gene-level (best sb-eQTL per gene by qval)
gene_lvl <- harm[, .SD[which.min(qval)],
                 by = ensembl_base][
  , .(ensembl_base, ensembl_gene_id, symbol, biotype_gencode,
      sb_variant_hg38 = variant_hg38, sb_rs_id = rs_id, sb_maf = maf,
      sb_pval = pval_sb, sb_slope = slope_sb, sb_qval = qval,
      sb_pval_f = pval_f, sb_slope_f = slope_f,
      sb_pval_m = pval_m, sb_slope_m = slope_m,
      sb_direction, sb_sig_q25, sb_sig_q10, sb_sig_q05)]

overlap <- gene_lvl
if (!is.null(cm)) overlap <- merge(overlap, cm, by = "ensembl_base", all.x = TRUE, sort = FALSE)
if (!is.null(cp)) overlap <- merge(overlap, cp, by = "ensembl_base", all.x = TRUE, sort = FALSE)
overlap <- merge(overlap, sx_sub, by = "ensembl_base", all.x = TRUE, sort = FALSE)

# Derived convergence flags
overlap[, is_coloc_05_polyfun  := !is.na(pp4_susie_pf) & pp4_susie_pf > 0.5]
overlap[, is_coloc_05_canonical := !is.na(pp4_susie_v1) & pp4_susie_v1 > 0.5]
overlap[, is_sex_biased_dge := sex_class %in% c("Female_biased","Male_biased","Divergent")]
overlap[, sb_eqtl_coloc_dge_triple :=
          sb_sig_q25 & (is_coloc_05_polyfun | is_coloc_05_canonical) & is_sex_biased_dge]

setorder(overlap, sb_qval)
fwrite(overlap, OUT_OVER)
cat(sprintf("  Overlap table -> %s  (rows=%d, cols=%d)\n",
            OUT_OVER, nrow(overlap), ncol(overlap)))

# ------------------------------------------------------------------
# 5. Summary
# ------------------------------------------------------------------
cat("\n[5/5] Summary\n")
cat(sprintf("  Liver sb-eQTL genes tested (with q-value): %d\n",
            sum(!is.na(overlap$sb_qval))))
cat(sprintf("  sb-eQTL sig q<=0.25 : %d\n",  sum(overlap$sb_sig_q25, na.rm = TRUE)))
cat(sprintf("  sb-eQTL sig q<=0.10 : %d\n",  sum(overlap$sb_sig_q10, na.rm = TRUE)))
cat(sprintf("  sb-eQTL sig q<=0.05 : %d\n",  sum(overlap$sb_sig_q05, na.rm = TRUE)))
cat(sprintf("  + COLOC PP4>0.5 (polyfun)        : %d\n",
            sum(overlap$sb_sig_q25 & overlap$is_coloc_05_polyfun, na.rm = TRUE)))
cat(sprintf("  + sex_class != Concordant        : %d\n",
            sum(overlap$sb_sig_q25 & overlap$is_sex_biased_dge, na.rm = TRUE)))
cat(sprintf("  TRIPLE convergent (sb-eQTL + COLOC + sex-biased DGE): %d\n",
            sum(overlap$sb_eqtl_coloc_dge_triple, na.rm = TRUE)))

cat("\nTop sb-eQTL hits (q<=0.25) with COLOC / sex-DGE annotation:\n")
print(overlap[sb_sig_q25 == TRUE,
              .(symbol, sb_qval, sb_direction,
                pp4_susie_pf, pp4_susie_v1,
                sex_class, sex_causal_signed)])

cat("\nDone.\n")
