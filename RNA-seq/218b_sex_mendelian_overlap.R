#!/usr/bin/env Rscript
# 218b_sex_mendelian_overlap.R — A12 — Sex-biased Mendelian liver disease overlap
#
# Curated sex-biased Mendelian liver disease genes (literature) crossed with:
#   - sex_class (from sex_deg_classification.csv)
#   - COLOC (from gene_level_coloc.csv)
#   - Okur 2026 monogenic panel (data/external/okur2026_monogenic/)
#
# Outputs:
#   RNA-seq/results/stratified_causal/sex_mendelian_overlap.csv

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n", sep = " ")
log_msg("218b sex × Mendelian × COLOC starting")

# ─── Curated sex-biased Mendelian liver disease genes ─────────────────────────
# Sources: OMIM, GeneReviews, primary literature. Effect direction = which sex
# has more severe disease penetrance / earlier onset.
sex_biased_mendelian <- data.table(
  gene = c(
    # Wilson disease — F worse: more severe, earlier hepatic presentation
    "ATP7B",
    # Hereditary hemochromatosis — M worse: M penetrance 28%, F 1% (menstrual
    # iron loss + estrogen hepcidin effect)
    "HFE", "HJV", "HAMP", "TFR2", "SLC40A1",
    # Alpha-1 antitrypsin deficiency — context-dependent; PiZZ males more cirrhosis
    "SERPINA1",
    # Glycogen storage disease type Ia — F has more hepatic adenomas
    "G6PC1",
    # Crigler-Najjar / Gilbert (UGT1A1) — M higher Gilbert prevalence
    "UGT1A1",
    # Citrin deficiency (NICCD/CTLN2) — M:F ~ 2:1 for adult-onset CTLN2
    "SLC25A13",
    # Wilson-like / X-linked MEDNIK (AP1S1) — X-linked
    "AP1S1",
    # ABCB4 (PFIC3) — F higher with intrahepatic cholestasis of pregnancy (ICP)
    "ABCB4",
    # ABCB11 (PFIC2/BRIC2) — F slight excess in BRIC2 episodes
    "ABCB11",
    # ATP8B1 (PFIC1/BRIC1)
    "ATP8B1",
    # PNPLA3 — common variant, included for cross-ref to MASLD genetics; M effect larger
    "PNPLA3"
  ),
  more_severe_sex = c(
    "F",                              # ATP7B
    "M", "M", "M", "M", "M",          # iron overload genes
    "M",                              # SERPINA1
    "F",                              # G6PC1
    "M",                              # UGT1A1
    "M",                              # SLC25A13
    "X-linked",                       # AP1S1
    "F",                              # ABCB4 (ICP)
    "F",                              # ABCB11 (BRIC2)
    "F",                              # ATP8B1
    "M"                               # PNPLA3 common variant
  ),
  mendelian_disease = c(
    "Wilson disease",
    "Hereditary hemochromatosis (HFE)",
    "Hemochromatosis type 2A (juvenile)",
    "Hemochromatosis type 2B (juvenile)",
    "Hemochromatosis type 3",
    "Ferroportin disease",
    "Alpha-1 antitrypsin deficiency",
    "GSD type Ia",
    "Crigler-Najjar / Gilbert syndrome",
    "Citrin deficiency / CTLN2",
    "MEDNIK syndrome",
    "PFIC3 / ICP",
    "PFIC2 / BRIC2",
    "PFIC1 / BRIC1",
    "MASLD (common variant)"
  ),
  source_evidence = c(
    "OMIM 277900; Litwin 2012 Hepatology (F more severe hepatic, M more severe neuro)",
    "GeneReviews HFE; Allen 2008 NEJM (M penetrance 28%, F 1%)",
    "OMIM 602390",
    "OMIM 613313",
    "OMIM 604250",
    "OMIM 606069",
    "GeneReviews AAT-D; M PiZZ cirrhosis risk higher",
    "OMIM 232200; Kishnani 2014 (F adenoma excess)",
    "OMIM 218800; M Gilbert prevalence 2-5× F",
    "OMIM 605814; Saheki 2008",
    "OMIM 609313",
    "OMIM 171060; Smith 2024 ICP review",
    "OMIM 605479",
    "OMIM 211600",
    "Romeo 2008 Nat Genet; sex-modifier effect"
  )
)

log_msg(sprintf("Curated Mendelian panel: %d genes", nrow(sex_biased_mendelian)))

# ─── Load Okur 2026 panel ─────────────────────────────────────────────────────
log_msg("Loading Okur 2026 monogenic panel")
okur_panel <- fread(file.path(BASE, "data/external/okur2026_monogenic",
                              "okur2026_panel_genes_partial.tsv"))
log_msg(sprintf("Okur panel: %d genes", nrow(okur_panel)))

# Cross-tab: which of our curated set is in Okur panel?
sex_biased_mendelian[, in_okur_panel := gene %in% okur_panel$gene]

# ─── Load sex_deg + COLOC ──────────────────────────────────────────────────────
log_msg("Loading sex-DEG classification (v3 mashr Bayesian preferred; v2 fallback)")
sex_v3_path <- file.path(BASE,
   "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/sex_deg_classification_v3.csv")
sex_v2_path <- file.path(BASE,
   "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
sex_deg_file <- if (file.exists(sex_v3_path)) sex_v3_path else sex_v2_path
log_msg(sprintf("  source: %s", basename(dirname(sex_deg_file))))
sex_deg <- fread(sex_deg_file)
# v3 column compatibility shim:
#  - v3 names interaction logFC as `beta_interaction`
#  - v3 lacks frequentist sig_M/sig_F flags; substitute Bayesian lfsr<0.05 (equivalent decision).
if (!"interaction_logFC" %in% names(sex_deg) && "beta_interaction" %in% names(sex_deg)) {
  sex_deg[, interaction_logFC := beta_interaction]
}
if (!"sig_M" %in% names(sex_deg) && "lfsr_M" %in% names(sex_deg)) {
  sex_deg[, sig_M := lfsr_M < 0.05]
}
if (!"sig_F" %in% names(sex_deg) && "lfsr_F" %in% names(sex_deg)) {
  sex_deg[, sig_F := lfsr_F < 0.05]
}
sex_deg[, ensembl_base := sub("\\..*", "", gene)]
gencode <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
sex_deg <- merge(sex_deg, gencode[, .(ensembl_base, gene_name)],
                 by = "ensembl_base", all.x = TRUE)
sex_deg[, symbol := gene_name]

log_msg("Loading COLOC")
coloc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc[, pp4_max := pmax(coloc_best_pp4, coloc_best_susie_pp4, na.rm = TRUE)]
coloc[is.infinite(pp4_max), pp4_max := NA_real_]

# ─── Merge ─────────────────────────────────────────────────────────────────────
ann <- merge(sex_biased_mendelian, sex_deg[, .(symbol, sex_class,
                                                logFC_M, logFC_F,
                                                interaction_logFC,
                                                interaction_padj,
                                                sig_M, sig_F)],
             by.x = "gene", by.y = "symbol", all.x = TRUE)
ann <- merge(ann, coloc[, .(symbol = gene, pp4_max, coloc_best_pp4,
                            coloc_best_susie_pp4, coloc_best_gwas)],
             by.x = "gene", by.y = "symbol", all.x = TRUE)

# Direction concordance: does sex_class match expected penetrance bias?
ann[, expected_rna_class := dplyr::case_when(
  more_severe_sex == "F" ~ "Female_biased",
  more_severe_sex == "M" ~ "Male_biased",
  TRUE                   ~ NA_character_
)]
ann[, direction_concordant := !is.na(sex_class) &
                              !is.na(expected_rna_class) &
                              sex_class == expected_rna_class]

ann[, has_coloc := !is.na(pp4_max) & pp4_max >= 0.5]
ann[, sex_dimorphic_rna := !is.na(sex_class) &
                            sex_class %in% c("Female_biased", "Male_biased",
                                             "Divergent")]

# ─── Fisher test: are curated Mendelian sex-biased genes enriched for sex_class? ─
# 2×2: rows = (Mendelian sex-biased curated set vs background sex_deg)
#      cols = (sex_dimorphic vs concordant)
bg <- sex_deg[!is.na(symbol)]
mend_set <- unique(sex_biased_mendelian$gene)

a <- sum(bg$symbol %in% mend_set &
         bg$sex_class %in% c("Female_biased","Male_biased","Divergent"), na.rm = TRUE)
b <- sum(bg$symbol %in% mend_set &
         !(bg$sex_class %in% c("Female_biased","Male_biased","Divergent")), na.rm = TRUE)
c_ <- sum(!(bg$symbol %in% mend_set) &
         bg$sex_class %in% c("Female_biased","Male_biased","Divergent"), na.rm = TRUE)
d <- sum(!(bg$symbol %in% mend_set) &
         !(bg$sex_class %in% c("Female_biased","Male_biased","Divergent")), na.rm = TRUE)
ft_dimorphic <- fisher.test(matrix(c(a,b,c_,d), nrow = 2))
log_msg(sprintf("Mendelian vs background — sex_dimorphic enrichment: OR=%.2f, p=%.3g (a=%d, b=%d, c=%d, d=%d)",
                ft_dimorphic$estimate, ft_dimorphic$p.value, a, b, c_, d))

# Direction-specific Fisher (Mendelian × expected direction)
fisher_results <- list()
for (xs in c("F", "M")) {
  exp_class <- if (xs == "F") "Female_biased" else "Male_biased"
  set <- sex_biased_mendelian[more_severe_sex == xs, gene]
  if (length(set) == 0) next
  a2 <- sum(bg$symbol %in% set & bg$sex_class == exp_class, na.rm = TRUE)
  b2 <- sum(bg$symbol %in% set & bg$sex_class != exp_class, na.rm = TRUE)
  c2 <- sum(!(bg$symbol %in% set) & bg$sex_class == exp_class, na.rm = TRUE)
  d2 <- sum(!(bg$symbol %in% set) & bg$sex_class != exp_class, na.rm = TRUE)
  ft <- fisher.test(matrix(c(a2,b2,c2,d2), nrow=2))
  fisher_results[[xs]] <- data.table(
    panel = sprintf("Mendelian_%s_severe", xs),
    expected_class = exp_class,
    n_panel = length(set),
    n_concordant = a2,
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    ci_low = ft$conf.int[1],
    ci_high = ft$conf.int[2]
  )
}
fisher_dt <- rbindlist(fisher_results)
fisher_dt[, padj := p.adjust(p_value, "BH")]
log_msg("Fisher tests by direction:"); print(fisher_dt)

# ─── COLOC × Mendelian fisher ─────────────────────────────────────────────────-
a3 <- sum(bg$symbol %in% mend_set &
          bg$symbol %in% coloc[pp4_max >= 0.5, gene], na.rm = TRUE)
b3 <- sum(bg$symbol %in% mend_set &
          !(bg$symbol %in% coloc[pp4_max >= 0.5, gene]), na.rm = TRUE)
c3 <- sum(!(bg$symbol %in% mend_set) &
          bg$symbol %in% coloc[pp4_max >= 0.5, gene], na.rm = TRUE)
d3 <- sum(!(bg$symbol %in% mend_set) &
          !(bg$symbol %in% coloc[pp4_max >= 0.5, gene]), na.rm = TRUE)
ft_coloc <- fisher.test(matrix(c(a3,b3,c3,d3), nrow=2))
log_msg(sprintf("Mendelian × COLOC enrichment: OR=%.2f, p=%.3g (a=%d b=%d c=%d d=%d)",
                ft_coloc$estimate, ft_coloc$p.value, a3, b3, c3, d3))

# Append meta-level fisher rows
fisher_dt <- rbind(fisher_dt,
  data.table(panel = "Mendelian_dimorphic_overall",
             expected_class = "any_dimorphic",
             n_panel = length(mend_set),
             n_concordant = a,
             odds_ratio = unname(ft_dimorphic$estimate),
             p_value = ft_dimorphic$p.value,
             ci_low = ft_dimorphic$conf.int[1],
             ci_high = ft_dimorphic$conf.int[2],
             padj = NA_real_),
  data.table(panel = "Mendelian_coloc_overall",
             expected_class = "coloc_pp4_0.5",
             n_panel = length(mend_set),
             n_concordant = a3,
             odds_ratio = unname(ft_coloc$estimate),
             p_value = ft_coloc$p.value,
             ci_low = ft_coloc$conf.int[1],
             ci_high = ft_coloc$conf.int[2],
             padj = NA_real_)
)

# ─── Write outputs ────────────────────────────────────────────────────────────-
out_path <- file.path(OUT_DIR, "sex_mendelian_overlap.csv")
fwrite(ann, out_path)
log_msg(sprintf("Wrote %s", out_path))

fisher_path <- file.path(OUT_DIR, "sex_mendelian_fisher.csv")
fwrite(fisher_dt, fisher_path)
log_msg(sprintf("Wrote %s", fisher_path))

log_msg("218b complete")
