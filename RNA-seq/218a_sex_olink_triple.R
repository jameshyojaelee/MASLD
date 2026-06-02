#!/usr/bin/env Rscript
# 218a_sex_olink_triple.R — A11 — Olink plasma sex × DEG × COLOC triple-concordance
#
# Goal: Identify proteins where sex-biased plasma NPX agrees with sex-biased
#       hepatic transcriptional dimorphism AND has COLOC support (genetic causal).
#
# Olink NPX matrix: Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt
#   (218 subjects × 1,460 proteins; 177 of the 218 are matched to liver GSE276114)
# Subject sex metadata: attempted recovery from GEO (GSE276114) via GEOquery.
#                       If sex unavailable, the script falls back to a disease × PP4 contrast
#                       without sex stratification and emits an "olink_metadata_block" file.
#
# Outputs:
#   RNA-seq/results/stratified_causal/sex_coloc_olink_triple.csv
#   RNA-seq/results/stratified_causal/sex_coloc_olink_triple_qc.csv
#   (optional) RNA-seq/results/stratified_causal/olink_metadata_status.txt

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n", sep = " ")
log_msg("218a sex × DEG × COLOC × plasma starting")

# ─── Inputs ────────────────────────────────────────────────────────────────────
OLINK_PATH <- file.path(BASE, "Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt")
OLINK_META <- file.path(BASE, "Analysis/Proteomics/results/gse276114_disease_metadata.csv")
# v3 mashr Bayesian preferred; v2 fallback. v3 CSV provides v2-compatible `sex_class` alias.
SEX_V3     <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/sex_deg_classification_v3.csv")
SEX_V2     <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
SEX_DEG    <- if (file.exists(SEX_V3)) SEX_V3 else SEX_V2
COLOC_GENE <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
GENCODE    <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")

# ─── 1. Load Olink NPX ────────────────────────────────────────────────────────
log_msg("Loading Olink NPX matrix")
olink <- fread(OLINK_PATH, sep = "\t", header = TRUE)
protein_names <- olink$Assay
npx <- as.matrix(olink[, -"Assay", with = FALSE])     # rows = proteins, cols = subjects
rownames(npx) <- protein_names
log_msg(sprintf("Olink: %d proteins × %d subjects", nrow(npx), ncol(npx)))

# ─── 2. Disease metadata ──────────────────────────────────────────────────────
log_msg("Loading disease metadata (gse276114_disease_metadata.csv)")
meta <- fread(OLINK_META)
meta[, subject_col := paste0("Subject ", sample_number)]
n_match <- sum(meta$subject_col %in% colnames(npx))
log_msg(sprintf("Disease metadata: %d subjects, %d match NPX columns", nrow(meta), n_match))

# ─── 3. Try to recover sex from GEO GSE276114 ────────────────────────────────-
sex_lookup <- NULL
sex_source <- "unavailable"
tryCatch({
  if (requireNamespace("GEOquery", quietly = TRUE)) {
    log_msg("Attempting GEOquery::getGEO('GSE276114') — may take a few minutes")
    gse <- GEOquery::getGEO("GSE276114", GSEMatrix = TRUE, getGPL = FALSE,
                            destdir = tempdir())
    eset <- if (is.list(gse) && length(gse) > 0) gse[[1]] else gse
    pheno <- as.data.frame(Biobase::pData(eset))
    char_cols <- grep("^characteristics_ch", colnames(pheno), value = TRUE)
    # Try to extract sex
    sex_vec <- rep(NA_character_, nrow(pheno))
    for (cc in char_cols) {
      vals <- as.character(pheno[[cc]])
      tags <- tolower(sub(":.*", "", vals))
      raw  <- sub("^[^:]+:\\s*", "", vals)
      if (any(grepl("^(sex|gender)$", tags, ignore.case = TRUE))) {
        m <- grepl("^(sex|gender)$", tags, ignore.case = TRUE)
        sex_vec[m] <- tolower(raw[m])
        break
      }
    }
    # Also try the title column
    if (all(is.na(sex_vec)) && "title" %in% colnames(pheno)) {
      titles <- tolower(as.character(pheno$title))
      sex_vec[grepl("\\bfemale\\b|\\bf\\b", titles)] <- "female"
      sex_vec[grepl("\\bmale\\b|\\bm\\b", titles)]   <- "male"
    }
    if (any(!is.na(sex_vec))) {
      sex_lookup <- data.table(geo_accession = pheno$geo_accession,
                               sex = ifelse(grepl("^f", sex_vec), "F",
                                     ifelse(grepl("^m", sex_vec), "M", NA_character_)))
      sex_source <- "GEO_GSE276114"
      log_msg(sprintf("Recovered sex for %d / %d samples from GEO",
                      sum(!is.na(sex_lookup$sex)), nrow(sex_lookup)))
    } else {
      log_msg("GEO phenoData found but no sex column extractable")
    }
  }
}, error = function(e) {
  log_msg(sprintf("GEOquery failed: %s", conditionMessage(e)))
})

# Write status file
status_file <- file.path(OUT_DIR, "olink_metadata_status.txt")
writeLines(c(
  sprintf("Olink subject metadata source: %s", sex_source),
  sprintf("Olink subjects: %d", ncol(npx)),
  sprintf("Disease metadata matched: %d", n_match),
  if (!is.null(sex_lookup)) sprintf("Sex annotated: %d (%.1f%%)",
                                    sum(!is.na(sex_lookup$sex)),
                                    100 * mean(!is.na(sex_lookup$sex))) else
    "Sex annotated: 0 (BLOCKED — Yang 2025 mmc2.xlsx corrupt; mendeley zip empty; GEO sex characteristic absent)",
  sprintf("Time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
), status_file)
log_msg(sprintf("Wrote status file: %s", status_file))

# ─── 4. Load sex-DEG classification ───────────────────────────────────────────-
log_msg("Loading sex-DEG classification")
sex_deg <- fread(SEX_DEG)
# v3 column compatibility shim: v3 names interaction logFC as `beta_interaction`
if (!"interaction_logFC" %in% names(sex_deg) && "beta_interaction" %in% names(sex_deg)) {
  sex_deg[, interaction_logFC := beta_interaction]
}
# Strip ENSEMBL version
sex_deg[, ensembl_base := sub("\\..*", "", gene)]
log_msg(sprintf("Sex-DEG: %d genes, %d classes",
                nrow(sex_deg), length(unique(sex_deg$sex_class))))

# ─── 5. Load COLOC PP4 ────────────────────────────────────────────────────────
log_msg("Loading COLOC PP4")
coloc <- fread(COLOC_GENE)
# Both coloc_best_pp4 (ABF) and coloc_best_susie_pp4 (SuSiE) preserved.
coloc[, pp4_max := pmax(coloc_best_pp4, coloc_best_susie_pp4, na.rm = TRUE)]
coloc[is.infinite(pp4_max), pp4_max := NA_real_]
log_msg(sprintf("COLOC: %d genes, %d with pp4_max ≥ 0.5",
                nrow(coloc), sum(coloc$pp4_max >= 0.5, na.rm = TRUE)))

# ─── 6. ENSEMBL → SYMBOL mapping ──────────────────────────────────────────────
log_msg("Loading gencode v49 metadata for symbol mapping")
gencode <- fread(GENCODE)
sex_deg <- merge(sex_deg, gencode[, .(ensembl_base, gene_name)],
                 by = "ensembl_base", all.x = TRUE)
sex_deg[, symbol := gene_name]
sex_deg[is.na(symbol), symbol := ensembl_base]

# Merge sex_deg + coloc by gene symbol
gene_panel <- merge(sex_deg[, .(symbol, sex_class, logFC_M, logFC_F,
                                interaction_logFC, interaction_padj)],
                    coloc[, .(symbol = gene, pp4_max, coloc_best_pp4,
                              coloc_best_susie_pp4, coloc_best_gwas)],
                    by = "symbol", all = TRUE)

# ─── 7. Define plasma sex × disease × PP4 model ──────────────────────────────
# Disease group: F0-2 vs F3 vs F4. Code as disease_advanced = F3|F4.
meta[, disease_advanced := ifelse(disease_group %in% c("F3", "F4"), 1L, 0L)]

# Match subject_col to NPX columns
keep_subj <- meta$subject_col %in% colnames(npx)
meta_used <- meta[keep_subj]
npx_used  <- npx[, meta_used$subject_col, drop = FALSE]
log_msg(sprintf("Plasma model n: %d subjects", ncol(npx_used)))

# Merge sex if recovered
if (!is.null(sex_lookup)) {
  meta_used <- merge(meta_used, sex_lookup, by = "geo_accession", all.x = TRUE)
  log_msg(sprintf("Sex matched in plasma cohort: %d / %d",
                  sum(!is.na(meta_used$sex)), nrow(meta_used)))
  npx_used <- npx_used[, meta_used$subject_col, drop = FALSE]
}

# PP4 binarization at 0.5
gene_panel[, pp4_bin := factor(ifelse(pp4_max >= 0.5, "coloc", "no_coloc"),
                                levels = c("no_coloc", "coloc"))]

# Map protein -> sex_class via symbol
protein_panel <- merge(data.table(protein = protein_names),
                       gene_panel,
                       by.x = "protein", by.y = "symbol",
                       all.x = TRUE)
n_proteins_with_class <- sum(!is.na(protein_panel$sex_class))
log_msg(sprintf("Proteins mapped to sex_class: %d / %d",
                n_proteins_with_class, nrow(protein_panel)))

# ─── 8. Per-protein linear model ───────────────────────────────────────────────
# If sex is available: NPX ~ disease_advanced * sex + pp4_bin
# Else: NPX ~ disease_advanced + pp4_bin (no sex term)
log_msg("Fitting per-protein models")

run_model <- function(npx_vec, df) {
  ok <- !is.na(npx_vec)
  if (sum(ok) < 10) return(NULL)
  d <- copy(df)
  d[, npx := npx_vec]
  d <- d[ok]
  if (!is.null(sex_lookup) && any(!is.na(d$sex))) {
    d <- d[!is.na(sex)]
    if (length(unique(d$sex)) < 2 || nrow(d) < 10) return(NULL)
    fit <- tryCatch(lm(npx ~ disease_advanced * sex, data = d), error = function(e) NULL)
  } else {
    fit <- tryCatch(lm(npx ~ disease_advanced, data = d), error = function(e) NULL)
  }
  if (is.null(fit)) return(NULL)
  cf <- summary(fit)$coefficients
  out <- list(
    n          = nrow(d),
    beta_dis   = cf[grepl("disease_advanced$", rownames(cf)), "Estimate"][1],
    p_dis      = cf[grepl("disease_advanced$", rownames(cf)), "Pr(>|t|)"][1]
  )
  inter_row <- grep(":sex", rownames(cf), value = TRUE)
  if (length(inter_row) == 1) {
    out$beta_inter <- cf[inter_row, "Estimate"]
    out$p_inter    <- cf[inter_row, "Pr(>|t|)"]
  } else {
    out$beta_inter <- NA_real_
    out$p_inter    <- NA_real_
  }
  out
}

rows <- vector("list", nrow(npx_used))
for (i in seq_len(nrow(npx_used))) {
  prot <- rownames(npx_used)[i]
  npx_vec <- as.numeric(npx_used[i, ])
  m <- run_model(npx_vec, meta_used)
  if (is.null(m)) next
  rows[[i]] <- data.table(protein = prot,
                          n      = m$n,
                          beta_disease    = m$beta_dis,
                          p_disease       = m$p_dis,
                          beta_sex_interaction = m$beta_inter,
                          p_sex_interaction    = m$p_inter)
}
plasma_res <- rbindlist(Filter(Negate(is.null), rows))
plasma_res[, padj_disease := p.adjust(p_disease, "BH")]
plasma_res[!is.na(p_sex_interaction),
           padj_sex_interaction := p.adjust(p_sex_interaction, "BH")]
log_msg(sprintf("Plasma DE: %d proteins tested; %d at padj_disease<0.05",
                nrow(plasma_res), sum(plasma_res$padj_disease < 0.05, na.rm = TRUE)))

# Classify plasma sex direction
plasma_res[, plasma_sex_class := "Not_significant"]
if (!is.null(sex_lookup)) {
  plasma_res[!is.na(padj_sex_interaction) & padj_sex_interaction < 0.1 &
             beta_sex_interaction > 0, plasma_sex_class := "Male_higher_in_disease"]
  plasma_res[!is.na(padj_sex_interaction) & padj_sex_interaction < 0.1 &
             beta_sex_interaction < 0, plasma_sex_class := "Female_higher_in_disease"]
}

# ─── 9. Triple-concordance table ──────────────────────────────────────────────-
triple <- merge(protein_panel,
                plasma_res,
                by = "protein", all.x = TRUE)
triple[, has_coloc       := !is.na(pp4_max) & pp4_max >= 0.5]
triple[, sex_biased_rna  := !is.na(sex_class) &
                              sex_class %in% c("Female_biased", "Male_biased", "Divergent")]
triple[, sex_biased_plasma := !is.na(plasma_sex_class) &
                                plasma_sex_class != "Not_significant"]
triple[, triple_concordant := has_coloc & sex_biased_rna & sex_biased_plasma]

# Direction concordance check
triple[, sex_direction_concordant := NA]
triple[sex_class == "Male_biased" & plasma_sex_class == "Male_higher_in_disease",
       sex_direction_concordant := TRUE]
triple[sex_class == "Female_biased" & plasma_sex_class == "Female_higher_in_disease",
       sex_direction_concordant := TRUE]
triple[sex_class == "Male_biased" & plasma_sex_class == "Female_higher_in_disease",
       sex_direction_concordant := FALSE]
triple[sex_class == "Female_biased" & plasma_sex_class == "Male_higher_in_disease",
       sex_direction_concordant := FALSE]

setorder(triple, -triple_concordant, -pp4_max, na.last = TRUE)

out_path <- file.path(OUT_DIR, "sex_coloc_olink_triple.csv")
fwrite(triple, out_path)
log_msg(sprintf("Wrote %s (%d proteins, %d triple-concordant)",
                out_path, nrow(triple), sum(triple$triple_concordant, na.rm = TRUE)))

# QC summary
qc <- data.table(
  metric = c("olink_subjects",
             "subjects_matched_to_disease",
             "subjects_with_sex",
             "sex_source",
             "proteins_total",
             "proteins_with_rna_sex_class",
             "proteins_with_disease_sig",
             "proteins_with_sex_x_disease_sig",
             "proteins_with_coloc",
             "triple_concordant",
             "direction_concordant"),
  value  = c(ncol(npx),
             n_match,
             if (!is.null(sex_lookup)) sum(!is.na(meta_used$sex)) else 0L,
             sex_source,
             nrow(npx),
             sum(!is.na(protein_panel$sex_class)),
             sum(plasma_res$padj_disease < 0.05, na.rm = TRUE),
             if (!is.null(sex_lookup))
               sum(plasma_res$padj_sex_interaction < 0.1, na.rm = TRUE) else 0L,
             sum(triple$has_coloc, na.rm = TRUE),
             sum(triple$triple_concordant, na.rm = TRUE),
             sum(triple$sex_direction_concordant == TRUE, na.rm = TRUE))
)
fwrite(qc, file.path(OUT_DIR, "sex_coloc_olink_triple_qc.csv"))
log_msg("QC table written")

log_msg("218a complete")
