#!/usr/bin/env Rscript
# 218c2_pathway_gsea_stratified.R — A10b — Pathway-level sex × COLOC fgsea
#                                  STRATIFIED BY GWAS FAMILY
#
# Motivation (from A10b audit):
#   The original A10 (218c) finding REACTOME_TRANSLATION (padj=0.011) was driven by
#   liver-enzyme GWAS (UKBB_ALT/AST/GGT) coloc hits, NOT disease-endpoint GWAS.
#   Enzyme GWAS (ALT/AST/GGT) inflate coloc for ANY liver-expressed gene because
#   elevated liver enzymes are a tissue-damage release marker — they tag
#   hepatocyte-expressed proteins, not MASLD-causal biology. So translation/
#   ribosomal coloc with enzymes likely reflects tissue-damage marker biology,
#   not MASLD causal architecture.
#
# This script recomputes the joint sex-interaction × max-PP4 score under three
# GWAS strata, then re-runs fgsea against Hallmark + Reactome:
#
#   stratum 1 — DISEASE-ONLY (NAFLD / NASH / HCC / Cirrhosis endpoints)
#   stratum 2 — ENZYME+PDFF (liver enzyme + hepatic fat imaging markers)
#   stratum 3 — ALL GWAS    (canonical reference; matches 218c)
#
# GWAS-FAMILY ASSIGNMENT (documented per task spec):
#
#   ── disease (n=13 GWAS) ──
#     2019_31311600_NAFLD_EUR, 2020_32298765_NAFLD_EUR,
#     2020_32514122_Cirrhosis_EAS, 2020_32514122_HCC_EAS,
#     2021_34841290_NAFLD_EUR,
#     2023_36280732_NAFLD_deCode_EUR,
#     2023_36280732_NAFLD_Intermountain_EUR,
#     2023_36280732_NAFLD_UKBB_EUR,
#     FinnGen_HCC, FinnGen_NAFLD, FinnGen_NASH,
#     Ghouse_Cirrhosis, Ghouse_HCC
#
#   ── enzyme_pdff (n=12 GWAS) ──   # task spec lumps PDFF with enzyme markers
#     BBJ_ALT, BBJ_AST, BBJ_GGT,
#     UKBB_ALT, UKBB_AST, UKBB_GGT,
#     PanUKBB_AFR_ALT, PanUKBB_AFR_AST, PanUKBB_AFR_GGT,
#     PanUKBB_CSA_ALT, PanUKBB_CSA_AST, PanUKBB_CSA_GGT,
#     2021_34128465_PDFF_EUR, 2021_34957434_PDFF_EUR,
#     2022_36402844_PDFF_EUR     # PDFF = MRI hepatic fat fraction → enzyme-class marker
#
# Outputs:
#   RNA-seq/results/stratified_causal/sex_coloc_pathway_gsea_disease_only.csv
#   RNA-seq/results/stratified_causal/sex_coloc_pathway_gsea_enzyme_only.csv
#   RNA-seq/results/stratified_causal/sex_coloc_pathway_gsea_comparison.csv
#   outputs/team_A/A10b_pathway_gsea_stratified.md

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR  <- file.path(BASE, "RNA-seq/results/stratified_causal")
REPORT_DIR <- file.path(BASE, "outputs/team_A")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n", sep = " ")
log_msg("218c2 stratified sex × COLOC pathway GSEA starting")

# ─── GWAS family classification ────────────────────────────────────────────────
GWAS_DISEASE <- c(
  "2019_31311600_NAFLD_EUR", "2020_32298765_NAFLD_EUR",
  "2020_32514122_Cirrhosis_EAS", "2020_32514122_HCC_EAS",
  "2021_34841290_NAFLD_EUR",
  "2023_36280732_NAFLD_deCode_EUR",
  "2023_36280732_NAFLD_Intermountain_EUR",
  "2023_36280732_NAFLD_UKBB_EUR",
  "FinnGen_HCC", "FinnGen_NAFLD", "FinnGen_NASH",
  "Ghouse_Cirrhosis", "Ghouse_HCC"
)
GWAS_ENZYME <- c(
  "BBJ_ALT", "BBJ_AST", "BBJ_GGT",
  "UKBB_ALT", "UKBB_AST", "UKBB_GGT",
  "PanUKBB_AFR_ALT", "PanUKBB_AFR_AST", "PanUKBB_AFR_GGT",
  "PanUKBB_CSA_ALT", "PanUKBB_CSA_AST", "PanUKBB_CSA_GGT",
  "2021_34128465_PDFF_EUR", "2021_34957434_PDFF_EUR",
  "2022_36402844_PDFF_EUR"
)

# ─── 1. Load inputs ────────────────────────────────────────────────────────────
log_msg("Loading sex_interaction_dream")
sex_interaction <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_interaction_dream.csv"))
sex_interaction[, ensembl_base := sub("\\..*", "", gene)]

gencode <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
sex_interaction <- merge(sex_interaction, gencode[, .(ensembl_base, gene_name)],
                         by = "ensembl_base", all.x = TRUE)
sex_interaction <- sex_interaction[!is.na(gene_name) & nchar(gene_name) > 0]
log_msg(sprintf("Sex interaction: %d genes with t and symbol", nrow(sex_interaction)))

log_msg("Loading susie_coloc_all_gwas (per-gene-per-GWAS)")
coloc_long <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
log_msg(sprintf("Long-form coloc: %d rows, %d unique GWAS",
                nrow(coloc_long), uniqueN(coloc_long$gwas_name)))

# Derive per-row best PP4: prefer SuSiE, fall back to ABF (matches 218c logic)
coloc_long[, pp4 := PP.H4.susie]
coloc_long[is.na(pp4), pp4 := PP.H4.abf]
coloc_long <- coloc_long[!is.na(pp4) & gene != ""]

# ─── 2. Sanity-check classification covers all observed GWAS ───────────────────
obs_gwas <- sort(unique(coloc_long$gwas_name))
unclassified <- setdiff(obs_gwas, c(GWAS_DISEASE, GWAS_ENZYME))
if (length(unclassified) > 0) {
  log_msg(sprintf("WARNING: %d unclassified GWAS: %s",
                  length(unclassified), paste(unclassified, collapse = ", ")))
}
log_msg(sprintf("Disease GWAS observed: %d / %d declared",
                sum(GWAS_DISEASE %in% obs_gwas), length(GWAS_DISEASE)))
log_msg(sprintf("Enzyme+PDFF GWAS observed: %d / %d declared",
                sum(GWAS_ENZYME %in% obs_gwas), length(GWAS_ENZYME)))

# ─── 3. MSigDB gene sets (Hallmark + Reactome) ─────────────────────────────────
log_msg("Loading MSigDB Hallmark + Reactome")
mh <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
               error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
mr <- tryCatch(msigdbr(species = "Homo sapiens", collection = "C2",
                       subcollection = "CP:REACTOME"),
               error = function(e) msigdbr::msigdbr(species = "Homo sapiens",
                                                    category = "C2",
                                                    subcategory = "CP:REACTOME"))
gs_name_col <- intersect(c("gs_name", "set_name"), colnames(mh))[1]
gs_gene_col <- intersect(c("gene_symbol", "ensembl_gene", "human_gene_symbol"),
                          colnames(mh))[1]
log_msg(sprintf("msigdbr cols: %s / %s", gs_name_col, gs_gene_col))

split_to_list <- function(df, name_col, gene_col) {
  setDT(df)
  split(df[[gene_col]], df[[name_col]])
}
hallmark <- split_to_list(mh, gs_name_col, gs_gene_col)
reactome <- split_to_list(mr, gs_name_col, gs_gene_col)
log_msg(sprintf("Hallmark sets: %d, Reactome sets: %d",
                length(hallmark), length(reactome)))

# ─── 4. Per-stratum: compute max_pp4 per gene → joint score → fgsea ────────────
run_stratum <- function(stratum_label, gwas_keep) {
  log_msg(sprintf("─── stratum: %s (n=%d GWAS) ───", stratum_label, length(gwas_keep)))
  sub <- coloc_long[gwas_name %in% gwas_keep]
  if (nrow(sub) == 0) {
    log_msg("  no rows — skipping")
    return(NULL)
  }
  per_gene <- sub[, .(pp4_max = max(pp4, na.rm = TRUE),
                      n_gwas_tested = .N,
                      best_gwas = gwas_name[which.max(pp4)]),
                  by = .(gene, ensembl)]
  log_msg(sprintf("  per-gene rows: %d (max PP4 across stratum)", nrow(per_gene)))

  merged <- merge(sex_interaction[, .(symbol = gene_name, t)],
                  per_gene[, .(symbol = gene, pp4_max)],
                  by = "symbol", all = TRUE)
  merged[is.na(pp4_max), pp4_max := 0]
  merged[is.na(t), t := 0]
  merged[, signed_score := t * pp4_max]
  # break ties via t (matches 218c)
  merged[signed_score == 0, signed_score := t * 1e-6]
  merged <- merged[, .(signed_score = signed_score[which.max(abs(signed_score))]),
                   by = symbol]
  ranks <- setNames(merged$signed_score, merged$symbol)
  ranks <- ranks[!is.na(ranks) & is.finite(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  log_msg(sprintf("  ranked %d genes; range [%.3g, %.3g]",
                  length(ranks), min(ranks), max(ranks)))

  filter_sets <- function(sets, min_overlap = 10) {
    ok <- vapply(sets, function(s) sum(s %in% names(ranks)) >= min_overlap, logical(1))
    sets[ok]
  }
  hk <- filter_sets(hallmark, 10)
  rk <- filter_sets(reactome, 10)

  run_fgsea <- function(sets, label) {
    if (length(sets) == 0) return(NULL)
    set.seed(42)
    fg <- fgsea(pathways = sets, stats = ranks,
                minSize = 3, maxSize = 500,
                eps = 0, nPermSimple = 10000)
    if (nrow(fg) == 0) return(NULL)
    fg[, collection := label]
    fg[, leading_edge := vapply(leadingEdge,
                                function(x) paste(head(x, 25), collapse = ";"),
                                character(1))]
    fg[, leadingEdge := NULL]
    fg[]
  }

  log_msg("  running fgsea Hallmark")
  res_h <- run_fgsea(hk, "MSigDB_Hallmark")
  log_msg("  running fgsea Reactome")
  res_r <- run_fgsea(rk, "MSigDB_C2_Reactome")
  res <- rbindlist(list(res_h, res_r), use.names = TRUE, fill = TRUE)
  res[, padj_collection := p.adjust(pval, "BH"), by = collection]
  res[, stratum := stratum_label]
  res[, n_gwas_in_stratum := length(gwas_keep)]
  res[, abs_NES := abs(NES)]
  setorder(res, padj_collection, -abs_NES)
  res[, abs_NES := NULL]
  res[]
}

res_disease  <- run_stratum("disease_only", GWAS_DISEASE)
res_enzyme   <- run_stratum("enzyme_only",  GWAS_ENZYME)
res_all      <- run_stratum("all_gwas",     c(GWAS_DISEASE, GWAS_ENZYME))

# ─── 5. Write per-stratum + comparison CSVs ────────────────────────────────────
write_stratum <- function(df, fn) {
  if (is.null(df) || nrow(df) == 0) {
    log_msg(sprintf("  WARNING: empty result for %s", fn)); return(invisible(NULL))
  }
  fwrite(df, file.path(OUT_DIR, fn))
  log_msg(sprintf("  wrote %s (%d rows)", fn, nrow(df)))
}
write_stratum(res_disease, "sex_coloc_pathway_gsea_disease_only.csv")
write_stratum(res_enzyme,  "sex_coloc_pathway_gsea_enzyme_only.csv")
write_stratum(res_all,     "sex_coloc_pathway_gsea_all_gwas.csv")

# Side-by-side comparison
build_compare <- function() {
  cols <- c("pathway", "collection", "NES", "pval", "padj_collection", "size",
            "leading_edge")
  tag <- function(df, suffix) {
    if (is.null(df)) return(NULL)
    out <- df[, ..cols]
    setnames(out, c("pathway", "collection",
                    paste0("NES_", suffix),
                    paste0("pval_", suffix),
                    paste0("padj_collection_", suffix),
                    paste0("size_", suffix),
                    paste0("leading_edge_", suffix)))
    out
  }
  d <- tag(res_disease, "disease")
  e <- tag(res_enzyme,  "enzyme")
  a <- tag(res_all,     "all")
  if (is.null(d) || is.null(e) || is.null(a)) {
    return(rbindlist(list(d, e, a), fill = TRUE))
  }
  m <- merge(d, e, by = c("pathway", "collection"), all = TRUE)
  m <- merge(m, a, by = c("pathway", "collection"), all = TRUE)
  m[, min_padj := pmin(padj_collection_disease, padj_collection_enzyme,
                        padj_collection_all, na.rm = TRUE)]
  setorder(m, min_padj)
  m[]
}
cmp <- build_compare()
fwrite(cmp, file.path(OUT_DIR, "sex_coloc_pathway_gsea_comparison.csv"))
log_msg(sprintf("Wrote comparison: %d pathways", nrow(cmp)))

# ─── 6. Markdown report ────────────────────────────────────────────────────────
md_lines <- c(
  "# A10b — Pathway GSEA Sensitivity Stratified by GWAS Family",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## GWAS family assignment",
  "",
  "**Disease GWAS** (NAFLD / NASH / HCC / Cirrhosis endpoints):",
  sprintf("- %s", paste(GWAS_DISEASE, collapse = ", ")),
  "",
  "**Enzyme + PDFF GWAS** (ALT / AST / GGT + MRI hepatic fat):",
  sprintf("- %s", paste(GWAS_ENZYME, collapse = ", ")),
  "",
  "## Stratum-level top-10 hits",
  ""
)
top10 <- function(df, hdr) {
  if (is.null(df) || nrow(df) == 0) return(c(sprintf("### %s (no results)", hdr), ""))
  d <- copy(df); d[, abs_NES := abs(NES)]
  setorder(d, padj_collection, -abs_NES)
  rows <- head(d, 10)
  lines <- c(sprintf("### %s", hdr), "",
             "| Pathway | Collection | NES | pval | padj | size |",
             "|---|---|---:|---:|---:|---:|")
  for (i in seq_len(nrow(rows))) {
    lines <- c(lines, sprintf("| %s | %s | %.3f | %.3g | %.3g | %d |",
                              rows$pathway[i], rows$collection[i],
                              rows$NES[i], rows$pval[i], rows$padj_collection[i],
                              rows$size[i]))
  }
  c(lines, "")
}
md_lines <- c(md_lines,
              top10(res_disease, "Disease-only stratum (NAFLD/NASH/HCC/Cirrhosis)"),
              top10(res_enzyme,  "Enzyme+PDFF stratum (ALT/AST/GGT/PDFF)"),
              top10(res_all,     "All-GWAS reference (matches canonical A10)"))

# REACTOME_TRANSLATION focused check
get_translation <- function(df, lab) {
  if (is.null(df)) return(sprintf("- %s: <NULL>", lab))
  hit <- df[pathway == "REACTOME_TRANSLATION"]
  if (nrow(hit) == 0) return(sprintf("- %s: NOT TESTED (no row)", lab))
  sprintf("- %s: NES=%.3f, pval=%.3g, padj_collection=%.3g, size=%d, leading_edge_n=%d",
          lab, hit$NES, hit$pval, hit$padj_collection, hit$size,
          length(strsplit(hit$leading_edge, ";")[[1]]))
}
md_lines <- c(md_lines,
              "## REACTOME_TRANSLATION across strata",
              "",
              get_translation(res_disease, "disease_only"),
              get_translation(res_enzyme,  "enzyme_only"),
              get_translation(res_all,     "all_gwas"),
              "")

# Verdict
verdict <- (function() {
  d <- res_disease[pathway == "REACTOME_TRANSLATION"]
  e <- res_enzyme[pathway  == "REACTOME_TRANSLATION"]
  if (nrow(d) == 0 || nrow(e) == 0) return("INDETERMINATE: pathway absent in one or both strata.")
  d_padj <- d$padj_collection; e_padj <- e$padj_collection
  survives_disease <- !is.na(d_padj) && d_padj < 0.05
  flagged_enzyme   <- !is.na(e_padj) && e_padj < 0.05
  if (survives_disease && flagged_enzyme)
    return(sprintf("BOTH STRATA SIGNIFICANT: disease padj=%.3g, enzyme padj=%.3g. Translation signal is supported by disease-endpoint GWAS, NOT exclusively driven by enzyme-marker leakage.", d_padj, e_padj))
  if (survives_disease && !flagged_enzyme)
    return(sprintf("DISEASE-ONLY SURVIVES, ENZYME N.S.: disease padj=%.3g, enzyme padj=%.3g. Translation signal is disease-anchored (not enzyme leakage).", d_padj, e_padj))
  if (!survives_disease && flagged_enzyme)
    return(sprintf("DISEASE FAILS, ENZYME SIGNIFICANT: disease padj=%.3g, enzyme padj=%.3g. The canonical A10 REACTOME_TRANSLATION hit was driven by enzyme-marker COLOC (tissue-damage release signal), NOT MASLD-causal architecture. Audit concern CONFIRMED.", d_padj, e_padj))
  return(sprintf("BOTH STRATA N.S.: disease padj=%.3g, enzyme padj=%.3g. Translation signal does not survive either restricted analysis; original A10 may reflect cross-stratum aggregation effect.", d_padj, e_padj))
})()
md_lines <- c(md_lines, "## Verdict", "", verdict, "")

writeLines(md_lines, file.path(REPORT_DIR, "A10b_pathway_gsea_stratified.md"))
log_msg(sprintf("Wrote report: %s",
                file.path(REPORT_DIR, "A10b_pathway_gsea_stratified.md")))

log_msg("218c2 complete")
