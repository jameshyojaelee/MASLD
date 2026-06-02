#!/usr/bin/env Rscript
# 218c_sex_coloc_pathway_gsea.R — A10 — Pathway-level sex × COLOC fgsea
#
# Joint signed score = signed_sex_interaction_tstat × pp4_susie
# Rank genes by score; run fgsea against MSigDB Hallmark + C2:CP:REACTOME + NMF k=6
# program gene sets.
#
# Outputs:
#   RNA-seq/results/stratified_causal/sex_coloc_pathway_gsea.csv

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n", sep = " ")
log_msg("218c sex × COLOC pathway GSEA starting")

# ─── 1. Load inputs ────────────────────────────────────────────────────────────
log_msg("Loading sex_interaction_dream + sex_deg_classification")
sex_interaction <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_interaction_dream.csv"))
sex_interaction[, ensembl_base := sub("\\..*", "", gene)]

# v3 mashr Bayesian preferred; v2 fallback. v3 CSV provides v2-compatible `sex_class` alias.
sex_v3_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/sex_deg_classification_v3.csv")
sex_v2_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
sex_class_file <- if (file.exists(sex_v3_path)) sex_v3_path else sex_v2_path
log_msg(sprintf("sex_class source: %s", basename(dirname(sex_class_file))))
sex_class <- fread(sex_class_file)
sex_class[, ensembl_base := sub("\\..*", "", gene)]

# ENSEMBL → SYMBOL
gencode <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
sex_interaction <- merge(sex_interaction, gencode[, .(ensembl_base, gene_name)],
                         by = "ensembl_base", all.x = TRUE)
sex_interaction <- sex_interaction[!is.na(gene_name) & nchar(gene_name) > 0]

log_msg(sprintf("Sex interaction: %d genes with t and symbol", nrow(sex_interaction)))

# ─── 2. Load COLOC PP4 ─────────────────────────────────────────────────────────
coloc <- fread(file.path(BASE,
   "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc[, pp4_susie := coloc_best_susie_pp4]
# Fall back to ABF where SuSiE missing
coloc[is.na(pp4_susie), pp4_susie := coloc_best_pp4]
coloc <- coloc[!is.na(pp4_susie) & gene != ""]
log_msg(sprintf("COLOC PP4: %d genes", nrow(coloc)))

# ─── 3. Joint signed score ─────────────────────────────────────────────────────
# signed_score = sign(t) * |t| * pp4   (uses t-statistic from sex_interaction)
merged <- merge(sex_interaction[, .(symbol = gene_name, t)],
                coloc[, .(symbol = gene, pp4_susie)],
                by = "symbol", all = TRUE)

# Genes without COLOC -> treat pp4 as 0 (no genetic weight; still rank by sex t)
merged[is.na(pp4_susie), pp4_susie := 0]
merged[is.na(t), t := 0]
merged[, signed_score := t * pp4_susie]

# For genes where both signed_score == 0 -> use t alone (small jitter to avoid ties)
merged[signed_score == 0, signed_score := t * 1e-6]

# Aggregate duplicates: take max abs
merged <- merged[, .(signed_score = signed_score[which.max(abs(signed_score))]),
                 by = symbol]

# Build ranking vector
ranks <- setNames(merged$signed_score, merged$symbol)
ranks <- ranks[!is.na(ranks) & is.finite(ranks)]
ranks <- sort(ranks, decreasing = TRUE)
log_msg(sprintf("Ranked %d genes", length(ranks)))
log_msg(sprintf("Range: [%.3g, %.3g]", min(ranks), max(ranks)))

# ─── 4. Build gene sets ────────────────────────────────────────────────────────-
log_msg("Loading MSigDB Hallmark + Reactome")
mh <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
               error = function(e) {
                 # Old msigdbr API
                 msigdbr::msigdbr(species = "Homo sapiens", category = "H")
               })
mr <- tryCatch(msigdbr(species = "Homo sapiens", collection = "C2",
                       subcollection = "CP:REACTOME"),
               error = function(e) {
                 msigdbr::msigdbr(species = "Homo sapiens", category = "C2",
                                  subcategory = "CP:REACTOME")
               })

# msigdbr v10 returns gs_name, gene_symbol; older versions same
gs_name_col <- intersect(c("gs_name", "set_name"), colnames(mh))[1]
gs_gene_col <- intersect(c("gene_symbol", "ensembl_gene", "human_gene_symbol"),
                          colnames(mh))[1]
log_msg(sprintf("msigdbr column names: %s / %s", gs_name_col, gs_gene_col))

split_to_list <- function(df, name_col, gene_col) {
  setDT(df)
  split(df[[gene_col]], df[[name_col]])
}
hallmark <- split_to_list(mh, gs_name_col, gs_gene_col)
reactome <- split_to_list(mr, gs_name_col, gs_gene_col)
log_msg(sprintf("Hallmark sets: %d, Reactome sets: %d",
                length(hallmark), length(reactome)))

# NMF k=6 program top loadings
log_msg("Loading NMF k=6 program gene sets")
nmf_atlas <- fread(file.path(BASE, "RNA-seq/results/subtypes/nmf_ksweep_program_atlas.csv"))
nmf6 <- nmf_atlas[k == 6]
nmf_sets <- setNames(
  lapply(nmf6$top10_genes, function(s) {
    g <- strsplit(s, ",")[[1]]
    g <- trimws(gsub("\"", "", g))
    g <- g[nchar(g) > 0 & !grepl("^ENSG", g)]
    g
  }),
  paste0("NMF_k6_", nmf6$program, "_", gsub("[^A-Za-z0-9]+", "_", nmf6$label))
)
nmf_sets <- nmf_sets[lengths(nmf_sets) >= 5]
log_msg(sprintf("NMF k=6 sets: %d (with ≥5 genes)", length(nmf_sets)))

# Filter to sets with ≥10 genes overlapping our ranking
ranked_genes <- names(ranks)
filter_sets <- function(sets, min_overlap = 10) {
  ok <- vapply(sets, function(s) sum(s %in% ranked_genes) >= min_overlap,
               logical(1))
  sets[ok]
}
hallmark <- filter_sets(hallmark, 10)
reactome <- filter_sets(reactome, 10)
nmf_sets <- filter_sets(nmf_sets, 3)   # programs only 15 genes, allow ≥3
log_msg(sprintf("After overlap filter: hallmark=%d, reactome=%d, nmf=%d",
                length(hallmark), length(reactome), length(nmf_sets)))

# ─── 5. fgsea ─────────────────────────────────────────────────────────────────-
run_fgsea <- function(sets, label) {
  if (length(sets) == 0) return(NULL)
  set.seed(42)
  fg <- fgsea(pathways = sets, stats = ranks,
              minSize = 3, maxSize = 500,
              eps = 0, nPermSimple = 10000)
  if (nrow(fg) == 0) return(NULL)
  fg[, collection := label]
  # Format leadingEdge as semicolon-separated
  fg[, leading_edge := vapply(leadingEdge,
                              function(x) paste(head(x, 25), collapse = ";"),
                              character(1))]
  fg[, leadingEdge := NULL]
  fg[]
}

log_msg("Running fgsea: Hallmark")
res_h <- run_fgsea(hallmark, "MSigDB_Hallmark")
log_msg("Running fgsea: Reactome")
res_r <- run_fgsea(reactome, "MSigDB_C2_Reactome")
log_msg("Running fgsea: NMF k=6")
res_n <- run_fgsea(nmf_sets, "NMF_k6_programs")

res <- rbindlist(list(res_h, res_r, res_n), use.names = TRUE, fill = TRUE)
# Recompute padj within each collection
res[, padj_collection := p.adjust(pval, "BH"), by = collection]
res[, abs_NES := abs(NES)]
setorder(res, padj_collection, -abs_NES)
res[, abs_NES := NULL]

out_path <- file.path(OUT_DIR, "sex_coloc_pathway_gsea.csv")
fwrite(res, out_path)
log_msg(sprintf("Wrote %s (%d pathways tested)", out_path, nrow(res)))

# Summary
sig <- res[padj_collection < 0.05]
log_msg(sprintf("Significant (padj_collection<0.05): %d", nrow(sig)))
if (nrow(sig) > 0) {
  log_msg("Top 10 by |NES|:")
  sig_print <- copy(sig)
  sig_print[, abs_NES := abs(NES)]
  print(head(sig_print[order(-abs_NES), .(pathway, collection, NES, padj_collection,
                                          size)], 10))
}

log_msg("218c complete")
