#!/usr/bin/env Rscript
# 202_geneset_enrichment.R — Competitive Gene-Set GWAS Enrichment
#
# Tests whether dream DEGs (and subsets) are enriched for GWAS gene-level
# associations vs genome-wide background. Produces stratified QQ plot data.
#
# Tests:
#   1. All DEGs, upregulated, downregulated, top 5%
#   2. Hepatocyte-specific DEGs (from deconvolution attribution)
#   3. Conserved genes (cross-species)
#   4. Sex-divergent DEGs
#   5. Per-cell-type DEGs (from scRNA pseudobulk)
#   6. Stratified QQ: COLOC PP.H4 by DEG significance strata
#
# Inputs:
#   - COLOC gene-level: GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#   - TWAS results: RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv
#   - Dream DEGs: RNA-seq/Human/.../results/integration/dream_results.csv
#   - Deconv attribution: RNA-seq/results/causal_inference/deconv_attribution_scores.csv
#   - Cross-species concordance: Analysis/Cross_Species_Concordance/ (Conserved)
#   - Sex DEGs: RNA-seq/Human/.../results/integration/sex_deg_classification.csv
#   - scRNA DE: Analysis/SingleCell/results_gpu_v2/
#
# Outputs:
#   - RNA-seq/results/gwas_rna_integration/geneset_enrichment_results.csv
#   - RNA-seq/results/gwas_rna_integration/stratified_qq_data.csv
#
# SLURM: cpu partition, 4 CPUs, 32GB, 48h

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

outdir <- file.path(BASE, "RNA-seq/results/gwas_rna_integration")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Load gene-level GWAS scores
# ===========================================================================
cat("Loading gene-level GWAS scores...\n")

coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]

twas <- fread(file.path(BASE,
  "RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv"))
twas_best <- twas[, .(twas_z = zscore[which.min(pvalue)], twas_p = min(pvalue)),
                  by = .(gene = gene_name)]

# Merge
gwas_genes <- merge(
  coloc[, .(gene, coloc_pp4 = coloc_best_pp4, coloc_n = coloc_n_gwas_h4_05)],
  twas_best, by = "gene", all = TRUE
)
cat("  Gene-level GWAS scores:", nrow(gwas_genes), "genes\n")

# Build ENSEMBL → symbol mapping from COLOC data
ensembl_to_symbol <- coloc[gene != "" & ensembl != "", .(ensembl, gene)]
ensembl_to_symbol[, ensembl_base := sub("\\.\\d+$", "", ensembl)]

map_ensembl_to_symbol <- function(ensembl_ids) {
  # Strip version suffix for matching
  base_ids <- sub("\\.\\d+$", "", ensembl_ids)
  mapped <- ensembl_to_symbol$gene[match(base_ids, ensembl_to_symbol$ensembl_base)]
  mapped[is.na(mapped)] <- ensembl_ids[is.na(mapped)]  # keep original if no map
  return(mapped)
}

# ===========================================================================
# 2. Load gene sets
# ===========================================================================
cat("\nLoading gene sets...\n")

# Dream DEGs (uses ENSEMBL IDs — map to symbols)
dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv"))
dream[, ensembl_base := sub("\\.\\d+$", "", gene)]
dream <- merge(dream, ensembl_to_symbol[, .(ensembl_base, symbol = gene)],
               by = "ensembl_base", all.x = TRUE)
dream[!is.na(symbol), gene_symbol := symbol]
dream[is.na(symbol), gene_symbol := gene]
cat("  Dream: mapped", sum(!is.na(dream$symbol)), "/", nrow(dream), "to symbols\n")

# NOTE: Uses padj<0.1 for competitive gene-set enrichment (more permissive for GSEA background).
# Primary DEGs defined at padj<0.05 + |logFC|>0.5 in Script 05b.
dream_sig <- dream[padj < 0.1]

gene_sets <- list(
  "All_DEGs" = dream_sig$gene_symbol,
  "DEGs_Up" = dream_sig[logFC > 0]$gene_symbol,
  "DEGs_Down" = dream_sig[logFC < 0]$gene_symbol,
  "DEGs_Top5pct" = dream[order(-abs(t))][1:round(nrow(dream) * 0.05)]$gene_symbol,
  "DEGs_LFC_gt_08" = dream_sig[abs(logFC) > 0.8]$gene_symbol
)

# Hepatocyte-intrinsic DEGs
deconv_file <- file.path(BASE,
  "RNA-seq/results/causal_inference/deconv_attribution_scores.csv")
if (file.exists(deconv_file)) {
  deconv <- fread(deconv_file)
  # Map ENSEMBL IDs to symbols
  gene_sets[["Hepatocyte_intrinsic"]] <- map_ensembl_to_symbol(
    deconv[category == "Hepatocyte_intrinsic"]$gene)
  gene_sets[["Composition_driven"]] <- map_ensembl_to_symbol(
    deconv[category == "Composition_driven"]$gene)
  cat("  Hepatocyte-intrinsic:", length(gene_sets[["Hepatocyte_intrinsic"]]), "genes\n")
}

# Conserved
cc_file <- file.path(BASE,
  "Analysis/Cross_Species_Concordance/results/gene_concordance_per_gene.csv")
if (file.exists(cc_file)) {
  cc <- fread(cc_file)
  if ("best_category" %in% names(cc)) {
    gene_sets[["Conserved"]] <- cc[best_category == "Conserved"]$human_symbol
  }
  cat("  Conserved:", length(gene_sets[["Conserved"]]), "genes\n")
} else {
  cat("  WARNING: Cross-species concordance file not found\n")
}

# Sex-divergent DEGs
sex_file <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
if (file.exists(sex_file)) {
  sex <- fread(sex_file)
  cls_col <- intersect(names(sex), c("sex_class", "class"))[1]
  if (!is.na(cls_col)) {
    # Map ENSEMBL IDs to symbols
    # Auto-detect v2 (interaction-based) vs v1 (stratified) class names
    sex_cls_202 <- unique(sex[[cls_col]])
    use_v2_202 <- "Female_biased" %in% sex_cls_202
    fem_lbl_202 <- if (use_v2_202) "Female_biased" else "Female_specific"
    mal_lbl_202 <- if (use_v2_202) "Male_biased"   else "Male_specific"
    div_lbl_202 <- if (use_v2_202) "Divergent"      else "Sex_divergent"
    gene_sets[["Female_specific"]] <- map_ensembl_to_symbol(
      sex[get(cls_col) == fem_lbl_202]$gene)
    gene_sets[["Male_specific"]] <- map_ensembl_to_symbol(
      sex[get(cls_col) == mal_lbl_202]$gene)
    gene_sets[["Sex_divergent"]] <- map_ensembl_to_symbol(
      sex[get(cls_col) == div_lbl_202]$gene)
  }
} else {
  cat("  WARNING: Sex DEG classification not found\n")
}

# scRNA cell-type DEGs
sc_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
sc_files <- list.files(sc_dir, pattern = "_de\\.csv$", full.names = TRUE)
for (f in sc_files) {
  dt <- fread(f)
  ct <- gsub("_de\\.csv$", "", basename(f))
  # scRNA pseudobulk uses lowercase `lfc` (not `logFC` as in the dream bulk CSV)
  lfc_col <- if ("logFC" %in% names(dt)) "logFC" else "lfc"
  ct_sig <- dt[padj < 0.05 & abs(get(lfc_col)) > 0.25]$gene
  if (length(ct_sig) > 10) {
    gene_sets[[paste0("sc_", ct)]] <- ct_sig
  }
}

cat("  Total gene sets:", length(gene_sets), "\n")
for (gs_name in names(gene_sets)) {
  cat("    ", gs_name, ":", length(gene_sets[[gs_name]]), "\n")
}

# ===========================================================================
# 3. Competitive enrichment test for each gene set
# ===========================================================================
cat("\n=== Running competitive enrichment tests ===\n")

run_enrichment <- function(gene_set_name, gene_set_genes, gwas_dt) {
  in_set <- gwas_dt$gene %in% gene_set_genes

  if (sum(in_set) < 5) {
    return(data.table(gene_set = gene_set_name, n_in_set = sum(in_set),
                      note = "Too few genes"))
  }

  # Wilcoxon rank-sum (one-sided: set > background)
  wt_coloc <- tryCatch(
    wilcox.test(gwas_dt$coloc_pp4[in_set], gwas_dt$coloc_pp4[!in_set],
                alternative = "greater"),
    error = function(e) list(p.value = NA, statistic = NA)
  )

  wt_twas <- tryCatch(
    wilcox.test(abs(gwas_dt$twas_z[in_set & !is.na(gwas_dt$twas_z)]),
                abs(gwas_dt$twas_z[!in_set & !is.na(gwas_dt$twas_z)]),
                alternative = "greater"),
    error = function(e) list(p.value = NA, statistic = NA)
  )

  # Effect sizes
  mean_in <- mean(gwas_dt$coloc_pp4[in_set], na.rm = TRUE)
  mean_out <- mean(gwas_dt$coloc_pp4[!in_set], na.rm = TRUE)

  # Fisher exact: proportion of COLOC genes
  n_coloc_in <- sum(gwas_dt$coloc_pp4[in_set] > 0.5, na.rm = TRUE)
  n_coloc_out <- sum(gwas_dt$coloc_pp4[!in_set] > 0.5, na.rm = TRUE)
  ft <- fisher.test(matrix(c(
    n_coloc_in, sum(in_set) - n_coloc_in,
    n_coloc_out, sum(!in_set) - n_coloc_out
  ), nrow = 2), alternative = "greater")

  data.table(
    gene_set = gene_set_name,
    n_in_set = sum(in_set),
    n_coloc_in_set = n_coloc_in,
    mean_coloc_in = mean_in,
    mean_coloc_out = mean_out,
    fold_enrichment = mean_in / max(mean_out, 1e-10),
    wilcox_p_coloc = wt_coloc$p.value,
    wilcox_p_twas = wt_twas$p.value,
    fisher_or = ft$estimate,
    fisher_p = ft$p.value
  )
}

enrichment_list <- lapply(names(gene_sets), function(gs_name) {
  run_enrichment(gs_name, gene_sets[[gs_name]], gwas_genes)
})

enrichment_dt <- rbindlist(enrichment_list, fill = TRUE)
# Remove rows with missing test results
enrichment_dt <- enrichment_dt[!is.na(wilcox_p_coloc)]

if (nrow(enrichment_dt) == 0) {
  cat("WARNING: No enrichment tests completed.\n")
} else {
  # FDR correction
  enrichment_dt[, fdr_coloc := p.adjust(wilcox_p_coloc, method = "BH")]
  enrichment_dt[, fdr_twas := p.adjust(wilcox_p_twas, method = "BH")]
}

setorder(enrichment_dt, wilcox_p_coloc)

cat("\n=== Gene-set enrichment results ===\n")
print(enrichment_dt[, .(gene_set, n_in_set, fold_enrichment,
  wilcox_p_coloc, fdr_coloc, fisher_or, fisher_p)])

# ===========================================================================
# 4. Stratified QQ plot data
# ===========================================================================
cat("\n=== Generating stratified QQ plot data ===\n")

# Stratify genes by DEG significance (merge using symbol)
dream_merged <- merge(dream[, .(gene = gene_symbol, t, logFC, padj)],
                      gwas_genes, by = "gene")

# Create 3 strata: top 10% DEGs, middle 80%, bottom 10% (non-DEGs)
dream_merged[, abs_t := abs(t)]
q90 <- quantile(dream_merged$abs_t, 0.90, na.rm = TRUE)
q50 <- quantile(dream_merged$abs_t, 0.50, na.rm = TRUE)

dream_merged[, stratum := fcase(
  abs_t >= q90, "Top 10% DEGs",
  abs_t >= q50, "Middle",
  default = "Bottom 50%"
)]

# For each stratum, compute expected vs observed COLOC PP.H4 distribution
qq_data <- lapply(unique(dream_merged$stratum), function(s) {
  sub <- dream_merged[stratum == s & !is.na(coloc_pp4)]
  n <- nrow(sub)
  if (n == 0) return(NULL)

  obs <- sort(sub$coloc_pp4, decreasing = TRUE)
  exp_rank <- (seq_len(n)) / (n + 1)

  data.table(
    stratum = s,
    observed_coloc = obs,
    expected_rank = exp_rank,
    rank = seq_len(n)
  )
})

qq_dt <- rbindlist(qq_data, fill = TRUE)
cat("  Stratified QQ data:", nrow(qq_dt), "rows\n")

# Lambda (genomic inflation factor) per stratum
lambda_dt <- dream_merged[!is.na(coloc_pp4), .(
  n_genes = .N,
  mean_coloc = mean(coloc_pp4),
  median_coloc = median(coloc_pp4),
  n_coloc_05 = sum(coloc_pp4 > 0.5),
  prop_coloc_05 = sum(coloc_pp4 > 0.5) / .N
), by = stratum]

cat("\nPer-stratum summary:\n")
print(lambda_dt)

# ===========================================================================
# 5. DEG direction concordance with COLOC gene direction
# ===========================================================================
cat("\n=== DEG-COLOC directional analysis ===\n")

# For COLOC genes that are also DEGs, are they directionally concordant?
# COLOC itself is direction-agnostic, but eQTL direction + GWAS direction
# gives expected expression direction
# We just check: are COLOC genes more likely to be DEGs? And which direction?
coloc_degs <- dream_merged[coloc_pp4 > 0.5]
cat("  COLOC genes that are DEGs (padj < 0.1):",
    sum(coloc_degs$padj < 0.1, na.rm = TRUE), "/", nrow(coloc_degs), "\n")
cat("  Direction of COLOC DEGs: Up =",
    sum(coloc_degs$logFC > 0 & coloc_degs$padj < 0.1, na.rm = TRUE),
    ", Down =",
    sum(coloc_degs$logFC < 0 & coloc_degs$padj < 0.1, na.rm = TRUE), "\n")

# ===========================================================================
# 6. Save results
# ===========================================================================
cat("\n=== Saving results ===\n")

fwrite(enrichment_dt, file.path(outdir, "geneset_enrichment_results.csv"))
cat("  Saved geneset_enrichment_results.csv:", nrow(enrichment_dt), "rows\n")

fwrite(qq_dt, file.path(outdir, "stratified_qq_data.csv"))
cat("  Saved stratified_qq_data.csv:", nrow(qq_dt), "rows\n")

fwrite(lambda_dt, file.path(outdir, "deg_stratum_coloc_summary.csv"))
cat("  Saved deg_stratum_coloc_summary.csv:", nrow(lambda_dt), "rows\n")

cat("\nDone.\n")
