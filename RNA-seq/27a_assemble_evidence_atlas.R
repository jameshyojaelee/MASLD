#!/usr/bin/env Rscript
# 27a_assemble_evidence_atlas.R
# Assemble Multi-Evidence Target Atlas
#
# Replaces the weighted composite scoring approach (old Script 27) with a
# unified evidence table containing raw values from all evidence layers.
# NOTE: Mouse (L2) columns are retained for supplementary use but excluded
# from layers_active count — mouse relegated to supplementary (6 main sources).
# No weights, no composite score — provides the multi-dimensional evidence
# database for downstream filtering (Streamlit app or 27b presets).
#
# Output: RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
#
# CLI flags (parsed at top; defaults preserve backwards compat):
#   --include-s8 TRUE|FALSE   Ingest S8 perturbation_modeling columns from
#                             Analysis/Perturbation/results/integration/
#                             d{1..5}_atlas_columns.csv. Adds 14 perturb_*
#                             columns + perturb_n_arms_applicable +
#                             perturb_n_arms_passing + perturb_validation_tier.
#                             Default FALSE. Also settable via INCLUDE_S8 env var.
#   --subset-n N              Limit the atlas to the first N rows before the
#                             post-assembly merges. Intended for dry-run /
#                             smoke tests; NEVER use this for production
#                             atlas rebuilds. Default unset (full atlas).
#                             Also settable via SUBSET_N env var.
#
# Example: Rscript 27a_assemble_evidence_atlas.R --include-s8 TRUE
#          INCLUDE_S8=TRUE SUBSET_N=100 Rscript 27a_assemble_evidence_atlas.R

library(data.table)
library(jsonlite)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUTDIR <- file.path(BASE, "RNA-seq/results/multi_evidence")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# CLI / env-var argument parsing (kept inline; 27a has no other arg surface)
# ---------------------------------------------------------------------------
.args_raw <- commandArgs(trailingOnly = TRUE)
.parse_truthy <- function(x) {
  if (is.null(x) || length(x) == 0) return(FALSE)
  toupper(as.character(x)[1]) %in% c("TRUE", "T", "1", "YES", "Y")
}
.parse_kv_flag <- function(args, names_to_match) {
  # Supports both `--include-s8 TRUE` and `--include-s8=TRUE`
  for (nm in names_to_match) {
    eq_match <- grep(paste0("^", nm, "="), args, value = TRUE)
    if (length(eq_match) > 0) return(sub(paste0("^", nm, "="), "", eq_match[1]))
    sp_idx <- which(args == nm)
    if (length(sp_idx) > 0 && sp_idx[1] < length(args)) {
      return(args[sp_idx[1] + 1])
    }
  }
  NULL
}
.s8_cli <- .parse_kv_flag(.args_raw, c("--include-s8", "--include_s8"))
INCLUDE_S8 <- if (!is.null(.s8_cli)) {
  .parse_truthy(.s8_cli)
} else {
  .parse_truthy(Sys.getenv("INCLUDE_S8", unset = "FALSE"))
}

.subset_cli <- .parse_kv_flag(.args_raw, c("--subset-n", "--subset_n"))
.subset_env <- Sys.getenv("SUBSET_N", unset = "")
.subset_src <- if (!is.null(.subset_cli)) .subset_cli else .subset_env
SUBSET_N <- suppressWarnings(as.integer(.subset_src))
if (length(SUBSET_N) == 0 || is.na(SUBSET_N) || SUBSET_N <= 0) SUBSET_N <- NA_integer_

cat("=== Multi-Evidence Target Atlas Assembly ===\n")
cat("Base directory:", BASE, "\n")
cat("Flags: INCLUDE_S8 =", INCLUDE_S8,
    "| SUBSET_N =", if (is.na(SUBSET_N)) "NA (full atlas)" else SUBSET_N, "\n\n")

# ================================================================
# Layer 1: Human DE evidence
# ================================================================
cat("Loading Layer 1: Human consensus DEGs...\n")
consensus <- fread(file.path(RDIR, "consensus_degs.csv"))
consensus[, ensembl_clean := sub("\\..*", "", gene)]

# Use the ortholog comparison to get HGNC symbols
ortho <- fread(file.path(RDIR, "human_mouse_ortholog_comparison.csv"))
ortho[, ensembl_clean := sub("\\..*", "", gene_base)]
symbol_map <- unique(ortho[!is.na(human_symbol) & human_symbol != "", .(ensembl_clean, human_symbol)])
symbol_map <- symbol_map[!duplicated(ensembl_clean)]

# Supplement with GENCODE symbols for genes missing from ortholog map
# (captures lncRNAs, miRNAs, and other non-protein-coding genes)
gencode_gtf <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
gencode_cmd <- paste0("zcat ", gencode_gtf,
  " | awk -F'\\t' '$3==\"gene\"'",
  " | sed 's/.*gene_id \"\\([^\"]*\\)\".*gene_name \"\\([^\"]*\\)\".*/\\1\\t\\2/'",
  " | sed 's/^\\(ENSG[0-9]*\\)\\.[0-9]*/\\1/'")
gencode_map <- fread(cmd = gencode_cmd, header = FALSE, col.names = c("ensembl_clean", "human_symbol"))
gencode_map <- gencode_map[!duplicated(ensembl_clean)]

# Add only genes NOT already in ortho-based symbol_map
n_ortho <- nrow(symbol_map)
missing_genes <- setdiff(gencode_map$ensembl_clean, symbol_map$ensembl_clean)
gencode_supplement <- gencode_map[ensembl_clean %in% missing_genes]
symbol_map <- rbind(symbol_map, gencode_supplement)
cat(sprintf("  Symbol map: %d genes (%d from orthologs + %d from GENCODE)\n",
            nrow(symbol_map), n_ortho, nrow(gencode_supplement)))

# Add biotype annotation from GENCODE
gencode_biotype_cmd <- paste0("zcat ", gencode_gtf,
  " | awk -F'\\t' '$3==\"gene\"'",
  " | sed 's/.*gene_id \"\\([^\"]*\\)\".*gene_type \"\\([^\"]*\\)\".*/\\1\\t\\2/'",
  " | sed 's/^\\(ENSG[0-9]*\\)\\.[0-9]*/\\1/'")
biotype_map <- fread(cmd = gencode_biotype_cmd, header = FALSE,
                     col.names = c("ensembl_clean", "gene_biotype"))
biotype_map <- biotype_map[!duplicated(ensembl_clean)]

consensus <- merge(consensus, symbol_map, by = "ensembl_clean", all.x = TRUE)

# Also pull in the human tier from consensus columns
# dream_sig + dream_dir already exist; consensus tier comes from the full consensus file
# If a separate tier column is needed, derive from the per_study overlap later
# For now, keep bulk_logFC, bulk_padj, t-statistic, dream_sig, dream_dir

cat("  Genes with HGNC symbol:", sum(!is.na(consensus$human_symbol)), "/", nrow(consensus), "\n")

# ================================================================
# Layer 1b: Human consensus tiers (if separate tier file exists)
# ================================================================
# The consensus_degs.csv only has dream results. The actual tier classification
# (Tier1_HighConfidence, Tier2_Moderate, Tier3_Exploratory) was computed in Script 07.
# Check if it's present as a column or needs separate loading.
tier_file <- file.path(RDIR, "consensus_degs.csv")
# dream_sig column already indicates significance; we need the tier mapping
# from the existing multi_evidence_scored_genes.csv (which has human_tier)
old_me_file <- file.path(OUTDIR, "multi_evidence_scored_genes.csv")
if (file.exists(old_me_file)) {
  old_me <- fread(old_me_file, select = c("ensembl_clean", "human_tier"))
  old_me <- old_me[!duplicated(ensembl_clean)]
} else {
  cat("  NOTE: multi_evidence_scored_genes.csv not found; human_consensus_tier will be NA\n")
  old_me <- data.table(ensembl_clean = character(0), human_tier = character(0))
}

# ================================================================
# Layer 2: Mouse DE evidence
# ================================================================
cat("Loading Layer 2: Mouse consensus DEGs...\n")
mouse_consensus <- fread(file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/mouse_consensus_degs.csv"))
mouse_consensus[, ensembl_clean := sub("\\..*", "", gene)]

# Map mouse Ensembl to human symbol via ortholog comparison
mouse_map <- unique(ortho[!is.na(human_symbol) & human_symbol != "",
                           .(mouse_gene_id, human_symbol)])
mouse_map <- mouse_map[!duplicated(mouse_gene_id)]

mouse_consensus <- merge(mouse_consensus, mouse_map,
                          by.x = "ensembl_clean", by.y = "mouse_gene_id", all.x = TRUE)

# Mouse meta-analysis (per-diet) for padj
mouse_meta <- fread(file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/meta_analysis/meta_per_diet.csv"))
mouse_meta[, ensembl_clean := sub("\\..*", "", gene)]
mouse_meta <- merge(mouse_meta, mouse_map,
                     by.x = "ensembl_clean", by.y = "mouse_gene_id", all.x = TRUE)
mouse_meta_layer <- mouse_meta[!is.na(human_symbol), .(
  human_symbol,
  mouse_meta_logFC = meta_logFC,
  mouse_meta_padj  = meta_padj,
  n_diets_sig      = n_diets
)]
mouse_meta_layer <- mouse_meta_layer[!duplicated(human_symbol)]

# Mouse consensus tier
mouse_tier_dt <- mouse_consensus[!is.na(human_symbol), .(
  human_symbol,
  mouse_consensus_tier = tier,
  mouse_ensembl = ensembl_clean
)]
mouse_tier_dt <- mouse_tier_dt[!duplicated(human_symbol)]

# Merge mouse meta + tier
mouse_layer <- merge(mouse_meta_layer, mouse_tier_dt, by = "human_symbol", all = TRUE)

cat("  Mouse genes with human symbol:", nrow(mouse_layer), "\n")

# ================================================================
# Layer 3: Cross-species concordance
# ================================================================
cat("Loading Layer 3: Cross-species concordance atlas...\n")
concordance <- fread(file.path(BASE, "Analysis/Cross_Species_Concordance/results/concordance_atlas_unified.csv"))
conc_layer <- concordance[, .(
  human_symbol,
  primary_category,
  translatability_score,
  n_concordant_diets = n_concordant,
  best_mouse_model   = best_signature,
  mouse_gene_id
)]
conc_layer[, is_conserved := primary_category == "Conserved"]
conc_layer <- conc_layer[!duplicated(human_symbol)]
cat("  Genes in concordance atlas:", nrow(conc_layer), "\n")
cat("  Conserved:", sum(conc_layer$is_conserved, na.rm = TRUE), "\n")

# ================================================================
# Layer 4: Causal inference (TWAS + COLOC + INTACT)
# ================================================================
# 2026-04-22: Mendelian Randomization (MR) permanently ditched from the paper.
# TWAS + COLOC + INTACT is the sole causal framework. The MR loader block
# (Ghodsian mr_results_all.csv → mr_beta/mr_pval/mr_sig) was removed here.
# Archived scripts: archive/mr_ditched_2026-04-22/.
cat("Loading Layer 4: Causal inference...\n")

# Ensembl → HGNC annotation cache
ann_cache <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")
if (file.exists(ann_cache)) {
  ann <- fread(ann_cache)
} else {
  ann <- data.table(gene_base = character(0), symbol = character(0))
}

# TWAS results (multi-GWAS combined)
twas_combined_file <- file.path(BASE, "RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv")
if (file.exists(twas_combined_file)) {
  twas_combined <- fread(twas_combined_file)
  twas_combined[, gene_clean := sub("\\..*", "", gene)]
  # Map to HGNC
  if ("gene_name" %in% names(twas_combined)) {
    sym_lookup <- unique(twas_combined[, .(gene_clean = sub("\\..*", "", gene), gene_name)])
    twas_layer <- twas_combined[!is.na(pvalue), .(
      twas_z    = zscore[which.min(pvalue)],
      twas_pval = min(pvalue, na.rm = TRUE)
    ), by = .(gene_clean)]
    twas_layer <- merge(twas_layer, sym_lookup, by = "gene_clean", all.x = TRUE)
    twas_layer[, human_symbol := gene_name]
  } else {
    twas_layer <- twas_combined[!is.na(pvalue), .(
      twas_z    = zscore[which.min(pvalue)],
      twas_pval = min(pvalue, na.rm = TRUE)
    ), by = .(gene_clean)]
    twas_layer <- merge(twas_layer, ann[, .(gene_base, symbol)],
                        by.x = "gene_clean", by.y = "gene_base", all.x = TRUE)
    twas_layer[, human_symbol := symbol]
  }
  twas_layer <- twas_layer[!is.na(human_symbol), .(human_symbol, twas_z, twas_pval)]
  twas_layer <- twas_layer[!duplicated(human_symbol)]
  cat("  TWAS genes:", nrow(twas_layer), "\n")
} else {
  cat("  Multi-GWAS TWAS file not found\n")
  twas_layer <- data.table(human_symbol = character(0),
                            twas_z = numeric(0), twas_pval = numeric(0))
}

# COLOC results
coloc_file <- file.path(BASE, "RNA-seq/results/causal_inference/ghodsian/coloc_results.csv")
if (file.exists(coloc_file)) {
  coloc_res <- fread(coloc_file)
  coloc_res[, gene_clean := sub("\\..*", "", gene)]
  coloc_res <- merge(coloc_res, ann[, .(gene_base, symbol)],
                     by.x = "gene_clean", by.y = "gene_base", all.x = TRUE)
  coloc_layer <- coloc_res[!is.na(symbol), .(
    coloc_pp4 = max(PP.H4, na.rm = TRUE)
  ), by = .(human_symbol = symbol)]
  coloc_layer <- coloc_layer[!duplicated(human_symbol)]
  cat("  COLOC genes:", nrow(coloc_layer), "\n")
} else {
  cat("  COLOC results not found\n")
  coloc_layer <- data.table(human_symbol = character(0), coloc_pp4 = numeric(0))
}

# COLOC gene-level summary (ABF + real SuSiE when available)
# NOTE: coloc_susie_best_pp4 is the canonical atlas column name (13 downstream scripts
# depend on it). When real SuSiE results exist in gene_level_coloc.csv, we source from
# coloc_best_susie_pp4; otherwise we fall back to ABF coloc_best_pp4.
susie_gene_file <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
if (file.exists(susie_gene_file)) {
  susie_gene <- fread(susie_gene_file)
  has_real_susie <- "coloc_best_susie_pp4" %in% names(susie_gene)

  # ABF columns (always present)
  susie_layer <- susie_gene[gene != "", .(
    human_symbol            = gene,
    coloc_abf_best_pp4      = coloc_best_pp4,
    coloc_abf_best_gwas     = coloc_best_gwas,
    coloc_abf_n_gwas_h4_05  = coloc_n_gwas_h4_05,
    coloc_abf_n_gwas_h4_08  = coloc_n_gwas_h4_08,
    coloc_n_gwas_tested     = coloc_n_gwas_tested,
    coloc_n_groups_h4_05    = coloc_n_groups_h4_05,
    coloc_n_groups_h4_08    = coloc_n_groups_h4_08,
    coloc_is_mhc            = is_mhc,
    coloc_ld_cluster_flag   = ld_cluster_flag
  )]
  susie_layer <- susie_layer[!duplicated(human_symbol)]

  # Real SuSiE columns (present after Script 07 re-run with SuSiE-COLOC results)
  if (has_real_susie) {
    susie_extra <- susie_gene[gene != "", .(
      human_symbol                = gene,
      coloc_real_susie_pp4        = coloc_best_susie_pp4,
      coloc_real_susie_gwas       = coloc_best_susie_gwas,
      coloc_susie_n_gwas_h4_05   = coloc_n_gwas_susie_h4_05,
      coloc_susie_n_gwas_h4_08   = coloc_n_gwas_susie_h4_08,
      coloc_susie_n_gwas_h4_09   = coloc_n_gwas_susie_h4_09,
      coloc_susie_n_pairs_total  = coloc_susie_n_pairs_total,
      coloc_susie_success_rate   = coloc_susie_success_rate,
      coloc_n_groups_susie_h4_05 = coloc_n_groups_susie_h4_05,
      coloc_n_groups_susie_h4_08 = coloc_n_groups_susie_h4_08
    )]
    susie_extra <- susie_extra[!duplicated(human_symbol)]
    susie_layer <- merge(susie_layer, susie_extra, by = "human_symbol", all.x = TRUE)
    cat("  Real SuSiE columns detected and loaded\n")
  }

  # coloc_susie_best_pp4: use real SuSiE when available, ABF fallback otherwise
  # Preserves column name for 13+ downstream scripts
  if (has_real_susie) {
    susie_layer[, coloc_susie_best_pp4 := coloc_real_susie_pp4]
    susie_layer[, coloc_susie_best_gwas := coloc_real_susie_gwas]
    # T0.4 (2026-04-22): susie_backed flag = TRUE iff gene has a non-NA SuSiE PP4.
    # Lets downstream figures/filters distinguish "genuinely SuSiE-colocalizing" genes
    # from those only supported by legacy ABF (18.3% of eGenes per CLAUDE.md).
    susie_layer[, susie_backed := !is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 > 0]
    # Fall back to ABF for genes where SuSiE failed/unavailable
    susie_layer[is.na(coloc_susie_best_pp4), coloc_susie_best_pp4 := coloc_abf_best_pp4]
    susie_layer[is.na(coloc_susie_best_gwas), coloc_susie_best_gwas := coloc_abf_best_gwas]
    susie_layer[, c("coloc_real_susie_pp4", "coloc_real_susie_gwas") := NULL]
  } else {
    susie_layer[, coloc_susie_best_pp4  := coloc_abf_best_pp4]
    susie_layer[, coloc_susie_best_gwas := coloc_abf_best_gwas]
    # T0.4: no SuSiE available -> susie_backed always FALSE
    susie_layer[, susie_backed := FALSE]
    cat("  NOTE: No real SuSiE columns — coloc_susie_best_pp4 sourced from ABF\n")
  }
  # Legacy column names for downstream compat
  susie_layer[, coloc_susie_n_signals := coloc_abf_n_gwas_h4_05]
  susie_layer[, coloc_n_gwas_h4_08 := coloc_abf_n_gwas_h4_08]

  cat("  COLOC gene-level:", nrow(susie_layer), "genes\n")
  cat("  ABF PP4 > 0.5:", sum(susie_layer$coloc_abf_best_pp4 > 0.5, na.rm = TRUE), "\n")
  cat("  ABF PP4 > 0.9:", sum(susie_layer$coloc_abf_best_pp4 > 0.9, na.rm = TRUE), "\n")
  if (has_real_susie) {
    cat("  SuSiE PP4 > 0.5:", sum(susie_layer$coloc_susie_best_pp4 > 0.5, na.rm = TRUE), "\n")
    cat("  SuSiE PP4 > 0.9:", sum(susie_layer$coloc_susie_best_pp4 > 0.9, na.rm = TRUE), "\n")
  }
} else {
  cat("  COLOC gene-level file not found\n")
  susie_layer <- data.table(human_symbol = character(0),
                             coloc_abf_best_pp4 = numeric(0),
                             coloc_abf_best_gwas = character(0),
                             coloc_susie_best_pp4 = numeric(0),
                             coloc_susie_n_signals = integer(0),
                             coloc_n_gwas_h4_08 = integer(0),
                             coloc_susie_best_gwas = character(0),
                             coloc_is_mhc = logical(0),
                             coloc_ld_cluster_flag = character(0),
                             susie_backed = logical(0))
}

# ================================================================
# Layer 4aa: SuSiEX cross-ancestry joint fine-mapping (EUR+EAS)
# ================================================================
cat("Loading Layer 4aa: SuSiEX cross-ancestry fine-mapping...\n")
susiex_gene_file <- file.path(BASE, "GWAS/finemapping/results/susiex/susiex_gene_summary.csv")
if (file.exists(susiex_gene_file)) {
  susiex_gene <- fread(susiex_gene_file)
  susiex_layer <- susiex_gene[, .(
    human_symbol           = GeneSymbol,
    susiex_max_pip         = susiex_max_pip,
    susiex_n_loci          = susiex_n_loci,
    susiex_n_cs            = susiex_n_cs,
    susiex_cs_size_joint   = susiex_cs_size_joint,
    susiex_trait_pairs     = trait_pairs
  )]
  susiex_layer <- susiex_layer[!duplicated(human_symbol)]
  cat("  SuSiEX genes:", nrow(susiex_layer), "\n")
  cat("    PIP > 0.5:", sum(susiex_layer$susiex_max_pip > 0.5, na.rm = TRUE), "\n")
  cat("    PIP > 0.8:", sum(susiex_layer$susiex_max_pip > 0.8, na.rm = TRUE), "\n")
} else {
  cat("  SuSiEX gene summary not found\n")
  susiex_layer <- data.table(human_symbol = character(0),
                              susiex_max_pip = numeric(0),
                              susiex_n_loci = integer(0),
                              susiex_n_cs = integer(0),
                              susiex_cs_size_joint = integer(0),
                              susiex_trait_pairs = character(0))
}

# ================================================================
# Layer 4ab: MESuSiE cross-ancestry fine-mapping (shared + ancestry-specific)
# ================================================================
cat("Loading Layer 4ab: MESuSiE cross-ancestry fine-mapping...\n")
mesusie_gene_file <- file.path(BASE, "GWAS/finemapping/results/mesusie/mesusie_gene_summary.csv")
if (file.exists(mesusie_gene_file)) {
  mesusie_gene <- fread(mesusie_gene_file)
  mesusie_layer <- mesusie_gene[, .(
    human_symbol             = GeneSymbol,
    mesusie_max_pip          = mesusie_max_pip,
    mesusie_max_pip_shared   = mesusie_max_pip_shared,
    mesusie_in_shared_cs     = mesusie_in_shared_cs,
    mesusie_in_eur_cs        = mesusie_in_eur_cs,
    mesusie_in_eas_cs        = mesusie_in_eas_cs,
    mesusie_n_loci           = mesusie_n_loci,
    mesusie_trait_pairs      = trait_pairs
  )]
  mesusie_layer <- mesusie_layer[!duplicated(human_symbol)]
  cat("  MESuSiE genes:", nrow(mesusie_layer), "\n")
  cat("    PIP > 0.5:", sum(mesusie_layer$mesusie_max_pip > 0.5, na.rm = TRUE), "\n")
  cat("    In shared CS:", sum(mesusie_layer$mesusie_in_shared_cs, na.rm = TRUE), "\n")
} else {
  cat("  MESuSiE gene summary not found\n")
  mesusie_layer <- data.table(human_symbol = character(0),
                               mesusie_max_pip = numeric(0),
                               mesusie_max_pip_shared = numeric(0),
                               mesusie_in_shared_cs = logical(0),
                               mesusie_in_eur_cs = logical(0),
                               mesusie_in_eas_cs = logical(0),
                               mesusie_n_loci = integer(0),
                               mesusie_trait_pairs = character(0))
}

# ================================================================
# Layer 4b: sc-eQTL COLOC (multi-cell-type, from Script 36)
# ================================================================
cat("Loading Layer 4b: sc-eQTL COLOC...\n")
sceqtl_coloc_file <- file.path(BASE, "RNA-seq/results/causal_inference/sceqtl/sceqtl_multi_celltype_summary.csv")
if (file.exists(sceqtl_coloc_file)) {
  sceqtl_coloc <- fread(sceqtl_coloc_file)
  sceqtl_layer <- sceqtl_coloc[, .(
    human_symbol         = gene,
    sceqtl_coloc_pp4_hep = fifelse(grepl("hepatocyte", cell_types_tested),
                                    best_PP.H4_all, NA_real_),
    sceqtl_coloc_best_pp4 = best_PP.H4_all,
    sceqtl_coloc_cell_type = best_cell_type,
    sceqtl_n_cell_types  = n_cell_types
  )]
  sceqtl_layer <- sceqtl_layer[!duplicated(human_symbol)]
  cat("  sc-eQTL COLOC genes:", nrow(sceqtl_layer), "\n")
  cat("  sc-eQTL PP.H4 > 0.5:", sum(sceqtl_layer$sceqtl_coloc_best_pp4 > 0.5, na.rm = TRUE), "\n")
} else {
  cat("  sc-eQTL COLOC results not found (run Script 36 first)\n")
  sceqtl_layer <- data.table(human_symbol = character(0),
                              sceqtl_coloc_pp4_hep = numeric(0),
                              sceqtl_coloc_best_pp4 = numeric(0),
                              sceqtl_coloc_cell_type = character(0),
                              sceqtl_n_cell_types = integer(0))
}

# ================================================================
# Layer 4c: Zenodo precomputed COLOC (from Script 37)
# ================================================================
cat("Loading Layer 4c: Zenodo precomputed COLOC...\n")
zenodo_coloc_file <- file.path(BASE, "RNA-seq/results/causal_inference/sceqtl/zenodo_coloc_integrated.csv")
if (file.exists(zenodo_coloc_file)) {
  zenodo_coloc <- fread(zenodo_coloc_file)
  # Get NAFLD-colocalizing genes (best across cell types)
  zenodo_nafld <- zenodo_coloc[nafld_coloc == TRUE]
  if (nrow(zenodo_nafld) > 0) {
    zenodo_layer <- zenodo_nafld[, .(
      zenodo_nafld_coloc = TRUE,
      zenodo_coloc_cell_types = paste(unique(cell_type), collapse = ";"),
      zenodo_n_traits_coloc = max(n_traits_coloc, na.rm = TRUE)
    ), by = .(human_symbol = gene)]
    zenodo_layer <- zenodo_layer[!duplicated(human_symbol)]
  } else {
    zenodo_layer <- data.table(human_symbol = character(0),
                                zenodo_nafld_coloc = logical(0),
                                zenodo_coloc_cell_types = character(0),
                                zenodo_n_traits_coloc = integer(0))
  }
  cat("  Zenodo NAFLD-colocalizing genes:", nrow(zenodo_layer), "\n")
} else {
  cat("  Zenodo COLOC results not found (run Script 37 first)\n")
  zenodo_layer <- data.table(human_symbol = character(0),
                              zenodo_nafld_coloc = logical(0),
                              zenodo_coloc_cell_types = character(0),
                              zenodo_n_traits_coloc = integer(0))
}

# ================================================================
# Layer 4d: Broadaway COLOC (from Script 35)
# ================================================================
cat("Loading Layer 4d: Broadaway COLOC...\n")
broadaway_coloc_file <- file.path(BASE, "RNA-seq/results/causal_inference/broadaway/coloc_results.csv")
if (file.exists(broadaway_coloc_file)) {
  broadaway_coloc <- fread(broadaway_coloc_file)
  if (nrow(broadaway_coloc) > 0 && "PP.H4" %in% names(broadaway_coloc)) {
    broadaway_layer <- broadaway_coloc[, .(
      human_symbol     = gene,
      broadaway_coloc_pp4 = PP.H4
    )]
    broadaway_layer <- broadaway_layer[!duplicated(human_symbol)]
    cat("  Broadaway COLOC genes:", nrow(broadaway_layer), "\n")
    cat("  Broadaway PP.H4 > 0.5:", sum(broadaway_layer$broadaway_coloc_pp4 > 0.5, na.rm = TRUE), "\n")
  } else {
    broadaway_layer <- data.table(human_symbol = character(0), broadaway_coloc_pp4 = numeric(0))
  }
} else {
  cat("  Broadaway COLOC results not found (run Script 35 first)\n")
  broadaway_layer <- data.table(human_symbol = character(0), broadaway_coloc_pp4 = numeric(0))
}

# ================================================================
# Layer 4e: UKBB ALT COLOC (from Script 36b)
# ================================================================
cat("Loading Layer 4e: UKBB ALT COLOC...\n")
ukbb_coloc_file <- file.path(BASE, "RNA-seq/results/causal_inference/sceqtl/sceqtl_coloc_ukbb_alt_summary.csv")
if (file.exists(ukbb_coloc_file)) {
  ukbb_coloc <- fread(ukbb_coloc_file)
  if (nrow(ukbb_coloc) > 0) {
    ukbb_layer <- ukbb_coloc[, .(
      human_symbol = gene,
      ukbb_alt_coloc_pp4 = best_PP.H4_all,
      ukbb_alt_coloc_cell_type = best_cell_type
    )]
    ukbb_layer <- ukbb_layer[!duplicated(human_symbol)]
    cat("  UKBB ALT COLOC genes:", nrow(ukbb_layer), "\n")
    cat("  UKBB ALT PP.H4 > 0.5:", sum(ukbb_layer$ukbb_alt_coloc_pp4 > 0.5, na.rm = TRUE), "\n")
  } else {
    ukbb_layer <- data.table(human_symbol = character(0),
                              ukbb_alt_coloc_pp4 = numeric(0),
                              ukbb_alt_coloc_cell_type = character(0))
  }
} else {
  cat("  UKBB ALT COLOC results not found (run Script 36b first)\n")
  ukbb_layer <- data.table(human_symbol = character(0),
                            ukbb_alt_coloc_pp4 = numeric(0),
                            ukbb_alt_coloc_cell_type = character(0))
}

# ================================================================
# Layer 4h: Broadaway x UKBB AST COLOC (from Script 35g)
# ================================================================
cat("Loading Layer 4h: Broadaway x UKBB AST COLOC...\n")
ast_coloc_file <- file.path(BASE, "RNA-seq/results/causal_inference/broadaway_ukbb_ast/coloc_results.csv")
if (file.exists(ast_coloc_file)) {
  ast_coloc <- fread(ast_coloc_file)
  if (nrow(ast_coloc) > 0 && "PP.H4" %in% names(ast_coloc)) {
    ast_layer <- ast_coloc[, .(
      human_symbol     = gene,
      ast_coloc_pp4    = PP.H4
    )]
    ast_layer <- ast_layer[!duplicated(human_symbol)]
    cat("  AST COLOC genes:", nrow(ast_layer), "\n")
    cat("  AST PP.H4 > 0.5:", sum(ast_layer$ast_coloc_pp4 > 0.5, na.rm = TRUE), "\n")
  } else {
    ast_layer <- data.table(human_symbol = character(0), ast_coloc_pp4 = numeric(0))
  }
} else {
  cat("  AST COLOC results not found (run Script 35g with GWAS_NAME=UKBB_AST first)\n")
  ast_layer <- data.table(human_symbol = character(0), ast_coloc_pp4 = numeric(0))
}

# ================================================================
# Layer 4i: Broadaway x UKBB GGT COLOC (from Script 35g)
# ================================================================
cat("Loading Layer 4i: Broadaway x UKBB GGT COLOC...\n")
ggt_coloc_file <- file.path(BASE, "RNA-seq/results/causal_inference/broadaway_ukbb_ggt/coloc_results.csv")
if (file.exists(ggt_coloc_file)) {
  ggt_coloc <- fread(ggt_coloc_file)
  if (nrow(ggt_coloc) > 0 && "PP.H4" %in% names(ggt_coloc)) {
    ggt_layer <- ggt_coloc[, .(
      human_symbol     = gene,
      ggt_coloc_pp4    = PP.H4
    )]
    ggt_layer <- ggt_layer[!duplicated(human_symbol)]
    cat("  GGT COLOC genes:", nrow(ggt_layer), "\n")
    cat("  GGT PP.H4 > 0.5:", sum(ggt_layer$ggt_coloc_pp4 > 0.5, na.rm = TRUE), "\n")
  } else {
    ggt_layer <- data.table(human_symbol = character(0), ggt_coloc_pp4 = numeric(0))
  }
} else {
  cat("  GGT COLOC results not found (run Script 35g with GWAS_NAME=UKBB_GGT first)\n")
  ggt_layer <- data.table(human_symbol = character(0), ggt_coloc_pp4 = numeric(0))
}

# ================================================================
# Layer 4j: Broadaway x PDFF COLOC (from Script 35h)
# ================================================================
cat("Loading Layer 4j: Broadaway x PDFF COLOC...\n")
pdff_coloc_file <- file.path(BASE, "RNA-seq/results/causal_inference/broadaway_pdff/coloc_results.csv")
if (file.exists(pdff_coloc_file)) {
  pdff_coloc <- fread(pdff_coloc_file)
  if (nrow(pdff_coloc) > 0 && "PP.H4" %in% names(pdff_coloc)) {
    pdff_layer <- pdff_coloc[, .(
      human_symbol     = gene,
      pdff_coloc_pp4   = PP.H4
    )]
    pdff_layer <- pdff_layer[!duplicated(human_symbol)]
    cat("  PDFF COLOC genes:", nrow(pdff_layer), "\n")
    cat("  PDFF PP.H4 > 0.5:", sum(pdff_layer$pdff_coloc_pp4 > 0.5, na.rm = TRUE), "\n")
  } else {
    pdff_layer <- data.table(human_symbol = character(0), pdff_coloc_pp4 = numeric(0))
  }
} else {
  cat("  PDFF COLOC results not found (run Script 35h first)\n")
  pdff_layer <- data.table(human_symbol = character(0), pdff_coloc_pp4 = numeric(0))
}

# ================================================================
# Layer 4k-m: FinnGen COLOC — NAFLD/NASH/HCC
# ARCHIVED 2026-04-08: FinnGen duplicate of Whitfield 2023 — columns retained as NA for backwards compatibility
# ================================================================
cat("Layer 4k-m: FinnGen COLOC — ARCHIVED (duplicate of Whitfield 2023). Columns retained as NA.\n")
fg_nafld_layer <- data.table(human_symbol = character(0), finngen_nafld_coloc_pp4 = numeric(0))
fg_nash_layer  <- data.table(human_symbol = character(0), finngen_nash_coloc_pp4 = numeric(0))
fg_hcc_layer   <- data.table(human_symbol = character(0), finngen_hcc_coloc_pp4 = numeric(0))

# ================================================================
# Layer 4n: BBJ COLOC — ALT (from Script 47, cross-ancestry)
# ================================================================
cat("Loading Layer 4n: BBJ ALT COLOC (cross-ancestry)...\n")
bbj_alt_file <- file.path(BASE, "RNA-seq/results/causal_inference/bbj_alt/coloc_results.csv")
if (file.exists(bbj_alt_file)) {
  bbj_alt <- fread(bbj_alt_file)
  if (nrow(bbj_alt) > 0 && "PP.H4" %in% names(bbj_alt)) {
    bbj_alt_layer <- bbj_alt[, .(
      human_symbol       = gene,
      bbj_alt_coloc_pp4  = PP.H4
    )]
    bbj_alt_layer <- bbj_alt_layer[!duplicated(human_symbol)]
    cat("  BBJ ALT COLOC genes:", nrow(bbj_alt_layer), "\n")
    cat("  BBJ ALT PP.H4 > 0.5:", sum(bbj_alt_layer$bbj_alt_coloc_pp4 > 0.5, na.rm = TRUE), "\n")
  } else {
    bbj_alt_layer <- data.table(human_symbol = character(0), bbj_alt_coloc_pp4 = numeric(0))
  }
} else {
  cat("  BBJ ALT COLOC not found (run Script 47 with GWAS_NAME=BBJ_ALT first)\n")
  bbj_alt_layer <- data.table(human_symbol = character(0), bbj_alt_coloc_pp4 = numeric(0))
}

# ================================================================
# Layer 4o: BBJ COLOC — AST (from Script 47, cross-ancestry)
# ================================================================
cat("Loading Layer 4o: BBJ AST COLOC (cross-ancestry)...\n")
bbj_ast_file <- file.path(BASE, "RNA-seq/results/causal_inference/bbj_ast/coloc_results.csv")
if (file.exists(bbj_ast_file)) {
  bbj_ast <- fread(bbj_ast_file)
  if (nrow(bbj_ast) > 0 && "PP.H4" %in% names(bbj_ast)) {
    bbj_ast_layer <- bbj_ast[, .(
      human_symbol       = gene,
      bbj_ast_coloc_pp4  = PP.H4
    )]
    bbj_ast_layer <- bbj_ast_layer[!duplicated(human_symbol)]
    cat("  BBJ AST COLOC genes:", nrow(bbj_ast_layer), "\n")
    cat("  BBJ AST PP.H4 > 0.5:", sum(bbj_ast_layer$bbj_ast_coloc_pp4 > 0.5, na.rm = TRUE), "\n")
  } else {
    bbj_ast_layer <- data.table(human_symbol = character(0), bbj_ast_coloc_pp4 = numeric(0))
  }
} else {
  cat("  BBJ AST COLOC not found (run Script 47 with GWAS_NAME=BBJ_AST first)\n")
  bbj_ast_layer <- data.table(human_symbol = character(0), bbj_ast_coloc_pp4 = numeric(0))
}

# ================================================================
# Layer 4p: BBJ COLOC — GGT (from Script 47, cross-ancestry)
# ================================================================
cat("Loading Layer 4p: BBJ GGT COLOC (cross-ancestry)...\n")
bbj_ggt_file <- file.path(BASE, "RNA-seq/results/causal_inference/bbj_ggt/coloc_results.csv")
if (file.exists(bbj_ggt_file)) {
  bbj_ggt <- fread(bbj_ggt_file)
  if (nrow(bbj_ggt) > 0 && "PP.H4" %in% names(bbj_ggt)) {
    bbj_ggt_layer <- bbj_ggt[, .(
      human_symbol       = gene,
      bbj_ggt_coloc_pp4  = PP.H4
    )]
    bbj_ggt_layer <- bbj_ggt_layer[!duplicated(human_symbol)]
    cat("  BBJ GGT COLOC genes:", nrow(bbj_ggt_layer), "\n")
    cat("  BBJ GGT PP.H4 > 0.5:", sum(bbj_ggt_layer$bbj_ggt_coloc_pp4 > 0.5, na.rm = TRUE), "\n")
  } else {
    bbj_ggt_layer <- data.table(human_symbol = character(0), bbj_ggt_coloc_pp4 = numeric(0))
  }
} else {
  cat("  BBJ GGT COLOC not found (run Script 47 with GWAS_NAME=BBJ_GGT first)\n")
  bbj_ggt_layer <- data.table(human_symbol = character(0), bbj_ggt_coloc_pp4 = numeric(0))
}

# ================================================================
# Layer 4q: Ghouse Cirrhosis COLOC (from Script 49, progression)
# ================================================================
# Layer 4q DROPPED 2026-06-06: Ghouse Cirrhosis is an etiology-mixed endpoint, not
# MASLD-specific (Broadaway excluded cirrhosis/HCC). Column retained as empty/NA so the
# 23-GWAS MASLD portfolio gives no genetic credit for cirrhosis/HCC colocalization.
cat("Layer 4q: Ghouse Cirrhosis COLOC — DROPPED 2026-06-06 (cirrhosis/HCC not MASLD-specific). Column NA.\n")
ghouse_layer <- data.table(human_symbol = character(0), ghouse_cirrhosis_coloc_pp4 = numeric(0))

# ================================================================
# Layer 4q2: Ghouse HCC COLOC (from Script 49, progression)
# ================================================================
# Layer 4q2 DROPPED 2026-06-06: Ghouse HCC is not MASLD-specific (Broadaway excluded
# cirrhosis/HCC). Column retained as empty/NA — no genetic credit for HCC colocalization.
cat("Layer 4q2: Ghouse HCC COLOC — DROPPED 2026-06-06 (cirrhosis/HCC not MASLD-specific). Column NA.\n")
ghouse_hcc_layer <- data.table(human_symbol = character(0), ghouse_hcc_coloc_pp4 = numeric(0))

# ================================================================
# Layer 4r: deCODE NAFL COLOC (from Script 49, Icelandic replication)
# ================================================================
cat("Loading Layer 4r: deCODE NAFL COLOC...\n")
decode_nafl_file <- file.path(BASE, "RNA-seq/results/causal_inference/decode_nafl/coloc_results.csv")
if (file.exists(decode_nafl_file)) {
  decode_nafl <- fread(decode_nafl_file)
  if (nrow(decode_nafl) > 0 && "PP.H4" %in% names(decode_nafl)) {
    decode_nafl_layer <- decode_nafl[, .(
      human_symbol           = gene,
      decode_nafl_coloc_pp4  = PP.H4
    )]
    decode_nafl_layer <- decode_nafl_layer[!duplicated(human_symbol)]
    cat("  deCODE NAFL COLOC genes:", nrow(decode_nafl_layer), "\n")
    cat("  deCODE NAFL PP.H4 > 0.5:", sum(decode_nafl_layer$decode_nafl_coloc_pp4 > 0.5, na.rm = TRUE), "\n")
  } else {
    decode_nafl_layer <- data.table(human_symbol = character(0), decode_nafl_coloc_pp4 = numeric(0))
  }
} else {
  cat("  deCODE NAFL COLOC not found (run Script 49 with GWAS_NAME=DECODE_NAFL first)\n")
  decode_nafl_layer <- data.table(human_symbol = character(0), decode_nafl_coloc_pp4 = numeric(0))
}

# ================================================================
# Layer 4s: deCODE Cirrhosis COLOC (from Script 49)
# ================================================================
# Layer 4s DROPPED 2026-06-06: deCODE (Sveinbjornsson) Cirrhosis is not MASLD-specific
# (Broadaway excluded cirrhosis/HCC). Column retained as empty/NA.
cat("Layer 4s: deCODE Cirrhosis COLOC — DROPPED 2026-06-06 (cirrhosis/HCC not MASLD-specific). Column NA.\n")
decode_cirrhosis_layer <- data.table(human_symbol = character(0), decode_cirrhosis_coloc_pp4 = numeric(0))

# ================================================================
# Layer 4t: deCODE HCC COLOC (from Script 49)
# ================================================================
# Layer 4t DROPPED 2026-06-06: deCODE (Sveinbjornsson) HCC is not MASLD-specific
# (Broadaway excluded cirrhosis/HCC). Column retained as empty/NA.
cat("Layer 4t: deCODE HCC COLOC — DROPPED 2026-06-06 (cirrhosis/HCC not MASLD-specific). Column NA.\n")
decode_hcc_layer <- data.table(human_symbol = character(0), decode_hcc_coloc_pp4 = numeric(0))

# ================================================================
# Layer 4u-z: Pan-UKBB AFR/CSA COLOC (from Script 50, multi-ancestry)
# ================================================================
panukbb_layers <- list()
for (pop in c("afr", "csa")) {
  for (trait in c("alt", "ast", "ggt")) {
    layer_name <- paste0("panukbb_", pop, "_", trait)
    col_name <- paste0(layer_name, "_coloc_pp4")
    pfile <- file.path(BASE, "RNA-seq/results/causal_inference", layer_name, "coloc_results.csv")
    cat("Loading Layer: Pan-UKBB", toupper(pop), toupper(trait), "COLOC...\n")
    if (file.exists(pfile)) {
      pdt <- fread(pfile)
      if (nrow(pdt) > 0 && "PP.H4" %in% names(pdt)) {
        player <- pdt[, .(human_symbol = gene, pp4 = PP.H4)]
        player <- player[!duplicated(human_symbol)]
        setnames(player, "pp4", col_name)
        cat("  ", layer_name, "genes:", nrow(player), " PP.H4>0.5:", sum(player[[col_name]] > 0.5, na.rm = TRUE), "\n")
        panukbb_layers[[layer_name]] <- player
      } else {
        panukbb_layers[[layer_name]] <- data.table(human_symbol = character(0))
        panukbb_layers[[layer_name]][[col_name]] <- numeric(0)
      }
    } else {
      cat("  ", layer_name, "not found\n")
      panukbb_layers[[layer_name]] <- data.table(human_symbol = character(0))
      panukbb_layers[[layer_name]][[col_name]] <- numeric(0)
    }
  }
}

# ================================================================
# Layer 4f: ieQTL disease interaction (from Script 38)
# ================================================================
cat("Loading Layer 4e: ieQTL disease interaction...\n")
ieqtl_file <- file.path(BASE, "RNA-seq/results/causal_inference/sceqtl/ieqtl_disease_genes.csv")
if (file.exists(ieqtl_file)) {
  ieqtl <- fread(ieqtl_file)
  if (nrow(ieqtl) > 0) {
    # Best ieQTL per gene (across cell types)
    ieqtl_layer <- ieqtl[, .(
      ieqtl_disease_interaction = TRUE,
      ieqtl_interaction_pval = min(interaction_pval, na.rm = TRUE),
      ieqtl_cell_type = cell_type[which.min(interaction_pval)]
    ), by = .(human_symbol = gene)]
    ieqtl_layer <- ieqtl_layer[!duplicated(human_symbol)]
    cat("  ieQTL disease genes:", nrow(ieqtl_layer), "\n")
  } else {
    ieqtl_layer <- data.table(human_symbol = character(0),
                               ieqtl_disease_interaction = logical(0),
                               ieqtl_interaction_pval = numeric(0),
                               ieqtl_cell_type = character(0))
  }
} else {
  cat("  ieQTL results not found (run Script 38 first)\n")
  ieqtl_layer <- data.table(human_symbol = character(0),
                             ieqtl_disease_interaction = logical(0),
                             ieqtl_interaction_pval = numeric(0),
                             ieqtl_cell_type = character(0))
}

# ================================================================
# Layer 4g: Pleiotropy profile (from Script 39)
# ================================================================
cat("Loading Layer 4g: Pleiotropy classification...\n")
pleio_file <- file.path(BASE, "RNA-seq/results/causal_inference/pleiotropy/pleiotropy_gene_summary.csv")
if (file.exists(pleio_file)) {
  pleio_dt <- fread(pleio_file)
  pleio_layer <- pleio_dt[, .(human_symbol, pleiotropy_class = best_pleiotropy_class,
                               pleiotropy_n_traits = max_traits_coloc,
                               pleiotropy_n_domains = domain_breadth,
                               pleiotropy_domains = all_domains,
                               is_nafld_specific)]
  pleio_layer <- pleio_layer[!duplicated(human_symbol)]
  cat("  Genes with pleiotropy data:", nrow(pleio_layer), "\n")
} else {
  cat("  Pleiotropy results not found (run Script 39 first)\n")
  pleio_layer <- data.table(human_symbol = character(0),
                             pleiotropy_class = character(0),
                             pleiotropy_n_traits = integer(0),
                             pleiotropy_n_domains = integer(0),
                             pleiotropy_domains = character(0),
                             is_nafld_specific = logical(0))
}

# Compute n_coloc_sources per gene
cat("Computing n_coloc_sources...\n")
coloc_sources <- data.table(human_symbol = character(0), source = character(0))
# GTEx COLOC
if (nrow(coloc_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    coloc_layer[coloc_pp4 > 0.5, .(human_symbol, source = "GTEx")])
}
# Broadaway COLOC
if (nrow(broadaway_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    broadaway_layer[broadaway_coloc_pp4 > 0.5, .(human_symbol, source = "Broadaway")],
    fill = TRUE)
}
# sc-eQTL COLOC
if (nrow(sceqtl_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    sceqtl_layer[sceqtl_coloc_best_pp4 > 0.5, .(human_symbol, source = "sc-eQTL")],
    fill = TRUE)
}
# UKBB ALT COLOC
if (nrow(ukbb_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    ukbb_layer[ukbb_alt_coloc_pp4 > 0.5, .(human_symbol, source = "UKBB_ALT")],
    fill = TRUE)
}
# AST COLOC
if (nrow(ast_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    ast_layer[ast_coloc_pp4 > 0.5, .(human_symbol, source = "AST")],
    fill = TRUE)
}
# GGT COLOC
if (nrow(ggt_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    ggt_layer[ggt_coloc_pp4 > 0.5, .(human_symbol, source = "GGT")],
    fill = TRUE)
}
# PDFF COLOC
if (nrow(pdff_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    pdff_layer[pdff_coloc_pp4 > 0.5, .(human_symbol, source = "PDFF")],
    fill = TRUE)
}
# Zenodo COLOC
if (nrow(zenodo_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    zenodo_layer[zenodo_nafld_coloc == TRUE, .(human_symbol, source = "Zenodo")],
    fill = TRUE)
}
# FinnGen COLOC (3 phenotypes) — ARCHIVED 2026-04-08: duplicate of Whitfield 2023
# (layers are empty; no rows will be added)
# BBJ COLOC (3 liver enzymes — cross-ancestry, East Asian)
if (nrow(bbj_alt_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    bbj_alt_layer[bbj_alt_coloc_pp4 > 0.5, .(human_symbol, source = "BBJ_ALT")],
    fill = TRUE)
}
if (nrow(bbj_ast_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    bbj_ast_layer[bbj_ast_coloc_pp4 > 0.5, .(human_symbol, source = "BBJ_AST")],
    fill = TRUE)
}
if (nrow(bbj_ggt_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    bbj_ggt_layer[bbj_ggt_coloc_pp4 > 0.5, .(human_symbol, source = "BBJ_GGT")],
    fill = TRUE)
}
# Ghouse Cirrhosis COLOC (progression)
if (nrow(ghouse_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    ghouse_layer[ghouse_cirrhosis_coloc_pp4 > 0.5, .(human_symbol, source = "Ghouse_Cirrhosis")],
    fill = TRUE)
}
# Ghouse HCC COLOC (progression)
if (nrow(ghouse_hcc_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    ghouse_hcc_layer[ghouse_hcc_coloc_pp4 > 0.5, .(human_symbol, source = "Ghouse_HCC")],
    fill = TRUE)
}
# deCODE COLOC (Icelandic European replication)
if (nrow(decode_nafl_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    decode_nafl_layer[decode_nafl_coloc_pp4 > 0.5, .(human_symbol, source = "deCODE_NAFL")],
    fill = TRUE)
}
if (nrow(decode_cirrhosis_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    decode_cirrhosis_layer[decode_cirrhosis_coloc_pp4 > 0.5, .(human_symbol, source = "deCODE_Cirrhosis")],
    fill = TRUE)
}
if (nrow(decode_hcc_layer) > 0) {
  coloc_sources <- rbind(coloc_sources,
    decode_hcc_layer[decode_hcc_coloc_pp4 > 0.5, .(human_symbol, source = "deCODE_HCC")],
    fill = TRUE)
}
# Pan-UKBB AFR/CSA COLOC (multi-ancestry)
for (lname in names(panukbb_layers)) {
  player <- panukbb_layers[[lname]]
  col_name <- paste0(lname, "_coloc_pp4")
  if (nrow(player) > 0 && col_name %in% names(player)) {
    hits <- player[player[[col_name]] > 0.5, .(human_symbol, source = toupper(lname))]
    if (nrow(hits) > 0) {
      coloc_sources <- rbind(coloc_sources, hits, fill = TRUE)
    }
  }
}
if (nrow(coloc_sources) > 0) {
  n_coloc_dt <- coloc_sources[, .(n_coloc_sources = uniqueN(source)), by = human_symbol]
  cat("  Genes with any COLOC signal:", nrow(n_coloc_dt), "\n")
  cat("  Genes with 2+ COLOC sources:", sum(n_coloc_dt$n_coloc_sources >= 2), "\n")
} else {
  n_coloc_dt <- data.table(human_symbol = character(0), n_coloc_sources = integer(0))
}

# ================================================================
# Layer 5: Sex stratification
# ================================================================
cat("Loading Layer 5: Sex-stratified analysis...\n")
# v3 fallback (added 2026-05-13): prefer mashr-based v3 classification when present.
# v3 file (sex_v3/sex_deg_classification_v3.csv) is the canonical output of the
# Module 06 (06_classify_v3.R) pipeline. v2 (sex_deg_classification.csv) is
# retained as the legacy fallback; this preserves backward compatibility while
# the v3 pipeline lands.
sex_v3_path <- file.path(RDIR, "sex_v3/sex_deg_classification_v3.csv")
sex_v2_path <- file.path(RDIR, "sex_deg_classification.csv")
sex_class_file <- if (file.exists(sex_v3_path)) sex_v3_path else sex_v2_path
message(sprintf("Atlas reading sex classification from: %s", basename(sex_class_file)))
if (file.exists(sex_class_file)) {
  sex_class_dt <- fread(sex_class_file)
  sex_class_dt[, ensembl_clean := sub("\\..*", "", gene)]
  sex_class_dt <- merge(sex_class_dt, symbol_map, by = "ensembl_clean", all.x = TRUE)

  # Handle both old column names (meta_logFC_M/F from Script 06) and

  # new column names (logFC_M/F from dream-based Script 26)
  m_col <- if ("logFC_M" %in% names(sex_class_dt)) "logFC_M" else "meta_logFC_M"
  f_col <- if ("logFC_F" %in% names(sex_class_dt)) "logFC_F" else "meta_logFC_F"

  # Pull interaction_padj from classification file if available (Script 26 v2)
  has_int_padj <- "interaction_padj" %in% names(sex_class_dt)

  sex_layer <- sex_class_dt[!is.na(human_symbol), .(
    human_symbol,
    sex_class,
    bulk_logFC_M = get(m_col),
    bulk_logFC_F = get(f_col),
    sex_interaction_padj_from_cls = if (has_int_padj) interaction_padj else NA_real_
  )]
  sex_layer <- sex_layer[!duplicated(human_symbol)]
  cat("  Genes with sex classification:", nrow(sex_layer), "\n")
} else {
  cat("  Sex classification not found\n")
  sex_layer <- data.table(human_symbol = character(0), sex_class = character(0),
                           bulk_logFC_M = numeric(0), bulk_logFC_F = numeric(0))
}

# Sex interaction p-value
# R3 Issue 2 fix (2026-05-13): prefer v3 sex_v3/sex_interaction_dream_v3.csv
# when available, falling back to the v2 file. Same pattern as the L942
# sex classification fallback. v3 csv schema is a strict superset (extra
# adj.P.Val + se_int columns) so the `padj` lookup on L988 works for both.
sex_int_v3_path <- file.path(RDIR, "sex_v3/sex_interaction_dream_v3.csv")
sex_int_v2_path <- file.path(RDIR, "sex_interaction_dream.csv")
sex_int_file <- if (file.exists(sex_int_v3_path)) sex_int_v3_path else sex_int_v2_path
message(sprintf("Atlas reading sex interaction from: %s", basename(sex_int_file)))
if (file.exists(sex_int_file)) {
  sex_int <- fread(sex_int_file)
  sex_int[, ensembl_clean := sub("\\..*", "", gene)]
  sex_int <- merge(sex_int, symbol_map, by = "ensembl_clean", all.x = TRUE)
  sex_int_layer <- sex_int[!is.na(human_symbol), .(
    human_symbol,
    sex_interaction_padj = padj
  )]
  sex_int_layer <- sex_int_layer[!duplicated(human_symbol)]
  sex_layer <- merge(sex_layer, sex_int_layer, by = "human_symbol", all.x = TRUE)
  cat("  Genes with sex interaction padj:", sum(!is.na(sex_layer$sex_interaction_padj)), "\n")
}

# ================================================================
# Layer 6: Pathway membership (GSEA leading edges)
# ================================================================
cat("Loading Layer 6: GSEA pathway membership...\n")
gsea <- fread(file.path(RDIR, "gsea_results.csv"))
sig_pathways <- gsea[padj < 0.05]
cat("  Significant pathways:", nrow(sig_pathways), "\n")

if (nrow(sig_pathways) > 0) {
  le_list <- strsplit(sig_pathways$leadingEdge, ";")
  le_dt <- data.table(
    ensembl_id = unlist(le_list),
    pathway    = rep(sig_pathways$pathway, sapply(le_list, length))
  )
  le_dt[, ensembl_clean := sub("\\..*", "", ensembl_id)]

  # Count pathways per gene
  pathway_counts <- le_dt[, .(n_leading_edge_pathways = uniqueN(pathway)), by = ensembl_clean]
  pathway_counts <- merge(pathway_counts, symbol_map, by = "ensembl_clean", all.x = TRUE)

  # Top 3 pathways per gene
  gene_pathways <- le_dt[, .(pathways = paste(head(unique(pathway), 3), collapse = ";")),
                          by = ensembl_clean]
  gene_pathways <- merge(gene_pathways, symbol_map, by = "ensembl_clean", all.x = TRUE)

  pathway_layer <- merge(
    pathway_counts[!is.na(human_symbol), .(human_symbol, n_leading_edge_pathways)],
    gene_pathways[!is.na(human_symbol), .(human_symbol, top_pathways = pathways)],
    by = "human_symbol", all = TRUE
  )
  pathway_layer <- pathway_layer[!duplicated(human_symbol)]
} else {
  pathway_layer <- data.table(human_symbol = character(0),
                               n_leading_edge_pathways = integer(0),
                               top_pathways = character(0))
}
cat("  Genes in leading edges:", nrow(pathway_layer), "\n")

# ================================================================
# Layer 7: Essentiality (DepMap Cas9 Chronos — liver cell lines)
# ================================================================
cat("Loading Layer 7: Essentiality (DepMap Cas9, liver-specific)...\n")
depmap_file <- file.path(BASE, "Analysis/downstream_analysis/essentiality/CRISPRGeneEffect.csv")
model_file  <- file.path(BASE, "Analysis/downstream_analysis/essentiality/Model.csv")
if (file.exists(depmap_file) && file.exists(model_file)) {
  # Identify liver cell lines
  models <- fread(model_file)
  liver_ids <- models[OncotreeLineage == "Liver", ModelID]

  # Read DepMap Chronos scores (rows = cell lines, cols = genes)
  depmap <- fread(depmap_file)
  setnames(depmap, names(depmap)[1], "ModelID")
  liver_ids <- intersect(liver_ids, depmap$ModelID)
  cat("  Liver cell lines with CRISPR data:", length(liver_ids), "\n")

  # Subset to liver lines and gene columns only
  depmap_liver <- depmap[ModelID %in% liver_ids]
  gene_cols <- setdiff(names(depmap_liver), "ModelID")

  # Parse gene symbols from "GENE (EntrezID)" format
  gene_symbols <- sub(" \\(.*", "", gene_cols)

  # Compute per-gene mean Chronos score across liver lines
  scores <- depmap_liver[, lapply(.SD, function(x) mean(x, na.rm = TRUE)),
                         .SDcols = gene_cols]
  n_lines <- depmap_liver[, lapply(.SD, function(x) sum(!is.na(x))),
                          .SDcols = gene_cols]

  ess_layer <- data.table(
    human_symbol         = gene_symbols,
    essentiality_chronos = as.numeric(scores[1, ]),
    n_liver_lines        = as.integer(n_lines[1, ])
  )
  ess_layer[, is_essential := essentiality_chronos < -0.5]
  ess_layer <- ess_layer[!duplicated(human_symbol)]

  cat("  Genes with essentiality data:", nrow(ess_layer), "\n")
  cat("  Essential genes (Chronos < -0.5):", sum(ess_layer$is_essential, na.rm = TRUE), "\n")
  cat("  Median liver lines per gene:", median(ess_layer$n_liver_lines), "\n")

  # Clean up large objects
  rm(depmap, depmap_liver, models, scores, n_lines); gc(verbose = FALSE)
} else {
  cat("  DepMap files not found\n")
  ess_layer <- data.table(human_symbol = character(0),
                           essentiality_chronos = numeric(0),
                           n_liver_lines = integer(0),
                           is_essential = logical(0))
}

# ================================================================
# Annotations: Deconvolution attribution
# ================================================================
cat("Loading annotations: Deconvolution attribution...\n")
deconv_file <- file.path(BASE, "RNA-seq/results/causal_inference/deconv_attribution_scores.csv")
if (file.exists(deconv_file)) {
  deconv <- fread(deconv_file)
  deconv[, ensembl_clean := sub("\\..*", "", gene)]
  deconv <- merge(deconv, symbol_map, by = "ensembl_clean", all.x = TRUE)
  deconv_layer <- deconv[!is.na(human_symbol), .(
    human_symbol,
    attribution_class = category
  )]
  deconv_layer <- deconv_layer[!duplicated(human_symbol)]
  cat("  Genes with attribution:", nrow(deconv_layer), "\n")
} else {
  cat("  Deconvolution attribution not found\n")
  deconv_layer <- data.table(human_symbol = character(0), attribution_class = character(0))
}

# ================================================================
# Annotations: Drug target convergence
# ================================================================
cat("Loading annotations: Drug target convergence...\n")
drug_file <- file.path(BASE, "RNA-seq/results/drug_repurposing/convergent_drug_targets.csv")
if (!file.exists(drug_file)) {
  drug_file <- file.path(BASE, "RNA-seq/results/drug_repurposing/multi_layer_drug_targets.csv")
}
if (file.exists(drug_file)) {
  drugs <- fread(drug_file)
  drug_layer <- drugs[, .(
    human_symbol    = symbol,
    opentargets_drug = as.logical(ot_has_drug),
    lincs_reversal  = as.logical(lincs_reversal)
  )]
  drug_layer <- drug_layer[!duplicated(human_symbol)]
  cat("  Drug target genes:", nrow(drug_layer), "\n")
} else {
  cat("  Drug targets not found\n")
  drug_layer <- data.table(human_symbol = character(0),
                            opentargets_drug = logical(0),
                            lincs_reversal = logical(0))
}

# FIX 2026-06-10: dgidb_druggable was previously sourced from the small
# convergent_drug_targets.csv (~300 genes; genome-wide targets THRB/EGFR/HMGCR
# came back blank — a broken annotation). Now sourced GENOME-WIDE from the DGIdb
# interactions-derived druggable set (data/dgidb/dgidb_druggable_genes.csv,
# 5,012 genes; built from dgidb.org/data/latest/interactions.tsv). Merged as its
# own layer below; non-matches filled FALSE after the merge.
dgidb_file <- file.path(BASE, "data/dgidb/dgidb_druggable_genes.csv")
if (file.exists(dgidb_file)) {
  dgidb_gw <- fread(dgidb_file)
  dgidb_layer <- dgidb_gw[, .(human_symbol = toupper(gene_name), dgidb_druggable = TRUE)]
  dgidb_layer <- dgidb_layer[!duplicated(human_symbol)]
  cat("  DGIdb genome-wide druggable genes:", nrow(dgidb_layer), "\n")
} else {
  cat("  WARNING: genome-wide DGIdb file not found; dgidb_druggable will be FALSE\n")
  dgidb_layer <- data.table(human_symbol = character(0), dgidb_druggable = logical(0))
}

# ================================================================
# Layer: Proteomics evidence
# ================================================================
cat("Loading proteomics concordance...\n")
prot_path <- file.path(BASE, "Analysis/Proteomics/results/protein_transcript_concordance_v3.csv")
if (file.exists(prot_path)) {
  prot_raw <- fread(prot_path)
  # Keep TRUE proteomics only. GSE276114 is liver RNA-seq (transcriptomics); its
  # presence in this concordance file was the prior RNA-vs-RNA contamination bug
  # (see mrna_protein_concordance.R header). Restrict to PRIDE DIA-MS datasets so
  # the proteomics evidence layer never counts an RNA-seq cohort as a protein dataset.
  prot_raw <- prot_raw[grepl("^PXD", dataset)]
  # Aggregate to gene level: n_prot_datasets, best_protein_logFC, best_protein_padj
  #
  # FIX 2026-06-01 (review G8-002/003/007):
  #  (1) The `dataset` column conflates a physical dataset with its contrast
  #      (PXD051911 / PXD051911_mash_vs_masl / PXD051911_nas_high_vs_low all =
  #      ONE study). Collapse to the physical accession so n_prot_datasets counts
  #      studies (<=3), not contrasts (was up to 6 -> ~2x inflation).
  #  (2) best_protein_logFC and best_protein_padj must come from the SAME row, and
  #      from the PRIMARY disease-vs-control contrast — not an uncorrected min-padj
  #      across up to 6 correlated contrasts (which inflated significance and could
  #      report a within-stage contrast's effect as the gene's disease result).
  prot_raw[, dataset_phys := sub("^(GSE[0-9]+|PXD[0-9]+).*$", "\\1", dataset)]
  n_ds <- prot_raw[, .(n_prot_datasets = uniqueN(dataset_phys)),
                   by = .(human_symbol = gene)]
  # primary contrast = disease_vs_control; fall back to any contrast only for
  # genes that have no disease_vs_control row (so genes are not dropped).
  prot_dvc      <- prot_raw[dream_comparator == "disease_vs_control" & !is.na(protein_padj)]
  prot_fallback <- prot_raw[!(gene %in% unique(prot_dvc$gene)) & !is.na(protein_padj)]
  prot_best_src <- rbind(prot_dvc, prot_fallback)
  best <- prot_best_src[order(protein_padj),
                        .(best_protein_logFC = protein_logFC[1],
                          best_protein_padj  = protein_padj[1]),
                        by = .(human_symbol = gene)]
  prot_layer <- merge(n_ds, best, by = "human_symbol", all.x = TRUE)
  prot_layer <- prot_layer[!duplicated(human_symbol)]
  cat(sprintf("  Proteomics: %d genes across %d physical datasets (max n_prot_datasets = %d)\n",
              nrow(prot_layer), uniqueN(prot_raw$dataset_phys), max(prot_layer$n_prot_datasets)))
} else {
  cat("  Proteomics concordance not found\n")
  prot_layer <- data.table(human_symbol = character(0),
                            n_prot_datasets = integer(0),
                            best_protein_logFC = numeric(0),
                            best_protein_padj = numeric(0))
}

# ================================================================
# Layer: Spatial transcriptomics evidence
# ================================================================
cat("Loading spatial validation results...\n")
spatial_dir <- file.path(BASE, "Analysis/Spatial/results/validation_bulk")
spatial_files <- list.files(spatial_dir, pattern = "^spatial_bulk_merged_.*\\.csv$", full.names = TRUE)
if (length(spatial_files) > 0) {
  spatial_all <- rbindlist(lapply(spatial_files, fread), fill = TRUE)
  # Compute per-gene: spatial_max_I (Moran's I), spatial_sig (any dataset padj < 0.05)
  spatial_layer <- spatial_all[, .(
    spatial_max_I  = max(I, na.rm = TRUE),
    spatial_sig    = any(padj_bh < 0.05, na.rm = TRUE),
    spatial_n_datasets = uniqueN(dataset)
  ), by = .(human_symbol = symbol)]
  spatial_layer <- spatial_layer[!duplicated(human_symbol)]
  spatial_layer[is.infinite(spatial_max_I), spatial_max_I := NA_real_]
  cat(sprintf("  Spatial: %d genes from %d datasets, %d spatially variable (padj<0.05)\n",
              nrow(spatial_layer), length(spatial_files),
              sum(spatial_layer$spatial_sig, na.rm = TRUE)))
} else {
  cat("  Spatial validation files not found\n")
  spatial_layer <- data.table(human_symbol = character(0),
                               spatial_max_I = numeric(0),
                               spatial_sig = logical(0),
                               spatial_n_datasets = integer(0))
}

# ================================================================
# MERGE: Assemble unified atlas
# ================================================================
cat("\n=== Assembling unified atlas ===\n")

# Start with all unique human symbols from consensus
# Sort by bulk_padj so the most significant isoform wins deduplication
atlas <- unique(consensus[!is.na(human_symbol) & human_symbol != "",
                           .(human_symbol, ensembl_id = ensembl_clean,
                             bulk_logFC, bulk_padj, bulk_tstat = t,
                             bulk_shrunk_logFC, bulk_lfsr)])
setorder(atlas, bulk_padj, na.last = TRUE)
atlas <- atlas[!duplicated(human_symbol)]
cat("Starting genes (from consensus):", nrow(atlas), "\n")

# Add gene biotype from GENCODE
atlas <- merge(atlas, biotype_map, by.x = "ensembl_id", by.y = "ensembl_clean", all.x = TRUE)
cat(sprintf("  Biotype annotated: %d / %d genes\n",
            sum(!is.na(atlas$gene_biotype)), nrow(atlas)))
cat("  Biotype distribution (top 5):\n")
bt_tab <- atlas[!is.na(gene_biotype), .N, by = gene_biotype][order(-N)]
for (i in seq_len(min(5, nrow(bt_tab)))) {
  cat(sprintf("    %s: %d\n", bt_tab$gene_biotype[i], bt_tab$N[i]))
}

# Add human tier from old multi_evidence (meta-analysis columns removed)
atlas <- merge(atlas, old_me[, .(ensembl_clean, human_consensus_tier = human_tier)],
               by.x = "ensembl_id", by.y = "ensembl_clean", all.x = TRUE)

# Mouse ortholog from concordance
mouse_ortho <- conc_layer[, .(human_symbol, mouse_ortholog = mouse_gene_id)]
atlas <- merge(atlas, mouse_ortho, by = "human_symbol", all.x = TRUE)

# Layer 2: Mouse DE
atlas <- merge(atlas, mouse_layer, by = "human_symbol", all.x = TRUE)

# Layer 3: Concordance
atlas <- merge(atlas, conc_layer[, .(human_symbol, primary_category,
                                      translatability_score, n_concordant_diets,
                                      is_conserved, best_mouse_model)],
               by = "human_symbol", all.x = TRUE)

# Layer 4: Causal inference (TWAS + GTEx COLOC; MR ditched 2026-04-22)
atlas <- merge(atlas, twas_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, coloc_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, susie_layer, by = "human_symbol", all.x = TRUE)

# Layer 4b-f: Extended causal evidence (sc-eQTL, Zenodo, Broadaway, UKBB ALT, ieQTL)
atlas <- merge(atlas, sceqtl_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, zenodo_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, broadaway_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, ukbb_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, ieqtl_layer, by = "human_symbol", all.x = TRUE)
# Layer 4h-j: Additional Broadaway COLOC (AST, GGT, PDFF)
atlas <- merge(atlas, ast_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, ggt_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, pdff_layer, by = "human_symbol", all.x = TRUE)
# Layer 4k-m: FinnGen COLOC (ARCHIVED — empty layers produce all-NA columns for backwards compat)
atlas <- merge(atlas, fg_nafld_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, fg_nash_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, fg_hcc_layer, by = "human_symbol", all.x = TRUE)
# Layer 4n-p: BBJ COLOC (ALT, AST, GGT — cross-ancestry)
atlas <- merge(atlas, bbj_alt_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, bbj_ast_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, bbj_ggt_layer, by = "human_symbol", all.x = TRUE)
# Layer 4q: Ghouse Cirrhosis + HCC (progression)
atlas <- merge(atlas, ghouse_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, ghouse_hcc_layer, by = "human_symbol", all.x = TRUE)
# Layer 4r-t: deCODE COLOC (Icelandic European replication)
atlas <- merge(atlas, decode_nafl_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, decode_cirrhosis_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, decode_hcc_layer, by = "human_symbol", all.x = TRUE)
# Layer 4u-z: Pan-UKBB AFR/CSA COLOC (multi-ancestry expansion)
for (lname in names(panukbb_layers)) {
  atlas <- merge(atlas, panukbb_layers[[lname]], by = "human_symbol", all.x = TRUE)
}
# Layer 4aa-ab: Cross-ancestry joint fine-mapping (SuSiEX + MESuSiE)
atlas <- merge(atlas, susiex_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, mesusie_layer, by = "human_symbol", all.x = TRUE)
# Layer 4g: Pleiotropy
atlas <- merge(atlas, pleio_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, n_coloc_dt, by = "human_symbol", all.x = TRUE)
atlas[is.na(n_coloc_sources), n_coloc_sources := 0L]
atlas[is.na(zenodo_nafld_coloc), zenodo_nafld_coloc := FALSE]
atlas[is.na(ieqtl_disease_interaction), ieqtl_disease_interaction := FALSE]

# Compute liver enzyme COLOC replication metrics (European: UKBB)
# n_liver_enzyme_coloc: count of liver enzyme GWAS with PP.H4>0.5 (ALT, AST, GGT)
atlas[, n_liver_enzyme_coloc := as.integer(
  (!is.na(broadaway_coloc_pp4) & broadaway_coloc_pp4 > 0.5) +
  (!is.na(ast_coloc_pp4) & ast_coloc_pp4 > 0.5) +
  (!is.na(ggt_coloc_pp4) & ggt_coloc_pp4 > 0.5)
)]
# best_liver_enzyme_pp4: max PP.H4 across ALT/AST/GGT
atlas[, best_liver_enzyme_pp4 := pmax(
  fifelse(is.na(broadaway_coloc_pp4), 0, broadaway_coloc_pp4),
  fifelse(is.na(ast_coloc_pp4), 0, ast_coloc_pp4),
  fifelse(is.na(ggt_coloc_pp4), 0, ggt_coloc_pp4)
)]
atlas[best_liver_enzyme_pp4 == 0, best_liver_enzyme_pp4 := NA_real_]

# Cross-ancestry metrics
# n_ancestry_gwas: count of ancestry groups with any causal COLOC signal
# European: any of UKBB ALT/AST/GGT, PDFF, Broadaway, GTEx, Ghouse, deCODE
# (FinnGen removed 2026-04-08: duplicate of Whitfield 2023)
# East Asian (BBJ): any of BBJ ALT/AST/GGT
atlas[, has_eur_coloc := (!is.na(broadaway_coloc_pp4) & broadaway_coloc_pp4 > 0.5) |
                          (!is.na(ast_coloc_pp4) & ast_coloc_pp4 > 0.5) |
                          (!is.na(ggt_coloc_pp4) & ggt_coloc_pp4 > 0.5) |
                          (!is.na(pdff_coloc_pp4) & pdff_coloc_pp4 > 0.5) |
                          (!is.na(coloc_pp4) & coloc_pp4 > 0.5) |
                          (!is.na(ghouse_cirrhosis_coloc_pp4) & ghouse_cirrhosis_coloc_pp4 > 0.5) |
                          (!is.na(ghouse_hcc_coloc_pp4) & ghouse_hcc_coloc_pp4 > 0.5) |
                          (!is.na(decode_nafl_coloc_pp4) & decode_nafl_coloc_pp4 > 0.5) |
                          (!is.na(decode_cirrhosis_coloc_pp4) & decode_cirrhosis_coloc_pp4 > 0.5) |
                          (!is.na(decode_hcc_coloc_pp4) & decode_hcc_coloc_pp4 > 0.5)]
atlas[, has_eas_coloc := (!is.na(bbj_alt_coloc_pp4) & bbj_alt_coloc_pp4 > 0.5) |
                          (!is.na(bbj_ast_coloc_pp4) & bbj_ast_coloc_pp4 > 0.5) |
                          (!is.na(bbj_ggt_coloc_pp4) & bbj_ggt_coloc_pp4 > 0.5)]
# Pan-UKBB African ancestry
atlas[, has_afr_coloc := FALSE]
for (trait in c("alt", "ast", "ggt")) {
  col <- paste0("panukbb_afr_", trait, "_coloc_pp4")
  if (col %in% names(atlas)) {
    atlas[, has_afr_coloc := has_afr_coloc | (!is.na(get(col)) & get(col) > 0.5)]
  }
}
# Pan-UKBB Central/South Asian ancestry
atlas[, has_csa_coloc := FALSE]
for (trait in c("alt", "ast", "ggt")) {
  col <- paste0("panukbb_csa_", trait, "_coloc_pp4")
  if (col %in% names(atlas)) {
    atlas[, has_csa_coloc := has_csa_coloc | (!is.na(get(col)) & get(col) > 0.5)]
  }
}
atlas[, n_ancestry_gwas := as.integer(has_eur_coloc) + as.integer(has_eas_coloc) +
                            as.integer(has_afr_coloc) + as.integer(has_csa_coloc)]
# Validated ancestry count: EUR+EAS only (AFR/CSA are exploratory — underpowered
# and prior-driven per 50c sensitivity analysis)
atlas[, n_validated_ancestry_gwas := as.integer(has_eur_coloc) + as.integer(has_eas_coloc)]
# cross_ancestry_coloc_replication: requires EUR+EAS concordance (validated ancestries)
atlas[, cross_ancestry_coloc_replication := n_validated_ancestry_gwas >= 2L]

# ---------------------------------------------------------------------------
# T1.10 sample-overlap dedup 2026-04-22
# ---------------------------------------------------------------------------
# The legacy n_ancestry_gwas / n_validated_ancestry_gwas counts are
# OR-aggregated within ancestry across many EUR GWAS that share samples:
#   - Ghodsian (coloc_pp4 = Ghodsian 2024 MASLD meta) ⊃ Namjou
#   - Ghodsian ↔ Sveinbjornsson share ~395K UKBB controls
#   - UKBB liver enzyme traits (ALT/AST/GGT/PDFF) derive from overlapping UKBB
#     samples (same cohort, different phenotypes)
# Treat these as a single "UKBB-EUR" sample-overlap group and take max(pp4)
# across members rather than OR-ing (which double-counts evidence).
#
# Non-overlapping independent EUR panels (separate cohorts):
#   - Broadaway (eQTL panel; cis-eQTLs from 1,183 non-UKBB donors)
#   - Ghouse (FinnGen — Finnish isolate)
#   - Decode (Iceland, independent cohort)
#
# Non-EUR ancestries are already each a single group (BBJ, Pan-UKBB AFR/CSA).
# ---------------------------------------------------------------------------
# Helper: max(pp4) across a set of columns, NA-safe
.t110_group_hit <- function(atlas_dt, cols, thr = 0.5) {
  cols_present <- intersect(cols, names(atlas_dt))
  if (length(cols_present) == 0L) return(rep(FALSE, nrow(atlas_dt)))
  mat <- as.matrix(atlas_dt[, ..cols_present])
  mat[is.na(mat)] <- 0
  apply(mat, 1L, function(row) any(row > thr))
}

# EUR sample-overlap groups
# Group A: UKBB-derived (liver enzymes + PDFF + MASLD meta overlap)
t110_ukbb_cols <- c("ukbb_alt_coloc_pp4", "broadaway_coloc_pp4",
                    "ast_coloc_pp4", "ggt_coloc_pp4", "pdff_coloc_pp4",
                    "coloc_pp4")  # coloc_pp4 = Ghodsian (UKBB-overlapping)
t110_ukbb_hit <- .t110_group_hit(atlas, t110_ukbb_cols)

# Group B: FinnGen-derived (Finnish isolate, independent of UKBB)
t110_finngen_cols <- c("finngen_nafld_coloc_pp4", "finngen_nash_coloc_pp4",
                       "finngen_hcc_coloc_pp4", "ghouse_cirrhosis_coloc_pp4",
                       "ghouse_hcc_coloc_pp4")
t110_finngen_hit <- .t110_group_hit(atlas, t110_finngen_cols)

# Group C: deCODE (Iceland, independent)
t110_decode_cols <- c("decode_nafl_coloc_pp4", "decode_cirrhosis_coloc_pp4",
                      "decode_hcc_coloc_pp4")
t110_decode_hit <- .t110_group_hit(atlas, t110_decode_cols)

# Non-EUR groups (each ancestry = one independent group)
t110_bbj_hit <- .t110_group_hit(atlas, c("bbj_alt_coloc_pp4",
                                         "bbj_ast_coloc_pp4",
                                         "bbj_ggt_coloc_pp4"))
t110_afr_hit <- .t110_group_hit(atlas, c("panukbb_afr_alt_coloc_pp4",
                                         "panukbb_afr_ast_coloc_pp4",
                                         "panukbb_afr_ggt_coloc_pp4"))
t110_csa_hit <- .t110_group_hit(atlas, c("panukbb_csa_alt_coloc_pp4",
                                         "panukbb_csa_ast_coloc_pp4",
                                         "panukbb_csa_ggt_coloc_pp4"))

# Dedup count: one hit per independent sample-overlap group
atlas[, n_ancestry_gwas_dedup := as.integer(t110_ukbb_hit) +
                                 as.integer(t110_finngen_hit) +
                                 as.integer(t110_decode_hit) +
                                 as.integer(t110_bbj_hit) +
                                 as.integer(t110_afr_hit) +
                                 as.integer(t110_csa_hit)]

# Validated dedup (drop underpowered/prior-fragile AFR+CSA, keep EUR groups + BBJ EAS)
atlas[, n_validated_ancestry_gwas_dedup := as.integer(t110_ukbb_hit) +
                                           as.integer(t110_finngen_hit) +
                                           as.integer(t110_decode_hit) +
                                           as.integer(t110_bbj_hit)]

# Cross-ancestry replication requires ≥1 non-EUR group + ≥1 EUR group
# (avoids declaring EUR-EUR overlap as cross-ancestry)
t110_any_eur_dedup <- t110_ukbb_hit | t110_finngen_hit | t110_decode_hit
t110_any_non_eur   <- t110_bbj_hit  # EAS only in validated set
atlas[, cross_ancestry_coloc_replication_dedup := t110_any_eur_dedup & t110_any_non_eur]

rm(.t110_group_hit, t110_ukbb_cols, t110_ukbb_hit, t110_finngen_cols, t110_finngen_hit,
   t110_decode_cols, t110_decode_hit, t110_bbj_hit, t110_afr_hit, t110_csa_hit,
   t110_any_eur_dedup, t110_any_non_eur)
# --- end T1.10 dedup block ---

# Disease progression COLOC metrics
# has_progression_coloc: any cirrhosis or HCC GWAS signal
# (FinnGen HCC removed 2026-04-08: duplicate of Whitfield 2023)
atlas[, has_progression_coloc :=
  (!is.na(ghouse_cirrhosis_coloc_pp4) & ghouse_cirrhosis_coloc_pp4 > 0.5) |
  (!is.na(ghouse_hcc_coloc_pp4) & ghouse_hcc_coloc_pp4 > 0.5) |
  (!is.na(decode_cirrhosis_coloc_pp4) & decode_cirrhosis_coloc_pp4 > 0.5) |
  (!is.na(decode_hcc_coloc_pp4) & decode_hcc_coloc_pp4 > 0.5)]

# Clean up intermediate columns
atlas[, c("has_eur_coloc", "has_eas_coloc", "has_afr_coloc", "has_csa_coloc") := NULL]

# Layer 5: Sex stratification
atlas <- merge(atlas, sex_layer, by = "human_symbol", all.x = TRUE)

# Layer 6: Pathways
atlas <- merge(atlas, pathway_layer, by = "human_symbol", all.x = TRUE)

# Layer 7: Essentiality
atlas <- merge(atlas, ess_layer, by = "human_symbol", all.x = TRUE)

# Annotations
atlas <- merge(atlas, deconv_layer, by = "human_symbol", all.x = TRUE)
atlas <- merge(atlas, drug_layer, by = "human_symbol", all.x = TRUE)
# Genome-wide DGIdb druggability (fix 2026-06-10) — own layer, fill non-matches FALSE
atlas <- merge(atlas, dgidb_layer, by = "human_symbol", all.x = TRUE)
atlas[is.na(dgidb_druggable), dgidb_druggable := FALSE]

# Proteomics validation
atlas <- merge(atlas, prot_layer, by = "human_symbol", all.x = TRUE)

# Spatial validation
atlas <- merge(atlas, spatial_layer, by = "human_symbol", all.x = TRUE)

cat("Atlas genes after all merges:", nrow(atlas), "\n")

# ================================================================
# Compute layers_active
# ================================================================
# Define what counts as "active" for each layer
atlas[, layers_active := 0L]

# L1: DEG significance — lfsr < 0.05 and |shrunk_logFC| > 0.5 (ashr canonical, 2026-06-02)
atlas[, l1_active := !is.na(bulk_lfsr) & bulk_lfsr < 0.05 & abs(bulk_shrunk_logFC) > 0.5]
# L2: Mouse — EXCLUDED from layers_active (supplementary only)
# L3: has concordance data
atlas[, l3_active := !is.na(primary_category) & primary_category != "" & primary_category != "Not_Significant"]
# L4: any causal signal (TWAS/COLOC/sc-eQTL/Zenodo/Broadaway/ieQTL/BBJ/Ghouse/deCODE)
# (FinnGen removed 2026-04-08: duplicate of Whitfield 2023)
# (MR removed 2026-04-22: ditched from paper; TWAS + COLOC + INTACT is the causal framework)
atlas[, l4_active := (!is.na(twas_pval) & twas_pval < 0.05) |
                      (!is.na(coloc_pp4) & coloc_pp4 > 0.5) |
                      (!is.na(sceqtl_coloc_best_pp4) & sceqtl_coloc_best_pp4 > 0.5) |
                      (!is.na(broadaway_coloc_pp4) & broadaway_coloc_pp4 > 0.5) |
                      (!is.na(ukbb_alt_coloc_pp4) & ukbb_alt_coloc_pp4 > 0.5) |
                      (!is.na(ast_coloc_pp4) & ast_coloc_pp4 > 0.5) |
                      (!is.na(ggt_coloc_pp4) & ggt_coloc_pp4 > 0.5) |
                      (!is.na(pdff_coloc_pp4) & pdff_coloc_pp4 > 0.5) |
                      (!is.na(bbj_alt_coloc_pp4) & bbj_alt_coloc_pp4 > 0.5) |
                      (!is.na(bbj_ast_coloc_pp4) & bbj_ast_coloc_pp4 > 0.5) |
                      (!is.na(bbj_ggt_coloc_pp4) & bbj_ggt_coloc_pp4 > 0.5) |
                      (!is.na(ghouse_cirrhosis_coloc_pp4) & ghouse_cirrhosis_coloc_pp4 > 0.5) |
                      (!is.na(ghouse_hcc_coloc_pp4) & ghouse_hcc_coloc_pp4 > 0.5) |
                      (!is.na(decode_nafl_coloc_pp4) & decode_nafl_coloc_pp4 > 0.5) |
                      (!is.na(decode_cirrhosis_coloc_pp4) & decode_cirrhosis_coloc_pp4 > 0.5) |
                      (!is.na(decode_hcc_coloc_pp4) & decode_hcc_coloc_pp4 > 0.5) |
                      (zenodo_nafld_coloc == TRUE) |
                      (ieqtl_disease_interaction == TRUE)]
# Phase 5 overhaul columns (added by Script 75)
# (mr_ivw_pval and mr_presso_pval removed 2026-04-22; MR ditched from paper)
for (.col in c("otters_broadaway_pval")) {
  if (.col %in% names(atlas)) {
    atlas[, l4_active := l4_active | (!is.na(get(.col)) & get(.col) < 0.05)]
  }
}
# cTWAS dropped (review A09#1). HyPrColoc l4 is set by Script 75's method-gated recalc
# (review A08#2) — excluded here to avoid a non-method-gated double-count.
for (.col in c("coloc_susie_best_pp4")) {
  if (.col %in% names(atlas)) {
    atlas[, l4_active := l4_active | (!is.na(get(.col)) & get(.col) > 0.5)]
  }
}
# Also include Pan-UKBB COLOC in l4_active
for (.pop in c("afr", "csa")) {
  for (.trait in c("alt", "ast", "ggt")) {
    .col <- paste0("panukbb_", .pop, "_", .trait, "_coloc_pp4")
    if (.col %in% names(atlas)) {
      atlas[, l4_active := l4_active | (!is.na(get(.col)) & get(.col) > 0.5)]
    }
  }
}
# L5: has sex-dimorphic classification (interaction-based)
# 2026-05-28 P0 review fix (P0-C): switched from a negative-exclusion list to a
# POSITIVE-inclusion list. The v3 sex classification emits "Uncertain" (~29,171
# genes) and "" (empty string, ~7,960 genes), neither of which was in the old
# exclusion set c("Not_significant","Concordant") — so ~29,295 non-dimorphic
# genes were wrongly flagged L5-active and inflated layers_active for 79% of genes.
# Only Female_biased / Male_biased / Divergent are genuinely sex-dimorphic (~41 genes).
atlas[, l5_active := sex_class %in% c("Female_biased", "Male_biased", "Divergent")]
# L6: gene is in >= 1 leading edge pathway
atlas[, l6_active := !is.na(n_leading_edge_pathways) & n_leading_edge_pathways >= 1]
# L7: has essentiality data
atlas[, l7_active := !is.na(essentiality_chronos)]

atlas[, layers_active := as.integer(l1_active) +
                          as.integer(l3_active) + as.integer(l4_active) +
                          as.integer(l5_active) + as.integer(l6_active) +
                          as.integer(l7_active)]

# Capture l4_active_count before cleanup (used in summary stats below)
l4_active_count <- sum(atlas$l4_active, na.rm = TRUE)

# Clean up intermediate active columns
atlas[, c("l1_active", "l3_active", "l4_active",
          "l5_active", "l6_active", "l7_active") := NULL]

# ================================================================
# Reorder columns to match schema
# ================================================================
# Priority columns go first; all remaining columns are appended alphabetically
# This avoids dropping new columns when the list isn't updated
priority_cols <- c(
  # Identifiers
  "human_symbol", "ensembl_id", "gene_biotype", "mouse_ortholog",
  # L1: Human DE
  "bulk_logFC", "bulk_padj", "bulk_tstat",
  "bulk_shrunk_logFC", "bulk_lfsr",
  "human_consensus_tier",
  # L2: Mouse DE
  "mouse_meta_logFC", "mouse_meta_padj", "n_diets_sig", "mouse_consensus_tier",
  # L3: Cross-species concordance
  "primary_category", "translatability_score", "n_concordant_diets",
  "is_conserved", "best_mouse_model",
  # L4: Causal inference — TWAS/GTEx (MR columns ditched 2026-04-22)
  "twas_z", "twas_pval", "coloc_pp4",
  # L4: ABF COLOC (multi-GWAS aggregated)
  "coloc_abf_best_pp4", "coloc_abf_best_gwas",
  "coloc_abf_n_gwas_h4_05", "coloc_abf_n_gwas_h4_08",
  "coloc_n_gwas_tested", "coloc_n_groups_h4_05", "coloc_n_groups_h4_08",
  # L4: SuSiE-COLOC (real SuSiE when available, ABF fallback)
  "coloc_susie_best_pp4", "coloc_susie_n_signals", "coloc_susie_best_gwas",
  "coloc_susie_n_gwas_h4_05", "coloc_susie_n_gwas_h4_08", "coloc_susie_n_gwas_h4_09",
  "coloc_susie_n_pairs_total", "coloc_susie_success_rate",
  "coloc_n_groups_susie_h4_05", "coloc_n_groups_susie_h4_08",
  # L4: Causal inference — sc-eQTL
  "sceqtl_coloc_pp4_hep", "sceqtl_coloc_best_pp4", "sceqtl_coloc_cell_type",
  "sceqtl_n_cell_types",
  # L4: Causal inference — Broadaway COLOC
  "broadaway_coloc_pp4",
  "ast_coloc_pp4", "ggt_coloc_pp4", "pdff_coloc_pp4",
  "ukbb_alt_coloc_pp4", "ukbb_alt_coloc_cell_type",
  "n_liver_enzyme_coloc", "best_liver_enzyme_pp4",
  # L4: Causal inference — FinnGen COLOC (ARCHIVED 2026-04-08: all-NA, kept for backwards compat)
  "finngen_nafld_coloc_pp4", "finngen_nash_coloc_pp4", "finngen_hcc_coloc_pp4",
  # L4: Causal inference — BBJ COLOC (cross-ancestry, EAS)
  "bbj_alt_coloc_pp4", "bbj_ast_coloc_pp4", "bbj_ggt_coloc_pp4",
  # L4: Causal inference — Ghouse progression COLOC
  "ghouse_cirrhosis_coloc_pp4", "ghouse_hcc_coloc_pp4",
  # L4: Causal inference — deCODE COLOC (Icelandic replication)
  "decode_nafl_coloc_pp4", "decode_cirrhosis_coloc_pp4", "decode_hcc_coloc_pp4",
  # L4: Causal inference — cross-ancestry joint fine-mapping
  "susiex_max_pip", "susiex_n_loci", "susiex_n_cs",
  "susiex_cs_size_joint", "susiex_trait_pairs",
  "mesusie_max_pip", "mesusie_max_pip_shared",
  "mesusie_in_shared_cs", "mesusie_in_eur_cs", "mesusie_in_eas_cs",
  "mesusie_n_loci", "mesusie_trait_pairs",
  # L4: Causal inference — derived COLOC metrics
  "n_coloc_sources", "n_ancestry_gwas", "cross_ancestry_coloc_replication",
  "has_progression_coloc",
  # L4: Causal inference — other
  "zenodo_nafld_coloc", "zenodo_coloc_cell_types", "zenodo_n_traits_coloc",
  "ieqtl_disease_interaction", "ieqtl_interaction_pval", "ieqtl_cell_type",
  "pleiotropy_class", "pleiotropy_n_traits", "pleiotropy_n_domains",
  "pleiotropy_domains", "is_nafld_specific",
  # L5: Sex stratification
  "sex_class", "bulk_logFC_M", "bulk_logFC_F", "sex_interaction_padj",
  # L6: Pathway biology
  "n_leading_edge_pathways", "top_pathways",
  # L7: Safety/essentiality
  "essentiality_chronos", "is_essential", "n_liver_lines",
  # Annotations
  "attribution_class", "dgidb_druggable", "opentargets_drug", "lincs_reversal",
  # Proteomics validation
  "n_prot_datasets", "best_protein_logFC", "best_protein_padj",
  # Spatial validation
  "spatial_max_I", "spatial_sig", "spatial_n_datasets",
  # Summary
  "layers_active"
)

# Keep priority columns that exist, then append all remaining columns
priority_present <- priority_cols[priority_cols %in% names(atlas)]
remaining_cols <- setdiff(names(atlas), priority_present)
remaining_cols <- sort(remaining_cols)  # alphabetical for consistency
col_order <- c(priority_present, remaining_cols)

cat("  Atlas columns:", length(col_order), "\n")
cat("  Priority columns present:", length(priority_present), "/", length(priority_cols), "\n")
cat("  Additional columns (auto-appended):", length(remaining_cols), "\n")
if (length(remaining_cols) > 0) {
  cat("    ", paste(head(remaining_cols, 20), collapse = ", "),
      if (length(remaining_cols) > 20) paste0("... (+", length(remaining_cols) - 20, " more)") else "", "\n")
}

atlas <- atlas[, ..col_order]

# Sort by bulk_padj ascending (most significant first)
atlas <- atlas[order(bulk_padj)]

# ================================================================
# Output
# ================================================================
out_file <- file.path(OUTDIR, "multi_evidence_atlas.csv")

# Re-run guard: preserve downstream columns added by Scripts 75/217 etc.
# T0.2 (2026-04-22): previously the left-join used bare column names, so columns that
# collided with upstream (e.g. re-derived in this 27a run AND present in the prior file)
# triggered data.table's auto-suffixing to `.x` / `.y`, which then leaked through 217's
# regex and caused duplicate `prog_*` columns. Fix: strip any collisions from the prior
# snapshot BEFORE merging, preferring the freshly computed atlas column.
if (file.exists(out_file)) {
  prior_header <- names(fread(out_file, nrows = 0))
  downstream_cols <- setdiff(prior_header, names(atlas))
  # Also drop any bare-name duplicates already masquerading under .x/.y in prior
  downstream_cols <- downstream_cols[!grepl("\\.[xy]$", downstream_cols)]
  # Drop renamed columns that the re-run guard would otherwise resurrect from
  # the prior atlas. 2026-04-27 LFC=0.5 migration also renamed
  # is_conserved_core -> is_conserved; the prior file still holds the old name.
  # 2026-06-08 dream->bulk hard cutover: 27a now writes these 7 columns FRESH as
  # bulk_* (S1 disease-vs-control + sex logFC_M/F). Blocklist their old dream_* names
  # so the re-run guard does NOT resurrect the stale dream-method values from the
  # prior atlas snapshot alongside the new bulk_* columns.
  renamed_legacy <- c("is_conserved_core",
                      "dream_logFC", "dream_padj", "dream_tstat",
                      "dream_shrunk_logFC", "dream_lfsr",
                      "dream_logFC_M", "dream_logFC_F", "dream_robustness_flag")
  if (any(renamed_legacy %in% downstream_cols)) {
    cat(sprintf("  Dropping legacy renamed columns from prior: %s\n",
                paste(intersect(renamed_legacy, downstream_cols), collapse = ", ")))
    downstream_cols <- setdiff(downstream_cols, renamed_legacy)
  }
  # 2026-05-28 P0 review fix (P0-G): MR (Mendelian Randomization) was permanently
  # retired 2026-04-22 and is no longer loaded above (see L206-209). The prior
  # atlas snapshot still carries stale mr_* columns (mr_beta/mr_pval/mr_sig/
  # mr_ivw_beta/mr_ivw_pval/mr_n_instruments/mr_n_gwas_sig/mr_bidirectional_pval)
  # with ~3,272 non-NA genes; without this blocklist the re-run guard would
  # resurrect them. grep("^mr_") catches all of them robustly (no hardcoding).
  mr_zombie_cols <- grep("^mr_", downstream_cols, value = TRUE)
  if (length(mr_zombie_cols) > 0) {
    cat(sprintf("  Dropping %d retired MR (mr_*) columns from prior: %s\n",
                length(mr_zombie_cols), paste(mr_zombie_cols, collapse = ", ")))
    downstream_cols <- setdiff(downstream_cols, mr_zombie_cols)
  }
  if (length(downstream_cols) > 0) {
    cat(sprintf("  Preserving %d downstream columns from prior atlas: %s\n",
                length(downstream_cols),
                paste(head(downstream_cols, 10), collapse = ", ")))
    prior_data <- fread(out_file, select = c("human_symbol", downstream_cols))
    # Defensive: if prior_data still has any columns that collide with freshly built
    # atlas names (shouldn't, since we setdiff'd, but re-assert), drop them to prevent
    # merge's .x/.y suffixing.
    collisions <- intersect(setdiff(names(prior_data), "human_symbol"), names(atlas))
    if (length(collisions) > 0) {
      cat(sprintf("  Dropping %d colliding columns from prior (upstream takes precedence): %s\n",
                  length(collisions), paste(head(collisions, 5), collapse = ", ")))
      prior_data[, (collisions) := NULL]
    }
    atlas <- merge(atlas, prior_data, by = "human_symbol", all.x = TRUE)
    # Re-apply column ordering
    priority_present <- priority_cols[priority_cols %in% names(atlas)]
    remaining_cols <- sort(setdiff(names(atlas), priority_present))
    atlas <- atlas[, c(priority_present, remaining_cols), with = FALSE]
  }
}

if (!is.na(SUBSET_N) && SUBSET_N > 0 && SUBSET_N < nrow(atlas)) {
  cat(sprintf("\n[SUBSET] --subset-n=%d active; truncating atlas to first %d rows for downstream merges (dry-run mode).\n",
              SUBSET_N, SUBSET_N))
  atlas <- atlas[seq_len(SUBSET_N)]
}

fwrite(atlas, out_file)
cat("\nSaved:", out_file, "\n")
cat("  Dimensions:", nrow(atlas), "genes x", ncol(atlas), "columns\n")

# ================================================================
# Summary statistics
# ================================================================
cat("\n=== Atlas Summary ===\n")
cat("Total genes:", nrow(atlas), "\n\n")

cat("Layer coverage (non-NA/non-zero genes):\n")
cat(sprintf("  L1 Human DE (lfsr<0.05, |shrunk_logFC|>0.5): %d (%0.1f%%)\n",
            sum(!is.na(atlas$bulk_lfsr) & atlas$bulk_lfsr < 0.05 & abs(atlas$bulk_shrunk_logFC) > 0.5),
            100 * sum(!is.na(atlas$bulk_lfsr) & atlas$bulk_lfsr < 0.05 & abs(atlas$bulk_shrunk_logFC) > 0.5) / nrow(atlas)))
cat(sprintf("  L2 Mouse DE (mouse_meta_padj<0.1): %d (%0.1f%%)\n",
            sum(!is.na(atlas$mouse_meta_padj) & atlas$mouse_meta_padj < 0.1),
            100 * sum(!is.na(atlas$mouse_meta_padj) & atlas$mouse_meta_padj < 0.1) / nrow(atlas)))
cat(sprintf("  L3 Concordance (any category):      %d (%0.1f%%)\n",
            sum(!is.na(atlas$primary_category)),
            100 * sum(!is.na(atlas$primary_category)) / nrow(atlas)))
# l4_active_count captured above (from l4_active column before cleanup)
cat(sprintf("  L4 Causal (any evidence):             %d (%0.1f%%)\n",
            l4_active_count, 100 * l4_active_count / nrow(atlas)))
cat(sprintf("    - sc-eQTL COLOC (PP.H4>0.5):       %d\n",
            sum(!is.na(atlas$sceqtl_coloc_best_pp4) & atlas$sceqtl_coloc_best_pp4 > 0.5)))
cat(sprintf("    - Zenodo NAFLD COLOC:               %d\n",
            sum(atlas$zenodo_nafld_coloc == TRUE, na.rm = TRUE)))
cat(sprintf("    - Broadaway ALT COLOC (PP.H4>0.5): %d\n",
            sum(!is.na(atlas$broadaway_coloc_pp4) & atlas$broadaway_coloc_pp4 > 0.5)))
cat(sprintf("    - Broadaway AST COLOC (PP.H4>0.5): %d\n",
            sum(!is.na(atlas$ast_coloc_pp4) & atlas$ast_coloc_pp4 > 0.5)))
cat(sprintf("    - Broadaway GGT COLOC (PP.H4>0.5): %d\n",
            sum(!is.na(atlas$ggt_coloc_pp4) & atlas$ggt_coloc_pp4 > 0.5)))
cat(sprintf("    - Broadaway PDFF COLOC (PP.H4>0.5):%d\n",
            sum(!is.na(atlas$pdff_coloc_pp4) & atlas$pdff_coloc_pp4 > 0.5)))
cat(sprintf("    - UKBB ALT COLOC (PP.H4>0.5):      %d\n",
            sum(!is.na(atlas$ukbb_alt_coloc_pp4) & atlas$ukbb_alt_coloc_pp4 > 0.5)))
# FinnGen summary prints removed 2026-04-08: GWAS archived (duplicate of Whitfield 2023)
cat(sprintf("    - BBJ ALT COLOC (PP.H4>0.5, EAS):  %d\n",
            sum(!is.na(atlas$bbj_alt_coloc_pp4) & atlas$bbj_alt_coloc_pp4 > 0.5)))
cat(sprintf("    - BBJ AST COLOC (PP.H4>0.5, EAS):  %d\n",
            sum(!is.na(atlas$bbj_ast_coloc_pp4) & atlas$bbj_ast_coloc_pp4 > 0.5)))
cat(sprintf("    - BBJ GGT COLOC (PP.H4>0.5, EAS):  %d\n",
            sum(!is.na(atlas$bbj_ggt_coloc_pp4) & atlas$bbj_ggt_coloc_pp4 > 0.5)))
cat(sprintf("    - ieQTL disease interaction:        %d\n",
            sum(atlas$ieqtl_disease_interaction == TRUE, na.rm = TRUE)))
cat(sprintf("    - Genes with 2+ COLOC sources:     %d\n",
            sum(atlas$n_coloc_sources >= 2, na.rm = TRUE)))
cat(sprintf("    - Liver enzyme COLOC 2+ GWAS:      %d\n",
            sum(atlas$n_liver_enzyme_coloc >= 2, na.rm = TRUE)))
cat(sprintf("    - Liver enzyme COLOC 3/3 GWAS:     %d\n",
            sum(atlas$n_liver_enzyme_coloc >= 3, na.rm = TRUE)))
cat(sprintf("    - Cross-ancestry replication:       %d\n",
            sum(atlas$cross_ancestry_coloc_replication, na.rm = TRUE)))
cat(sprintf("    - Genes with 2 ancestry groups:    %d\n",
            sum(atlas$n_ancestry_gwas >= 2, na.rm = TRUE)))
# 2026-05-28 P0 review fix (P0-C): report the same positive-inclusion L5 definition
# used by l5_active above (Female_biased/Male_biased/Divergent only), so the printed
# coverage no longer counts "Uncertain"/"" genes as sex-dimorphic.
cat(sprintf("  L5 Sex-dimorphic (interaction):      %d (%0.1f%%)\n",
            sum(atlas$sex_class %in% c("Female_biased", "Male_biased", "Divergent")),
            100 * sum(atlas$sex_class %in% c("Female_biased", "Male_biased", "Divergent")) / nrow(atlas)))
cat(sprintf("  L6 Pathway (>=1 leading edge):      %d (%0.1f%%)\n",
            sum(!is.na(atlas$n_leading_edge_pathways) & atlas$n_leading_edge_pathways >= 1),
            100 * sum(!is.na(atlas$n_leading_edge_pathways) & atlas$n_leading_edge_pathways >= 1) / nrow(atlas)))
cat(sprintf("  L7 Essentiality (has data):         %d (%0.1f%%)\n",
            sum(!is.na(atlas$essentiality_chronos)),
            100 * sum(!is.na(atlas$essentiality_chronos)) / nrow(atlas)))
cat(sprintf("  Proteomics (any dataset):           %d (%0.1f%%)\n",
            sum(!is.na(atlas$n_prot_datasets) & atlas$n_prot_datasets > 0),
            100 * sum(!is.na(atlas$n_prot_datasets) & atlas$n_prot_datasets > 0) / nrow(atlas)))
cat(sprintf("  Spatial (any SVG padj<0.05):        %d (%0.1f%%)\n",
            sum(atlas$spatial_sig == TRUE, na.rm = TRUE),
            100 * sum(atlas$spatial_sig == TRUE, na.rm = TRUE) / nrow(atlas)))

cat("\nLayers active distribution:\n")
print(atlas[, .N, by = layers_active][order(layers_active)])

cat("\nConserved:", sum(atlas$is_conserved, na.rm = TRUE), "genes\n")

# lncRNA verification
if ("gene_biotype" %in% names(atlas)) {
  lnc <- atlas[gene_biotype == "lncRNA"]
  cat(sprintf("\nlncRNA verification:\n"))
  cat(sprintf("  lncRNAs in atlas: %d\n", nrow(lnc)))
  cat(sprintf("  with dream DEG (lfsr<0.05, |shrunk_logFC|>0.5): %d\n",
              sum(!is.na(lnc$bulk_lfsr) & lnc$bulk_lfsr < 0.05 & abs(lnc$bulk_shrunk_logFC) > 0.5)))
  cat(sprintf("  with Broadaway COLOC: %d\n",
              sum(!is.na(lnc$broadaway_coloc_pp4))))
  cat(sprintf("  with Broadaway PP.H4 > 0.5: %d\n",
              sum(lnc$broadaway_coloc_pp4 > 0.5, na.rm = TRUE)))
  cat(sprintf("  with TWAS signal (p<0.05): %d\n",
              sum(!is.na(lnc$twas_pval) & lnc$twas_pval < 0.05)))
  cat(sprintf("  with any L4 causal signal: %d\n",
              sum((!is.na(lnc$twas_pval) & lnc$twas_pval < 0.05) |
                  (!is.na(lnc$coloc_pp4) & lnc$coloc_pp4 > 0.5) |
                  (!is.na(lnc$broadaway_coloc_pp4) & lnc$broadaway_coloc_pp4 > 0.5) |
                  (!is.na(lnc$sceqtl_coloc_best_pp4) & lnc$sceqtl_coloc_best_pp4 > 0.5))))
  # Check specific lncRNAs
  known_lnc <- c("NEAT1", "MALAT1", "MEG3", "H19", "XIST", "HOTAIR")
  in_atlas <- known_lnc[known_lnc %in% atlas$human_symbol]
  cat(sprintf("  Known lncRNAs in atlas: %s\n",
              if (length(in_atlas) > 0) paste(in_atlas, collapse = ", ") else "none"))
}

# Validation: check known MASLD genes
known_genes <- c("COL1A1", "ACTA2", "TGFB1", "CIDEC", "FAP", "LGALS3", "FASN",
                 "THRB", "PPARG", "ACACB")
cat("\nValidation — known MASLD genes in atlas:\n")
for (g in known_genes) {
  row <- atlas[human_symbol == g]
  if (nrow(row) > 0) {
    cat(sprintf("  %s: bulk_padj=%.2e, layers_active=%d, conserved=%s\n",
                g, row$bulk_padj[1], row$layers_active[1],
                as.character(row$is_conserved[1])))
  } else {
    cat(sprintf("  %s: not in atlas\n", g))
  }
}

# ================================================================
# Layer correlation matrix (for downstream visualization)
# ================================================================
cat("\nComputing layer correlations...\n")

# Create continuous proxies for correlation
cor_data <- atlas[, .(
  L1 = -log10(pmax(bulk_lfsr, 1e-300, na.rm = TRUE)) * sign(bulk_shrunk_logFC),
  L2 = ifelse(is.na(mouse_meta_logFC), 0, -log10(pmax(mouse_meta_padj, 1e-300, na.rm = TRUE)) * sign(mouse_meta_logFC)),
  L3 = ifelse(is.na(n_concordant_diets), 0, n_concordant_diets),
  L4 = ifelse(is.na(twas_z), 0, abs(twas_z)) +
       # MR term removed 2026-04-22 (MR ditched from paper)
       ifelse(!is.na(sceqtl_coloc_best_pp4), sceqtl_coloc_best_pp4 * 2, 0) +
       ifelse(!is.na(broadaway_coloc_pp4), broadaway_coloc_pp4 * 2, 0) +
       ifelse(!is.na(ast_coloc_pp4), ast_coloc_pp4 * 2, 0) +
       ifelse(!is.na(ggt_coloc_pp4), ggt_coloc_pp4 * 2, 0) +
       ifelse(!is.na(pdff_coloc_pp4), pdff_coloc_pp4 * 2, 0) +
       ifelse(!is.na(ukbb_alt_coloc_pp4), ukbb_alt_coloc_pp4 * 2, 0) +
       # FinnGen removed 2026-04-08: duplicate of Whitfield 2023 (columns all-NA)
       ifelse(!is.na(bbj_alt_coloc_pp4), bbj_alt_coloc_pp4 * 2, 0) +
       ifelse(!is.na(bbj_ast_coloc_pp4), bbj_ast_coloc_pp4 * 2, 0) +
       ifelse(!is.na(bbj_ggt_coloc_pp4), bbj_ggt_coloc_pp4 * 2, 0) +
       ifelse(zenodo_nafld_coloc == TRUE, 2, 0) +
       ifelse(ieqtl_disease_interaction == TRUE, 1, 0),
  L5 = ifelse(is.na(sex_class) | sex_class %in% c("Not_significant", "Concordant"), 0,
              ifelse(sex_class %in% c("Male_biased", "Female_biased",
                                       "Male_specific", "Female_specific"), 0.5,
                     ifelse(sex_class == "Divergent", 1.0, 0))),
  L6 = ifelse(is.na(n_leading_edge_pathways), 0, n_leading_edge_pathways),
  L7 = ifelse(is.na(essentiality_chronos), 0, essentiality_chronos)
)]

# Cap extreme values
cor_data[L1 > 50, L1 := 50]
cor_data[L1 < -50, L1 := -50]
cor_data[L2 > 50, L2 := 50]
cor_data[L2 < -50, L2 := -50]

cor_mat <- cor(as.matrix(cor_data), use = "pairwise.complete.obs")
colnames(cor_mat) <- rownames(cor_mat) <- c("L1_Human_DE", "L2_Mouse_DE", "L3_Concordance",
                                             "L4_Causal", "L5_Sex", "L6_Pathway", "L7_Safety")
cat("Layer correlation matrix:\n")
print(round(cor_mat, 3))

# Save correlation matrix
cor_out <- as.data.table(as.data.frame(as.table(cor_mat)))
setnames(cor_out, c("Layer1", "Layer2", "Correlation"))
fwrite(cor_out, file.path(OUTDIR, "atlas_layer_correlations.csv"))

# ================================================================
# Post-assembly: Merge ferroptosis + zonation annotations (if available)
# ================================================================
cat("\n--- Post-assembly: merging ferroptosis + zonation annotations ---\n")

ferr_path <- file.path(dirname(OUTDIR), "pathway_programs/ferroptosis_degs.csv")
zon_path <- file.path(dirname(OUTDIR), "zonation/deg_zonation_classification.csv")

if (file.exists(ferr_path)) {
  ferr <- fread(ferr_path)
  # Keep only ferroptosis_class for genes with non-None class
  ferr_slim <- ferr[, .(human_symbol, ferroptosis_class)]
  ferr_slim <- ferr_slim[ferroptosis_class != "None"]
  if (nrow(ferr_slim) > 0) {
    if ("ferroptosis_class" %in% names(atlas)) atlas[, ferroptosis_class := NULL]
    atlas <- merge(atlas, ferr_slim, by = "human_symbol", all.x = TRUE)
    cat(sprintf("  Ferroptosis: %d genes annotated\n",
                sum(!is.na(atlas$ferroptosis_class))))
  }
} else {
  cat("  Ferroptosis annotations not found (run 41_ferroptosis_lipotoxicity.R first)\n")
}

if (file.exists(zon_path)) {
  zon <- fread(zon_path)
  if ("zonation_class" %in% names(zon)) {
    zon_slim <- zon[zonation_class != "Non-zoned", .(human_symbol, zonation_class)]
    if (nrow(zon_slim) > 0) {
      if ("zonation_class" %in% names(atlas)) atlas[, zonation_class := NULL]
      atlas <- merge(atlas, zon_slim, by = "human_symbol", all.x = TRUE)
      cat(sprintf("  Zonation: %d genes annotated\n",
                  sum(!is.na(atlas$zonation_class))))
    }
  }
} else {
  cat("  Zonation annotations not found (run 43_zonation_classification.R first)\n")
}

# ---------------------------------------------------------------------------
# Mega-analysis validation columns (added 2026-05-19)
# Source: results/mega_validation/concordance/atlas_columns.tsv
# Spec:   docs/superpowers/specs/2026-05-19-mega-validation-design.md
# Adds 13 columns: edgeR-QL, voomLmFit, metafor (logFC/SE/tau2/I2/HKSJ_padj),
#                  mashr (n_sig_cohorts/sharing_class/pan_cohort),
#                  dream_robustness_flag.
# ---------------------------------------------------------------------------
mv_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/mega_validation/concordance/atlas_columns.tsv")
if (file.exists(mv_path)) {
  mv <- fread(mv_path)
  # Arm outputs use versioned Ensembl IDs (ENSG00000000003.17); atlas uses
  # stripped form via ensembl_clean. Translate before merge.
  mv[, ensembl_id := sub("\\.\\d+$", "", gene)]
  mv[, gene := NULL]
  mv_cols <- setdiff(names(mv), "ensembl_id")
  for (col in intersect(mv_cols, names(atlas))) atlas[, (col) := NULL]
  atlas <- merge(atlas, mv, by = "ensembl_id", all.x = TRUE)
  cat(sprintf("  + Merged %d mega_validation columns (%d genes covered)\n",
              length(mv_cols),
              sum(!is.na(atlas$metafor_HKSJ_padj))))
} else {
  cat("  ! mega_validation atlas_columns.tsv not found;",
      "skipping (run mega_validation/05_concordance.R first)\n")
}

# ---------------------------------------------------------------------------
# S8: Perturbation modeling evidence (Phase 4, 2026-05-21)
# Source: Analysis/Perturbation/results/integration/d{1..5}_atlas_columns.csv
# Spec:   Analysis/Perturbation/scripts/shared/atlas_writer.py (Python reference)
#
# Adds 14 perturb_* columns + perturb_n_arms_applicable + perturb_n_arms_passing
# + perturb_validation_tier (composite "—"/partial/validated/strong tier
# computed from APPLICABLE arms only — see compute_validation_tier docstring
# in atlas_writer.py).
#
# Gated by --include-s8 flag (default FALSE) so existing rebuilds don't
# pick up perturb_* cols until the perturbation campaign is wired in.
# When the input CSV is missing, we log + leave that arm's columns NA;
# tier still computes correctly from the applicability rules.
# ---------------------------------------------------------------------------
if (INCLUDE_S8) {
  cat("\n--- S8: Perturbation modeling ingestion (INCLUDE_S8=TRUE) ---\n")
  s8_dir <- file.path(BASE, "Analysis/Perturbation/results/integration")
  s8_files <- list(
    d1 = file.path(s8_dir, "d1_atlas_columns.csv"),
    d2 = file.path(s8_dir, "d2_atlas_columns.csv"),
    d3 = file.path(s8_dir, "d3_atlas_columns.csv"),
    d4 = file.path(s8_dir, "d4_atlas_columns.csv"),
    d5 = file.path(s8_dir, "d5_atlas_columns.csv")
  )
  # 14 columns S8 contributes (mirrors atlas_writer.S8_COLUMNS minus the
  # composite tier, which we compute here ourselves).
  S8_COLUMNS_DATA <- c(
    # D1 (mechanism / downstream robustness)
    "perturb_mechanism_pathway", "perturb_mechanism_coa",
    "perturb_downstream_robust_n", "perturb_downstream_robust_genes",
    # D2 (LINCS / signature reversal)
    "perturb_reversal_score", "perturb_reversal_rank",
    "perturb_reversal_stage_specific",
    # D3 (combinatorial synergy)
    "perturb_synergy_partners", "perturb_synergy_top_n", "perturb_synergy_class",
    # D4 (CCC ligand prediction)
    "perturb_ccc_ligand_role", "perturb_ccc_predicted_receivers",
    # D5 (cross-species concordance)
    "perturb_crossspecies_concordance"
  )

  .read_arm_csv <- function(path, arm) {
    if (!file.exists(path)) {
      cat(sprintf("  WARNING: missing %s; %s columns will be NA\n", path, arm))
      return(NULL)
    }
    # Empty file (0-byte) or header-only -> skip
    if (file.info(path)$size == 0L) {
      cat(sprintf("  WARNING: %s is empty; %s columns will be NA\n", path, arm))
      return(NULL)
    }
    df <- tryCatch(fread(path), error = function(e) {
      cat(sprintf("  WARNING: failed to read %s (%s); %s columns will be NA\n",
                  path, conditionMessage(e), arm))
      return(NULL)
    })
    if (is.null(df) || nrow(df) == 0L) {
      cat(sprintf("  NOTE: %s has 0 data rows; %s columns will be NA\n", path, arm))
      return(NULL)
    }
    # Standardize join column to `human_symbol` (atlas convention).
    if (!"human_symbol" %in% names(df)) {
      if ("gene" %in% names(df)) {
        setnames(df, "gene", "human_symbol")
      } else {
        cat(sprintf("  WARNING: %s has neither human_symbol nor gene column (cols: %s); skipping\n",
                    path, paste(head(names(df), 8), collapse = ", ")))
        return(NULL)
      }
    }
    df <- df[!is.na(human_symbol) & human_symbol != ""]
    df <- df[!duplicated(human_symbol)]
    cat(sprintf("  %s: %d rows from %s (cols: %s)\n",
                arm, nrow(df), basename(path),
                paste(setdiff(names(df), "human_symbol"), collapse = ", ")))
    df
  }

  # Drop any pre-existing perturb_* columns so re-runs are idempotent
  # (atlas_writer.py uses overwrite=True semantics; we mirror that here).
  existing_perturb <- grep("^perturb_", names(atlas), value = TRUE)
  if (length(existing_perturb) > 0) {
    cat(sprintf("  Overwrite: dropping %d existing perturb_* cols (%s)\n",
                length(existing_perturb),
                paste(head(existing_perturb, 5), collapse = ", ")))
    atlas[, (existing_perturb) := NULL]
  }

  # Per-arm merge.
  for (arm in names(s8_files)) {
    arm_df <- .read_arm_csv(s8_files[[arm]], arm)
    if (is.null(arm_df)) next
    arm_cols <- setdiff(names(arm_df), "human_symbol")
    # Only carry forward known S8 schema columns (defensive against extra cols)
    keep <- intersect(arm_cols, S8_COLUMNS_DATA)
    if (length(keep) == 0L) {
      cat(sprintf("  NOTE: %s has no S8 schema columns (got %s); skipping\n",
                  arm, paste(arm_cols, collapse = ", ")))
      next
    }
    # Drop collisions before merge to avoid data.table .x/.y suffixing
    collisions <- intersect(keep, names(atlas))
    if (length(collisions) > 0) atlas[, (collisions) := NULL]
    atlas <- merge(atlas, arm_df[, c("human_symbol", keep), with = FALSE],
                   by = "human_symbol", all.x = TRUE)
  }

  # Pre-populate any S8 schema column that no arm wrote (so tier logic
  # downstream is uniform regardless of which arms produced output).
  for (col in S8_COLUMNS_DATA) {
    if (!col %in% names(atlas)) atlas[, (col) := NA]
  }

  # ---- Composite validation tier (mirror of atlas_writer.compute_validation_tier)
  # Applicability rules:
  #   D1 / D2 / D3: apply to every gene
  #   D4: applies only if perturb_ccc_ligand_role is non-NA (hepatocyte ligands)
  #   D5: applies only if perturb_crossspecies_concordance is non-NA
  # Pass rules:
  #   D1: perturb_downstream_robust_n >= 1
  #   D2: |perturb_reversal_score| > 0
  #   D3: perturb_synergy_top_n >= 1
  #   D4: perturb_ccc_ligand_role == "active" (case-insensitive, non-NA)
  #   D5: perturb_crossspecies_concordance > 0.5
  d1_app <- rep(TRUE, nrow(atlas))
  d2_app <- rep(TRUE, nrow(atlas))
  d3_app <- rep(TRUE, nrow(atlas))
  d4_app <- !is.na(atlas[["perturb_ccc_ligand_role"]])
  d5_app <- !is.na(atlas[["perturb_crossspecies_concordance"]])

  .num <- function(x) suppressWarnings(as.numeric(x))
  .int <- function(x) suppressWarnings(as.integer(x))
  d1_pass <- !is.na(atlas[["perturb_downstream_robust_n"]]) &
             .int(atlas[["perturb_downstream_robust_n"]]) >= 1L
  d2_pass <- !is.na(atlas[["perturb_reversal_score"]]) &
             abs(.num(atlas[["perturb_reversal_score"]])) > 0
  d3_pass <- !is.na(atlas[["perturb_synergy_top_n"]]) &
             .int(atlas[["perturb_synergy_top_n"]]) >= 1L
  d4_pass <- d4_app &
             tolower(as.character(atlas[["perturb_ccc_ligand_role"]])) == "active"
  d5_pass <- d5_app &
             !is.na(atlas[["perturb_crossspecies_concordance"]]) &
             .num(atlas[["perturb_crossspecies_concordance"]]) > 0.5

  # NA-safe coalesce (TRUE iff pass AND applicable AND non-NA)
  d1_pass[is.na(d1_pass)] <- FALSE
  d2_pass[is.na(d2_pass)] <- FALSE
  d3_pass[is.na(d3_pass)] <- FALSE
  d4_pass[is.na(d4_pass)] <- FALSE
  d5_pass[is.na(d5_pass)] <- FALSE

  n_app  <- as.integer(d1_app) + as.integer(d2_app) + as.integer(d3_app) +
            as.integer(d4_app) + as.integer(d5_app)
  n_pass <- as.integer(d1_pass) + as.integer(d2_pass) + as.integer(d3_pass) +
            as.integer(d4_pass) + as.integer(d5_pass)
  frac <- n_pass / pmax(n_app, 1L)

  tier <- rep("—", nrow(atlas))
  tier[n_pass >= 1L & frac > 0 & frac < 0.6] <- "partial"
  tier[frac >= 0.6 & frac < 1.0] <- "validated"
  tier[frac == 1.0 & n_app >= 3L] <- "strong"

  atlas[, perturb_n_arms_applicable := n_app]
  atlas[, perturb_n_arms_passing := n_pass]
  atlas[, perturb_validation_tier := tier]

  cat(sprintf("  S8 ingestion complete: %d genes; tier dist — strong=%d, validated=%d, partial=%d, none=%d\n",
              nrow(atlas),
              sum(tier == "strong"),
              sum(tier == "validated"),
              sum(tier == "partial"),
              sum(tier == "—")))
  any_perturb <- rowSums(!is.na(atlas[, S8_COLUMNS_DATA, with = FALSE])) > 0L
  cat(sprintf("  Genes with >=1 perturb_* value populated: %d\n", sum(any_perturb)))
} else {
  cat("\n--- S8: skipped (INCLUDE_S8=FALSE; pass --include-s8 TRUE to enable)\n")
}

# ---------------------------------------------------------------------------
# T1.1 Druggability / translation bundle (2026-06-01)
# Source: RNA-seq/results/multi_evidence/druggability_atlas_columns.tsv (Script 76)
# Adds Pharos TDL, gnomAD LOEUF/mis_z, curated MASH clinical pipeline, plus the
# derived genetic_support flag (COLOC PP4>0.5 OR TWAS p<0.05) and minikel_approval_or
# (Minikel 2024 Nature: genetic support -> 2.6x drug-approval odds). Reframes the
# drug-validation arm as a transparent, externally anchored druggability ladder.
# ---------------------------------------------------------------------------
drug_path <- file.path(OUTDIR, "druggability_atlas_columns.tsv")
if (file.exists(drug_path)) {
  drug <- fread(drug_path)
  drug_cols <- setdiff(names(drug), "human_symbol")
  # Defensive drop-then-merge: the re-run guard may have resurrected these from
  # the prior atlas snapshot; the fresh table takes precedence (prevents .x/.y).
  for (col in intersect(drug_cols, names(atlas))) atlas[, (col) := NULL]
  atlas <- merge(atlas, drug, by = "human_symbol", all.x = TRUE)
  atlas[is.na(clintrial_active), clintrial_active := FALSE]

  # Derived genetic-support flag from COLOC/TWAS columns already present.
  coloc_gs_cols <- intersect(c("coloc_best_susie_pp4_polyfun", "coloc_best_pp4_polyfun",
                               "coloc_susie_best_pp4", "coloc_best_pp4", "coloc_abf_best_pp4"),
                             names(atlas))
  gs <- rep(FALSE, nrow(atlas))
  for (cc in coloc_gs_cols) gs <- gs | (!is.na(atlas[[cc]]) & atlas[[cc]] > 0.5)
  if ("twas_pval" %in% names(atlas)) gs <- gs | (!is.na(atlas$twas_pval) & atlas$twas_pval < 0.05)
  for (col in c("genetic_support", "minikel_approval_or"))
    if (col %in% names(atlas)) atlas[, (col) := NULL]
  atlas[, genetic_support := gs]
  atlas[, minikel_approval_or := fifelse(genetic_support, 2.6, 1.0)]
  cat(sprintf(paste0("  + Druggability: Pharos-TDL=%d gnomAD-LOEUF=%d clinical-trial-active=%d ",
                     "genetic-support=%d (coloc cols used: %s)\n"),
              sum(!is.na(atlas$pharos_tdl)), sum(!is.na(atlas$gnomad_loeuf)),
              sum(atlas$clintrial_active %in% TRUE), sum(atlas$genetic_support),
              paste(coloc_gs_cols, collapse = ",")))
} else {
  cat("  ! druggability_atlas_columns.tsv not found;",
      "skipping (run 76_druggability_translation.R first)\n")
}

# ---------------------------------------------------------------------------
# T1.7 Human Protein Atlas target characterization (2026-06-01)
# Source: RNA-seq/results/multi_evidence/hpa_atlas_columns.tsv (Script 77)
# Adds hpa_rna_tissue_specificity, hpa_liver_elevated, hpa_secretome_location,
# hpa_subcellular_main, hpa_protein_class, hpa_liver_hcc_prognostic.
# ---------------------------------------------------------------------------
hpa_path <- file.path(OUTDIR, "hpa_atlas_columns.tsv")
if (file.exists(hpa_path)) {
  hpa_cols_dt <- fread(hpa_path)
  hpa_new <- setdiff(names(hpa_cols_dt), "human_symbol")
  for (col in intersect(hpa_new, names(atlas))) atlas[, (col) := NULL]
  atlas <- merge(atlas, hpa_cols_dt, by = "human_symbol", all.x = TRUE)
  cat(sprintf("  + HPA: %d genes annotated (%d liver-elevated)\n",
              sum(!is.na(atlas$hpa_rna_tissue_specificity)),
              sum(atlas$hpa_liver_elevated %in% TRUE)))
} else {
  cat("  ! hpa_atlas_columns.tsv not found; skipping (run 77_hpa_annotation.R first)\n")
}

# ---------------------------------------------------------------------------
# T1.4 GTEx liver sQTL credible-set overlap (2026-06-01)
# Source: RNA-seq/results/multi_evidence/sqtl_atlas_columns.tsv (Script 58b)
# Orthogonal genetic axis: a MASLD credible-set variant is a significant liver
# sQTL for the gene. Adds sqtl_credset_hit, sqtl_n_introns, sqtl_min_pval,
# sqtl_max_credset_pip, sqtl_studies. (HSD17B13 = top hit, as expected.)
# ---------------------------------------------------------------------------
sqtl_path <- file.path(OUTDIR, "sqtl_atlas_columns.tsv")
if (file.exists(sqtl_path)) {
  sqtl_dt <- fread(sqtl_path)
  sqtl_new <- setdiff(names(sqtl_dt), "human_symbol")
  for (col in intersect(sqtl_new, names(atlas))) atlas[, (col) := NULL]
  atlas <- merge(atlas, sqtl_dt, by = "human_symbol", all.x = TRUE)
  atlas[is.na(sqtl_credset_hit), sqtl_credset_hit := FALSE]
  cat(sprintf("  + sQTL: %d genes with a liver-sQTL credible-set hit\n",
              sum(atlas$sqtl_credset_hit %in% TRUE)))
} else {
  cat("  ! sqtl_atlas_columns.tsv not found; skipping (run 58b_sqtl_credset_overlap.R first)\n")
}

# ---------------------------------------------------------------------------
# Population / clinical-genetics cluster (2026-06-01; Script 78)
# Source: RNA-seq/results/multi_evidence/popgen_atlas_columns.tsv
# Future-proof: merges ALL columns present (ClinVar MASLD P/LP now; gnomAD AF,
# OMIM monogenic, MASLD rare-variant burden, AlphaMissense added as data lands).
# clinvar_* are MASLD-condition-filtered (relevance gate); see plan v2.
# ---------------------------------------------------------------------------
popgen_path <- file.path(OUTDIR, "popgen_atlas_columns.tsv")
if (file.exists(popgen_path)) {
  pg_dt <- fread(popgen_path)
  pg_new <- setdiff(names(pg_dt), "human_symbol")
  for (col in intersect(pg_new, names(atlas))) atlas[, (col) := NULL]
  atlas <- merge(atlas, pg_dt, by = "human_symbol", all.x = TRUE)
  if ("clinvar_n_plp_masld" %in% names(atlas))
    atlas[is.na(clinvar_n_plp_masld), clinvar_n_plp_masld := 0L]
  cat(sprintf("  + Pop-genetics: %d columns merged (%s)\n",
              length(pg_new), paste(head(pg_new, 8), collapse = ", ")))
} else {
  cat("  ! popgen_atlas_columns.tsv not found; skipping (run 78_population_genetics.R first)\n")
}

# ---------------------------------------------------------------------------
# T2 regulatory variant-to-gene CROSS-VALIDATION overlay (2026-06-01; 55b/c/d/e)
# Source: RNA-seq/results/multi_evidence/regulatory_atlas_columns.tsv
# ABC enhancer->gene + Currin liver caQTL + SCREEN cCRE, gated on MASLD credible
# sets. Tier-2 (validates the atlas's own MASLD scATAC; HepG2 activity = QC only,
# NOT MASLD evidence). abc_/caqtl_/ccre_/regulatory_ columns.
# ---------------------------------------------------------------------------
reg_path <- file.path(OUTDIR, "regulatory_atlas_columns.tsv")
if (file.exists(reg_path)) {
  reg_dt <- fread(reg_path)
  reg_new <- setdiff(names(reg_dt), "human_symbol")
  for (col in intersect(reg_new, names(atlas))) atlas[, (col) := NULL]
  atlas <- merge(atlas, reg_dt, by = "human_symbol", all.x = TRUE)
  for (hc in intersect(c("abc_v2g_hit", "caqtl_credset_hit", "ccre_credset_hit",
                         "regulatory_masld_gwas_driven"), names(atlas)))
    atlas[is.na(get(hc)), (hc) := FALSE]
  if ("n_regulatory_layers" %in% names(atlas))
    atlas[is.na(n_regulatory_layers), n_regulatory_layers := 0L]
  cat(sprintf("  + Regulatory V2G (T2): %d genes (%d MASLD-GWAS-driven)\n",
              sum(atlas$n_regulatory_layers %in% 1:3),
              sum(atlas$regulatory_masld_gwas_driven %in% TRUE)))
} else {
  cat("  ! regulatory_atlas_columns.tsv not found; skipping (run 55b/c/d/e first)\n")
}

# ---------------------------------------------------------------------------
# Secreted-protein / non-invasive-biomarker axis (2026-06-01; Script 79)
# Source: RNA-seq/results/multi_evidence/secreted_protein_atlas_columns.tsv
# Chains tissue DE -> HPA secretion route -> measured MASLD plasma proteomics;
# distinguishes regulated-secretion biomarkers from cell-death leakage markers.
# (pQTL genetic arm omitted — UKB-PPP gated.) secreted_* columns.
# ---------------------------------------------------------------------------
secp_path <- file.path(OUTDIR, "secreted_protein_atlas_columns.tsv")
if (file.exists(secp_path)) {
  secp_dt <- fread(secp_path)
  secp_new <- setdiff(names(secp_dt), "human_symbol")
  for (col in intersect(secp_new, names(atlas))) atlas[, (col) := NULL]
  atlas <- merge(atlas, secp_dt, by = "human_symbol", all.x = TRUE)
  for (lc in intersect(c("secreted_to_blood", "plasma_measured_change",
                         "secreted_plasma_concordant"), names(atlas)))
    atlas[is.na(get(lc)), (lc) := FALSE]
  cat(sprintf("  + Secreted-protein axis: %d measured-secreted biomarkers\n",
              sum(atlas$secreted_biomarker_tier == "measured_secreted_biomarker", na.rm = TRUE)))
} else {
  cat("  ! secreted_protein_atlas_columns.tsv not found; skipping (run 79_secreted_protein_axis.R first)\n")
}

# ---------------------------------------------------------------------------
# T1.5 Metabolite/lipid-QTL COLOC (2026-06-01; 60b/61b)
# Source: RNA-seq/results/multi_evidence/mqtl_atlas_columns.tsv
# coloc.abf of MASLD-disease GWAS (90 non-UKBB credible-set loci) vs Chen-2023
# plasma metabolites + Ottensmann lipid species. Supporting metabolic-causal
# layer (pathway-mediated, NOT primary S3). mqtl_* columns.
# ---------------------------------------------------------------------------
mqtl_path <- file.path(OUTDIR, "mqtl_atlas_columns.tsv")
if (file.exists(mqtl_path)) {
  mqtl_dt <- fread(mqtl_path)
  mqtl_new <- setdiff(names(mqtl_dt), "human_symbol")
  for (col in intersect(mqtl_new, names(atlas))) atlas[, (col) := NULL]
  atlas <- merge(atlas, mqtl_dt, by = "human_symbol", all.x = TRUE)
  if ("mqtl_known_masld_biomarker" %in% names(atlas))
    atlas[is.na(mqtl_known_masld_biomarker), mqtl_known_masld_biomarker := FALSE]
  if ("mqtl_n_colocalizing_metabolites" %in% names(atlas))
    atlas[is.na(mqtl_n_colocalizing_metabolites), mqtl_n_colocalizing_metabolites := 0L]
  cat(sprintf("  + Metabolite/lipid COLOC: %d genes PP4>0.5 (%d PP4>0.8)\n",
              sum(atlas$mqtl_best_pp4 > 0.5, na.rm = TRUE),
              sum(atlas$mqtl_best_pp4 > 0.8, na.rm = TRUE)))
} else {
  cat("  ! mqtl_atlas_columns.tsv not found; skipping (run 61b_combine_metabolite_coloc.R first)\n")
}

# Re-save atlas with new columns
fwrite(atlas, file.path(OUTDIR, "multi_evidence_atlas.csv"))
cat(sprintf("  Atlas re-saved: %d genes x %d columns\n", nrow(atlas), ncol(atlas)))

cat("\nDone. Atlas assembly complete.\n")
