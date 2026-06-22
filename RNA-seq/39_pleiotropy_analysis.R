#!/usr/bin/env Rscript
# 39_pleiotropy_analysis.R
# ---------------------------------------------------------------------------
# Pleiotropy Classification + Enrichment + FinnGen Cross-Phenotype COLOC
#
# Classifies Zenodo eGenes by pleiotropic architecture across 25 GWAS traits,
# tests enrichments against DEGs/drug targets, and optionally runs FinnGen
# cross-phenotype COLOC for hepatocyte eGenes.
#
# Inputs:
#   - zenodo_coloc_integrated.csv (from Script 37)
#   - Raw Zenodo COLOC files (for per-trait Boolean matrix)
#   - Dream results, concordance atlas, drug targets
#   - (Optional) FinnGen GWAS files for cross-phenotype COLOC
#
# Outputs (in RNA-seq/results/causal_inference/pleiotropy/):
#   - pleiotropy_classification.csv: per gene-celltype
#   - pleiotropy_gene_summary.csv: gene-level summary
#   - pleiotropy_deg_enrichment.csv: Fisher's exact tests
#   - pleiotropy_class_characteristics.csv: per-class summary
#   - pleiotropy_nafld_specific_targets.csv: ranked therapeutic targets
#   - pleiotropy_summary.csv: integrated summary
#   - pleiotropy_heatmap.pdf, trait_cooccurrence_heatmap.pdf,
#     pleiotropy_class_barplot.pdf
#   - (Optional) finngen_coloc_results.csv, finngen_coloc_heatmap.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
})

cat("=== Script 39: Pleiotropy Analysis ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR    <- Sys.getenv("MASLD_PROJECT_ROOT",
                          unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ZENODO_DIR  <- file.path(BASE_DIR, "RNA-seq/References/eQTL/MASLD_sc_eQTL_2025/Zenodo")
COLOC_DIR   <- file.path(ZENODO_DIR, "coloc")
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/pleiotropy")
SCEQTL_DIR  <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/sceqtl")
FINNGEN_DIR <- file.path(BASE_DIR, "GWAS/MR_Data/FinnGen")

# Input files
ZENODO_INTEGRATED <- file.path(SCEQTL_DIR, "zenodo_coloc_integrated.csv")
DREAM_FILE <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
GENE_CACHE <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")
CONCORDANCE_FILE <- file.path(BASE_DIR,
  "Analysis/Cross_Species_Concordance/results/concordance_atlas_unified.csv")
DRUG_FILE <- file.path(BASE_DIR, "RNA-seq/results/drug_repurposing/convergent_drug_targets.csv")
ATLAS_FILE <- file.path(BASE_DIR, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

# Hepatocyte eQTL file (for FinnGen COLOC)
HEP_EQTL_FILE <- file.path(BASE_DIR,
  "RNA-seq/References/eQTL/MASLD_sc_eQTL_2025/hepatocyte_twas_formatted.txt.gz")

# Cell types
CELL_TYPES <- c("hepatocyte", "endothelial_cell", "cholangiocyte", "stellate_cell")
NAFLD_COL  <- "Nat22_cALT_NAFLD_GWAS"
PADJ_THR   <- 0.1

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Trait domain classification
trait_domains <- list(
  nafld        = "Nat22_cALT_NAFLD_GWAS",
  liver_enzyme = c("GCST90019492_UKBB_ALT", "GCST90019497_UKBB_AST",
                   "GCST90019507_UKBB_GGT", "GCST90019494_UKBB_ALP"),
  lipid        = c("GCST90019512_UKBB_LDL", "GCST90019510_UKBB_HDL",
                   "GCST90019523_UKBB_Triglyceride", "GCST90019501_UKBB_TotalCholesterol",
                   "GCST90019496_UKBB_ApoB", "GCST90019495_UKBB_ApoA1"),
  metabolic    = c("GCST90019509_UKBB_HbA1C", "GCST90019506_UKBB_eGFR",
                   "GCST90019499_UKBB_CRP", "GCST90019524_UKBB_Urate",
                   "GCST90019504_UKBB_CystatinC"),
  hepatobiliary = c("GCST90019521_UKBB_TotalBilirubin", "GCST90019500_UKBB_Ca",
                    "GCST90019516_UKBB_Phosphate"),
  other        = c("GCST90019519_UKBB_UrineSodium", "GCST90019517_UKBB_UrineK",
                   "GCST90019511_UKBB_IGF1", "GCST90019526_UKBB_VitD")
)

# Non-NAFLD domain names (for domain counting)
non_nafld_domains <- setdiff(names(trait_domains), "nafld")

# ==============================================================================
# Part A: Load Data and Compute Trait Domain Counts
# ==============================================================================
cat("--- Part A: Trait Domain Classification + Pleiotropy Assignment ---\n\n")

# Load zenodo_coloc_integrated.csv
cat("Loading zenodo_coloc_integrated.csv...\n")
stopifnot(file.exists(ZENODO_INTEGRATED))
zenodo_int <- fread(ZENODO_INTEGRATED)
cat(sprintf("  Integrated data: %d rows, %d unique genes\n",
            nrow(zenodo_int), uniqueN(zenodo_int$gene)))

# Re-load raw Zenodo COLOC files to get per-trait Boolean matrix
# (zenodo_int has trait_list as text, but we need individual trait columns)
cat("Re-loading raw Zenodo COLOC files for per-trait breakdown...\n")
coloc_list <- list()
for (ct in CELL_TYPES) {
  f <- file.path(COLOC_DIR, sprintf("eGenes_coloc_list_%s.txt.gz", ct))
  if (!file.exists(f)) {
    cat(sprintf("  WARNING: Missing COLOC file for %s\n", ct))
    next
  }
  dt <- fread(f, header = TRUE, sep = "\t")
  dt[, cell_type := ct]
  coloc_list[[ct]] <- dt
}
coloc_raw <- rbindlist(coloc_list, use.names = TRUE, fill = TRUE)
cat(sprintf("  Raw COLOC: %d gene-celltype records\n", nrow(coloc_raw)))

# Identify trait columns (excluding gene, cell_type, is_liver_eGene_*)
egene_cols <- grep("^is_liver_eGene_", names(coloc_raw), value = TRUE)
meta_cols  <- c("gene", "cell_type", egene_cols)
all_trait_cols <- setdiff(names(coloc_raw), meta_cols)

# Verify all expected traits are present
expected_traits <- unique(unlist(trait_domains))
missing_traits <- setdiff(expected_traits, all_trait_cols)
if (length(missing_traits) > 0) {
  cat(sprintf("  WARNING: Missing traits: %s\n", paste(missing_traits, collapse = ", ")))
}

# Compute per-domain counts for each gene-celltype row
for (domain in names(trait_domains)) {
  domain_traits <- intersect(trait_domains[[domain]], all_trait_cols)
  col_name <- paste0("n_", domain)
  if (length(domain_traits) > 0) {
    coloc_raw[, (col_name) := rowSums(.SD == TRUE, na.rm = TRUE),
              .SDcols = domain_traits]
  } else {
    coloc_raw[, (col_name) := 0L]
  }
}

# Total colocalizations (including NAFLD)
coloc_raw[, n_total_coloc := rowSums(.SD == TRUE, na.rm = TRUE),
          .SDcols = all_trait_cols]

# NAFLD colocalization flag
coloc_raw[, nafld_coloc := get(NAFLD_COL) == TRUE]

# Count non-NAFLD domains with >= 1 TRUE trait
coloc_raw[, n_domains_coloc := 0L]
for (domain in non_nafld_domains) {
  col_name <- paste0("n_", domain)
  coloc_raw[, n_domains_coloc := n_domains_coloc + as.integer(get(col_name) >= 1L)]
}

# Domain profile: semicolon-separated list of domains with >= 1 coloc
coloc_raw[, domain_profile := {
  profiles <- character(.N)
  for (i in seq_len(.N)) {
    active_domains <- character(0)
    for (domain in non_nafld_domains) {
      if (get(paste0("n_", domain))[i] >= 1L) {
        active_domains <- c(active_domains, domain)
      }
    }
    profiles[i] <- if (length(active_domains) > 0) {
      paste(active_domains, collapse = ";")
    } else {
      NA_character_
    }
  }
  profiles
}]

# Pleiotropy classification (priority order)
classify_pleiotropy <- function(n_total, n_domains, n_liver, n_lipid, n_metab,
                                 nafld_col, n_total_coloc_val) {
  if (n_total >= 5 || n_domains >= 4) return("Broadly_pleiotropic")
  if (n_domains >= 3) return("Multi_domain")
  if (n_liver >= 2 && n_domains <= 2) return("Liver_enzyme_pleiotropic")
  if (n_lipid >= 2 && n_domains <= 2) return("Lipid_pleiotropic")
  if (n_metab >= 2 && n_domains <= 2) return("Metabolic_pleiotropic")
  if (nafld_col == TRUE && n_total_coloc_val <= 2) return("NAFLD_specific")
  if (n_total_coloc_val == 1 && nafld_col == FALSE) return("Single_trait")
  if (n_total_coloc_val == 0) return("Not_colocalizing")
  return("Other_pleiotropic")
}

coloc_raw[, pleiotropy_class := mapply(
  classify_pleiotropy,
  n_total_coloc, n_domains_coloc, n_liver_enzyme, n_lipid, n_metabolic,
  nafld_coloc, n_total_coloc
)]

cat("\nPleiotropy class distribution (gene-celltype rows):\n")
print(coloc_raw[, .N, by = pleiotropy_class][order(-N)])

# Merge dream DEG info from the integrated file
coloc_raw <- merge(
  coloc_raw,
  zenodo_int[, .(gene, cell_type, is_deg, bulk_logFC, bulk_padj,
                 consensus_tier, layers_active, in_atlas)],
  by = c("gene", "cell_type"), all.x = TRUE
)

# Build output: pleiotropy_classification.csv
domain_count_cols <- paste0("n_", names(trait_domains))
out_class <- coloc_raw[, c("gene", "cell_type", "nafld_coloc",
                            "n_total_coloc", domain_count_cols,
                            "n_domains_coloc", "domain_profile",
                            "pleiotropy_class",
                            "is_deg", "bulk_logFC", "bulk_padj",
                            "consensus_tier"), with = FALSE]

out_class_file <- file.path(RESULTS_DIR, "pleiotropy_classification.csv")
fwrite(out_class, out_class_file)
cat(sprintf("\nWrote: %s (%d rows)\n", out_class_file, nrow(out_class)))

# Gene-level summary (collapse across cell types)
cat("\nBuilding gene-level summary...\n")

# Priority ordering for best class
class_priority <- c("Broadly_pleiotropic", "Multi_domain",
                     "Liver_enzyme_pleiotropic", "Lipid_pleiotropic",
                     "Metabolic_pleiotropic", "Other_pleiotropic",
                     "NAFLD_specific", "Single_trait", "Not_colocalizing")

gene_summary <- coloc_raw[, {
  # Best class: highest priority across cell types
  classes <- pleiotropy_class
  best_idx <- min(match(classes, class_priority), na.rm = TRUE)
  best_class <- class_priority[best_idx]

  # Cell types with any colocalization
  ct_with_coloc <- unique(cell_type[n_total_coloc >= 1])

  # Union of domains across cell types
  all_doms <- unique(unlist(strsplit(
    domain_profile[!is.na(domain_profile)], ";")))
  if (length(all_doms) == 0) all_doms <- NA_character_

  list(
    best_pleiotropy_class = best_class,
    max_traits_coloc = max(n_total_coloc, na.rm = TRUE),
    n_cell_types_with_coloc = length(ct_with_coloc),
    domain_breadth = max(n_domains_coloc, na.rm = TRUE),
    all_cell_types = paste(ct_with_coloc, collapse = ";"),
    all_domains = paste(all_doms, collapse = ";"),
    is_nafld_specific = best_class == "NAFLD_specific",
    is_deg = any(is_deg == TRUE, na.rm = TRUE),
    bulk_logFC = bulk_logFC[1],
    bulk_padj = bulk_padj[1]
  )
}, by = gene]

# Rename for atlas consistency
setnames(gene_summary, "gene", "human_symbol")

out_summary_file <- file.path(RESULTS_DIR, "pleiotropy_gene_summary.csv")
fwrite(gene_summary, out_summary_file)
cat(sprintf("Wrote: %s (%d unique genes)\n", out_summary_file, nrow(gene_summary)))

cat("\nGene-level pleiotropy class distribution:\n")
print(gene_summary[, .N, by = best_pleiotropy_class][order(-N)])

# Sanity checks
cat("\nSanity checks:\n")
for (g in c("PCCB", "HSD17B13", "EFHD1", "ETS2")) {
  row <- gene_summary[human_symbol == g]
  if (nrow(row) > 0) {
    cat(sprintf("  %s: class=%s, max_traits=%d, domains=%d, is_deg=%s\n",
                g, row$best_pleiotropy_class, row$max_traits_coloc,
                row$domain_breadth, row$is_deg))
  }
}

# ==============================================================================
# Part B: Visualizations
# ==============================================================================
cat("\n--- Part B: Visualizations ---\n\n")

# B1: Gene x trait Boolean heatmap (genes with >= 2 colocalizations)
cat("Generating gene x trait heatmap...\n")
heatmap_genes <- coloc_raw[n_total_coloc >= 2, unique(gene)]
cat(sprintf("  Genes with >= 2 colocalizations: %d\n", length(heatmap_genes)))

if (length(heatmap_genes) > 0) {
  # Build Boolean matrix for heatmap
  hm_data <- coloc_raw[gene %in% heatmap_genes]

  # For multi-cell-type genes, use the cell type with most colocs
  hm_best <- hm_data[, .SD[which.max(n_total_coloc)], by = gene]

  # Create matrix: genes x traits
  trait_mat <- as.matrix(hm_best[, ..all_trait_cols])
  rownames(trait_mat) <- hm_best$gene
  # Convert TRUE/FALSE to 1/0
  trait_mat <- ifelse(trait_mat == TRUE, 1, 0)

  # Column annotation: trait domains
  trait_domain_map <- character(length(all_trait_cols))
  names(trait_domain_map) <- all_trait_cols
  for (domain in names(trait_domains)) {
    for (tr in trait_domains[[domain]]) {
      if (tr %in% names(trait_domain_map)) {
        trait_domain_map[tr] <- domain
      }
    }
  }

  col_annotation <- data.frame(
    Domain = trait_domain_map[all_trait_cols],
    row.names = all_trait_cols
  )

  # Row annotation: pleiotropy class + cell type
  row_annotation <- data.frame(
    Class = hm_best$pleiotropy_class,
    Cell_Type = hm_best$cell_type,
    row.names = hm_best$gene
  )

  # Domain colors
  domain_colors <- list(
    Domain = c(nafld = "#E41A1C", liver_enzyme = "#FF7F00",
               lipid = "#377EB8", metabolic = "#4DAF4A",
               hepatobiliary = "#984EA3", other = "#999999"),
    Class = c(Broadly_pleiotropic = "#E41A1C", Multi_domain = "#FF7F00",
              Liver_enzyme_pleiotropic = "#FF7F00",
              Lipid_pleiotropic = "#377EB8",
              Metabolic_pleiotropic = "#4DAF4A",
              Other_pleiotropic = "#FFFF33",
              NAFLD_specific = "#984EA3",
              Single_trait = "#A65628",
              Not_colocalizing = "#F0F0F0")
  )

  # Simplify column names for display
  col_labels <- gsub("GCST90019\\d+_UKBB_", "", all_trait_cols)
  col_labels <- gsub("Nat22_cALT_NAFLD_GWAS", "NAFLD", col_labels)
  colnames(trait_mat) <- col_labels

  # Update annotation rownames to match
  rownames(col_annotation) <- col_labels

  # Order columns by domain
  domain_order <- c("nafld", "liver_enzyme", "lipid", "metabolic",
                    "hepatobiliary", "other")
  col_order <- order(match(col_annotation$Domain, domain_order))
  trait_mat <- trait_mat[, col_order]
  col_annotation <- col_annotation[col_order, , drop = FALSE]

  # Limit to top 80 genes by n_total_coloc for readability
  if (nrow(trait_mat) > 80) {
    top_genes <- hm_best[order(-n_total_coloc)][1:80, gene]
    trait_mat <- trait_mat[top_genes, ]
    row_annotation <- row_annotation[top_genes, , drop = FALSE]
  }

  tryCatch({
    pdf(file.path(RESULTS_DIR, "pleiotropy_heatmap.pdf"),
        width = 14, height = max(10, nrow(trait_mat) * 0.18))
    pheatmap(trait_mat,
             color = c("white", "#2166AC"),
             cluster_cols = FALSE,
             clustering_method = "ward.D2",
             annotation_col = col_annotation,
             annotation_row = row_annotation,
             annotation_colors = domain_colors,
             show_colnames = TRUE,
             fontsize_row = 6,
             fontsize_col = 8,
             main = "COLOC Pleiotropy: Gene x GWAS Trait (Boolean)",
             legend = FALSE)
    dev.off()
    cat("  Wrote: pleiotropy_heatmap.pdf\n")
  }, error = function(e) {
    cat(sprintf("  WARNING: Heatmap generation failed: %s\n", e$message))
  })
}

# B2: Trait co-occurrence matrix (hepatocyte)
cat("Generating trait co-occurrence heatmap...\n")
hep_raw <- coloc_raw[cell_type == "hepatocyte"]
if (nrow(hep_raw) > 0) {
  hep_trait_mat <- as.matrix(hep_raw[, ..all_trait_cols])
  hep_trait_mat <- ifelse(hep_trait_mat == TRUE, 1, 0)

  # Co-occurrence: t(mat) %*% mat gives count of genes sharing each trait pair
  cooccur <- t(hep_trait_mat) %*% hep_trait_mat

  # Simplify names
  simple_names <- gsub("GCST90019\\d+_UKBB_", "", rownames(cooccur))
  simple_names <- gsub("Nat22_cALT_NAFLD_GWAS", "NAFLD", simple_names)
  rownames(cooccur) <- colnames(cooccur) <- simple_names

  tryCatch({
    pdf(file.path(RESULTS_DIR, "trait_cooccurrence_heatmap.pdf"),
        width = 12, height = 10)
    pheatmap(cooccur,
             color = colorRampPalette(c("white", "#F7FBFF", "#6BAED6", "#08306B"))(50),
             display_numbers = TRUE,
             number_format = "%d",
             fontsize_number = 7,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             main = "Trait Co-occurrence (Hepatocyte eGenes)")
    dev.off()
    cat("  Wrote: trait_cooccurrence_heatmap.pdf\n")
  }, error = function(e) {
    cat(sprintf("  WARNING: Co-occurrence heatmap failed: %s\n", e$message))
  })
}

# B3: Pleiotropy class barplot by cell type
cat("Generating pleiotropy class barplot...\n")
class_ct <- coloc_raw[, .N, by = .(cell_type, pleiotropy_class)]

tryCatch({
  pdf(file.path(RESULTS_DIR, "pleiotropy_class_barplot.pdf"),
      width = 10, height = 6)

  # Reshape for stacked barplot
  class_wide <- dcast(class_ct, cell_type ~ pleiotropy_class, value.var = "N",
                      fill = 0)
  ct_names <- class_wide$cell_type
  class_wide[, cell_type := NULL]
  bar_mat <- as.matrix(class_wide)
  rownames(bar_mat) <- ct_names

  class_colors <- c(
    Broadly_pleiotropic = "#E41A1C",
    Multi_domain = "#FF7F00",
    Liver_enzyme_pleiotropic = "#FF9933",
    Lipid_pleiotropic = "#377EB8",
    Metabolic_pleiotropic = "#4DAF4A",
    Other_pleiotropic = "#FFFF33",
    NAFLD_specific = "#984EA3",
    Single_trait = "#A65628",
    Not_colocalizing = "#CCCCCC"
  )
  cols_present <- intersect(colnames(bar_mat), names(class_colors))
  bar_mat <- bar_mat[, cols_present, drop = FALSE]

  barplot(t(bar_mat),
          beside = FALSE,
          col = class_colors[cols_present],
          main = "Pleiotropy Classification by Cell Type",
          ylab = "Number of eGenes",
          las = 2,
          cex.names = 0.8)
  legend("topright", legend = cols_present, fill = class_colors[cols_present],
         cex = 0.6, ncol = 2)

  dev.off()
  cat("  Wrote: pleiotropy_class_barplot.pdf\n")
}, error = function(e) {
  cat(sprintf("  WARNING: Barplot generation failed: %s\n", e$message))
})

# ==============================================================================
# Part C: DEG-Pleiotropy Enrichment
# ==============================================================================
cat("\n--- Part C: DEG-Pleiotropy Enrichment ---\n\n")

# Load Conserved genes
conserved <- character(0)
if (file.exists(CONCORDANCE_FILE)) {
  conc <- fread(CONCORDANCE_FILE)
  conserved <- conc[primary_category == "Conserved", human_symbol]
  cat(sprintf("  Conserved genes: %d\n", length(conserved)))
}

# Load drug targets
drug_symbols <- character(0)
if (file.exists(DRUG_FILE)) {
  drugs <- fread(DRUG_FILE)
  drug_symbols <- unique(drugs$symbol)
  cat(sprintf("  Drug target genes: %d\n", length(drug_symbols)))
}

# Load atlas for layers_active
if (file.exists(ATLAS_FILE)) {
  atlas_layers <- fread(ATLAS_FILE, select = c("human_symbol", "layers_active"))
  setkey(atlas_layers, human_symbol)
} else {
  atlas_layers <- data.table(human_symbol = character(0), layers_active = integer(0))
  setkey(atlas_layers, human_symbol)
}

# Universe: all Zenodo eGenes (gene-level)
universe_genes <- unique(gene_summary$human_symbol)
n_universe <- length(universe_genes)
cat(sprintf("  Enrichment universe (Zenodo eGenes): %d genes\n", n_universe))

# Annotations on gene_summary
gene_summary[, is_conserved := human_symbol %in% conserved]
gene_summary[, is_drug_target := human_symbol %in% drug_symbols]
gene_summary[, layers_active_atlas := atlas_layers[human_symbol, layers_active]]

# --- Enrichment Test 1: Broadly pleiotropic enriched among DEGs ---
cat("\nTest 1: Broadly pleiotropic among DEGs vs non-DEGs\n")
broad_genes <- gene_summary[best_pleiotropy_class == "Broadly_pleiotropic", human_symbol]
a <- sum(gene_summary[human_symbol %in% broad_genes, is_deg == TRUE])
b <- sum(gene_summary[human_symbol %in% broad_genes, is_deg == FALSE])
c_val <- sum(gene_summary[!human_symbol %in% broad_genes, is_deg == TRUE])
d <- sum(gene_summary[!human_symbol %in% broad_genes, is_deg == FALSE])
ft1 <- fisher.test(matrix(c(a, b, c_val, d), nrow = 2))
cat(sprintf("  Broadly_pleiotropic DEGs: %d/%d, non-Broad DEGs: %d/%d\n",
            a, a + b, c_val, c_val + d))
cat(sprintf("  Fisher OR = %.2f, p = %.3e\n", ft1$estimate, ft1$p.value))

# --- Enrichment Test 2: NAFLD-specific among Conserved ---
cat("\nTest 2: NAFLD-specific among Conserved\n")
nafld_spec_genes <- gene_summary[best_pleiotropy_class == "NAFLD_specific", human_symbol]
# Universe: NAFLD-colocalizing eGenes
nafld_coloc_genes <- gene_summary[grepl("nafld|NAFLD", best_pleiotropy_class) |
                                   max_traits_coloc >= 1, human_symbol]
# More precisely: genes with nafld_coloc TRUE in any cell type
nafld_coloc_genes <- unique(coloc_raw[nafld_coloc == TRUE, gene])
n_nafld_universe <- length(nafld_coloc_genes)

a2 <- sum(nafld_spec_genes %in% conserved)
b2 <- length(nafld_spec_genes) - a2
c2 <- sum(setdiff(nafld_coloc_genes, nafld_spec_genes) %in% conserved)
d2 <- length(setdiff(nafld_coloc_genes, nafld_spec_genes)) - c2
if (n_nafld_universe > 0) {
  ft2 <- fisher.test(matrix(c(a2, b2, c2, d2), nrow = 2))
  cat(sprintf("  NAFLD_specific in Core: %d/%d, other NAFLD-coloc in Core: %d/%d\n",
              a2, a2 + b2, c2, c2 + d2))
  cat(sprintf("  Fisher OR = %.2f, p = %.3e\n", ft2$estimate, ft2$p.value))
} else {
  ft2 <- list(estimate = NA, p.value = NA)
  cat("  No NAFLD-colocalizing genes found\n")
}

# --- Enrichment Test 3: Drug targets among pleiotropic genes ---
cat("\nTest 3: Drug targets among pleiotropic vs non-pleiotropic\n")
pleiotropic_genes <- gene_summary[
  best_pleiotropy_class %in% c("Broadly_pleiotropic", "Multi_domain",
                                "Liver_enzyme_pleiotropic", "Lipid_pleiotropic",
                                "Metabolic_pleiotropic", "Other_pleiotropic"),
  human_symbol]
# Restrict universe to genes with drug target annotation possible
drug_universe <- gene_summary[human_symbol %in% drug_symbols | is_drug_target == TRUE, human_symbol]
# Actually use full Zenodo eGene universe
a3 <- sum(pleiotropic_genes %in% drug_symbols)
b3 <- length(pleiotropic_genes) - a3
c3 <- sum(setdiff(universe_genes, pleiotropic_genes) %in% drug_symbols)
d3 <- length(setdiff(universe_genes, pleiotropic_genes)) - c3
ft3 <- fisher.test(matrix(c(a3, b3, c3, d3), nrow = 2))
cat(sprintf("  Pleiotropic drug targets: %d/%d, non-pleiotropic drug targets: %d/%d\n",
            a3, a3 + b3, c3, c3 + d3))
cat(sprintf("  Fisher OR = %.2f, p = %.3e\n", ft3$estimate, ft3$p.value))

# Save enrichment results
enrichment_dt <- data.table(
  test = c("Broadly_pleiotropic_among_DEGs",
           "NAFLD_specific_in_Conserved",
           "Drug_targets_in_pleiotropic"),
  group_a = c(a, a2, a3),
  group_b = c(b, b2, b3),
  group_c = c(c_val, c2, c3),
  group_d = c(d, d2, d3),
  odds_ratio = c(ft1$estimate, ft2$estimate, ft3$estimate),
  p_value = c(ft1$p.value, ft2$p.value, ft3$p.value),
  universe_size = c(n_universe, n_nafld_universe, n_universe)
)
fwrite(enrichment_dt, file.path(RESULTS_DIR, "pleiotropy_deg_enrichment.csv"))
cat("Wrote: pleiotropy_deg_enrichment.csv\n")

# Per-class characteristics table
cat("\nComputing per-class characteristics...\n")
class_chars <- gene_summary[, .(
  n_genes = .N,
  median_bulk_logFC = median(abs(bulk_logFC), na.rm = TRUE),
  frac_degs = mean(is_deg, na.rm = TRUE),
  frac_conserved = mean(is_conserved, na.rm = TRUE),
  frac_drug_targets = mean(is_drug_target, na.rm = TRUE),
  median_layers_active = as.double(median(layers_active_atlas, na.rm = TRUE)),
  median_traits = as.double(median(max_traits_coloc, na.rm = TRUE)),
  median_domains = as.double(median(domain_breadth, na.rm = TRUE))
), by = best_pleiotropy_class]

setorder(class_chars, -n_genes)
fwrite(class_chars, file.path(RESULTS_DIR, "pleiotropy_class_characteristics.csv"))
cat("Wrote: pleiotropy_class_characteristics.csv\n")
print(class_chars)

# ==============================================================================
# Part D: FinnGen Cross-Phenotype COLOC (Tier 2, conditional)
# ==============================================================================
cat("\n--- Part D: FinnGen Cross-Phenotype COLOC ---\n\n")

finngen_traits <- list(
  NAFLD     = list(file = "finngen_R12_NAFLD.gz",
                   cases = 3504, controls = 496844),
  Cirrhosis = list(file = "finngen_R12_CHIRHEP_NAS.gz",
                   cases = 1703, controls = 494803),
  HCC       = list(file = "finngen_R12_C3_HEPATOCELLU_CARC_EXALLC.gz",
                   cases = 947, controls = 378749),
  Obesity   = list(file = "finngen_R12_E4_OBESITY.gz",
                   cases = 31499, controls = 468693)
)

# Check if FinnGen files exist
finngen_files_exist <- all(sapply(finngen_traits, function(x) {
  file.exists(file.path(FINNGEN_DIR, x$file))
}))

if (!finngen_files_exist) {
  cat("  FinnGen GWAS files not found. Skipping Part D.\n")
  cat("  To download, run:\n")
  cat("    bash GWAS/download_finngen_gwas.sh\n")
} else if (!file.exists(HEP_EQTL_FILE)) {
  cat("  Hepatocyte eQTL file not found. Skipping Part D.\n")
} else {
  cat("  FinnGen files found. Running cross-phenotype COLOC for hepatocyte eGenes.\n\n")

  # Check coloc package
  if (!requireNamespace("coloc", quietly = TRUE)) {
    cat("  WARNING: coloc package not available. Skipping Part D.\n")
  } else {
    library(coloc)

    # Load hepatocyte eQTL data
    cat("  Loading hepatocyte eQTL data...\n")
    hep_eqtl <- fread(HEP_EQTL_FILE)
    cat(sprintf("    eQTL rows: %d, eGenes: %d\n",
                nrow(hep_eqtl), uniqueN(hep_eqtl$gene)))

    # Determine eQTL columns
    cat("    eQTL columns:", paste(names(hep_eqtl), collapse = ", "), "\n")

    # Standardize eQTL column names to: gene, chr, pos, beta, se, pval
    eqtl_cols <- names(hep_eqtl)
    # gene
    if ("gene_id" %in% eqtl_cols && !"gene" %in% eqtl_cols) setnames(hep_eqtl, "gene_id", "gene")
    # position
    if ("position" %in% eqtl_cols && !"pos" %in% eqtl_cols) setnames(hep_eqtl, "position", "pos")
    if ("variant_pos" %in% eqtl_cols && !"pos" %in% eqtl_cols) setnames(hep_eqtl, "variant_pos", "pos")
    # chromosome
    if ("chromosome" %in% eqtl_cols && !"chr" %in% eqtl_cols) setnames(hep_eqtl, "chromosome", "chr")
    if ("variant_chr" %in% eqtl_cols && !"chr" %in% eqtl_cols) setnames(hep_eqtl, "variant_chr", "chr")
    # beta/se/pval — use _all condition (all samples) as primary
    if ("beta_all" %in% eqtl_cols && !"beta" %in% eqtl_cols) setnames(hep_eqtl, "beta_all", "beta")
    if ("se_all" %in% eqtl_cols && !"se" %in% eqtl_cols) setnames(hep_eqtl, "se_all", "se")
    if ("pval_all" %in% eqtl_cols && !"pval" %in% eqtl_cols) setnames(hep_eqtl, "pval_all", "pval")
    if ("slope" %in% eqtl_cols && !"beta" %in% eqtl_cols) setnames(hep_eqtl, "slope", "beta")
    if ("slope_se" %in% eqtl_cols && !"se" %in% eqtl_cols) setnames(hep_eqtl, "slope_se", "se")
    if ("pval_nominal" %in% eqtl_cols && !"pval" %in% eqtl_cols) setnames(hep_eqtl, "pval_nominal", "pval")

    # Ensure chr is integer
    if (is.character(hep_eqtl$chr)) {
      hep_eqtl[, chr := as.integer(sub("chr", "", chr))]
    }

    # Get unique eGenes from hepatocyte
    hep_egenes <- unique(hep_eqtl$gene)
    cat(sprintf("    Hepatocyte eGenes to test: %d\n", length(hep_egenes)))

    # Pre-split eQTL by chromosome and key by pos for fast merge
    cat("  Pre-splitting eQTL by chromosome...\n")
    setkey(hep_eqtl, chr, pos)
    eqtl_by_chr <- split(hep_eqtl, by = "chr")
    # Gene-to-chromosome lookup
    gene_chr_map <- unique(hep_eqtl[, .(gene, chr)])

    # Process each FinnGen trait
    finngen_results <- list()

    for (trait_name in names(finngen_traits)) {
      trait_info <- finngen_traits[[trait_name]]
      gwas_file <- file.path(FINNGEN_DIR, trait_info$file)
      s_param <- trait_info$cases / (trait_info$cases + trait_info$controls)

      cat(sprintf("\n  Processing FinnGen %s (s=%.4f)...\n", trait_name, s_param))

      # Load FinnGen GWAS (select only needed columns for speed)
      gwas <- fread(gwas_file,
                    select = c("#chrom", "pos", "beta", "sebeta"))
      setnames(gwas, c("#chrom", "sebeta"), c("chrom", "se"))
      gwas[, chr := as.integer(sub("chr", "", chrom))]
      gwas[, chrom := NULL]

      cat(sprintf("    GWAS rows: %d, chr range: %s\n",
                  nrow(gwas), paste(range(gwas$chr, na.rm=TRUE), collapse="-")))

      # Pre-split GWAS by chromosome and key by pos
      setkey(gwas, chr, pos)
      gwas_by_chr <- split(gwas, by = "chr")

      cat(sprintf("    eQTL chromosomes: %d, GWAS chromosomes: %d\n",
                  length(eqtl_by_chr), length(gwas_by_chr)))

      # Run COLOC for each eGene
      trait_results <- list()
      n_tested <- 0
      n_skipped <- 0

      n_total_genes <- length(hep_egenes)
      for (i_eg in seq_along(hep_egenes)) {
        eg <- hep_egenes[i_eg]
        if (i_eg %% 1000 == 0) cat(sprintf("      Progress: %d/%d genes (tested=%d)\n", i_eg, n_total_genes, n_tested))
        eg_chr <- gene_chr_map[gene == eg, chr[1]]
        chr_key <- as.character(eg_chr)

        eqtl_chr_dt <- eqtl_by_chr[[chr_key]]
        if (is.null(eqtl_chr_dt)) { n_skipped <- n_skipped + 1; next }
        eqtl_sub <- eqtl_chr_dt[gene == eg]
        if (nrow(eqtl_sub) < 10) { n_skipped <- n_skipped + 1; next }

        gwas_chr_dt <- gwas_by_chr[[chr_key]]
        if (is.null(gwas_chr_dt)) { n_skipped <- n_skipped + 1; next }

        # Fast keyed merge on pos, deduplicate multiallelic sites
        merged <- merge(eqtl_sub, gwas_chr_dt, by = "pos",
                        suffixes = c("_eqtl", "_gwas"))
        merged <- merged[!duplicated(pos)]
        if (nrow(merged) < 10) { n_skipped <- n_skipped + 1; next }

        # Prepare coloc datasets
        eqtl_dataset <- list(
          beta     = merged$beta_eqtl,
          varbeta  = merged$se_eqtl^2,
          snp      = paste0("chr", eg_chr, ":", merged$pos),
          position = merged$pos,
          type     = "quant",
          sdY      = 1,
          N        = 312
        )

        gwas_dataset <- list(
          beta     = merged$beta_gwas,
          varbeta  = merged$se_gwas^2,
          snp      = paste0("chr", eg_chr, ":", merged$pos),
          position = merged$pos,
          type     = "cc",
          s        = s_param,
          N        = trait_info$cases + trait_info$controls
        )

        res <- tryCatch({
          suppressMessages(coloc.abf(eqtl_dataset, gwas_dataset,
                    p1 = 1e-4, p2 = 1e-4, p12 = 5e-6))
        }, error = function(e) NULL)

        if (!is.null(res)) {
          trait_results[[eg]] <- data.table(
            gene = eg,
            finngen_trait = trait_name,
            n_snps = nrow(merged),
            PP.H0 = res$summary["PP.H0.abf"],
            PP.H1 = res$summary["PP.H1.abf"],
            PP.H2 = res$summary["PP.H2.abf"],
            PP.H3 = res$summary["PP.H3.abf"],
            PP.H4 = res$summary["PP.H4.abf"]
          )
          n_tested <- n_tested + 1
        } else {
          n_skipped <- n_skipped + 1
        }
      }

      cat(sprintf("    Tested: %d, Skipped: %d (of %d eGenes)\n",
                  n_tested, n_skipped, length(hep_egenes)))

      if (length(trait_results) > 0) {
        trait_dt <- rbindlist(trait_results)
        cat(sprintf("    PP.H4 > 0.5: %d genes\n",
                    sum(trait_dt$PP.H4 > 0.5, na.rm = TRUE)))
        cat(sprintf("    PP.H4 > 0.3 (suggestive): %d genes\n",
                    sum(trait_dt$PP.H4 > 0.3, na.rm = TRUE)))
        finngen_results[[trait_name]] <- trait_dt
      }
    }

    # Combine all FinnGen results
    if (length(finngen_results) > 0) {
      finngen_all <- rbindlist(finngen_results)
      finngen_out <- file.path(RESULTS_DIR, "finngen_coloc_results.csv")
      fwrite(finngen_all, finngen_out)
      cat(sprintf("\nWrote: %s (%d rows)\n", finngen_out, nrow(finngen_all)))

      # Replication: genes with Zenodo NAFLD=TRUE and FinnGen NAFLD PP.H4
      zenodo_nafld_genes <- coloc_raw[cell_type == "hepatocyte" & nafld_coloc == TRUE, gene]
      finngen_nafld <- finngen_all[finngen_trait == "NAFLD"]
      replication <- finngen_nafld[gene %in% zenodo_nafld_genes]
      if (nrow(replication) > 0) {
        cat("\nFinnGen NAFLD replication of Zenodo hits:\n")
        cat(sprintf("  Zenodo NAFLD hep eGenes tested in FinnGen: %d\n",
                    nrow(replication)))
        cat(sprintf("  FinnGen PP.H4 > 0.3: %d\n",
                    sum(replication$PP.H4 > 0.3, na.rm = TRUE)))
        cat(sprintf("  FinnGen PP.H4 > 0.5: %d\n",
                    sum(replication$PP.H4 > 0.5, na.rm = TRUE)))
      }

      # FinnGen COLOC heatmap: gene x trait PP.H4
      cat("\nGenerating FinnGen COLOC heatmap...\n")
      finngen_wide <- dcast(finngen_all, gene ~ finngen_trait, value.var = "PP.H4")
      genes_with_signal <- finngen_all[PP.H4 > 0.3, unique(gene)]
      if (length(genes_with_signal) > 0) {
        hm_genes <- intersect(genes_with_signal, finngen_wide$gene)
        finngen_mat <- as.matrix(finngen_wide[gene %in% hm_genes, -"gene"])
        rownames(finngen_mat) <- finngen_wide[gene %in% hm_genes, gene]
        finngen_mat[is.na(finngen_mat)] <- 0

        tryCatch({
          pdf(file.path(RESULTS_DIR, "finngen_coloc_heatmap.pdf"),
              width = 8, height = max(6, nrow(finngen_mat) * 0.3))
          pheatmap(finngen_mat,
                   color = colorRampPalette(c("white", "#FEE0D2", "#FC9272",
                                              "#DE2D26"))(50),
                   cluster_cols = FALSE,
                   main = "FinnGen Cross-Phenotype COLOC (Hepatocyte eGenes)",
                   fontsize_row = 8)
          dev.off()
          cat("  Wrote: finngen_coloc_heatmap.pdf\n")
        }, error = function(e) {
          cat(sprintf("  WARNING: FinnGen heatmap failed: %s\n", e$message))
        })
      } else {
        cat("  No genes with PP.H4 > 0.3 — skipping FinnGen heatmap\n")
      }
    }
  }
}

# ==============================================================================
# Part F: Integrated Summary
# ==============================================================================
cat("\n--- Part F: Integrated Summary ---\n\n")

# Merge classification with enrichment context
summary_out <- gene_summary[, .(
  human_symbol,
  best_pleiotropy_class,
  max_traits_coloc,
  domain_breadth,
  n_cell_types_with_coloc,
  all_cell_types,
  all_domains,
  is_nafld_specific,
  is_deg,
  bulk_logFC,
  bulk_padj,
  is_conserved,
  is_drug_target,
  layers_active = layers_active_atlas
)]

fwrite(summary_out, file.path(RESULTS_DIR, "pleiotropy_summary.csv"))
cat("Wrote: pleiotropy_summary.csv\n")

# NAFLD-specific therapeutic targets
nafld_targets <- gene_summary[
  is_nafld_specific == TRUE & is_deg == TRUE & !is.na(bulk_padj) & bulk_padj < PADJ_THR
]
nafld_targets[, abs_logFC := abs(bulk_logFC)]
setorder(nafld_targets, -abs_logFC)
nafld_targets[, abs_logFC := NULL]

nafld_targets_out <- nafld_targets[, .(
  human_symbol,
  bulk_logFC,
  bulk_padj,
  max_traits_coloc,
  is_conserved,
  is_drug_target,
  all_cell_types,
  layers_active = layers_active_atlas
)]

fwrite(nafld_targets_out,
       file.path(RESULTS_DIR, "pleiotropy_nafld_specific_targets.csv"))
cat(sprintf("Wrote: pleiotropy_nafld_specific_targets.csv (%d targets)\n",
            nrow(nafld_targets_out)))

if (nrow(nafld_targets_out) > 0) {
  cat("\nTop NAFLD-specific therapeutic targets (DEGs, ranked by |LFC|):\n")
  print(head(nafld_targets_out, 15))
}

# ==============================================================================
# Final Summary
# ==============================================================================
cat("\n=== Script 39 Complete ===\n")
cat("End time:", format(Sys.time()), "\n")
cat(sprintf("\nOutputs in: %s\n", RESULTS_DIR))
cat("  pleiotropy_classification.csv\n")
cat("  pleiotropy_gene_summary.csv\n")
cat("  pleiotropy_deg_enrichment.csv\n")
cat("  pleiotropy_class_characteristics.csv\n")
cat("  pleiotropy_nafld_specific_targets.csv\n")
cat("  pleiotropy_summary.csv\n")
cat("  pleiotropy_heatmap.pdf\n")
cat("  trait_cooccurrence_heatmap.pdf\n")
cat("  pleiotropy_class_barplot.pdf\n")
if (finngen_files_exist && file.exists(HEP_EQTL_FILE)) {
  cat("  finngen_coloc_results.csv\n")
  cat("  finngen_coloc_heatmap.pdf\n")
}
