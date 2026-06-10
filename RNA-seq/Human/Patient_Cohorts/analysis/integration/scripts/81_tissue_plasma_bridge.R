#!/usr/bin/env Rscript
# 81_tissue_plasma_bridge.R
# Phase 2B: Tissue-to-Plasma Bridge Analysis
#
# Cross-references 16,333 dream DEGs with 1,460 Olink plasma proteins to
# identify which tissue-level transcriptomic changes are detectable in
# circulating plasma proteins. Builds the bridge between tissue-based
# staging classifiers and non-invasive plasma biomarkers.
#
# Inputs:
#   - Multi-evidence atlas (dream DEGs with gene symbols)
#   - Olink plasma NPX data (1,460 proteins x 218 subjects)
#   - GSE276114 liver metadata (fibrosis staging for 177 samples)
#   - Staging classifier panel genes
#
# Outputs:
#   - tissue_plasma_bridge.csv: Full bridge table (gene, tissue DE, plasma available)
#   - bridge_classification.csv: Gene classification (concordant/discordant/tissue-only/plasma-only)
#   - panel_plasma_availability.csv: Which classifier panel genes have plasma equivalents
#   - plasma_degs.csv: Differential plasma proteins between fibrosis stages
#   - tissue_plasma_correlation.csv: Tissue LFC vs plasma LFC for paired genes
#
# Usage: Rscript 81_tissue_plasma_bridge.R
# SLURM: cpu, 8 CPUs, 64GB RAM, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUTDIR <- file.path(INT, "results/staging_classifier")
OLINK_DIR <- file.path(BASE, "Analysis/Proteomics/data/olink_plasma")
PROT_DIR <- file.path(BASE, "Analysis/Proteomics/results")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 81: Tissue-to-Plasma Bridge ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ============================================================
# STEP 1: Load all data sources
# ============================================================
cat("=== STEP 1: Loading Data ===\n")

# 1a. Multi-evidence atlas (dream DEGs with gene symbols)
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
cat("Atlas:", nrow(atlas), "genes\n")

# Dream DEGs
dream_degs <- atlas[!is.na(bulk_padj) & bulk_padj < 0.1]
cat("Dream DEGs (padj<0.1):", nrow(dream_degs), "\n")

# 1b. Olink plasma data
olink_file <- file.path(OLINK_DIR, "olink.qc.finished.mendeley.data.txt")
stopifnot(file.exists(olink_file))
olink <- fread(olink_file, header = TRUE)
olink_proteins <- olink$Assay
olink_mat <- as.matrix(olink[, -1, with = FALSE])
rownames(olink_mat) <- olink_proteins
cat("Olink plasma:", nrow(olink_mat), "proteins x", ncol(olink_mat), "subjects\n")

# 1c. GSE276114 liver metadata (fibrosis stage + disease)
liver_meta_file <- file.path(PROT_DIR, "gse276114_disease_metadata.csv")
liver_meta <- NULL
if (file.exists(liver_meta_file)) {
  liver_meta <- fread(liver_meta_file)
  cat("Liver metadata:", nrow(liver_meta), "samples\n")
  cat("  Disease groups:", paste(liver_meta[, .N, by = disease_group][, paste0(disease_group, "=", N)], collapse = ", "), "\n")
  cat("  Diseases:", paste(liver_meta[, .N, by = disease][, paste0(disease, "=", N)], collapse = ", "), "\n")
}

# 1d. Staging classifier panel genes (if available)
panel_files <- list(
  panel_10 = file.path(OUTDIR, "minimal_panel_10.csv"),
  panel_15 = file.path(OUTDIR, "minimal_panel_15.csv"),
  panel_25 = file.path(OUTDIR, "minimal_panel_25.csv"),
  panel_50 = file.path(OUTDIR, "minimal_panel_50.csv")
)

panel_genes <- list()
for (pname in names(panel_files)) {
  if (file.exists(panel_files[[pname]])) {
    pg <- fread(panel_files[[pname]])
    # Map Ensembl IDs to gene symbols via atlas
    pg[, ensembl_base := sub("[.][0-9]+$", "", gene)]
    pg_mapped <- merge(pg, atlas[, .(ensembl_id, human_symbol)],
                        by.x = "ensembl_base", by.y = "ensembl_id", all.x = TRUE)
    panel_genes[[pname]] <- pg_mapped$human_symbol[!is.na(pg_mapped$human_symbol)]
    cat("  Panel", pname, ":", length(panel_genes[[pname]]), "genes with symbols\n")
  }
}

cat("\n")

# ============================================================
# STEP 2: Build the tissue-plasma bridge
# ============================================================
cat("=== STEP 2: Tissue-Plasma Bridge ===\n")

olink_set <- data.table(human_symbol = olink_proteins, in_plasma = TRUE)
dream_set <- atlas[, .(human_symbol, bulk_logFC, bulk_padj, bulk_tstat,
                        ensembl_id, gene_biotype)]
dream_set[, is_deg := !is.na(bulk_padj) & bulk_padj < 0.1]

# Merge: all atlas genes x Olink availability
bridge <- merge(dream_set, olink_set, by = "human_symbol", all.x = TRUE)
bridge[is.na(in_plasma), in_plasma := FALSE]

# Classify
bridge[, classification := fifelse(
  is_deg & in_plasma, "tissue_DEG_plasma_detectable",
  fifelse(is_deg & !in_plasma, "tissue_DEG_plasma_undetectable",
  fifelse(!is_deg & in_plasma, "not_DEG_but_plasma_detectable",
  "neither"))
)]

class_counts <- bridge[, .N, by = classification]
cat("Bridge classification:\n")
print(class_counts)
cat("\n")

# Key metric: what fraction of DEGs are plasma-detectable?
n_deg_plasma <- bridge[classification == "tissue_DEG_plasma_detectable", .N]
n_deg_total <- bridge[is_deg == TRUE, .N]
cat("DEGs detectable in plasma:", n_deg_plasma, "/", n_deg_total,
    "(", round(n_deg_plasma / n_deg_total * 100, 1), "%)\n")

# What fraction of Olink proteins are DEGs?
n_olink_deg <- bridge[in_plasma == TRUE & is_deg == TRUE, .N]
n_olink_total <- bridge[in_plasma == TRUE, .N]
cat("Olink proteins that are DEGs:", n_olink_deg, "/", n_olink_total,
    "(", round(n_olink_deg / n_olink_total * 100, 1), "%)\n\n")

# Save bridge
fwrite(bridge, file.path(OUTDIR, "tissue_plasma_bridge.csv"))
cat("Saved tissue_plasma_bridge.csv:", nrow(bridge), "rows\n")

# ============================================================
# STEP 3: Plasma differential analysis (fibrosis staging)
# ============================================================
cat("\n=== STEP 3: Plasma Differential Analysis ===\n")

# Map Olink subjects to liver metadata
# Paper: 218 discovery subjects, 177 have liver biopsies (samples 1-177?)
# Olink columns are "Subject 1", "Subject 2", etc.
# Liver metadata has sample_number 1-178 (some may be missing)
# Attempt mapping by number

if (!is.null(liver_meta)) {
  # Extract subject numbers from Olink column names
  subject_cols <- colnames(olink_mat)
  subject_nums <- as.integer(gsub("Subject ", "", subject_cols))

  # Map to liver metadata
  liver_meta[, subject_col := paste0("Subject ", sample_number)]
  matched <- liver_meta[subject_col %in% subject_cols]
  cat("Subjects matched to liver metadata:", nrow(matched), "of", nrow(liver_meta), "\n")

  if (nrow(matched) >= 30) {
    # Subset Olink to matched subjects
    olink_matched <- olink_mat[, matched$subject_col, drop = FALSE]

    # Binary: advanced fibrosis (F3+F4) vs early (F0-2)
    matched[, fib_binary := fifelse(disease_group %in% c("F3", "F4"), "advanced", "early")]
    cat("  Advanced fibrosis:", sum(matched$fib_binary == "advanced"),
        "  Early:", sum(matched$fib_binary == "early"), "\n")

    # Limma DE on plasma proteins
    design <- model.matrix(~ fib_binary + disease, data = matched)

    # Check for rank deficiency
    if (qr(design)$rank == ncol(design)) {
      fit <- lmFit(olink_matched, design)
      fit <- eBayes(fit)

      plasma_de <- topTable(fit, coef = "fib_binaryearly", number = Inf, sort.by = "none")
      plasma_de$protein <- rownames(plasma_de)
      plasma_de <- as.data.table(plasma_de)
      setnames(plasma_de, c("logFC", "adj.P.Val"), c("plasma_logFC", "plasma_padj"),
               skip_absent = TRUE)

      n_plasma_sig <- sum(plasma_de$plasma_padj < 0.05, na.rm = TRUE)
      cat("  Plasma DE proteins (padj<0.05):", n_plasma_sig, "\n")

      fwrite(plasma_de, file.path(OUTDIR, "plasma_degs.csv"))
      cat("  Saved plasma_degs.csv:", nrow(plasma_de), "rows\n")
    } else {
      cat("  WARNING: Design matrix rank-deficient, skipping DE\n")
      plasma_de <- NULL
    }
  } else {
    cat("  Too few matched subjects for DE analysis\n")
    plasma_de <- NULL
  }

  # Also do MASLD-only analysis
  masld_matched <- matched[disease == "MASLD"]
  cat("\n  MASLD-only subset:", nrow(masld_matched), "samples\n")
  if (nrow(masld_matched) >= 20) {
    olink_masld <- olink_mat[, masld_matched$subject_col, drop = FALSE]
    design_masld <- model.matrix(~ fib_binary, data = masld_matched)

    if (qr(design_masld)$rank == ncol(design_masld)) {
      fit_masld <- lmFit(olink_masld, design_masld)
      fit_masld <- eBayes(fit_masld)
      plasma_de_masld <- topTable(fit_masld, coef = "fib_binaryearly", number = Inf, sort.by = "none")
      plasma_de_masld$protein <- rownames(plasma_de_masld)
      plasma_de_masld <- as.data.table(plasma_de_masld)
      setnames(plasma_de_masld, c("logFC", "adj.P.Val"), c("plasma_logFC_masld", "plasma_padj_masld"),
               skip_absent = TRUE)

      n_masld_sig <- sum(plasma_de_masld$plasma_padj_masld < 0.05, na.rm = TRUE)
      cat("  MASLD-only plasma DE (padj<0.05):", n_masld_sig, "\n")

      fwrite(plasma_de_masld, file.path(OUTDIR, "plasma_degs_masld_only.csv"))
    }
  }
} else {
  cat("  No liver metadata available, skipping plasma DE\n")
  plasma_de <- NULL
}

# ============================================================
# STEP 4: Tissue LFC vs Plasma LFC correlation
# ============================================================
cat("\n=== STEP 4: Tissue-Plasma Concordance ===\n")

if (!is.null(plasma_de)) {
  # Merge tissue dream LFC with plasma LFC
  concordance <- merge(
    bridge[is_deg == TRUE & in_plasma == TRUE, .(human_symbol, bulk_logFC, bulk_padj)],
    plasma_de[, .(protein, plasma_logFC, plasma_padj)],
    by.x = "human_symbol", by.y = "protein"
  )

  if (nrow(concordance) > 10) {
    rho <- cor(concordance$bulk_logFC, concordance$plasma_logFC, method = "spearman",
               use = "pairwise.complete")
    r <- cor(concordance$bulk_logFC, concordance$plasma_logFC, method = "pearson",
             use = "pairwise.complete")

    # Direction concordance (among significant in both)
    both_sig <- concordance[bulk_padj < 0.1 & plasma_padj < 0.05]
    if (nrow(both_sig) > 0) {
      dir_conc <- mean(sign(both_sig$bulk_logFC) == sign(both_sig$plasma_logFC))
    } else {
      dir_conc <- NA
    }

    # Classify each gene
    concordance[, direction_match := sign(bulk_logFC) == sign(plasma_logFC)]
    concordance[, tissue_plasma_class := fifelse(
      plasma_padj < 0.05 & direction_match, "concordant",
      fifelse(plasma_padj < 0.05 & !direction_match, "discordant",
      "tissue_only_significant")
    )]

    cat("  Genes in both tissue & plasma:", nrow(concordance), "\n")
    cat("  Spearman rho (tissue vs plasma LFC):", round(rho, 3), "\n")
    cat("  Pearson r:", round(r, 3), "\n")
    cat("  Significant in BOTH (padj tissue<0.1 & plasma<0.05):", nrow(both_sig), "\n")
    if (!is.na(dir_conc)) {
      cat("  Direction concordance (among dual-sig):", round(dir_conc * 100, 1), "%\n")
    }
    cat("  Concordant (same direction, plasma sig):", sum(concordance$tissue_plasma_class == "concordant"), "\n")
    cat("  Discordant (opposite direction, plasma sig):", sum(concordance$tissue_plasma_class == "discordant"), "\n")

    fwrite(concordance, file.path(OUTDIR, "tissue_plasma_correlation.csv"))
    cat("  Saved tissue_plasma_correlation.csv\n")
  }
}

# ============================================================
# STEP 5: Panel gene plasma availability
# ============================================================
cat("\n=== STEP 5: Classifier Panel Plasma Availability ===\n")

panel_avail_list <- list()
for (pname in names(panel_genes)) {
  genes <- panel_genes[[pname]]
  if (length(genes) == 0) next

  in_olink <- genes %in% olink_proteins
  cat("  ", pname, ":", sum(in_olink), "/", length(genes), "genes in Olink",
      "(", round(sum(in_olink)/length(genes)*100, 1), "%)\n")

  if (sum(in_olink) > 0) {
    cat("    In Olink:", paste(genes[in_olink], collapse = ", "), "\n")
  }
  if (sum(!in_olink) > 0) {
    cat("    Missing:", paste(genes[!in_olink], collapse = ", "), "\n")
  }

  panel_avail_list[[pname]] <- data.table(
    panel = pname,
    gene = genes,
    in_olink = in_olink
  )
}

if (length(panel_avail_list) > 0) {
  panel_avail <- rbindlist(panel_avail_list)
  fwrite(panel_avail, file.path(OUTDIR, "panel_plasma_availability.csv"))
  cat("  Saved panel_plasma_availability.csv\n")
}

# ============================================================
# STEP 6: Top plasma-translatable DEGs
# ============================================================
cat("\n=== STEP 6: Top Plasma-Translatable DEGs ===\n")

# Rank bridge genes by tissue significance among plasma-detectable
translatable <- bridge[classification == "tissue_DEG_plasma_detectable"]
translatable <- translatable[order(bulk_padj)]

cat("Top 30 plasma-translatable DEGs (by tissue significance):\n")
print(translatable[1:30, .(human_symbol, bulk_logFC = round(bulk_logFC, 2),
                            bulk_padj = formatC(bulk_padj, format = "e", digits = 1),
                            gene_biotype)])

# Add plasma DE info if available
if (!is.null(plasma_de)) {
  translatable <- merge(translatable, plasma_de[, .(protein, plasma_logFC, plasma_padj)],
                         by.x = "human_symbol", by.y = "protein", all.x = TRUE)
  translatable[, dual_significant := !is.na(plasma_padj) & plasma_padj < 0.05 & bulk_padj < 0.1]
  n_dual <- sum(translatable$dual_significant, na.rm = TRUE)
  cat("\nDual-significant (tissue padj<0.1 AND plasma padj<0.05):", n_dual, "genes\n")
}

fwrite(translatable, file.path(OUTDIR, "plasma_translatable_degs.csv"))
cat("Saved plasma_translatable_degs.csv:", nrow(translatable), "genes\n")

# ============================================================
# Summary
# ============================================================
cat("\n=== SUMMARY ===\n")
cat("Total atlas genes:", nrow(atlas), "\n")
cat("Dream DEGs:", nrow(dream_degs), "\n")
cat("Olink plasma proteins:", length(olink_proteins), "\n")
cat("DEGs in plasma:", n_deg_plasma, "(", round(n_deg_plasma/n_deg_total*100, 1), "% of DEGs)\n")
cat("Olink proteins that are DEGs:", n_olink_deg, "(", round(n_olink_deg/n_olink_total*100, 1), "% of Olink)\n")
if (!is.null(plasma_de)) {
  cat("Plasma DE proteins (adv vs early fib, padj<0.05):", sum(plasma_de$plasma_padj < 0.05, na.rm = TRUE), "\n")
}

cat("\n=== 81_tissue_plasma_bridge.R completed:", as.character(Sys.time()), "===\n")
