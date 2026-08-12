#!/usr/bin/env Rscript
# differential_proteomics.R — Multi-contrast proteomics validation (v3)
#
# Each proteomics dataset runs multiple contrasts matched to specific
# dream comparators. Covariates (BMI, age) added where available.
#
#   Section B:  PXD052937 — (1) MASLD vs Normal, (2) MASH vs MASL
#   Section B2: PXD051911 — (1) MASLD vs No_MASLD, (2) MASH vs MASL, (3) High NAS vs Low NAS
#   Section C:  Matched-contrast concordance (each contrast to its dream comparator)
#   Section D:  Ranked enrichment (fgsea)
#   Section E:  Effect-size stratified detection
#   Section F:  Positive control + drug target validation
#   QC:         PCA + missing data summaries per dataset
#
# Outputs (all to Analysis/Proteomics/results/):
#   protein_differential_results_v3.csv
#   protein_transcript_concordance_v3.csv
#   protein_ranked_enrichment.csv
#   protein_effectsize_detection.csv
#   protein_validation_summary.csv
#   qc/pca_*.csv, qc/missingness_*.csv
#
# NOTE: GSE276114 was removed. GEO confirms it is bulk RNA-seq
# (Expression profiling by high throughput sequencing), not SomaScan proteomics.

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(fgsea)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PROTEO_DIR  <- file.path(BASE, "Analysis/Proteomics")
RESULTS_DIR <- file.path(PROTEO_DIR, "results")
QC_DIR      <- file.path(RESULTS_DIR, "qc")
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

DISEASE_SIG_DIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                              "results/disease_signatures")

cat("=== Differential Proteomics Analysis (v3 — multi-contrast) ===\n")
cat("Time:", as.character(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# Helper: symbol mapping via the MASLD Gene Catalog
# ---------------------------------------------------------------------------
atlas_f <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
gene_map <- NULL
if (file.exists(atlas_f)) {
  atlas_raw <- fread(atlas_f, select = c("ensembl_id", "human_symbol"))
  atlas_raw[, ensembl_clean := sub("\\..*", "", ensembl_id)]
  gene_map <- atlas_raw[!is.na(human_symbol) & human_symbol != "",
                        .(ensembl_clean, symbol = human_symbol)]
  gene_map <- unique(gene_map, by = "ensembl_clean")
  cat("Loaded gene map:", nrow(gene_map), "Ensembl->symbol mappings\n\n")
}

# ---------------------------------------------------------------------------
# Helper: run limma DE on expression matrix
# ---------------------------------------------------------------------------
run_limma_de <- function(expr, meta, group_col, contrast_str, dataset_name,
                         covariate_cols = NULL) {
  # Build design formula
  meta[, .group := get(group_col)]
  if (!is.null(covariate_cols) && length(covariate_cols) > 0) {
    form_str <- paste0("~ 0 + .group + ", paste(covariate_cols, collapse = " + "))
  } else {
    form_str <- "~ 0 + .group"
  }
  design <- model.matrix(as.formula(form_str), data = meta)
  # Clean column names: remove ".group" prefix
  colnames(design) <- gsub("^\\.group", "", colnames(design))

  fit <- lmFit(expr[, meta$sample_id, drop = FALSE], design)
  contrast_mat <- makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- contrasts.fit(fit, contrast_mat)
  fit2 <- eBayes(fit2, robust = TRUE)

  de <- data.table(topTable(fit2, number = Inf, sort.by = "none"))
  de[, gene := rownames(expr)]
  de[, dataset := dataset_name]

  setnames(de, "adj.P.Val", "padj", skip_absent = TRUE)
  setnames(de, "P.Value", "pvalue", skip_absent = TRUE)

  cat("    Significant (padj<0.05):", sum(de$padj < 0.05, na.rm = TRUE), "\n")
  cat("    Significant (padj<0.1):", sum(de$padj < 0.1, na.rm = TRUE), "\n")

  # Remove temp column

  meta[, .group := NULL]
  return(de)
}

# ---------------------------------------------------------------------------
# Helper: QC — PCA + missingness per dataset
# ---------------------------------------------------------------------------
run_qc <- function(expr, meta, dataset_name, group_col = "group") {
  cat("  QC for", dataset_name, "...\n")

  # PCA (on complete cases, impute NA with row median for PCA only)
  expr_imp <- expr
  row_medians <- apply(expr_imp, 1, median, na.rm = TRUE)
  for (i in seq_len(nrow(expr_imp))) {
    expr_imp[i, is.na(expr_imp[i, ])] <- row_medians[i]
  }
  pca <- tryCatch(prcomp(t(expr_imp), center = TRUE, scale. = TRUE),
                  error = function(e) NULL)
  if (!is.null(pca)) {
    pca_df <- data.frame(
      sample_id = colnames(expr),
      PC1 = pca$x[, 1], PC2 = pca$x[, 2],
      PC3 = if (ncol(pca$x) >= 3) pca$x[, 3] else NA,
      var_PC1 = round(summary(pca)$importance[2, 1] * 100, 1),
      var_PC2 = round(summary(pca)$importance[2, 2] * 100, 1)
    )
    pca_df$group <- meta[[group_col]][match(pca_df$sample_id, meta$sample_id)]
    fwrite(pca_df, file.path(QC_DIR, paste0("pca_", dataset_name, ".csv")))
    cat("    PCA saved (var explained: PC1=", pca_df$var_PC1[1],
        "%, PC2=", pca_df$var_PC2[1], "%)\n")
  }

  # Missingness per group
  na_counts <- colSums(is.na(expr))
  na_df <- data.frame(
    sample_id = colnames(expr),
    n_missing = na_counts,
    pct_missing = round(100 * na_counts / nrow(expr), 1)
  )
  na_df$group <- meta[[group_col]][match(na_df$sample_id, meta$sample_id)]
  fwrite(na_df, file.path(QC_DIR, paste0("missingness_", dataset_name, ".csv")))

  for (g in unique(na_df$group)) {
    sub <- na_df[na_df$group == g, ]
    cat("    Group", g, ": median missing =", median(sub$pct_missing), "%\n")
  }
}

all_results <- list()

# =========================================================================
# Section B: PXD052937 — Plasma DIA-MS (two contrasts)
#   Contrast 1: MASLD (MASL+MASH, excl cirrhosis) vs Normal → dream_results.csv
#   Contrast 2: MASH vs MASL → nafl_vs_nash_dream.csv
# =========================================================================
cat("\n--- Section B: PXD052937 (plasma DIA-MS) ---\n")

pxd_matrix_f <- file.path(RESULTS_DIR, "pxd052937_protein_matrix.csv")
pxd_meta_f   <- file.path(RESULTS_DIR, "pxd052937_disease_metadata.csv")

uniprot_map <- data.table()

if (file.exists(pxd_matrix_f) && file.exists(pxd_meta_f)) {
  pxd_mat  <- fread(pxd_matrix_f)
  pxd_meta <- fread(pxd_meta_f)

  gene_col    <- names(pxd_mat)[1]
  genes       <- pxd_mat[[gene_col]]
  sample_cols <- setdiff(names(pxd_mat), gene_col)

  meta_all <- pxd_meta[, .(sample_id = file_name, condition = condition_label)]
  matched <- intersect(meta_all$sample_id, sample_cols)
  meta_all <- meta_all[sample_id %in% matched]
  cat("  Total samples matched:", length(matched), "\n")

  # Prepare expression matrix (shared across contrasts)
  expr_pxd <- as.matrix(pxd_mat[, ..matched])
  rownames(expr_pxd) <- genes
  na_frac <- rowSums(is.na(expr_pxd)) / ncol(expr_pxd)
  keep <- na_frac < 0.5
  expr_pxd <- expr_pxd[keep, ]
  cat("  Retained", nrow(expr_pxd), "proteins\n")

  # --- Contrast 1: MASLD (excl cirrhosis) vs Normal ---
  cat("\n  [B1] MASLD (MASL+MASH) vs Normal (excl cirrhosis):\n")
  meta_b1 <- copy(meta_all)
  meta_b1[, group := fifelse(condition == "Normal", "Control",
                    fifelse(condition %in% c("MASL", "MASH"), "MASLD", NA_character_))]
  meta_b1 <- meta_b1[!is.na(group)]
  meta_b1 <- meta_b1[match(intersect(meta_b1$sample_id, colnames(expr_pxd)), sample_id)]
  cat("    Control:", sum(meta_b1$group == "Control"),
      ", MASLD:", sum(meta_b1$group == "MASLD"),
      " (", sum(meta_all$condition == "Cirrhosis"), " cirrhosis excluded)\n")

  if (sum(meta_b1$group == "Control") < 10) {
    cat("    WARNING: <10 controls — individual gene-level DE underpowered.\n")
    cat("    Ranked enrichment (fgsea) is the appropriate validation metric.\n")
  }

  run_qc(expr_pxd[, meta_b1$sample_id], meta_b1, "PXD052937_masld_vs_normal")

  if (length(unique(meta_b1$group)) == 2) {
    all_results[["PXD052937"]] <- run_limma_de(
      expr_pxd, meta_b1, "group", "MASLD - Control", "PXD052937"
    )
  }

  # --- Contrast 2: MASH vs MASL ---
  cat("\n  [B2] MASH vs MASL:\n")
  meta_b2 <- meta_all[condition %in% c("MASL", "MASH")]
  meta_b2[, nash_group := condition]
  meta_b2 <- meta_b2[match(intersect(meta_b2$sample_id, colnames(expr_pxd)), sample_id)]
  cat("    MASH:", sum(meta_b2$nash_group == "MASH"),
      ", MASL:", sum(meta_b2$nash_group == "MASL"), "\n")

  if (length(unique(meta_b2$nash_group)) == 2) {
    all_results[["PXD052937_mash_vs_masl"]] <- run_limma_de(
      expr_pxd, meta_b2, "nash_group", "MASH - MASL", "PXD052937_mash_vs_masl"
    )
  }

  # Build UniProt→gene symbol mapping from Candidates.tsv
  candidates_f <- file.path(BASE, "data/PXD052937/20220805_163627_Plasma_liver2/Candidates.tsv")
  if (file.exists(candidates_f)) {
    cat("\n  Building UniProt→gene mapping from Candidates.tsv...\n")
    cand <- fread(candidates_f, sep = "\t")
    map_raw <- unique(cand[, .(ProteinGroups, Genes)])
    map_raw[, protein_id := sub(";.*", "", ProteinGroups)]
    map_raw[, gene_name  := sub(";.*", "", Genes)]
    uniprot_map <- map_raw[gene_name != "" & !is.na(gene_name),
                           .(protein_id, gene_name)]
    uniprot_map <- unique(uniprot_map, by = "protein_id")
    cat("  Mapped", nrow(uniprot_map), "UniProt IDs to gene symbols\n")
  }
} else {
  cat("  SKIPPED: missing protein matrix or metadata\n")
}

# =========================================================================
# Section B2: PXD051911 — Liver DIA-MS (three contrasts)
#   Boel et al. 2025, Communications Medicine
#   Contrast 1: MASLD vs No_MASLD (+ BMI, age) → dream_results.csv
#   Contrast 2: MASH vs MASL (+ BMI, age) → nafl_vs_nash_dream.csv
#   Contrast 3: High NAS vs Low NAS (+ BMI, age) → nas_score_dream.csv
# =========================================================================
cat("\n--- Section B2: PXD051911 (liver DIA-MS, multi-contrast) ---\n")

pxd051_matrix_f <- file.path(BASE, "data/PXD051911/liver_protein_quant.txt")
pxd051_meta_f   <- file.path(BASE, "data/PXD051911/meta_data.txt")

if (file.exists(pxd051_matrix_f) && file.exists(pxd051_meta_f)) {
  pxd051_mat  <- fread(pxd051_matrix_f, sep = "\t")
  pxd051_meta <- fread(pxd051_meta_f, sep = "\t")

  gene_symbols <- pxd051_mat[["Genes"]]
  sample_cols  <- setdiff(names(pxd051_mat), c("ProteinAccessions", "Genes", "ProteinDescriptions"))

  liver_meta <- pxd051_meta[liver_proteomics_filename != "" &
                            liver_proteomics_filename != "NA" &
                            !is.na(liver_proteomics_filename)]
  cat("  Samples with liver proteomics:", nrow(liver_meta), "\n")

  matched <- intersect(liver_meta$liver_proteomics_filename, sample_cols)
  cat("  Matched", length(matched), "of", nrow(liver_meta), "samples to matrix\n")

  if (length(matched) >= 10) {
    # Build full metadata with covariates
    meta_full <- liver_meta[liver_proteomics_filename %in% matched,
                            .(sample_id = liver_proteomics_filename,
                              saf_diagnosis, kleiner_fibrosis_grade,
                              nafld_activity_score,
                              age = as.numeric(alder),
                              bmi = as.numeric(bmi),
                              sex = factor(gender))]
    # Scale BMI for model stability
    meta_full[, bmi_scaled := as.numeric(scale(bmi))]

    cat("  Covariates: age range [", min(meta_full$age), "-", max(meta_full$age),
        "], BMI range [", round(min(meta_full$bmi), 1), "-", round(max(meta_full$bmi), 1), "]\n")

    # Prepare expression matrix
    expr_051 <- as.matrix(pxd051_mat[, ..matched])
    rownames(expr_051) <- gene_symbols

    # Log2-transform raw intensities (DIA-MS reports linear-scale abundances)
    expr_051 <- log2(expr_051)

    na_frac <- rowSums(is.na(expr_051)) / ncol(expr_051)
    keep <- na_frac < 0.5
    expr_051 <- expr_051[keep, ]
    cat("  Retained", nrow(expr_051), "proteins (", sum(!keep), "dropped for >50% missing)\n")

    meta_full <- meta_full[match(colnames(expr_051), sample_id)]

    # --- Contrast 1: MASLD vs No_MASLD (with covariates) ---
    cat("\n  [B2-1] MASLD (MASL+MASH) vs No_MASLD (adjusted for BMI, age):\n")
    meta_c1 <- copy(meta_full)
    meta_c1[, group := fifelse(saf_diagnosis == "No_MASLD", "Control", "MASLD")]
    cat("    Control:", sum(meta_c1$group == "Control"),
        ", MASLD:", sum(meta_c1$group == "MASLD"), "\n")

    run_qc(expr_051, meta_c1, "PXD051911_masld_vs_ctrl")

    if (length(unique(meta_c1$group)) == 2 && sum(meta_c1$group == "Control") >= 5) {
      all_results[["PXD051911"]] <- run_limma_de(
        expr_051, meta_c1, "group", "MASLD - Control", "PXD051911",
        covariate_cols = c("bmi_scaled", "age")
      )
    }

    # --- Contrast 2: MASH vs MASL (exclude No_MASLD) ---
    cat("\n  [B2-2] MASH vs MASL (adjusted for BMI, age):\n")
    meta_c2 <- meta_full[saf_diagnosis %in% c("MASH", "MASL")]
    meta_c2[, nash_group := saf_diagnosis]
    cat("    MASH:", sum(meta_c2$nash_group == "MASH"),
        ", MASL:", sum(meta_c2$nash_group == "MASL"), "\n")

    if (length(unique(meta_c2$nash_group)) == 2) {
      all_results[["PXD051911_mash_vs_masl"]] <- run_limma_de(
        expr_051, meta_c2, "nash_group", "MASH - MASL", "PXD051911_mash_vs_masl",
        covariate_cols = c("bmi_scaled", "age")
      )
    }

    # --- Contrast 3: High NAS (>=4) vs Low NAS (<4) ---
    cat("\n  [B2-3] High NAS (>=4) vs Low NAS (<4) (adjusted for BMI, age):\n")
    meta_c3 <- copy(meta_full)
    meta_c3[, nas_group := fifelse(as.integer(nafld_activity_score) >= 4,
                                    "HighNAS", "LowNAS")]
    cat("    HighNAS:", sum(meta_c3$nas_group == "HighNAS"),
        ", LowNAS:", sum(meta_c3$nas_group == "LowNAS"), "\n")

    if (length(unique(meta_c3$nas_group)) == 2 &&
        min(table(meta_c3$nas_group)) >= 10) {
      all_results[["PXD051911_nas_high_vs_low"]] <- run_limma_de(
        expr_051, meta_c3, "nas_group", "HighNAS - LowNAS", "PXD051911_nas_high_vs_low",
        covariate_cols = c("bmi_scaled", "age")
      )
    } else {
      cat("    SKIPPED: insufficient samples in one group\n")
    }
  }
} else {
  cat("  SKIPPED: missing liver_protein_quant.txt or meta_data.txt\n")
}

# =========================================================================
# Save combined DE results
# =========================================================================
if (length(all_results) > 0) {
  combined <- rbindlist(all_results, fill = TRUE)
  out_f <- file.path(RESULTS_DIR, "protein_differential_results_v3.csv")
  fwrite(combined, out_f)
  cat("\n=== Protein differential results saved ===\n")
  cat("  File:", out_f, "\n")
  cat("  Contrasts:", paste(names(all_results), collapse = ", "), "\n")
  cat("  Total rows:", nrow(combined), "\n")
} else {
  cat("\nERROR: No datasets produced DE results.\n")
  quit(status = 1)
}

# =========================================================================
# Section C: Matched-contrast concordance
#   Each proteomics contrast is matched to its specific dream comparator
# =========================================================================
cat("\n--- Section C: Protein-transcript concordance (matched contrasts) ---\n")

# Pre-load transcript comparators — ALL on the C2 axis (limma-voom quality-weighted).
# disease_vs_control: canonical_deg_results.csv (C2; repointed 2026-06-08).
# nafl_vs_nash: nafl_vs_nash_lvqw.csv (C2/LVQW; swapped off the retired dream file
#   2026-06-20 — same columns, 27,638 genes). The `dream_*` OUTPUT column names are
#   kept only for downstream-consumer compatibility (27a keys on dream_comparator);
#   the DATA they carry is C2, not dream.
dream_f     <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                         "results/integration/canonical_deg_results.csv")
fib_dream_f <- file.path(DISEASE_SIG_DIR, "adv_vs_early_fibrosis_dream.csv")  # unused (no contrast maps to it)
nn_dream_f  <- file.path(DISEASE_SIG_DIR, "nafl_vs_nash_lvqw.csv")

# Helper: load a dream file and map Ensembl→symbol
load_dream_slim <- function(dream_path, lfc_col = "logFC", padj_col = NULL) {
  if (!file.exists(dream_path)) return(NULL)
  d <- fread(dream_path)
  # Detect padj column name
  if (is.null(padj_col)) {
    padj_col <- intersect(c("padj", "adj.P.Val"), names(d))[1]
  }
  if ("symbol" %in% names(d)) {
    slim <- d[, .(symbol, bulk_logFC = get(lfc_col), bulk_padj = get(padj_col))]
  } else if (!is.null(gene_map)) {
    d[, ensembl_clean := sub("\\..*", "", gene)]
    slim <- merge(d[, .(ensembl_clean, bulk_logFC = get(lfc_col),
                         bulk_padj = get(padj_col))],
                  gene_map, by = "ensembl_clean")
    slim[, ensembl_clean := NULL]
  } else {
    return(NULL)
  }
  slim[!duplicated(symbol)]
}

dream_slim     <- load_dream_slim(dream_f)
fib_dream_slim <- load_dream_slim(fib_dream_f)
nn_dream_slim  <- load_dream_slim(nn_dream_f)

# Map proteomics contrasts → dream comparators
contrast_dream_map <- list(
  PXD052937                = list(dream = dream_slim,     label = "disease_vs_control"),
  PXD052937_mash_vs_masl   = list(dream = nn_dream_slim,  label = "nafl_vs_nash"),
  PXD051911                = list(dream = dream_slim,     label = "disease_vs_control"),
  PXD051911_mash_vs_masl   = list(dream = nn_dream_slim,  label = "nafl_vs_nash"),
  # NOTE (review G2-005): no NAS-high-vs-low dream contrast exists, so the protein
  # NAS-high-vs-low contrast is compared against the disease-vs-control dream as a
  # proxy. The `label` names the PROTEIN contrast; the dream side is disease_vs_control.
  # These rows are intentionally EXCLUDED from the atlas best_protein_* selection
  # (27a keys on dream_comparator == "disease_vs_control").
  PXD051911_nas_high_vs_low = list(dream = dream_slim,    label = "nas_high_vs_low")
)

conc_parts <- list()

for (ds_name in names(all_results)) {
  cdm <- contrast_dream_map[[ds_name]]
  if (is.null(cdm) || is.null(cdm$dream)) next

  ds <- all_results[[ds_name]]
  ds_part <- ds[, .(gene, protein_logFC = logFC, protein_padj = padj, protein_t = t,
                     dataset = ds_name)]

  # Map UniProt → gene symbol for PXD052937 variants
  if (grepl("^PXD052937", ds_name) && nrow(uniprot_map) > 0) {
    ds_part <- merge(ds_part, uniprot_map, by.x = "gene", by.y = "protein_id", all.x = FALSE)
    ds_part[, gene := gene_name]
    ds_part[, gene_name := NULL]
  }

  ds_conc <- merge(ds_part, cdm$dream, by.x = "gene", by.y = "symbol", all = FALSE)
  ds_conc[, dream_comparator := cdm$label]

  if (nrow(ds_conc) > 0) {
    conc_parts[[ds_name]] <- ds_conc
    rho <- cor(ds_conc$protein_logFC, ds_conc$bulk_logFC,
               method = "spearman", use = "complete.obs")
    dir_conc <- sum(sign(ds_conc$protein_logFC) == sign(ds_conc$bulk_logFC), na.rm = TRUE)
    # Filtered concordance: only genes with |LFC| > 0.5 in both
    filt_mask <- abs(ds_conc$protein_logFC) > 0.5 & abs(ds_conc$bulk_logFC) > 0.5
    filt_conc <- if (sum(filt_mask) > 0)
      sum(sign(ds_conc$protein_logFC[filt_mask]) == sign(ds_conc$bulk_logFC[filt_mask]))
    else NA
    filt_n <- sum(filt_mask)

    cat("  ", ds_name, "→", cdm$label, ":", nrow(ds_conc), "genes overlap\n")
    cat("    Spearman rho:", round(rho, 3), "\n")
    cat("    Direction concordance:", dir_conc, "/", nrow(ds_conc),
        "(", round(100 * dir_conc / nrow(ds_conc), 1), "%)\n")
    if (!is.na(filt_conc)) {
      cat("    Filtered (|LFC|>0.5 both):", filt_conc, "/", filt_n,
          "(", round(100 * filt_conc / filt_n, 1), "%)\n")
    }
  }
}

if (length(conc_parts) > 0) {
  concordance <- rbindlist(conc_parts, fill = TRUE)
  concordance[, direction_concordant := sign(protein_logFC) == sign(bulk_logFC)]
  concordance[, direction_concordant_filtered :=
    sign(protein_logFC) == sign(bulk_logFC) &
    abs(protein_logFC) > 0.5 & abs(bulk_logFC) > 0.5]
  conc_f <- file.path(RESULTS_DIR, "protein_transcript_concordance_v3.csv")
  fwrite(concordance, conc_f)
  cat("  Saved:", conc_f, "\n")
} else {
  cat("  WARNING: No concordance data generated\n")
}

# =========================================================================
# Section D: Ranked enrichment (fgsea)
# =========================================================================
cat("\n--- Section D: Ranked enrichment (fgsea) ---\n")

concordance_f <- file.path(BASE,
  "Analysis/Cross_Species_Concordance/results/concordance_atlas_unified.csv")
drug_targets_f <- file.path(BASE,
  "RNA-seq/results/drug_repurposing/convergent_drug_targets.csv")

gene_sets <- list()

if (file.exists(dream_f)) {
  dream <- fread(dream_f)
  # C2 canonical_deg_results.csv already carries a `symbol` column; only fall back
  # to the Ensembl->symbol gene_map merge for the legacy dream layout (no `symbol`).
  if (!("symbol" %in% names(dream)) && !is.null(gene_map)) {
    dream[, ensembl_clean := sub("\\..*", "", gene)]
    dream <- merge(dream, gene_map, by = "ensembl_clean", all.x = FALSE)
  }
  bulk_sig <- dream[padj < 0.1]
  if (nrow(bulk_sig) > 0) {
    gene_sets[["dream_DEG_up"]] <- bulk_sig[logFC > 0, symbol]
    gene_sets[["dream_DEG_down"]] <- bulk_sig[logFC < 0, symbol]
    gene_sets[["dream_DEG_all"]] <- bulk_sig[, symbol]
  }
}

if (file.exists(concordance_f)) {
  conc_atlas <- fread(concordance_f)
  cc_genes <- conc_atlas[primary_category == "Conserved", human_symbol]
  # Gene-set key "Conserved" matches the atlas primary_category value and the stratum
  # label emitted by mrna_protein_concordance.R; figures display it as "Conserved Core".
  if (length(cc_genes) > 0) gene_sets[["Conserved"]] <- cc_genes
}

if (file.exists(drug_targets_f)) {
  drug_dt <- fread(drug_targets_f)
  if ("symbol" %in% names(drug_dt)) {
    gene_sets[["drug_targets"]] <- drug_dt[, symbol]
  }
}

cat("  Gene sets loaded:", paste(names(gene_sets), collapse = ", "), "\n")
cat("  Sizes:", paste(sapply(gene_sets, length), collapse = ", "), "\n")

fgsea_results <- list()

for (ds_name in names(all_results)) {
  ds <- all_results[[ds_name]]
  if (!"t" %in% names(ds)) next

  ds_ranked <- ds[!is.na(t), .(gene, t)]

  # Map UniProt → symbol for PXD052937 variants
  if (grepl("^PXD052937", ds_name) && nrow(uniprot_map) > 0) {
    ds_ranked <- merge(ds_ranked, uniprot_map, by.x = "gene", by.y = "protein_id", all.x = FALSE)
    ds_ranked[, gene := gene_name]
    ds_ranked[, gene_name := NULL]
  }

  ds_ranked <- ds_ranked[!duplicated(gene)]
  ds_ranked <- ds_ranked[!is.na(gene) & gene != ""]
  ranks <- setNames(ds_ranked$t, ds_ranked$gene)
  ranks <- sort(ranks)

  if (length(ranks) < 100) next

  gs_filtered <- lapply(gene_sets, function(gs) intersect(gs, names(ranks)))
  gs_filtered <- gs_filtered[sapply(gs_filtered, length) >= 5]

  if (length(gs_filtered) == 0) next

  cat("  Running fgsea on", ds_name, "(", length(ranks), "genes ranked)...\n")
  # Diagnostic: gene set overlap
  for (gs_name in names(gs_filtered)) {
    cat("    Gene set", gs_name, ":", length(gs_filtered[[gs_name]]), "genes in ranked list\n")
  }

  set.seed(42)  # reproducible fgsea permutation p-values (review G2-009)
  res <- tryCatch({
    fgsea(pathways = gs_filtered, stats = ranks, minSize = 5, maxSize = 15000,
          nPermSimple = 10000)
  }, error = function(e) {
    cat("    fgsea ERROR:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(res)) {
    res <- as.data.table(res)
    res[, dataset := ds_name]
    fgsea_results[[ds_name]] <- res

    for (i in seq_len(nrow(res))) {
      cat("    ", res$pathway[i], ": NES=", round(res$NES[i], 2),
          ", padj=", format(res$padj[i], digits = 3), "\n")
    }
  }
}

if (length(fgsea_results) > 0) {
  fgsea_combined <- rbindlist(fgsea_results, fill = TRUE)
  fgsea_combined[, leadingEdge := sapply(leadingEdge, function(x) paste(x, collapse = ";"))]
  fgsea_f <- file.path(RESULTS_DIR, "protein_ranked_enrichment.csv")
  fwrite(fgsea_combined, fgsea_f)
  cat("  Saved:", fgsea_f, "\n")
} else {
  cat("  WARNING: No fgsea results generated\n")
}

# =========================================================================
# Section E: Effect-size stratified detection
# =========================================================================
cat("\n--- Section E: Effect-size stratified detection ---\n")

if (file.exists(dream_f)) {
  dream <- fread(dream_f)
  # C2 canonical_deg_results.csv already carries `symbol`; only map via gene_map
  # for the legacy dream layout (Ensembl-only `gene`, no `symbol`).
  if (!("symbol" %in% names(dream)) && !is.null(gene_map)) {
    dream[, ensembl_clean := sub("\\..*", "", gene)]
    dream <- merge(dream, gene_map, by = "ensembl_clean", all.x = FALSE)
  }
  bulk_sig <- dream[padj < 0.1]

  bins <- c(0, 0.5, 1, 1.5, 2, Inf)
  bin_labels <- c("0-0.5", "0.5-1", "1-1.5", "1.5-2", "2+")
  bulk_sig[, lfc_bin := cut(abs(logFC), breaks = bins, labels = bin_labels,
                              include.lowest = TRUE, right = FALSE)]

  eff_parts <- list()

  # Only primary contrasts for effect-size detection
  primary_ds <- intersect(c("PXD052937", "PXD051911"),
                           names(all_results))

  for (ds_name in primary_ds) {
    ds <- all_results[[ds_name]]
    if (!"logFC" %in% names(ds)) next

    if (grepl("^PXD052937", ds_name) && nrow(uniprot_map) > 0) {
      ds <- merge(ds, uniprot_map, by.x = "gene", by.y = "protein_id", all.x = FALSE)
      ds[, gene := gene_name]
      ds[, gene_name := NULL]
    }

    for (b in bin_labels) {
      bin_genes <- bulk_sig[lfc_bin == b, symbol]
      n_dream <- length(bin_genes)
      if (n_dream < 5) next

      detected <- ds[gene %in% bin_genes]
      n_detected <- nrow(detected)

      bulk_dirs <- bulk_sig[symbol %in% detected$gene, .(symbol, bulk_dir = sign(logFC))]
      if (nrow(bulk_dirs) > 0 && nrow(detected) > 0) {
        merged <- merge(detected[, .(gene, protein_dir = sign(logFC))],
                        bulk_dirs, by.x = "gene", by.y = "symbol")
        n_concordant <- sum(merged$protein_dir == merged$bulk_dir, na.rm = TRUE)
      } else {
        n_concordant <- NA_integer_
      }

      eff_parts[[paste(ds_name, b)]] <- data.table(
        dataset      = ds_name,
        lfc_bin      = b,
        n_dream_degs = n_dream,
        n_detected   = n_detected,
        detection_rate = n_detected / n_dream,
        n_concordant = n_concordant,
        concordance_rate = ifelse(!is.na(n_concordant) & n_detected > 0,
                                  n_concordant / n_detected, NA_real_)
      )
    }
  }

  if (length(eff_parts) > 0) {
    eff_dt <- rbindlist(eff_parts)
    eff_f <- file.path(RESULTS_DIR, "protein_effectsize_detection.csv")
    fwrite(eff_dt, eff_f)
    cat("  Saved:", eff_f, "\n")
    print(eff_dt)
  }
}

# =========================================================================
# Section F: Positive control + drug target validation
# =========================================================================
cat("\n--- Section F: Positive control + drug target validation ---\n")

validation_parts <- list()

posctrl_f <- file.path(BASE, "RNA-seq/results/validation/positive_control_validation.csv")
if (file.exists(posctrl_f)) {
  posctrl <- fread(posctrl_f)
  pc_symbols <- if ("gene" %in% names(posctrl) && !any(grepl("^ENSG", posctrl$gene)))
    posctrl$gene else if ("symbol" %in% names(posctrl)) posctrl$symbol else character(0)
  cat("  Positive controls:", length(pc_symbols), "\n")
}

if (file.exists(drug_targets_f)) {
  drug_dt <- fread(drug_targets_f)
  dt_symbols <- drug_dt$symbol
  cat("  Drug targets:", length(dt_symbols), "\n")
}

for (ds_name in names(all_results)) {
  ds <- all_results[[ds_name]]
  ds_syms <- ds[, .(gene, logFC, padj)]

  if (grepl("^PXD052937", ds_name) && nrow(uniprot_map) > 0) {
    ds_syms <- merge(ds_syms, uniprot_map, by.x = "gene", by.y = "protein_id", all.x = FALSE)
    ds_syms[, gene := gene_name]
    ds_syms[, gene_name := NULL]
  }

  all_detected <- unique(ds_syms$gene)

  if (exists("pc_symbols") && length(pc_symbols) > 0) {
    pc_detected <- intersect(pc_symbols, all_detected)
    pc_sig <- ds_syms[gene %in% pc_symbols & padj < 0.05]
    validation_parts[[paste0(ds_name, "_posctrl")]] <- data.table(
      dataset        = ds_name,
      validation_set = "positive_controls",
      n_total        = length(pc_symbols),
      n_detected     = length(pc_detected),
      detection_rate = length(pc_detected) / length(pc_symbols),
      n_significant  = nrow(pc_sig),
      sig_rate       = nrow(pc_sig) / max(length(pc_detected), 1)
    )
  }

  if (exists("dt_symbols") && length(dt_symbols) > 0) {
    dt_detected <- intersect(dt_symbols, all_detected)
    dt_sig <- ds_syms[gene %in% dt_symbols & padj < 0.05]
    validation_parts[[paste0(ds_name, "_drug")]] <- data.table(
      dataset        = ds_name,
      validation_set = "drug_targets",
      n_total        = length(dt_symbols),
      n_detected     = length(dt_detected),
      detection_rate = length(dt_detected) / length(dt_symbols),
      n_significant  = nrow(dt_sig),
      sig_rate       = nrow(dt_sig) / max(length(dt_detected), 1)
    )
  }
}

if (length(validation_parts) > 0) {
  val_dt <- rbindlist(validation_parts, fill = TRUE)
  val_f <- file.path(RESULTS_DIR, "protein_validation_summary.csv")
  fwrite(val_dt, val_f)
  cat("  Saved:", val_f, "\n")
  print(val_dt)
}

cat("\n=== Differential Proteomics v3 Complete ===\n")
cat("Time:", as.character(Sys.time()), "\n")
