#!/usr/bin/env Rscript
# 144_subtype_progression_interaction.R
# ---------------------------------------------------------------------------
# Test whether NMF subtypes (S1/S2) associate with disease progression.
#
# Analyses:
#   a. Chi-squared: S1/S2 enrichment in NAFL vs NASH
#   b. Chi-squared: S1/S2 enrichment in F3-F4 vs F0-F2
#   c. Ordinal regression: subtype as predictor of fibrosis stage (POLR)
#   d. Fisher exact: S1/S2 overlap with progression DEGs (C2, C3)
#   e. Subtype-specific DEG comparison: which progression genes are S2-specific?
#
# Output: results/progression/
#   subtype_progression_interaction.csv  — statistical test results
#   subtype_stage_enrichment.csv         — S1/S2 per stage distribution
#   subtype_progression_overlap.csv      — overlap between S2 and progression DEGs
#
# Usage: Rscript 144_subtype_progression_interaction.R
# Compute: login node OK (~2 min)
# Requires: data.table, MASS (for polr)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(MASS)  # polr for ordinal regression
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
ODIR <- file.path(INT, "results/progression")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 144: Subtype x Progression Interaction ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ── 1. Load NMF subtype assignments ────────────────────────────────────
cat("--- Loading NMF subtype assignments ---\n")

nmf_path <- file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv")
if (!file.exists(nmf_path)) {
  stop("NMF assignments not found: ", nmf_path)
}
nmf <- fread(nmf_path)
cat("  NMF assignments:", nrow(nmf), "samples\n")
cat("  Subtypes:", paste(names(table(nmf$nmf_subtype)), collapse = ", "), "\n")
cat("  Distribution:\n")
print(table(nmf$nmf_subtype))

# ── 2. Load modeling metadata ──────────────────────────────────────────
cat("\n--- Loading modeling metadata ---\n")

mm_path <- file.path(INT, "results/staging_classifier/modeling_metadata.csv")
if (!file.exists(mm_path)) {
  stop("modeling_metadata.csv not found: ", mm_path)
}
modeling_meta <- fread(mm_path)
cat("  Modeling metadata:", nrow(modeling_meta), "samples\n")

# Merge NMF subtypes with metadata
merged <- merge(nmf, modeling_meta, by = "sample_id", all = FALSE)
cat("  Merged (NMF + metadata):", nrow(merged), "samples\n")

# ── 3. Define stage groups ──────────────────────────────────────────────
cat("\n--- Defining stage groups ---\n")

# Fibrosis binary: Advanced (F3-F4) vs Early (F0-F2)
merged[, fib_binary := fifelse(fibrosis_stage %in% c(3, 4), "Advanced",
                       fifelse(fibrosis_stage %in% c(0, 1, 2), "Early", NA_character_))]

# NAS binary: high (>=5) vs low (<5)
merged[, nas_binary := fifelse(!is.na(nas_score) & nas_score >= 5, "NAS_high",
                      fifelse(!is.na(nas_score) & nas_score < 5, "NAS_low", NA_character_))]

# Diagnosis
merged[, diagnosis := diagnosis_harmonized]

cat("  Fibrosis binary:", sum(!is.na(merged$fib_binary)), "samples\n")
cat("  NAS binary:", sum(!is.na(merged$nas_binary)), "samples\n")
cat("  Diagnosis:", sum(!is.na(merged$diagnosis) & merged$diagnosis != ""), "samples\n")

# Store test results
test_results <- list()

# ── 4a. Chi-squared: S1/S2 in NAFL vs NASH ──────────────────────────────
cat("\n--- Test A: S1/S2 enrichment in NAFL vs NASH ---\n")

sub_a <- merged[!is.na(diagnosis) & diagnosis %in% c("NAFL", "NASH") &
                !is.na(nmf_subtype)]
if (nrow(sub_a) >= 10) {
  tbl_a <- table(sub_a$nmf_subtype, sub_a$diagnosis)
  cat("  Contingency table:\n")
  print(tbl_a)

  chi_a <- chisq.test(tbl_a)
  cat(sprintf("  Chi-squared: X2=%.3f, df=%d, p=%.2e\n",
    chi_a$statistic, chi_a$parameter, chi_a$p.value))

  # Effect size (Cramer's V)
  n_total <- sum(tbl_a)
  k <- min(nrow(tbl_a), ncol(tbl_a))
  cramers_v <- sqrt(chi_a$statistic / (n_total * (k - 1)))
  cat(sprintf("  Cramer's V: %.3f\n", cramers_v))

  # Proportions
  prop_tbl <- prop.table(tbl_a, margin = 2)
  cat("  Column proportions:\n")
  print(round(prop_tbl, 3))

  test_results[["A_NAFL_vs_NASH"]] <- data.table(
    test_id = "A",
    description = "S1/S2 enrichment: NAFL vs NASH",
    test_type = "chi_squared",
    statistic = as.numeric(chi_a$statistic),
    df = chi_a$parameter,
    p_value = chi_a$p.value,
    effect_size = as.numeric(cramers_v),
    effect_label = "Cramers_V",
    n_samples = nrow(sub_a),
    n_group1 = sum(sub_a$diagnosis == "NAFL"),
    n_group2 = sum(sub_a$diagnosis == "NASH"),
    group1_label = "NAFL",
    group2_label = "NASH",
    s1_in_group1 = tbl_a["S1", "NAFL"],
    s2_in_group1 = tbl_a["S2", "NAFL"],
    s1_in_group2 = tbl_a["S1", "NASH"],
    s2_in_group2 = tbl_a["S2", "NASH"]
  )
} else {
  cat("  Too few samples for NAFL vs NASH test\n")
}

# ── 4b. Chi-squared: S1/S2 in F3-F4 vs F0-F2 ────────────────────────────
cat("\n--- Test B: S1/S2 enrichment in Adv vs Early Fibrosis ---\n")

sub_b <- merged[!is.na(fib_binary) & !is.na(nmf_subtype)]
if (nrow(sub_b) >= 10) {
  tbl_b <- table(sub_b$nmf_subtype, sub_b$fib_binary)
  cat("  Contingency table:\n")
  print(tbl_b)

  chi_b <- chisq.test(tbl_b)
  cat(sprintf("  Chi-squared: X2=%.3f, df=%d, p=%.2e\n",
    chi_b$statistic, chi_b$parameter, chi_b$p.value))

  cramers_v_b <- sqrt(chi_b$statistic / (sum(tbl_b) * (min(nrow(tbl_b), ncol(tbl_b)) - 1)))
  cat(sprintf("  Cramer's V: %.3f\n", cramers_v_b))

  prop_tbl_b <- prop.table(tbl_b, margin = 2)
  cat("  Column proportions:\n")
  print(round(prop_tbl_b, 3))

  test_results[["B_Fib_binary"]] <- data.table(
    test_id = "B",
    description = "S1/S2 enrichment: F3-F4 vs F0-F2",
    test_type = "chi_squared",
    statistic = as.numeric(chi_b$statistic),
    df = chi_b$parameter,
    p_value = chi_b$p.value,
    effect_size = as.numeric(cramers_v_b),
    effect_label = "Cramers_V",
    n_samples = nrow(sub_b),
    n_group1 = sum(sub_b$fib_binary == "Early"),
    n_group2 = sum(sub_b$fib_binary == "Advanced"),
    group1_label = "F0-F2",
    group2_label = "F3-F4",
    s1_in_group1 = tbl_b["S1", "Early"],
    s2_in_group1 = tbl_b["S2", "Early"],
    s1_in_group2 = tbl_b["S1", "Advanced"],
    s2_in_group2 = tbl_b["S2", "Advanced"]
  )
} else {
  cat("  Too few samples for fibrosis binary test\n")
}

# ── 4c. Ordinal regression: subtype predicting fibrosis stage ────────────
cat("\n--- Test C: Ordinal regression (subtype -> fibrosis stage) ---\n")

sub_c <- merged[!is.na(fibrosis_stage) & fibrosis_stage >= 0 & !is.na(nmf_subtype)]
if (nrow(sub_c) >= 20 && length(unique(sub_c$fibrosis_stage)) >= 3) {
  sub_c[, fib_ord := ordered(fibrosis_stage)]
  sub_c[, subtype_factor := factor(nmf_subtype, levels = c("S1", "S2"))]

  # Also include dataset as covariate if available
  if ("dataset" %in% names(sub_c) && length(unique(sub_c$dataset)) > 1) {
    sub_c[, dataset_factor := factor(dataset)]
    polr_fit <- tryCatch(
      polr(fib_ord ~ subtype_factor + dataset_factor, data = sub_c, Hess = TRUE),
      error = function(e) {
        cat(sprintf("  polr with dataset failed: %s\n  Retrying without dataset.\n", e$message))
        NULL
      }
    )
    if (is.null(polr_fit)) {
      polr_fit <- tryCatch(
        polr(fib_ord ~ subtype_factor, data = sub_c, Hess = TRUE),
        error = function(e) { cat(sprintf("  polr failed: %s\n", e$message)); NULL }
      )
    }
  } else {
    polr_fit <- tryCatch(
      polr(fib_ord ~ subtype_factor, data = sub_c, Hess = TRUE),
      error = function(e) { cat(sprintf("  polr failed: %s\n", e$message)); NULL }
    )
  }

  if (!is.null(polr_fit)) {
    coef_tbl <- coef(summary(polr_fit))
    # Compute p-value from t-value (approximation)
    s2_row <- coef_tbl["subtype_factorS2", , drop = FALSE]
    s2_coef <- s2_row[1, "Value"]
    s2_se <- s2_row[1, "Std. Error"]
    s2_tval <- s2_row[1, "t value"]
    s2_pval <- 2 * pt(abs(s2_tval), df = nrow(sub_c) - length(coef(polr_fit)), lower.tail = FALSE)
    s2_or <- exp(s2_coef)

    cat(sprintf("  S2 coefficient: %.3f (SE=%.3f, t=%.3f, p=%.2e)\n",
      s2_coef, s2_se, s2_tval, s2_pval))
    cat(sprintf("  S2 odds ratio (cumulative): %.3f\n", s2_or))
    cat(sprintf("  Interpretation: S2 has %.1fx %s odds of higher fibrosis vs S1\n",
      ifelse(s2_or > 1, s2_or, 1/s2_or),
      ifelse(s2_or > 1, "higher", "lower")))

    test_results[["C_ordinal_regression"]] <- data.table(
      test_id = "C",
      description = "Ordinal regression: S2 predicting fibrosis stage",
      test_type = "polr_ordinal",
      statistic = s2_tval,
      df = nrow(sub_c) - length(coef(polr_fit)),
      p_value = s2_pval,
      effect_size = s2_or,
      effect_label = "odds_ratio_S2_vs_S1",
      n_samples = nrow(sub_c),
      n_group1 = sum(sub_c$subtype_factor == "S1"),
      n_group2 = sum(sub_c$subtype_factor == "S2"),
      group1_label = "S1",
      group2_label = "S2",
      s1_in_group1 = NA_integer_,
      s2_in_group1 = NA_integer_,
      s1_in_group2 = NA_integer_,
      s2_in_group2 = NA_integer_
    )
  }
} else {
  cat("  Too few samples or fibrosis levels for ordinal regression\n")
}

# ── 5. Load progression DEGs ────────────────────────────────────────────
cat("\n--- Loading progression DEGs ---\n")

# C2: NAFL vs NASH
c2_path <- file.path(INT, "results/disease_signatures/nafl_vs_nash_dream.csv")
c2_deg <- NULL
if (file.exists(c2_path)) {
  c2_dream <- fread(c2_path)
  c2_deg <- c2_dream[adj.P.Val < 0.1, symbol]
  c2_deg <- c2_deg[!is.na(c2_deg) & c2_deg != ""]
  cat(sprintf("  C2 NAFL-vs-NASH DEGs: %d (padj<0.1)\n", length(c2_deg)))
} else {
  cat("  WARNING: nafl_vs_nash_dream.csv not found\n")
}

# C3: Advanced vs Early Fibrosis
c3_path <- file.path(INT, "results/disease_signatures/adv_vs_early_fibrosis_dream.csv")
c3_deg <- NULL
if (file.exists(c3_path)) {
  c3_dream <- fread(c3_path)
  c3_deg <- c3_dream[adj.P.Val < 0.1, symbol]
  c3_deg <- c3_deg[!is.na(c3_deg) & c3_deg != ""]
  cat(sprintf("  C3 Adv-vs-Early Fibrosis DEGs: %d (padj<0.1)\n", length(c3_deg)))
} else {
  cat("  WARNING: adv_vs_early_fibrosis_dream.csv not found\n")
}

# ── 6. Load NMF basis (subtype-defining genes) ──────────────────────────
cat("\n--- Loading NMF subtype-defining genes ---\n")

# Check for NMF basis matrix, gene weights, or subtype markers
nmf_basis_path <- file.path(BASE, "RNA-seq/results/subtypes/nmf_basis.csv")
nmf_genes_path <- file.path(BASE, "RNA-seq/results/subtypes/nmf_gene_weights.csv")
nmf_markers_path <- file.path(BASE, "RNA-seq/results/subtypes/subtype_markers.csv")

s2_genes <- character(0)
s1_genes <- character(0)

for (np in c(nmf_basis_path, nmf_genes_path, nmf_markers_path)) {
  if (file.exists(np)) {
    nmf_genes <- fread(np)
    cat(sprintf("  NMF gene data loaded from %s: %d rows\n", basename(np), nrow(nmf_genes)))

    # Wide format: columns named S1/S2 (basis matrix or gene weights)
    if (all(c("S1", "S2") %in% names(nmf_genes))) {
      gene_col <- intersect(c("gene", "symbol", "gene_symbol", "human_symbol"), names(nmf_genes))[1]
      if (!is.na(gene_col)) {
        nmf_genes[, s2_specificity := S2 - S1]
        s2_genes <- nmf_genes[order(-s2_specificity)][1:min(500, nrow(nmf_genes)), get(gene_col)]
        s1_genes <- nmf_genes[order(s2_specificity)][1:min(500, nrow(nmf_genes)), get(gene_col)]
        s2_genes <- s2_genes[!is.na(s2_genes) & s2_genes != ""]
        s1_genes <- s1_genes[!is.na(s1_genes) & s1_genes != ""]
        cat(sprintf("  S2-specific genes (top 500): %d\n", length(s2_genes)))
        cat(sprintf("  S1-specific genes (top 500): %d\n", length(s1_genes)))
      }
    # Long format: subtype_markers.csv with columns gene, subtype, direction, t
    } else if (all(c("subtype", "gene") %in% names(nmf_genes))) {
      cat("  Detected long-format subtype markers\n")
      s2_up <- nmf_genes[subtype == "S2" & direction == "up"][order(-t)]
      s1_up <- nmf_genes[subtype == "S1" & direction == "up"][order(-t)]
      s2_genes <- head(s2_up$gene, 500)
      s1_genes <- head(s1_up$gene, 500)
      s2_genes <- s2_genes[!is.na(s2_genes) & s2_genes != ""]
      s1_genes <- s1_genes[!is.na(s1_genes) & s1_genes != ""]
      cat(sprintf("  S2-specific genes (top 500 up): %d\n", length(s2_genes)))
      cat(sprintf("  S1-specific genes (top 500 up): %d\n", length(s1_genes)))
    }
    break
  }
}

if (length(s2_genes) == 0) {
  cat("  WARNING: NMF gene weights not found — skipping gene-level overlap\n")
} else {
  # If s2_genes are Ensembl IDs, map to symbols for overlap with DEGs
  if (any(grepl("^ENSG", s2_genes))) {
    annot_file <- file.path(INT, "results/disease_signatures/unified_disease_signatures.csv")
    if (file.exists(annot_file)) {
      id_map <- fread(annot_file, select = c("gene", "symbol"))
      id_map <- unique(id_map[!is.na(symbol) & symbol != ""])
      s2_genes <- id_map[gene %in% s2_genes, symbol]
      s1_genes <- id_map[gene %in% s1_genes, symbol]
      cat(sprintf("  Mapped to symbols: S2=%d, S1=%d\n", length(s2_genes), length(s1_genes)))
    }
  }
}

# ── 4d. Fisher exact: S2 overlap with progression DEGs ──────────────────
cat("\n--- Test D: S2-specific gene overlap with progression DEGs ---\n")

if (length(s2_genes) > 0) {
  # Use unified disease signatures as universe (has symbols)
  annot_file <- file.path(INT, "results/disease_signatures/unified_disease_signatures.csv")
  if (file.exists(annot_file)) {
    gene_universe <- fread(annot_file, select = "symbol")$symbol
    gene_universe <- gene_universe[!is.na(gene_universe) & gene_universe != ""]
    gene_universe <- unique(gene_universe)
  } else {
    gene_universe <- unique(c(s2_genes, s1_genes, c2_deg, c3_deg))
  }

  for (contrast_name in c("C2", "C3")) {
    deg_set <- if (contrast_name == "C2") c2_deg else c3_deg
    if (is.null(deg_set) || length(deg_set) == 0) next

    # S2 x progression DEG overlap
    overlap <- intersect(s2_genes, deg_set)
    n_overlap <- length(overlap)
    n_s2 <- sum(s2_genes %in% gene_universe)
    n_deg <- sum(deg_set %in% gene_universe)
    n_univ <- length(gene_universe)

    ft <- fisher.test(matrix(c(
      n_overlap,
      n_s2 - n_overlap,
      n_deg - n_overlap,
      n_univ - n_s2 - n_deg + n_overlap
    ), nrow = 2), alternative = "greater")

    cat(sprintf("  %s: S2 x progression DEG overlap = %d (OR=%.2f, p=%.2e)\n",
      contrast_name, n_overlap, ft$estimate, ft$p.value))

    test_results[[paste0("D_S2_overlap_", contrast_name)]] <- data.table(
      test_id = paste0("D_", contrast_name),
      description = sprintf("Fisher: S2-specific x %s progression DEG overlap", contrast_name),
      test_type = "fisher_exact",
      statistic = ft$estimate,
      df = NA_integer_,
      p_value = ft$p.value,
      effect_size = ft$estimate,
      effect_label = "odds_ratio",
      n_samples = n_univ,
      n_group1 = n_s2,
      n_group2 = n_deg,
      group1_label = "S2_genes",
      group2_label = paste0(contrast_name, "_DEGs"),
      s1_in_group1 = NA_integer_,
      s2_in_group1 = NA_integer_,
      s1_in_group2 = NA_integer_,
      s2_in_group2 = n_overlap
    )
  }
} else {
  cat("  Skipped: no S2 gene set available\n")
}

# ── 4e. Subtype-specific progression DEGs ────────────────────────────────
cat("\n--- Test E: Subtype-specific progression DEG characterization ---\n")

overlap_results <- list()

if (length(s2_genes) > 0 && (!is.null(c2_deg) || !is.null(c3_deg))) {
  for (contrast_name in c("C2", "C3")) {
    deg_set <- if (contrast_name == "C2") c2_deg else c3_deg
    dream_dt <- if (contrast_name == "C2") {
      if (exists("c2_dream")) c2_dream else NULL
    } else {
      if (exists("c3_dream")) c3_dream else NULL
    }

    if (is.null(deg_set) || length(deg_set) == 0) next

    # Classify each DEG
    s2_overlap <- intersect(deg_set, s2_genes)
    s1_overlap <- intersect(deg_set, s1_genes)
    both_overlap <- intersect(s2_overlap, s1_overlap)
    s2_only <- setdiff(s2_overlap, s1_overlap)
    s1_only <- setdiff(s1_overlap, s2_overlap)
    neither <- setdiff(deg_set, union(s2_genes, s1_genes))

    cat(sprintf("  %s progression DEGs: %d total\n", contrast_name, length(deg_set)))
    cat(sprintf("    S2-specific: %d (%.1f%%)\n", length(s2_only),
      100 * length(s2_only) / length(deg_set)))
    cat(sprintf("    S1-specific: %d (%.1f%%)\n", length(s1_only),
      100 * length(s1_only) / length(deg_set)))
    cat(sprintf("    Both: %d\n", length(both_overlap)))
    cat(sprintf("    Neither: %d\n", length(neither)))

    # Build overlap table with LFC direction from dream
    if (!is.null(dream_dt)) {
      for (gene in c(s2_only, s1_only, both_overlap)) {
        gene_row <- dream_dt[symbol == gene]
        if (nrow(gene_row) > 0) {
          overlap_results[[length(overlap_results) + 1]] <- data.table(
            gene = gene,
            contrast = contrast_name,
            logFC = gene_row$logFC[1],
            padj = gene_row$adj.P.Val[1],
            subtype_class = fifelse(
              gene %in% s2_only, "S2_specific",
              fifelse(gene %in% s1_only, "S1_specific", "shared"))
          )
        }
      }
    }
  }
}

if (length(overlap_results) > 0) {
  overlap_dt <- rbindlist(overlap_results)
  fwrite(overlap_dt, file.path(ODIR, "subtype_progression_overlap.csv"))
  cat(sprintf("\n  Overlap table written: %d entries\n", nrow(overlap_dt)))

  # Summary by subtype class
  cat("  Summary by subtype class:\n")
  print(overlap_dt[, .(
    n_genes = .N,
    mean_absLFC = round(mean(abs(logFC), na.rm = TRUE), 3),
    median_padj = signif(median(padj, na.rm = TRUE), 3)
  ), by = .(contrast, subtype_class)])
} else {
  cat("  No overlap data to write\n")
}

# ── 7. Build stage enrichment table ─────────────────────────────────────
cat("\n--- Building stage enrichment table ---\n")

# S1/S2 distribution per fibrosis stage
stage_enrichment <- list()

if (any(!is.na(merged$fibrosis_stage))) {
  fib_dist <- merged[!is.na(fibrosis_stage) & !is.na(nmf_subtype),
    .N, by = .(nmf_subtype, fibrosis_stage)]
  fib_dist[, total_in_stage := sum(N), by = fibrosis_stage]
  fib_dist[, fraction := N / total_in_stage]
  fib_dist[, stage_type := "fibrosis"]
  setnames(fib_dist, "fibrosis_stage", "stage_value")
  stage_enrichment[["fibrosis"]] <- fib_dist
  cat("  Fibrosis stage distribution:\n")
  print(dcast(fib_dist, stage_value ~ nmf_subtype, value.var = "fraction",
    fill = 0)[, lapply(.SD, round, 3), .SDcols = is.numeric, by = stage_value])
}

if (any(!is.na(merged$nas_score))) {
  nas_dist <- merged[!is.na(nas_score) & !is.na(nmf_subtype),
    .N, by = .(nmf_subtype, nas_score)]
  nas_dist[, total_in_stage := sum(N), by = nas_score]
  nas_dist[, fraction := N / total_in_stage]
  nas_dist[, stage_type := "nas"]
  setnames(nas_dist, "nas_score", "stage_value")
  stage_enrichment[["nas"]] <- nas_dist
  cat("\n  NAS score distribution:\n")
  print(dcast(nas_dist, stage_value ~ nmf_subtype, value.var = "fraction",
    fill = 0)[, lapply(.SD, round, 3), .SDcols = is.numeric, by = stage_value])
}

if (any(!is.na(merged$diagnosis) & merged$diagnosis != "")) {
  diag_dist <- merged[!is.na(diagnosis) & diagnosis != "" & !is.na(nmf_subtype),
    .N, by = .(nmf_subtype, diagnosis)]
  diag_dist[, total_in_stage := sum(N), by = diagnosis]
  diag_dist[, fraction := N / total_in_stage]
  diag_dist[, stage_type := "diagnosis"]
  setnames(diag_dist, "diagnosis", "stage_value")
  # Ensure stage_value is character for rbind
  diag_dist[, stage_value := as.character(stage_value)]
  stage_enrichment[["diagnosis"]] <- diag_dist
}

if (length(stage_enrichment) > 0) {
  # Ensure all stage_value columns are character before binding
  for (i in seq_along(stage_enrichment)) {
    stage_enrichment[[i]][, stage_value := as.character(stage_value)]
  }
  stage_dt <- rbindlist(stage_enrichment, fill = TRUE)
  fwrite(stage_dt, file.path(ODIR, "subtype_stage_enrichment.csv"))
  cat(sprintf("\n  Stage enrichment table written: %d rows\n", nrow(stage_dt)))
}

# ── 8. Write test results ──────────────────────────────────────────────
cat("\n--- Writing test results ---\n")

if (length(test_results) > 0) {
  tests_dt <- rbindlist(test_results, fill = TRUE)
  fwrite(tests_dt, file.path(ODIR, "subtype_progression_interaction.csv"))
  cat(sprintf("  Test results written: %d tests\n", nrow(tests_dt)))

  cat("\n  Test summary:\n")
  print(tests_dt[, .(test_id, description,
    statistic = round(statistic, 3),
    p_value = signif(p_value, 3),
    effect_size = round(effect_size, 3))])
} else {
  cat("  WARNING: No test results to write\n")
}

# ── Summary ──────────────────────────────────────────────────────────────
cat("\n=== RESULTS SUMMARY ===\n")
cat(sprintf("  NMF samples merged: %d\n", nrow(merged)))

# S2 = female-enriched progressor subtype (from MEMORY.md)
if (nrow(merged[!is.na(nmf_subtype)]) > 0) {
  s2_frac_adv <- merged[fib_binary == "Advanced" & nmf_subtype == "S2", .N] /
                 max(merged[fib_binary == "Advanced", .N], 1)
  s2_frac_early <- merged[fib_binary == "Early" & nmf_subtype == "S2", .N] /
                   max(merged[fib_binary == "Early", .N], 1)
  cat(sprintf("  S2 fraction in Advanced fibrosis: %.1f%%\n", 100 * s2_frac_adv))
  cat(sprintf("  S2 fraction in Early fibrosis: %.1f%%\n", 100 * s2_frac_early))
}

if (length(test_results) > 0) {
  sig_tests <- sum(sapply(test_results, function(x) x$p_value < 0.05))
  cat(sprintf("  Significant tests (p<0.05): %d / %d\n", sig_tests, length(test_results)))
}

cat("\nOutputs written to:", ODIR, "\n")
cat("  subtype_progression_interaction.csv\n")
cat("  subtype_stage_enrichment.csv\n")
cat("  subtype_progression_overlap.csv\n")
cat("\nDone:", format(Sys.time()), "\n")
