#!/usr/bin/env Rscript
# 73_minimal_panel_extraction.R
# Extract minimal gene panels via stability ranking + RFE and benchmark
# against published MASLD staging signatures.
#
# Phase 5 (Shared Outputs) — aggregates feature selections across all 3 Plans.
#
# Inputs (from Plans 1-3):
#   - selected_features_per_fold.csv  (Script 63, Plan 1: elastic net coefficients)
#   - plan1_feature_importances.csv   (Script 64, Plan 1: RF/XGB/LGBM importance)
#   - tier1_top_features.csv          (Script 65, Plan 2: tier-1 elastic net)
#   - attention_weights_all.csv       (Script 71, Plan 3: transformer attention)
#   - rank_expression_matrix.rds      (expression data)
#   - modeling_metadata.csv           (sample metadata with LOCO folds)
#
# Outputs:
#   - stability_ranking_all_genes.csv
#   - minimal_panel_10.csv, minimal_panel_15.csv, minimal_panel_25.csv, minimal_panel_50.csv
#   - panel_performance_curve.csv
#   - published_panel_benchmark.csv
#   - panel_overlap_analysis.csv
#
# Usage: Rscript 73_minimal_panel_extraction.R
# SLURM: cpu, 8 CPUs, 64GB RAM, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pROC)
})

set.seed(42)

# =============================================================================
# PATHS
# =============================================================================
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results")
OUTDIR <- file.path(RDIR, "staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 73: Minimal Panel Extraction & Published Benchmark ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# HELPER: QWK (manual, consistent with Script 63)
# =============================================================================
compute_qwk <- function(y_true, y_pred, n_classes = NULL) {
  labels <- sort(unique(c(y_true, y_pred)))
  n <- length(labels)
  if (n < 2) return(NA_real_)

  O <- matrix(0, nrow = n, ncol = n)
  for (i in seq_along(y_true)) {
    r <- which(labels == y_true[i])
    c <- which(labels == y_pred[i])
    if (length(r) == 1 && length(c) == 1) O[r, c] <- O[r, c] + 1
  }
  W <- outer(seq_len(n), seq_len(n), function(i, j) (i - j)^2 / (n - 1)^2)
  row_sums <- rowSums(O); col_sums <- colSums(O); total <- sum(O)
  if (total == 0) return(NA_real_)
  E <- outer(row_sums, col_sums) / total
  num <- sum(W * O); den <- sum(W * E)
  if (den == 0) return(1.0)
  return(1 - num / den)
}

# =============================================================================
# HELPER: Gene symbol → Ensembl ID mapping
# =============================================================================
# The expression matrix uses versioned Ensembl IDs (ENSG00000136709.13),
# but published panels use gene symbols (MT1F, CTGF, etc.)
# Build a lookup from the multi-evidence atlas which has both.

ATLAS_PATH <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
symbol_to_ensembl <- list()  # will be populated after loading expression data

# Known HGNC alias mappings (published name -> GENCODE v49 name)
ALIAS_MAP <- c(
  "CTGF"   = "CCN2",
  "LCNL1"  = "LCN15",
  "GDNF"   = "GDNF",
  "A2M"    = "A2M",
  "FSTL1"  = "FSTL1",
  "MFAP4"  = "MFAP4"
)

build_symbol_lookup <- function(atlas_path, available_ensembl_ids) {
  # Build symbol -> versioned_ensembl_id mapping
  lookup <- list()
  if (!file.exists(atlas_path)) {
    cat("  WARNING: Atlas not found, symbol lookup will be empty\n")
    return(lookup)
  }
  atlas <- fread(atlas_path, select = c("human_symbol", "ensembl_id"))
  # Strip version from available IDs for matching
  available_base <- sub("\\.[0-9]+$", "", available_ensembl_ids)
  base_to_versioned <- setNames(available_ensembl_ids, available_base)

  for (i in seq_len(nrow(atlas))) {
    sym <- atlas$human_symbol[i]
    eid <- atlas$ensembl_id[i]
    if (!is.na(sym) && !is.na(eid) && nchar(sym) > 0) {
      # Match unversioned atlas ID to versioned matrix ID
      if (eid %in% names(base_to_versioned)) {
        lookup[[sym]] <- base_to_versioned[[eid]]
      }
    }
  }
  cat("  Symbol lookup built:", length(lookup), "mappings\n")
  return(lookup)
}

resolve_alias <- function(gene, available_genes, symbol_lookup = symbol_to_ensembl) {
  # Direct match (gene is already an Ensembl ID)
  if (gene %in% available_genes) return(gene)

  # Try symbol -> Ensembl lookup
  if (gene %in% names(symbol_lookup)) {
    eid <- symbol_lookup[[gene]]
    if (eid %in% available_genes) return(eid)
  }

  # Try HGNC alias -> symbol -> Ensembl
  if (gene %in% names(ALIAS_MAP)) {
    alias <- ALIAS_MAP[[gene]]
    if (alias %in% available_genes) return(alias)
    if (alias %in% names(symbol_lookup)) {
      eid <- symbol_lookup[[alias]]
      if (eid %in% available_genes) return(eid)
    }
  }

  return(NA_character_)
}

# =============================================================================
# LOAD DATA
# =============================================================================
cat("Loading expression data...\n")
rank_mat <- readRDS(file.path(OUTDIR, "rank_expression_matrix.rds"))
cat("  Expression:", nrow(rank_mat), "genes x", ncol(rank_mat), "samples\n")

cat("Loading metadata...\n")
meta <- fread(file.path(OUTDIR, "modeling_metadata.csv"))
cat("  Metadata:", nrow(meta), "samples\n")

# Align
common_samples <- intersect(colnames(rank_mat), meta$sample_id)
meta <- meta[match(common_samples, meta$sample_id)]
X_full <- t(rank_mat[, common_samples])  # samples x genes
stopifnot(all(rownames(X_full) == meta$sample_id))
all_genes <- colnames(X_full)
cat("  Aligned:", nrow(X_full), "samples x", ncol(X_full), "genes\n")

# Build symbol -> Ensembl lookup for published panel benchmarking
cat("Building symbol-to-Ensembl lookup...\n")
symbol_to_ensembl <- build_symbol_lookup(ATLAS_PATH, all_genes)
cat("\n")

# =============================================================================
# SECTION 1: STABILITY RANKING ACROSS ALL PLANS
# =============================================================================
cat("=== SECTION 1: Stability Ranking ===\n")

# Track (gene -> list of selection events)
gene_selections <- list()
total_model_runs <- 0

# --- Plan 1a: Elastic net selected features per fold (Script 63) ---
f63 <- file.path(OUTDIR, "selected_features_per_fold.csv")
if (file.exists(f63)) {
  cat("  Loading Plan 1 elastic net features...\n")
  sel63 <- fread(f63)
  cat("    Rows:", nrow(sel63), " Cols:", paste(names(sel63), collapse = ", "), "\n")
  # Expect columns: gene, fold (and possibly target, coef)
  if ("gene" %in% names(sel63) && "fold" %in% names(sel63)) {
    folds <- unique(sel63$fold)
    for (fld in folds) {
      genes_in_fold <- unique(sel63[fold == fld, gene])
      for (g in genes_in_fold) {
        gene_selections[[g]] <- c(gene_selections[[g]], paste0("P1_enet_", fld))
      }
      total_model_runs <- total_model_runs + 1
    }
    cat("    Folds:", length(folds),
        " Unique genes:", length(unique(sel63$gene)), "\n")
  } else {
    cat("    WARNING: unexpected columns; skipping\n")
  }
} else {
  cat("  Plan 1 elastic net features not found; skipping\n")
}

# --- Plan 1b: RF/XGB/LGBM feature importances (Script 64) ---
f64 <- file.path(OUTDIR, "plan1_feature_importances.csv")
if (file.exists(f64)) {
  cat("  Loading Plan 1 tree-based importances...\n")
  imp64 <- fread(f64)
  cat("    Rows:", nrow(imp64), " Cols:", paste(names(imp64), collapse = ", "), "\n")
  # Expect: gene, model, fold, importance (or similar)
  # Select genes with importance > 0
  if (all(c("gene", "model", "fold") %in% names(imp64))) {
    imp_col <- intersect(c("importance", "gain", "mean_importance", "feature_importance"),
                         names(imp64))
    if (length(imp_col) == 0) imp_col <- "importance"
    imp_col <- imp_col[1]

    for (mf in unique(paste(imp64$model, imp64$fold, sep = "_"))) {
      parts <- strsplit(mf, "_")[[1]]
      model_name <- parts[1]
      fold_name <- paste(parts[-1], collapse = "_")
      sub <- imp64[model == model_name & fold == fold_name]
      if (imp_col %in% names(sub)) {
        top_genes <- sub[get(imp_col) > 0, gene]
      } else {
        top_genes <- sub$gene  # all genes in this fold
      }
      for (g in top_genes) {
        gene_selections[[g]] <- c(gene_selections[[g]], paste0("P1_tree_", mf))
      }
      total_model_runs <- total_model_runs + 1
    }
    cat("    Unique model-fold combos:", length(unique(paste(imp64$model, imp64$fold))),
        " Unique genes:", length(unique(imp64$gene)), "\n")
  } else if ("gene" %in% names(imp64)) {
    # Simpler format: just gene + importance, one row per gene (aggregated)
    top_genes <- imp64[imp64[[2]] > 0, gene]
    for (g in top_genes) {
      gene_selections[[g]] <- c(gene_selections[[g]], "P1_tree_aggregated")
    }
    total_model_runs <- total_model_runs + 1
    cat("    Aggregated format, genes with importance > 0:", length(top_genes), "\n")
  } else {
    cat("    WARNING: unexpected columns; skipping\n")
  }
} else {
  cat("  Plan 1 tree importances not found; skipping\n")
}

# --- Plan 2: Tier-1 top features (Script 65) ---
f65 <- file.path(OUTDIR, "tier1_top_features.csv")
if (file.exists(f65)) {
  cat("  Loading Plan 2 tier-1 features...\n")
  tier1 <- fread(f65)
  cat("    Rows:", nrow(tier1), "\n")
  if ("gene" %in% names(tier1)) {
    for (g in tier1$gene) {
      gene_selections[[g]] <- c(gene_selections[[g]], "P2_tier1")
    }
    total_model_runs <- total_model_runs + 1
  }
} else {
  cat("  Plan 2 tier-1 features not found; skipping\n")
}

# --- Plan 3: Attention weights (Script 71) ---
f71 <- file.path(OUTDIR, "attention_weights_all.csv")
if (file.exists(f71)) {
  cat("  Loading Plan 3 attention weights...\n")
  attn <- fread(f71)
  cat("    Rows:", nrow(attn), " Cols:", paste(names(attn), collapse = ", "), "\n")
  # Expect: gene, mean_attention, std_attention (or per-seed columns)
  if ("gene" %in% names(attn)) {
    # Use genes with above-median attention as "selected"
    attn_col <- intersect(c("mean_attention", "attention_mean", "mean_weight"),
                          names(attn))
    if (length(attn_col) > 0) {
      attn_col <- attn_col[1]
      med_attn <- median(attn[[attn_col]], na.rm = TRUE)
      selected <- attn[get(attn_col) > med_attn, gene]
    } else {
      # Check for per-seed columns like seed_42, seed_123, ...
      seed_cols <- grep("^seed_|^s[0-9]", names(attn), value = TRUE)
      if (length(seed_cols) > 0) {
        for (sc in seed_cols) {
          med_val <- median(attn[[sc]], na.rm = TRUE)
          sel_genes <- attn[get(sc) > med_val, gene]
          for (g in sel_genes) {
            gene_selections[[g]] <- c(gene_selections[[g]], paste0("P3_attn_", sc))
          }
          total_model_runs <- total_model_runs + 1
        }
        selected <- character(0)  # already handled per seed
      } else {
        # fallback: use all genes as weakly selected
        selected <- attn$gene
      }
    }
    if (length(selected) > 0) {
      for (g in selected) {
        gene_selections[[g]] <- c(gene_selections[[g]], "P3_attention")
      }
      total_model_runs <- total_model_runs + 1
    }
    cat("    Attention-selected genes:", length(selected), "\n")
  }
} else {
  cat("  Plan 3 attention weights not found; skipping\n")
}

cat("\n  Total model runs aggregated:", total_model_runs, "\n")
cat("  Genes with any selection:", length(gene_selections), "\n")

# Compute selection frequency for each gene
stability_dt <- data.table(
  gene = names(gene_selections),
  n_selections = sapply(gene_selections, length),
  selection_sources = sapply(gene_selections, function(x) paste(unique(x), collapse = ";"))
)
stability_dt[, selection_frequency := n_selections / max(1, total_model_runs)]
setorder(stability_dt, -selection_frequency, -n_selections)
stability_dt[, stability_rank := .I]

cat("  Top 10 genes by stability:\n")
print(head(stability_dt[, .(gene, stability_rank, selection_frequency, n_selections)], 10))

fwrite(stability_dt, file.path(OUTDIR, "stability_ranking_all_genes.csv"))
cat("  Saved stability_ranking_all_genes.csv\n\n")

# =============================================================================
# SECTION 2: RECURSIVE FEATURE ELIMINATION (RFE) WITH LOCO-CV
# =============================================================================
cat("=== SECTION 2: Recursive Feature Elimination ===\n")

# Binary target: F >= 3 vs F < 3 (most clinically relevant)
fib_mask <- meta$fib_stage >= 0
meta_fib <- meta[fib_mask]
X_fib <- X_full[fib_mask, ]
y_binary <- as.integer(meta_fib$fib_stage >= 3)

cat("  Binary F>=3 samples:", nrow(meta_fib),
    " (pos:", sum(y_binary), " neg:", sum(y_binary == 0), ")\n")

# LOCO folds
loco_col <- if ("loco_fold_fibrosis" %in% names(meta_fib)) {
  "loco_fold_fibrosis"
} else {
  "dataset"
}
fold_ids <- meta_fib[[loco_col]]
unique_folds <- sort(unique(fold_ids[fold_ids != "excluded" & fold_ids != "NA"]))
cat("  LOCO folds:", length(unique_folds), "\n")

# Top 500 stability-ranked genes (intersected with expression matrix genes)
top_stability_genes <- stability_dt$gene[stability_dt$gene %in% all_genes]
if (length(top_stability_genes) < 500) {
  cat("  WARNING: only", length(top_stability_genes),
      "stability-ranked genes in matrix; using all\n")
  rfe_gene_pool <- top_stability_genes
} else {
  rfe_gene_pool <- top_stability_genes[1:500]
}
cat("  RFE starting pool:", length(rfe_gene_pool), "genes\n")

# Panel sizes to evaluate
panel_sizes <- c(500, 200, 100, 50, 30, 25, 20, 15, 10)
panel_sizes <- panel_sizes[panel_sizes <= length(rfe_gene_pool)]

# Helper: evaluate panel on LOCO with elastic net
evaluate_panel <- function(gene_panel, X, y, folds, fold_labels, alpha = 0.5) {
  gene_panel <- gene_panel[gene_panel %in% colnames(X)]
  if (length(gene_panel) < 2) return(data.table(auroc = NA_real_, n_genes = length(gene_panel)))

  fold_aurocs <- numeric()
  for (fld in fold_labels) {
    test_idx  <- which(folds == fld)
    train_idx <- which(folds != fld & folds != "excluded" & folds != "NA")
    if (length(test_idx) < 5 || length(train_idx) < 20) next
    if (length(unique(y[test_idx])) < 2) next

    X_tr <- X[train_idx, gene_panel, drop = FALSE]
    X_te <- X[test_idx,  gene_panel, drop = FALSE]
    y_tr <- y[train_idx]
    y_te <- y[test_idx]

    # Class weights
    tab <- table(y_tr)
    w <- length(y_tr) / (2 * tab)
    wt <- as.numeric(w[as.character(y_tr)])

    fit <- tryCatch({
      cv.glmnet(
        x = X_tr, y = factor(y_tr),
        family = "binomial", alpha = alpha,
        weights = wt, nfolds = 5, type.measure = "auc"
      )
    }, error = function(e) NULL)

    if (is.null(fit)) next

    prob <- tryCatch({
      as.numeric(predict(fit, newx = X_te, s = "lambda.min", type = "response"))
    }, error = function(e) rep(0.5, nrow(X_te)))

    auroc <- tryCatch({
      as.numeric(pROC::auc(pROC::roc(y_te, prob, quiet = TRUE)))
    }, error = function(e) NA_real_)

    if (!is.na(auroc)) fold_aurocs <- c(fold_aurocs, auroc)
  }

  data.table(
    mean_auroc = if (length(fold_aurocs) > 0) mean(fold_aurocs) else NA_real_,
    sd_auroc   = if (length(fold_aurocs) > 1) sd(fold_aurocs)   else NA_real_,
    n_folds    = length(fold_aurocs),
    n_genes    = length(gene_panel)
  )
}

# Run RFE: at each panel size, pick top-N genes by stability and evaluate
perf_curve <- list()
current_genes <- rfe_gene_pool

for (ps in panel_sizes) {
  panel <- current_genes[1:min(ps, length(current_genes))]
  cat("  Evaluating panel size:", length(panel), "...")
  res <- evaluate_panel(panel, X_fib, y_binary, fold_ids, unique_folds)
  res[, panel_size := length(panel)]
  res[, genes := paste(panel, collapse = ";")]
  perf_curve[[length(perf_curve) + 1]] <- res
  cat(" AUROC:", sprintf("%.3f", res$mean_auroc), "\n")
}

perf_dt <- rbindlist(perf_curve)
setcolorder(perf_dt, c("panel_size", "n_genes", "mean_auroc", "sd_auroc", "n_folds", "genes"))
fwrite(perf_dt, file.path(OUTDIR, "panel_performance_curve.csv"))
cat("  Saved panel_performance_curve.csv\n")

# Save individual panel files
for (ps in c(10, 15, 25, 50)) {
  row <- perf_dt[panel_size == ps]
  if (nrow(row) == 0) next
  panel_genes <- strsplit(row$genes, ";")[[1]]
  panel_info <- stability_dt[gene %in% panel_genes]
  panel_info[, panel_size := ps]
  panel_info[, panel_auroc := row$mean_auroc]
  fwrite(panel_info, file.path(OUTDIR, paste0("minimal_panel_", ps, ".csv")))
  cat("  Saved minimal_panel_", ps, ".csv (", nrow(panel_info), " genes)\n", sep = "")
}

# =============================================================================
# SECTION 3: PUBLISHED PANEL BENCHMARK
# =============================================================================
cat("\n=== SECTION 3: Published Panel Benchmark ===\n")

# Define published signatures
published_panels <- list(
  SteatoSITE_15 = c("MT1F", "KCNH7", "COL25A1", "RASD2", "CTGF", "STC1",
                     "GDNF", "PRRX1", "FGF7", "LCNL1", "DPEP1", "CHRDL2",
                     "LHX6", "POU4F1", "CDH16"),
  Govaere_25 = c("AKR1B10", "DUSP6", "GDF15", "THBS2", "A2M", "CDH2",
                  "COL1A1", "COL3A1", "COL4A1", "COL4A2", "COL6A3", "DCN",
                  "FBN1", "FSTL1", "IGFBP7", "LUM", "MFAP4", "MMP2",
                  "POSTN", "SPARC", "SPP1", "TAGLN", "THY1", "TIMP1", "VCAN")
)

# Also include our best panels for comparison
our_panels <- list()
for (ps in c(15, 25)) {
  row <- perf_dt[panel_size == ps]
  if (nrow(row) > 0) {
    our_panels[[paste0("Our_", ps)]] <- strsplit(row$genes, ";")[[1]]
  }
}

all_panels <- c(published_panels, our_panels)

benchmark_results <- list()

for (pname in names(all_panels)) {
  raw_genes <- all_panels[[pname]]

  # Resolve aliases
  resolved <- vapply(raw_genes, resolve_alias, character(1),
                     available_genes = all_genes)
  found <- resolved[!is.na(resolved)]
  missing <- raw_genes[is.na(resolved)]

  cat("  Panel:", pname, "| Requested:", length(raw_genes),
      "| Found:", length(found))
  if (length(missing) > 0) cat(" | Missing:", paste(missing, collapse = ", "))
  cat("\n")

  if (length(found) < 3) {
    cat("    Too few genes; skipping\n")
    benchmark_results[[length(benchmark_results) + 1]] <- data.table(
      panel_name = pname, n_genes_requested = length(raw_genes),
      n_genes_found = length(found), missing_genes = paste(missing, collapse = ";"),
      mean_auroc = NA_real_, sd_auroc = NA_real_, per_fold_auroc = NA_character_
    )
    next
  }

  res <- evaluate_panel(found, X_fib, y_binary, fold_ids, unique_folds)

  # Also compute per-fold AUROCs for detailed comparison
  fold_aurocs <- numeric()
  for (fld in unique_folds) {
    test_idx  <- which(fold_ids == fld)
    train_idx <- which(fold_ids != fld & fold_ids != "excluded" & fold_ids != "NA")
    if (length(test_idx) < 5 || length(train_idx) < 20) next
    if (length(unique(y_binary[test_idx])) < 2) next

    fit <- tryCatch({
      cv.glmnet(
        x = X_fib[train_idx, found, drop = FALSE],
        y = factor(y_binary[train_idx]),
        family = "binomial", alpha = 0.5, nfolds = 5, type.measure = "auc"
      )
    }, error = function(e) NULL)
    if (is.null(fit)) next

    prob <- tryCatch({
      as.numeric(predict(fit, newx = X_fib[test_idx, found, drop = FALSE],
                         s = "lambda.min", type = "response"))
    }, error = function(e) rep(0.5, length(test_idx)))

    auroc <- tryCatch({
      as.numeric(pROC::auc(pROC::roc(y_binary[test_idx], prob, quiet = TRUE)))
    }, error = function(e) NA_real_)
    fold_aurocs <- c(fold_aurocs, auroc)
  }

  benchmark_results[[length(benchmark_results) + 1]] <- data.table(
    panel_name = pname,
    n_genes_requested = length(raw_genes),
    n_genes_found = length(found),
    missing_genes = paste(missing, collapse = ";"),
    mean_auroc = if (length(fold_aurocs) > 0) mean(fold_aurocs, na.rm = TRUE) else NA_real_,
    sd_auroc   = if (length(fold_aurocs) > 1) sd(fold_aurocs, na.rm = TRUE)   else NA_real_,
    per_fold_auroc = paste(round(fold_aurocs, 4), collapse = ";")
  )
}

bench_dt <- rbindlist(benchmark_results, fill = TRUE)
fwrite(bench_dt, file.path(OUTDIR, "published_panel_benchmark.csv"))
cat("  Saved published_panel_benchmark.csv\n\n")

# =============================================================================
# SECTION 4: CROSS-REFERENCE WITH ATLAS AND CONSERVED CORE
# =============================================================================
cat("=== SECTION 4: Cross-Reference Analysis ===\n")

# Load multi-evidence atlas (bayesian posterior for top 200)
atlas_top200 <- character(0)
atlas_path <- file.path(BASE, "RNA-seq/results/multi_evidence/bayesian_posterior.csv")
if (file.exists(atlas_path)) {
  atlas <- fread(atlas_path)
  cat("  Loaded bayesian_posterior.csv:", nrow(atlas), "rows\n")
  if ("posterior_rank" %in% names(atlas) && "human_symbol" %in% names(atlas)) {
    atlas_top200 <- atlas[posterior_rank <= 200, human_symbol]
  } else if ("human_symbol" %in% names(atlas)) {
    # Sort by first numeric column desc and take top 200
    num_cols <- names(atlas)[sapply(atlas, is.numeric)]
    if (length(num_cols) > 0) {
      setorderv(atlas, num_cols[1], order = -1)
      atlas_top200 <- atlas$human_symbol[1:min(200, nrow(atlas))]
    }
  }
  cat("  Atlas top 200 genes:", length(atlas_top200), "\n")
} else {
  cat("  Convergence score not found; checking multi_evidence_atlas.csv...\n")
  alt_path <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
  if (file.exists(alt_path)) {
    atlas <- fread(alt_path, select = c("human_symbol", "sources_active"))
    setorder(atlas, -sources_active)
    atlas_top200 <- atlas$human_symbol[1:min(200, nrow(atlas))]
    cat("  Atlas top 200 by sources_active:", length(atlas_top200), "\n")
  }
}

# Load Conserved (723 genes)
cc_genes <- character(0)
cc_path <- file.path(BASE,
  "Analysis/Cross_Species_Concordance/results/concordance_atlas_unified.csv")
if (file.exists(cc_path)) {
  cc <- fread(cc_path)
  if ("primary_category" %in% names(cc)) {
    cc_genes <- cc[primary_category == "Conserved", human_symbol]
  }
  cat("  Conserved genes:", length(cc_genes), "\n")
} else {
  cat("  Concordance atlas not found; skipping Conserved\n")
}

# Build overlap table
# For each gene in our best panel (25 or 15), check membership
best_panel_size <- if (25 %in% perf_dt$panel_size) 25 else {
  if (15 %in% perf_dt$panel_size) 15 else max(perf_dt$panel_size)
}
best_row <- perf_dt[panel_size == best_panel_size]
our_best_genes <- if (nrow(best_row) > 0) strsplit(best_row$genes[1], ";")[[1]] else character(0)

# Union of all gene sets for overlap analysis
steatosite_resolved <- vapply(published_panels$SteatoSITE_15, resolve_alias,
                              character(1), available_genes = all_genes)
steatosite_found <- steatosite_resolved[!is.na(steatosite_resolved)]

govaere_resolved <- vapply(published_panels$Govaere_25, resolve_alias,
                           character(1), available_genes = all_genes)
govaere_found <- govaere_resolved[!is.na(govaere_resolved)]

all_panel_genes <- unique(c(our_best_genes, steatosite_found, govaere_found))

overlap_dt <- data.table(
  gene = all_panel_genes,
  in_our_panel = all_panel_genes %in% our_best_genes,
  in_steatosite = all_panel_genes %in% steatosite_found,
  in_govaere = all_panel_genes %in% govaere_found,
  in_conserved = all_panel_genes %in% cc_genes,
  in_atlas_top200 = all_panel_genes %in% atlas_top200
)

# Add stability info where available
overlap_dt <- merge(overlap_dt,
                    stability_dt[, .(gene, stability_rank, selection_frequency)],
                    by = "gene", all.x = TRUE)
setorder(overlap_dt, stability_rank)

fwrite(overlap_dt, file.path(OUTDIR, "panel_overlap_analysis.csv"))
cat("  Saved panel_overlap_analysis.csv (", nrow(overlap_dt), " genes)\n")

cat("\n  Overlap summary:\n")
cat("    Our panel & SteatoSITE:", sum(overlap_dt$in_our_panel & overlap_dt$in_steatosite), "\n")
cat("    Our panel & Govaere:", sum(overlap_dt$in_our_panel & overlap_dt$in_govaere), "\n")
cat("    Our panel & Conserved:", sum(overlap_dt$in_our_panel & overlap_dt$in_conserved), "\n")
cat("    Our panel & Atlas top 200:", sum(overlap_dt$in_our_panel & overlap_dt$in_atlas_top200), "\n")

# =============================================================================
# SECTION 5: NANOSTRING / qPCR FEASIBILITY
# =============================================================================
cat("\n=== SECTION 5: NanoString Feasibility ===\n")

# Load raw logCPM for expression-level analysis (if available)
raw_logcpm_path <- file.path(OUTDIR, "raw_logcpm_matrix.rds")
if (file.exists(raw_logcpm_path)) {
  cat("  Loading raw logCPM for dynamic range analysis...\n")
  raw_logcpm <- readRDS(raw_logcpm_path)

  # For each panel gene, compute expression statistics
  for (ps in c(10, 15, 25)) {
    panel_file <- file.path(OUTDIR, paste0("minimal_panel_", ps, ".csv"))
    if (!file.exists(panel_file)) next

    panel_info <- fread(panel_file)
    panel_genes <- panel_info$gene[panel_info$gene %in% rownames(raw_logcpm)]

    if (length(panel_genes) == 0) next

    expr_stats <- data.table(
      gene = panel_genes,
      mean_logcpm = rowMeans(raw_logcpm[panel_genes, , drop = FALSE], na.rm = TRUE),
      sd_logcpm   = apply(raw_logcpm[panel_genes, , drop = FALSE], 1, sd, na.rm = TRUE),
      min_logcpm  = apply(raw_logcpm[panel_genes, , drop = FALSE], 1, min, na.rm = TRUE),
      max_logcpm  = apply(raw_logcpm[panel_genes, , drop = FALSE], 1, max, na.rm = TRUE)
    )
    expr_stats[, dynamic_range := max_logcpm - min_logcpm]
    expr_stats[, nanostring_feasible := mean_logcpm > 0 & dynamic_range > 1]

    panel_info <- merge(panel_info, expr_stats, by = "gene", all.x = TRUE)
    fwrite(panel_info, panel_file)
    cat("    Updated minimal_panel_", ps, ".csv with expression stats\n", sep = "")

    n_feasible <- sum(expr_stats$nanostring_feasible, na.rm = TRUE)
    cat("    Panel", ps, ": ", n_feasible, "/", nrow(expr_stats),
        " genes NanoString-feasible (mean logCPM > 0, dynamic range > 1)\n")
  }
} else {
  cat("  Raw logCPM not found; skipping NanoString feasibility\n")
}

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=== Summary ===\n")
cat("Stability ranking: ", nrow(stability_dt), " genes ranked from ",
    total_model_runs, " model runs\n")
cat("RFE performance curve: ", nrow(perf_dt), " panel sizes evaluated\n")
if (nrow(perf_dt) > 0) {
  best_size <- perf_dt[which.max(mean_auroc), panel_size]
  best_auroc <- perf_dt[which.max(mean_auroc), mean_auroc]
  cat("  Best panel size:", best_size, "genes, AUROC:", sprintf("%.3f", best_auroc), "\n")
  cat("  Panel at 15 genes: AUROC =",
      sprintf("%.3f", perf_dt[panel_size == 15, mean_auroc]), "\n")
}
cat("Published benchmark:", nrow(bench_dt), " panels evaluated\n")
for (i in seq_len(nrow(bench_dt))) {
  cat("  ", bench_dt$panel_name[i], ": AUROC =",
      sprintf("%.3f", bench_dt$mean_auroc[i]),
      " (", bench_dt$n_genes_found[i], "/", bench_dt$n_genes_requested[i], " genes found)\n")
}

cat("\nFinished:", as.character(Sys.time()), "\n")
