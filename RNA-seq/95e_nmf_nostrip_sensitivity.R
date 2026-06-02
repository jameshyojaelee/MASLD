#!/usr/bin/env Rscript
# 95e_nmf_nostrip_sensitivity.R — NMF confounder-strip sensitivity (R1 §9.2)
#
# Fits NMF at k=4 and k=6 WITHOUT the confounder strip applied in Script 95,
# then compares program loadings (W matrix) against the stripped versions via
# Hungarian-matched cosine similarity, Pearson correlation, and top-50 gene
# Jaccard overlap.
#
# If cosine > 0.9 across all programs → strip is conservative but not
# distorting; programs reflect genuine biology, not strip artefact.
#
# Inputs:
#   - merged_dge.rds (same as Script 95)
#   - unified_metadata.csv + sample_qc_report.csv (same QC filter)
#   - nmf_results_cache_clean.rds (canonical stripped cache from Script 95)
#   - multi_evidence_atlas.csv (Ensembl → symbol map)
#
# Outputs:
#   RNA-seq/results/subtypes/nmf_strip_sensitivity.csv
#     columns: k, stripped_program, unstripped_program, cosine, pearson,
#              jaccard_top50, stripped_top10, unstripped_top10
#
# SLURM: bigmem, 16 CPU, 500G, 48h, NMF-nostrip

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(edgeR)
  library(limma)
  library(NMF)
  library(doParallel)
})

cat("=== 95e NMF no-strip sensitivity (R1 §9.2) ===\n")
cat("Start:", format(Sys.time()), "\n\n")

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

dge_path    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
meta_path   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
qc_path     <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv")
atlas_path  <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
cache_path  <- file.path(BASE, "RNA-seq/results/subtypes/nmf_results_cache_clean.rds")
out_dir     <- file.path(BASE, "RNA-seq/results/subtypes")
out_csv     <- file.path(out_dir, "nmf_strip_sensitivity.csv")

N_TOP_GENES <- 5000
NMF_RUNS    <- 50
K_VALUES    <- c(4L, 6L)

# .pbackend capped at 2 workers (known NMF package limitation — CLAUDE.md)
NMF_WORKERS <- 2L
cl <- makeCluster(NMF_WORKERS)
registerDoParallel(cl)
on.exit(try(stopCluster(cl), silent = TRUE))

# ============================================================
# 1. Build UNSTRIPPED matrix (same samples/QC as Script 95)
# ============================================================
cat("=== Step 1: Build unstripped input matrix ===\n")

dge  <- readRDS(dge_path)
meta <- fread(meta_path)
qc   <- fread(qc_path)
meta <- meta %>% left_join(qc %>% dplyr::select(sample_id, pass_technical),
                           by = "sample_id")
disease_samples <- meta %>%
  dplyr::filter(group_binary == "Disease" & pass_technical == TRUE) %>%
  dplyr::pull(sample_id)
available <- intersect(disease_samples, colnames(dge))
dge_sub   <- dge[, available]
meta_sub  <- meta %>% dplyr::filter(sample_id %in% available)

cat(sprintf("  Disease QC-passing samples: %d\n", length(available)))

dge_sub <- calcNormFactors(dge_sub, method = "TMM")
v       <- voom(dge_sub, design = NULL, plot = FALSE)
logcpm  <- removeBatchEffect(v$E, batch = meta_sub$dataset)

cat(sprintf("  Full log-CPM (NO strip): %d genes x %d samples\n",
            nrow(logcpm), ncol(logcpm)))

# Symbol map (same as Script 95)
atlas  <- fread(atlas_path, select = c("ensembl_id", "human_symbol"))
atlas[, ens_base := sub("\\..*", "", ensembl_id)]
id2sym <- setNames(atlas$human_symbol, atlas$ens_base)

# Top-IQR gene selection (no strip — applied directly to full matrix)
gene_iqr   <- apply(logcpm, 1, IQR)
top_genes  <- names(sort(gene_iqr, decreasing = TRUE))[1:min(N_TOP_GENES, length(gene_iqr))]
mat_nostrip <- logcpm[top_genes, ]
mat_nn_nostrip <- mat_nostrip - apply(mat_nostrip, 1, min)
mat_nn_nostrip[mat_nn_nostrip < 0] <- 0

nostrip_sample_ids <- colnames(mat_nn_nostrip)
nostrip_gene_ids   <- rownames(mat_nn_nostrip)

cat(sprintf("  Unstripped top-IQR: %d genes x %d samples\n",
            nrow(mat_nn_nostrip), ncol(mat_nn_nostrip)))

# ============================================================
# 2. Load canonical STRIPPED cache
# ============================================================
cat("\n=== Step 2: Load stripped cache ===\n")
stopifnot(file.exists(cache_path))
cached       <- readRDS(cache_path)
mat_nn_strip <- cached$mat_nn
strip_gene_ids   <- rownames(mat_nn_strip)
strip_sample_ids <- colnames(mat_nn_strip)

cat(sprintf("  Stripped matrix: %d genes x %d samples\n",
            nrow(mat_nn_strip), ncol(mat_nn_strip)))

# Ensure sample ordering is identical
stopifnot(all(strip_sample_ids %in% nostrip_sample_ids))
# Reorder unstripped to match stripped sample order
mat_nn_nostrip <- mat_nn_nostrip[, strip_sample_ids]
nostrip_sample_ids <- strip_sample_ids

# ============================================================
# 3. Fit NMF on unstripped matrix at k=4, k=6
# ============================================================
cat("\n=== Step 3: Fit NMF on unstripped matrix ===\n")
nostrip_fits <- list()
for (k in K_VALUES) {
  cat(sprintf("  k=%d: fitting (nrun=%d, workers=%d)\n", k, NMF_RUNS, NMF_WORKERS))
  t0 <- Sys.time()
  set.seed(42 + k)
  fit <- nmf(mat_nn_nostrip, rank = k, nrun = NMF_RUNS,
             method = "brunet", seed = "random",
             .pbackend = "par", .options = list(verbose = FALSE))
  nostrip_fits[[as.character(k)]] <- fit
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("    done in %.1f min (cophenetic=%.3f)\n", el, cophcor(fit)))
}

# ============================================================
# 4. Compare stripped vs unstripped via Hungarian matching
# ============================================================
cat("\n=== Step 4: Program comparison ===\n")

# Helper: cosine similarity between two vectors
cosine_sim <- function(a, b) {
  sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)) + 1e-30)
}

# Helper: Hungarian assignment (minimise cost = -cosine) via R's built-in
# solve_LSAP from the clue package (bundled with R base in many installs).
# Fallback: greedy matching if clue unavailable.
hungarian_match <- function(cost_matrix) {
  # cost_matrix: n_stripped x n_unstripped (values = cosine similarity)
  # We want max-weight matching, so use -cost for LSAP (minimisation).
  if (requireNamespace("clue", quietly = TRUE)) {
    assignment <- clue::solve_LSAP(-cost_matrix)
    return(as.integer(assignment))
  }
  # Greedy fallback
  n <- nrow(cost_matrix)
  used <- logical(ncol(cost_matrix))
  result <- integer(n)
  for (i in seq_len(n)) {
    avail <- which(!used)
    best <- avail[which.max(cost_matrix[i, avail])]
    result[i] <- best
    used[best] <- TRUE
  }
  result
}

results <- list()

for (k in K_VALUES) {
  cat(sprintf("\n  -- k=%d comparison --\n", k))

  # Get stripped NMF fit from cache
  strip_fit <- cached$nmf_results[[as.character(k)]]
  if (is.null(strip_fit)) {
    cat(sprintf("    WARNING: k=%d not found in stripped cache, skipping\n", k))
    next
  }
  nostrip_fit <- nostrip_fits[[as.character(k)]]

  # Extract H matrices (programs x samples) — compare on sample space
  H_strip   <- coef(strip_fit)
  H_nostrip <- coef(nostrip_fit)
  colnames(H_strip)   <- strip_sample_ids
  colnames(H_nostrip) <- nostrip_sample_ids

  # Ensure same sample order
  common_samples <- intersect(colnames(H_strip), colnames(H_nostrip))
  H_strip   <- H_strip[, common_samples, drop = FALSE]
  H_nostrip <- H_nostrip[, common_samples, drop = FALSE]

  # Build cosine similarity matrix (k x k) on H (sample loadings)
  cos_mat <- matrix(0, nrow = k, ncol = k)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      cos_mat[i, j] <- cosine_sim(H_strip[i, ], H_nostrip[j, ])
    }
  }

  # Hungarian matching: row = stripped program, col = unstripped program
  assignment <- hungarian_match(cos_mat)

  # Also compute W-based comparisons on common genes
  W_strip   <- basis(strip_fit)
  W_nostrip <- basis(nostrip_fit)
  rownames(W_strip)   <- strip_gene_ids
  rownames(W_nostrip) <- nostrip_gene_ids

  common_genes <- intersect(rownames(W_strip), rownames(W_nostrip))
  cat(sprintf("    Common genes in W: %d\n", length(common_genes)))

  W_strip_c   <- W_strip[common_genes, , drop = FALSE]
  W_nostrip_c <- W_nostrip[common_genes, , drop = FALSE]

  # Symbol map for top-gene reporting
  strip_ver <- function(x) sub("\\.[0-9]+$", "", x)
  sym_map <- id2sym[strip_ver(common_genes)]
  names(sym_map) <- common_genes

  for (i in seq_len(k)) {
    j <- assignment[i]

    # H-based cosine (already computed)
    h_cosine <- cos_mat[i, j]

    # W-based cosine on common genes
    w_cosine <- cosine_sim(W_strip_c[, i], W_nostrip_c[, j])

    # Pearson on H
    h_pearson <- cor(H_strip[i, ], H_nostrip[j, ], method = "pearson")

    # Pearson on W (common genes)
    w_pearson <- cor(W_strip_c[, i], W_nostrip_c[, j], method = "pearson")

    # Top-50 Jaccard on W loadings
    # For stripped: rank by W column, take top 50
    ord_strip   <- order(W_strip[, i], decreasing = TRUE)[1:50]
    ord_nostrip <- order(W_nostrip[, j], decreasing = TRUE)[1:50]
    top50_strip_ids   <- rownames(W_strip)[ord_strip]
    top50_nostrip_ids <- rownames(W_nostrip)[ord_nostrip]
    # Convert to symbols for Jaccard (handles Ensembl version mismatches)
    top50_strip_sym   <- id2sym[strip_ver(top50_strip_ids)]
    top50_nostrip_sym <- id2sym[strip_ver(top50_nostrip_ids)]
    top50_strip_sym   <- top50_strip_sym[!is.na(top50_strip_sym)]
    top50_nostrip_sym <- top50_nostrip_sym[!is.na(top50_nostrip_sym)]

    inter <- length(intersect(top50_strip_sym, top50_nostrip_sym))
    union <- length(union(top50_strip_sym, top50_nostrip_sym))
    jaccard <- if (union > 0) inter / union else 0

    # Top-10 gene labels for reporting
    top10_strip   <- paste(head(top50_strip_sym, 10), collapse = ",")
    top10_nostrip <- paste(head(top50_nostrip_sym, 10), collapse = ",")

    cat(sprintf("    P%d(strip) <-> P%d(nostrip): H_cos=%.3f W_cos=%.3f H_r=%.3f W_r=%.3f J50=%.3f\n",
                i, j, h_cosine, w_cosine, h_pearson, w_pearson, jaccard))

    results[[length(results) + 1]] <- data.table(
      k = k,
      stripped_program   = paste0("P", i),
      unstripped_program = paste0("P", j),
      cosine_H           = round(h_cosine, 4),
      cosine_W           = round(w_cosine, 4),
      pearson_H          = round(h_pearson, 4),
      pearson_W          = round(w_pearson, 4),
      jaccard_top50      = round(jaccard, 4),
      n_common_genes     = length(common_genes),
      stripped_top10     = top10_strip,
      unstripped_top10   = top10_nostrip
    )
  }
}

# ============================================================
# 5. Write output
# ============================================================
out_dt <- rbindlist(results)
fwrite(out_dt, out_csv)
cat(sprintf("\n=== Output written: %s ===\n", out_csv))
cat(sprintf("    %d rows (k x programs)\n", nrow(out_dt)))

# Summary verdict
mean_h_cos <- mean(out_dt$cosine_H)
min_h_cos  <- min(out_dt$cosine_H)
mean_w_cos <- mean(out_dt$cosine_W)
min_w_cos  <- min(out_dt$cosine_W)
mean_jac   <- mean(out_dt$jaccard_top50)

cat("\n=== Summary ===\n")
cat(sprintf("  H cosine:  mean=%.3f  min=%.3f\n", mean_h_cos, min_h_cos))
cat(sprintf("  W cosine:  mean=%.3f  min=%.3f\n", mean_w_cos, min_w_cos))
cat(sprintf("  Jaccard50: mean=%.3f\n", mean_jac))
if (min_h_cos > 0.9) {
  cat("  VERDICT: All H cosine > 0.9 -> strip is conservative, not distorting.\n")
} else if (min_h_cos > 0.8) {
  cat("  VERDICT: Some H cosine in [0.8, 0.9] -> strip mildly reshapes some programs.\n")
} else {
  cat("  VERDICT: H cosine < 0.8 for some programs -> strip materially alters decomposition.\n")
}

cat("\nDone:", format(Sys.time()), "\n")
