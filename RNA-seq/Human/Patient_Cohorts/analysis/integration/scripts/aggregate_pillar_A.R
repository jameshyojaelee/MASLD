#!/usr/bin/env Rscript
# aggregate_pillar_A.R
# ---------------------------------------------------------------------------
# Aggregate Pillar A (sample-level stability) per-iter outputs into per-gene
# stability metrics + per-iter summary statistics. Run after CPSS, bootstrap,
# and k-fold arrays have all completed.
#
# Inputs:
#   robustness/cpss/iter_{001..100}_half{A,B}.csv    (200 fits)
#   robustness/bootstrap/iter_{0001..1000}.csv       (1000 fits)
#   robustness/kfold/iter_{01..60}.csv               (60 fits)
# Outputs:
#   audit_sensitivity/pillar_A_stability.csv         (per-gene)
#   audit_sensitivity/pillar_A_iter_summary.csv      (per-iter rho/Jaccard)
#   audit_sensitivity/pillar_A_kfold_summary.csv     (per-fold metrics)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(yaml)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CPSS_DIR  <- file.path(RDIR, "robustness/cpss")
BOOT_DIR  <- file.path(RDIR, "robustness/bootstrap")
KFOLD_DIR <- file.path(RDIR, "robustness/kfold")

# Threshold parameterization. Set LFC_THR via env var to compute Pillar A at a
# different LFC cutoff. Output filenames carry a tag derived from the threshold
# (e.g. "lfc05" for >0.5, "nolfc" for 0). Two canonical settings:
#   LFC_THR=0.5  -> Tier 1 primary DEG list (~1,551 genes, kallisto canonical; CLAUDE.md primary)
#   LFC_THR=0    -> padj-only sig DEG list (13,615 genes; sensitivity sweep)
PADJ_THR <- 0.05
LFC_THR  <- as.numeric(Sys.getenv("LFC_THR", "0.5"))
.thr_tag_env <- Sys.getenv("THR_TAG", "")
THR_TAG <- if (nzchar(.thr_tag_env)) {
  .thr_tag_env
} else if (LFC_THR == 0) {
  "nolfc"
} else {
  paste0("lfc", sub("^0\\.", "", as.character(LFC_THR)))
}
cat(sprintf("Aggregating at padj<%.2f & |logFC|>%.2f  (tag = %s)\n", PADJ_THR, LFC_THR, THR_TAG))

# --- Load full dream results (canonical ~4,370 DEG list at primary threshold, kallisto canonical) ---
full <- fread(file.path(RDIR, "dream_results.csv"))
deg_full <- full[padj < PADJ_THR & abs(logFC) > LFC_THR, gene]
cat(sprintf("Canonical DEGs at padj<%.2f & |LFC|>%.2f: %d\n",
            PADJ_THR, LFC_THR, length(deg_full)))

# Helper: per-iter selected genes
sel_from_dt <- function(dt) dt[padj < PADJ_THR & abs(logFC) > LFC_THR, gene]

# ============================================================================
# CPSS — per-gene π̂ (max selection probability over the two halves)
# ============================================================================
cat("\n── CPSS aggregation ──\n")
cpss_files <- list.files(CPSS_DIR, pattern = "^iter_\\d{3}_half[AB]\\.csv$", full.names = TRUE)
cat("Files found:", length(cpss_files), "\n")

if (length(cpss_files) > 0) {
  # Tag each file with its pair index and half label so we can compute π̂ correctly
  meta <- data.table(
    file = cpss_files,
    pair = as.integer(sub("^iter_(\\d{3}).*$", "\\1", basename(cpss_files))),
    half = sub("^iter_\\d{3}_(half[AB])\\.csv$", "\\1", basename(cpss_files)))

  read_sel <- function(f) {
    dt <- fread(f, select = c("gene", "logFC", "padj"))
    dt[padj < PADJ_THR & abs(logFC) > LFC_THR, gene]
  }
  sel_list <- lapply(meta$file, read_sel)
  meta[, n_sel := vapply(sel_list, length, integer(1))]

  # Selection per (pair, half) → per-gene presence matrix
  all_genes <- full$gene
  n_pairs <- length(unique(meta$pair))
  cat(sprintf("Pairs found: %d (expected B=100; will compute on what's present)\n", n_pairs))

  # For each pair, compute s^A_g = 1{g selected in half A}, s^B_g = 1{g in B}
  # CPSS selection probability per pair: max(s^A, s^B)
  # π̂_g = mean over pairs of max(s^A, s^B)
  # Shah-Samworth PFER bound: V <= q^2 / ((2*pi - 1)^2 * p)  (r-concave version)
  pair_max <- matrix(0L, nrow = length(all_genes), ncol = n_pairs,
                     dimnames = list(all_genes, NULL))
  pair_index <- sort(unique(meta$pair))
  for (j in seq_along(pair_index)) {
    pi <- pair_index[j]
    rows <- meta[pair == pi]
    sels <- unique(unlist(sel_list[as.integer(rownames(meta)[meta$pair == pi])]))
    # Equivalent: sels = union across both halves of the pair
    pair_max[match(sels, all_genes), j] <- 1L
  }
  pi_hat <- rowMeans(pair_max)

  # Average size of selected set across all 2*B fits
  q_hat <- mean(meta$n_sel)
  p <- nrow(full)

  # Shah-Samworth PFER bound (r-concave variant): E[V] <= q^2 / ((2*pi - 1)^2 * p)
  # For each gene's pi_hat, this gives an INDIVIDUAL upper bound — interpret as
  # bound at threshold = pi_hat. We compute a bound at a fixed pi_thr=0.7 too.
  pfer_at_pi <- function(pi, q, p) {
    out <- numeric(length(pi))
    out[] <- NA_real_
    valid <- pi > 0.5
    out[valid] <- q^2 / ((2 * pi[valid] - 1)^2 * p)
    out
  }
  pi_thr_default <- 0.7
  pfer_global <- q_hat^2 / ((2 * pi_thr_default - 1)^2 * p)
  cat(sprintf("Shah-Samworth PFER bound at π_thr=%.2f: %.3f (q=%.0f, p=%d)\n",
              pi_thr_default, pfer_global, q_hat, p))

  cpss_dt <- data.table(
    gene = all_genes,
    cpss_pi_hat = round(pi_hat, 4),
    cpss_pfer_at_pi_hat = round(pfer_at_pi(pi_hat, q_hat, p), 3))
} else {
  cat("No CPSS files; skipping.\n")
  cpss_dt <- data.table(gene = full$gene,
                        cpss_pi_hat = NA_real_,
                        cpss_pfer_at_pi_hat = NA_real_)
}

# ============================================================================
# Bootstrap — per-gene selection frequency + bootstrap CI on logFC
# ============================================================================
cat("\n── Bootstrap aggregation ──\n")
boot_files <- list.files(BOOT_DIR, pattern = "^iter_\\d{4}\\.csv$", full.names = TRUE)
cat("Files found:", length(boot_files), "\n")

if (length(boot_files) > 0) {
  # Stack into one long table with iter id; vectorized aggregation
  boot_long <- rbindlist(lapply(boot_files, function(f) {
    dt <- fread(f, select = c("gene", "logFC", "padj"))
    dt[, iter := as.integer(sub("^iter_(\\d{4})\\.csv$", "\\1", basename(f)))]
    dt
  }))
  cat("Bootstrap rows:", nrow(boot_long), "\n")
  boot_long[, sel := padj < PADJ_THR & abs(logFC) > LFC_THR]

  boot_dt <- boot_long[, .(
    bootstrap_freq = mean(sel),
    bootstrap_logFC_mean = mean(logFC),
    bootstrap_logFC_sd   = sd(logFC),
    bootstrap_logFC_q025 = quantile(logFC, 0.025, na.rm = TRUE),
    bootstrap_logFC_q975 = quantile(logFC, 0.975, na.rm = TRUE),
    bootstrap_n_iter     = .N), by = gene]
  boot_dt[, bootstrap_freq := round(bootstrap_freq, 4)]
  boot_dt[, bootstrap_logFC_mean := round(bootstrap_logFC_mean, 4)]
  boot_dt[, bootstrap_logFC_sd   := round(bootstrap_logFC_sd, 4)]
  boot_dt[, bootstrap_logFC_q025 := round(bootstrap_logFC_q025, 4)]
  boot_dt[, bootstrap_logFC_q975 := round(bootstrap_logFC_q975, 4)]
} else {
  cat("No bootstrap files; skipping.\n")
  boot_dt <- data.table(gene = full$gene,
                        bootstrap_freq = NA_real_,
                        bootstrap_logFC_mean = NA_real_,
                        bootstrap_logFC_sd = NA_real_,
                        bootstrap_logFC_q025 = NA_real_,
                        bootstrap_logFC_q975 = NA_real_,
                        bootstrap_n_iter = 0L)
}

# ============================================================================
# K-fold — per-gene fold recurrence + per-fold ρ / Jaccard / recovery%
# ============================================================================
cat("\n── K-fold aggregation ──\n")
kfold_files <- list.files(KFOLD_DIR, pattern = "^iter_\\d{2}\\.csv$", full.names = TRUE)
cat("Files found:", length(kfold_files), "\n")

kfold_iter_summary <- data.table()
kfold_per_gene <- data.table(gene = full$gene, kfold10_recur = 0L, kfold50_recur = 0L)

if (length(kfold_files) > 0) {
  kf_long <- rbindlist(lapply(kfold_files, function(f) {
    dt <- fread(f)
    dt[, iter := as.integer(sub("^iter_(\\d{2})\\.csv$", "\\1", basename(f)))]
    dt
  }))
  cat("K-fold rows:", nrow(kf_long), "\n")

  # Per-gene recurrence (count of folds where gene was in DEG list)
  kf_long[, sel := padj < PADJ_THR & abs(logFC) > LFC_THR]
  rec10 <- kf_long[K == 10 & sel == TRUE, .N, by = gene]; setnames(rec10, "N", "kfold10_recur")
  rec50 <- kf_long[K == 50 & sel == TRUE, .N, by = gene]; setnames(rec50, "N", "kfold50_recur")
  kfold_per_gene <- merge(kfold_per_gene[, .(gene)], rec10, by = "gene", all.x = TRUE)
  kfold_per_gene <- merge(kfold_per_gene, rec50, by = "gene", all.x = TRUE)
  kfold_per_gene[is.na(kfold10_recur), kfold10_recur := 0L]
  kfold_per_gene[is.na(kfold50_recur), kfold50_recur := 0L]

  # Per-fold metrics: ρ, Jaccard, % canonical DEG recovery
  full_thin <- full[, .(gene, full_logFC = logFC)]
  iter_meta <- unique(kf_long[, .(iter, K, fold)])
  kfold_iter_summary <- rbindlist(lapply(seq_len(nrow(iter_meta)), function(j) {
    it <- iter_meta[j, iter]
    dt <- kf_long[iter == it]
    if (nrow(dt) == 0) return(NULL)
    sel_genes <- dt[sel == TRUE, gene]
    m <- merge(full_thin, dt[, .(gene, fold_logFC = logFC)], by = "gene")
    rho <- suppressWarnings(cor(m$full_logFC, m$fold_logFC, method = "spearman"))
    inter <- length(intersect(deg_full, sel_genes))
    union_n <- length(union(deg_full, sel_genes))
    jacc <- if (union_n > 0) inter / union_n else 0
    rec_pct <- if (length(deg_full) > 0) inter / length(deg_full) else 0
    data.table(iter = it, K = unique(dt$K), fold = unique(dt$fold),
               n_sel = length(sel_genes),
               rho = round(rho, 4),
               jaccard = round(jacc, 4),
               recovery_pct = round(rec_pct * 100, 2))
  }))
}

# ============================================================================
# Per-iter summary for CPSS + bootstrap (rho with full, Jaccard, recovery%)
# ============================================================================
cat("\n── Per-iter summary (CPSS + bootstrap) ──\n")
full_thin <- full[, .(gene, full_logFC = logFC)]

iter_summary_one <- function(file_list, iter_pattern) {
  if (length(file_list) == 0) return(data.table())
  rbindlist(lapply(file_list, function(f) {
    dt <- fread(f, select = c("gene", "logFC", "padj"))
    sel_genes <- dt[padj < PADJ_THR & abs(logFC) > LFC_THR, gene]
    m <- merge(full_thin, dt[, .(gene, fold_logFC = logFC)], by = "gene")
    rho <- suppressWarnings(cor(m$full_logFC, m$fold_logFC, method = "spearman"))
    inter <- length(intersect(deg_full, sel_genes))
    union_n <- length(union(deg_full, sel_genes))
    data.table(file = basename(f),
               n_sel = length(sel_genes),
               rho = round(rho, 4),
               jaccard = round(if (union_n > 0) inter / union_n else 0, 4),
               recovery_pct = round(if (length(deg_full) > 0) inter / length(deg_full) * 100 else 0, 2))
  }))
}
cpss_iter_summary <- iter_summary_one(cpss_files, "cpss")
boot_iter_summary <- iter_summary_one(boot_files, "boot")
if (nrow(cpss_iter_summary) > 0) cpss_iter_summary[, source := "cpss"]
if (nrow(boot_iter_summary) > 0) boot_iter_summary[, source := "bootstrap"]
iter_summary <- rbind(cpss_iter_summary, boot_iter_summary, fill = TRUE)
if (nrow(kfold_iter_summary) > 0) {
  kfold_iter_summary[, file := paste0("kfold_iter_", sprintf("%02d", iter), ".csv")]
  kfold_iter_summary[, source := paste0("kfold_K", K)]
  iter_summary <- rbind(iter_summary,
                         kfold_iter_summary[, .(file, n_sel, rho, jaccard, recovery_pct, source)],
                         fill = TRUE)
}

fwrite(iter_summary, file.path(OUT_DIR, paste0("pillar_A_iter_summary_", THR_TAG, ".csv")))
cat("Saved iter summary:", nrow(iter_summary), "rows\n")

if (nrow(kfold_iter_summary) > 0) {
  fwrite(kfold_iter_summary, file.path(OUT_DIR, paste0("pillar_A_kfold_summary_", THR_TAG, ".csv")))
  cat("Saved kfold summary\n")
}

# ============================================================================
# Master per-gene Pillar A table
# ============================================================================
pillar_A <- merge(cpss_dt, boot_dt, by = "gene", all = TRUE)
pillar_A <- merge(pillar_A, kfold_per_gene, by = "gene", all = TRUE)
pillar_A[is.na(kfold10_recur), kfold10_recur := 0L]
pillar_A[is.na(kfold50_recur), kfold50_recur := 0L]

# Annotate canonical DEG status from full dream
pillar_A <- merge(pillar_A,
                   full[, .(gene, dream_logFC = logFC, dream_padj = padj,
                            is_canonical_DEG = (padj < PADJ_THR & abs(logFC) > LFC_THR))],
                   by = "gene", all = TRUE)

fwrite(pillar_A, file.path(OUT_DIR, paste0("pillar_A_stability_", THR_TAG, ".csv")))
# Also write the canonical (lfc05) version under the legacy filename so the
# master atlas keeps loading without changes.
if (THR_TAG == "lfc05") fwrite(pillar_A, file.path(OUT_DIR, "pillar_A_stability.csv"))
cat("\nSaved per-gene Pillar A:", nrow(pillar_A), "rows ×", ncol(pillar_A), "columns\n")

# Quick reporting
deg <- pillar_A[is_canonical_DEG == TRUE]
if (nrow(deg) > 0) {
  cat("\n=== Headline numbers (canonical DEGs) ===\n")
  cat(sprintf("Mean CPSS π̂: %.3f\n", mean(deg$cpss_pi_hat, na.rm = TRUE)))
  cat(sprintf("CPSS π̂ ≥ 0.7: %d / %d (%.1f%%)\n",
              sum(deg$cpss_pi_hat >= 0.7, na.rm = TRUE), nrow(deg),
              100 * mean(deg$cpss_pi_hat >= 0.7, na.rm = TRUE)))
  cat(sprintf("Mean bootstrap freq: %.3f\n", mean(deg$bootstrap_freq, na.rm = TRUE)))
  cat(sprintf("Bootstrap freq ≥ 0.9: %d / %d (%.1f%%)\n",
              sum(deg$bootstrap_freq >= 0.9, na.rm = TRUE), nrow(deg),
              100 * mean(deg$bootstrap_freq >= 0.9, na.rm = TRUE)))
  cat(sprintf("Mean kfold10 recur: %.2f / 10\n", mean(deg$kfold10_recur, na.rm = TRUE)))
  cat(sprintf("Mean kfold50 recur: %.2f / 50\n", mean(deg$kfold50_recur, na.rm = TRUE)))
}

cat("\nDone Pillar A aggregation.\n")
