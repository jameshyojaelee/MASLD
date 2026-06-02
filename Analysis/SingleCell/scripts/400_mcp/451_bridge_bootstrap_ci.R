#!/usr/bin/env Rscript
# 451_bridge_bootstrap_ci.R — bootstrap CI on bulk↔sc bridge correlations.
# Addresses adversarial reviewer Item 2: are the max|r| anchors statistically
# distinct from the scRNA-unique threshold r=0.30?
#
# Two bootstrap modes (BOOTSTRAP_LEVEL env var):
#
#   gene  (default, original): resample genes with replacement over the
#         intersecting gene set; recompute Pearson r between cNMF spectra
#         and bulk W. Tests "is r gene-stable?" but NOT donor-generalizable.
#
#   donor (added 2026-04-25, adversarial-#3 closure): resample bulk donors
#         (samples) with replacement; for each resampled donor set, recompute
#         a per-gene "bulk program signature" as the donor-weighted mean
#         expression weighted by H[k, donor]. This treats the donor cohort
#         as the unit of observation and asks "does r generalize across
#         patients?". Implementation:
#           sig[g, k] = sum_d (mat[g, d] * H[k, d]) / sum_d H[k, d]
#         which is a weighted mean of expression where the weight is the
#         sample's program loading (a standard "soft cohort centroid").
#         Bootstrapping donors with replacement perturbs both the weights
#         and the per-gene signal in a biologically meaningful way.
#
# Output: bulk_sc_bridge_bootstrap_ci.tsv (gene-mode, default)
#         bulk_sc_bridge_bootstrap_ci_donor.tsv (donor-mode)
#
# Env:
#   NMF_BOOT_N        — number of bootstrap iterations (default 1000)
#   NMF_BOOT_SEED     — RNG seed (default 42)
#   BOOTSTRAP_LEVEL   — "gene" (default) or "donor"

suppressPackageStartupMessages({
  library(data.table); library(NMF)
})

root <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
out_dir <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/integration")
B <- as.integer(Sys.getenv("NMF_BOOT_N", "1000"))
SEED <- as.integer(Sys.getenv("NMF_BOOT_SEED", "42"))
SCRNA_K <- 16
THRESH <- 0.30
BOOT_LEVEL <- Sys.getenv("BOOTSTRAP_LEVEL", "gene")
stopifnot(BOOT_LEVEL %in% c("gene", "donor"))
cat(sprintf("[451] BOOTSTRAP_LEVEL = %s, B = %d, SEED = %d\n", BOOT_LEVEL, B, SEED))

# cNMF gene-spectra score for scRNA k=16 (same file used by 450_bulk_sc_bridge.R)
spectra_f <- file.path(root,
  sprintf("Analysis/SingleCell/results_gpu_v2/mcp/cnmf_runs/global/global.gene_spectra_score.k_%d.dt_0_03.txt",
          SCRNA_K))
stopifnot("spectra file missing" = file.exists(spectra_f))
spectra <- as.matrix(fread(spectra_f), rownames = 1)  # programs x genes (symbols)

# Ensembl mapping
gm <- fread(file.path(root, "data/gencode_v49_gene_metadata.tsv.gz"))
sym2ens <- setNames(gm$ensembl_base, gm$gene_name)
keep <- colnames(spectra) %in% names(sym2ens) & !is.na(sym2ens[colnames(spectra)])
spectra <- spectra[, keep]
colnames(spectra) <- sym2ens[colnames(spectra)]

# Bulk-cache variants to bootstrap
variants <- list(
  list(tag = "canonical",
       cache = file.path(root, "RNA-seq/results/subtypes/nmf_results_cache_clean.rds"),
       labels = file.path(root, "RNA-seq/results/subtypes/program_labels.csv"),
       k = 6L),
  list(tag = "protonly_k6",
       cache = file.path(root, "RNA-seq/results/subtypes/nmf_results_cache_clean_protonly.rds"),
       labels = file.path(root, "RNA-seq/results/subtypes/program_labels_protonly.csv"),
       k = 6L),
  list(tag = "protonly_k8",
       cache = file.path(root, "RNA-seq/results/subtypes/nmf_results_cache_clean_protonly.rds"),
       labels = file.path(root, "RNA-seq/results/subtypes/program_labels_protonly_k8.csv"),
       k = 8L),
  list(tag = "nonprotonly_k6",
       cache = file.path(root, "RNA-seq/results/subtypes/nmf_results_cache_clean_nonprotonly.rds"),
       labels = file.path(root, "RNA-seq/results/subtypes/program_labels_nonprotonly.csv"),
       k = 6L)
)

load_labels <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path); lines <- lines[!grepl("^#", lines)]
  fread(text = paste(lines, collapse = "\n"))
}

boot_one_variant <- function(v) {
  cat(sprintf("\n--- variant %s (k=%d) ---\n", v$tag, v$k))
  bulk <- readRDS(v$cache)
  fit <- bulk$nmf_results[[as.character(v$k)]]
  stopifnot(!is.null(fit))
  W_bulk <- basis(fit)
  H_bulk <- coef(fit)  # programs x donors
  rownames(W_bulk) <- sub("\\..*$", "", rownames(W_bulk))
  labels <- load_labels(v$labels)
  lbl_col <- if (!is.null(labels) && "biological_label" %in% names(labels))
               "biological_label" else NULL
  if (!is.null(lbl_col) && nrow(labels) == ncol(W_bulk)) {
    prog_names <- paste0(LETTERS[seq_len(ncol(W_bulk))], "=",
                         labels[[lbl_col]][seq_len(ncol(W_bulk))])
    colnames(W_bulk) <- prog_names
    rownames(H_bulk) <- prog_names
  } else {
    colnames(W_bulk) <- LETTERS[seq_len(ncol(W_bulk))]
    rownames(H_bulk) <- LETTERS[seq_len(ncol(W_bulk))]
  }

  common <- intersect(colnames(spectra), rownames(W_bulk))
  cat(sprintf("  common genes: %d\n", length(common)))
  Wc <- W_bulk[common, , drop = FALSE]
  Sc <- spectra[, common, drop = FALSE]  # programs x genes

  # Point estimate using the canonical W-vs-spectra correlation
  R_point <- cor(t(Sc), Wc, method = "pearson", use = "pairwise.complete.obs")
  point_max_abs <- apply(R_point, 1, function(r) max(abs(r)))
  point_anchor_idx <- apply(R_point, 1, function(r) which.max(abs(r)))
  point_anchor_name <- colnames(Wc)[point_anchor_idx]

  set.seed(SEED)

  if (BOOT_LEVEL == "gene") {
    # Gene-resampling bootstrap (original)
    G <- length(common)
    max_abs_boot <- matrix(NA_real_, nrow = nrow(Sc), ncol = B)
    for (b in seq_len(B)) {
      idx <- sample.int(G, G, replace = TRUE)
      Rb <- cor(t(Sc[, idx, drop = FALSE]), Wc[idx, , drop = FALSE],
                method = "pearson", use = "pairwise.complete.obs")
      max_abs_boot[, b] <- apply(Rb, 1, function(r) max(abs(r)))
    }
    boot_method_str <- "gene"
    n_resampled_unit <- length(common)
  } else {
    # Donor-resampling bootstrap. Because W (NMF basis, genes x prog) is fixed
    # across donors, we cannot directly bootstrap W. Instead we bootstrap the
    # _donor-weighted-mean expression signature_: for each program k, define
    #   sig[g, k] = sum_d mat[g, d] * H[k, d] / sum_d H[k, d]
    # which is the donor-cohort mean expression weighted by the sample's
    # program loading (a soft cohort centroid). We then correlate this
    # signature against the cNMF spectra. Resampling donors with replacement
    # perturbs sig in a way that asks "does the bridge correlation generalize
    # across patients?"
    #
    # NOTE: signature-mode point estimates are systematically lower than the
    # canonical W-mode point estimates (W is concentrated, sig is diffuse;
    # see `point_max_abs_r_signature` for the donor-mode reference statistic).
    # The two regimes test different but complementary stability questions.
    mat <- bulk$mat
    if (is.null(mat)) stop("bulk cache missing $mat for donor bootstrap")
    rownames(mat) <- sub("\\..*$", "", rownames(mat))
    mat_c <- mat[common, , drop = FALSE]
    H_c <- H_bulk  # programs x donors
    stopifnot(all(colnames(H_c) == colnames(mat_c)) || ncol(H_c) == ncol(mat_c))
    if (is.null(colnames(H_c))) colnames(H_c) <- colnames(mat_c)
    n_donors <- ncol(mat_c)
    cat(sprintf("  bulk donors: %d  (donor-resampling)\n", n_donors))

    # Signature-mode point estimate (full cohort) — different from W-mode
    sig_full <- mat_c %*% t(H_c)
    H_sum_full <- rowSums(H_c)
    H_sum_full[H_sum_full == 0] <- NA_real_
    sig_full <- sweep(sig_full, 2, H_sum_full, "/")
    colnames(sig_full) <- rownames(H_c)
    R_sig_point <- cor(t(Sc), sig_full, method = "pearson",
                       use = "pairwise.complete.obs")
    point_max_abs_sig <- apply(R_sig_point, 1, function(r) max(abs(r), na.rm = TRUE))

    max_abs_boot <- matrix(NA_real_, nrow = nrow(Sc), ncol = B)
    for (b in seq_len(B)) {
      idx <- sample.int(n_donors, n_donors, replace = TRUE)
      mat_b <- mat_c[, idx, drop = FALSE]
      H_b <- H_c[, idx, drop = FALSE]
      sig_b <- (mat_b %*% t(H_b))
      H_sum <- rowSums(H_b)
      H_sum[H_sum == 0] <- NA_real_
      sig_b <- sweep(sig_b, 2, H_sum, "/")
      colnames(sig_b) <- rownames(H_c)
      Rb <- cor(t(Sc), sig_b, method = "pearson",
                use = "pairwise.complete.obs")
      max_abs_boot[, b] <- apply(Rb, 1, function(r) max(abs(r), na.rm = TRUE))
    }
    boot_method_str <- "donor"
    n_resampled_unit <- n_donors
  }

  cat(sprintf("  cor any-NA=%s  max|r| range=[%.3f,%.3f]\n",
              any(is.na(R_point)), min(point_max_abs), max(point_max_abs)))
  out <- data.table(
    variant = v$tag,
    bulk_k = v$k,
    cnmf_program = seq_len(nrow(Sc)),
    point_max_abs_r = round(point_max_abs, 4),
    point_anchor = point_anchor_name,
    boot_ci_lo  = round(apply(max_abs_boot, 1, quantile, 0.025, na.rm = TRUE), 4),
    boot_ci_hi  = round(apply(max_abs_boot, 1, quantile, 0.975, na.rm = TRUE), 4),
    boot_mean   = round(rowMeans(max_abs_boot, na.rm = TRUE), 4),
    p_above_thresh = round(rowMeans(max_abs_boot >= THRESH, na.rm = TRUE), 4),
    scRNA_unique_ci = apply(max_abs_boot, 1,
      function(x) quantile(x, 0.975, na.rm = TRUE) < THRESH),
    n_common_genes = length(common),
    bootstrap_method = boot_method_str,
    n_resampled_unit = n_resampled_unit
  )
  if (BOOT_LEVEL == "donor") {
    # Add the signature-mode point estimate so reviewers can compare apples-to-apples
    out[, point_max_abs_r_signature := round(point_max_abs_sig, 4)]
  }
  out
}

results <- rbindlist(lapply(variants, boot_one_variant))
out_basename <- if (BOOT_LEVEL == "donor")
                  "bulk_sc_bridge_bootstrap_ci_donor.tsv" else
                  "bulk_sc_bridge_bootstrap_ci.tsv"
out_path <- file.path(out_dir, out_basename)
fwrite(results, out_path, sep = "\t")
cat("\nWrote:", out_path, "\n\n")

# Summary per variant
cat("=== Bootstrap summary (B =", B, ") ===\n")
for (v in unique(results$variant)) {
  sub <- results[variant == v]
  cat(sprintf("  %s: %d/%d scRNA-unique by point estimate (max|r| < %.2f); %d/%d by CI upper < %.2f\n",
              v,
              sum(sub$point_max_abs_r < THRESH), nrow(sub),
              THRESH,
              sum(sub$scRNA_unique_ci), nrow(sub),
              THRESH))
}

# Rotten-anchor focus: cnmf P1 and P4
cat("\n=== Rotten-anchor bootstrap CIs (cNMF programs 1, 4) ===\n")
print(results[cnmf_program %in% c(1L, 4L)], row.names = FALSE)
