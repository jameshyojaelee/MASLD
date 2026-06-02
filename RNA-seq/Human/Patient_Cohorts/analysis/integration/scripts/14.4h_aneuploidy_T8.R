#!/usr/bin/env Rscript
# 14.4h — Sex chromosome aneuploidy detection (T8)
# Author: Team A, agent A9
# Date: 2026-05-11
#
# Goal: Compute Mahalanobis distance per sample from XIST + Y-marker ensemble
#       (DDX3Y + EIF1AY + RPS4Y1) log1p-CPM. Flag samples beyond chi-sq 2df
#       p<0.001 tail (d^2 > 13.816) and classify aneuploidy candidates by
#       quadrant (XXY, XO, low-quality).
#
# C2 revisions applied:
#   * Use POOLED within-cluster covariance (LDA convention) — not per-cluster
#     covariance — because the M cluster (52 samples in 5-cohort controls) is
#     too small for stable per-cluster cov estimation. Ledoit-Wolf shrinkage
#     applied as a safety belt.
#   * Y-marker = mean of log1p-CPM of DDX3Y + EIF1AY + RPS4Y1 (ensemble; reduces
#     single-gene dropout in 3' poly-A capture libraries).
#
# Inputs:
#   * RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_counts_raw.rds
#   * RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv
#   * RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv
#
# Outputs (RNA-seq/results/audit_sensitivity/sex_aneuploidy/):
#   * aneuploidy_flags.csv  (one row per sample with all metrics + flag)
#   * aneuploidy_summary.tsv (counts by candidate karyotype)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INTEG <- file.path(PROJ, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUT   <- file.path(PROJ, "RNA-seq/results/audit_sensitivity/sex_aneuploidy")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------- 1. Load counts and metadata ----------
counts_path <- file.path(INTEG, "results/integration/merged_counts_raw.rds")
meta_path   <- file.path(INTEG, "metadata/unified_metadata.csv")
qc_path     <- file.path(INTEG, "qc/sample_qc_report.csv")

stopifnot(file.exists(counts_path))
stopifnot(file.exists(meta_path))
stopifnot(file.exists(qc_path))

counts <- readRDS(counts_path)
if (is.list(counts) && !is.matrix(counts) && "counts" %in% names(counts)) {
  counts <- counts$counts
}
cat(sprintf("Count matrix loaded: %d genes x %d samples\n",
            nrow(counts), ncol(counts)))

meta <- fread(meta_path)
qc   <- fread(qc_path)

# Use only samples present in counts (technical pass not required — we want all
# samples to flag aneuploidy candidates including failed-QC ones)
shared <- intersect(colnames(counts), meta$sample_id)
counts <- counts[, shared]
meta   <- meta[match(shared, meta$sample_id)]
qc     <- qc[match(shared, qc$sample_id)]
cat(sprintf("Working set: %d samples after meta+counts intersect\n", length(shared)))

# ---------- 2. Locate sex-marker genes ----------
find_gene <- function(ensembl_base, all_genes) {
  hits <- grep(paste0("^", ensembl_base), all_genes, value = TRUE)
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

XIST_ENSG   <- "ENSG00000229807"
DDX3Y_ENSG  <- "ENSG00000067048"   # corrected DDX3Y (was 00000012817 = TSIX in some refs)
EIF1AY_ENSG <- "ENSG00000198692"
RPS4Y1_ENSG <- "ENSG00000129824"

all_genes <- rownames(counts)
xist_id   <- find_gene(XIST_ENSG,   all_genes)
ddx3y_id  <- find_gene(DDX3Y_ENSG,  all_genes)
eif1ay_id <- find_gene(EIF1AY_ENSG, all_genes)
rps4y1_id <- find_gene(RPS4Y1_ENSG, all_genes)

cat(sprintf("XIST=%s | DDX3Y=%s | EIF1AY=%s | RPS4Y1=%s\n",
            xist_id, ddx3y_id, eif1ay_id, rps4y1_id))

# Fallback DDX3Y id if first ENSG not found (different references use 00000012817)
if (is.na(ddx3y_id)) {
  ddx3y_id <- find_gene("ENSG00000012817", all_genes)
  cat(sprintf("DDX3Y fallback id: %s\n", ddx3y_id))
}

# ---------- 3. CPM-normalize ----------
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")
cpm_mat <- cpm(dge, log = FALSE)

get_lcpm <- function(gid) {
  if (is.na(gid)) return(rep(NA_real_, ncol(cpm_mat)))
  log1p(cpm_mat[gid, ])
}
xist_lcpm   <- get_lcpm(xist_id)
ddx3y_lcpm  <- get_lcpm(ddx3y_id)
eif1ay_lcpm <- get_lcpm(eif1ay_id)
rps4y1_lcpm <- get_lcpm(rps4y1_id)

# Y-marker ensemble: rowMeans of three Y genes, ignoring NAs
y_mat <- cbind(ddx3y_lcpm, eif1ay_lcpm, rps4y1_lcpm)
y_logcpm_mean <- rowMeans(y_mat, na.rm = TRUE)

# ---------- 4. k-means on (xist, y_mean) -> cluster identity ----------
sex_mat <- cbind(x = xist_lcpm, y = y_logcpm_mean)
sex_mat <- sex_mat[complete.cases(sex_mat), , drop = FALSE]
stopifnot(nrow(sex_mat) >= 100)

set.seed(42)
km <- kmeans(sex_mat, centers = 2, nstart = 25)
cluster_xist_means <- tapply(sex_mat[, "x"], km$cluster, mean)
female_cluster <- which.max(cluster_xist_means)
male_cluster   <- if (female_cluster == 1) 2 else 1

# Build per-sample assignment
sex_df <- data.table(
  sample_id = colnames(counts),
  xist_lcpm = xist_lcpm,
  y_logcpm_mean = y_logcpm_mean
)
sex_df[, complete := !is.na(xist_lcpm) & !is.na(y_logcpm_mean)]

assigned_sex <- rep(NA_character_, nrow(sex_df))
assigned_sex[sex_df$complete] <- ifelse(km$cluster == female_cluster, "F", "M")
sex_df[, assigned_sex := assigned_sex]

# ---------- 5. POOLED covariance (C2 fix) ----------
# Build per-cluster mu and pooled Sigma (LDA convention).
mu_F <- colMeans(sex_mat[km$cluster == female_cluster, , drop = FALSE])
mu_M <- colMeans(sex_mat[km$cluster == male_cluster,   , drop = FALSE])
n_F  <- sum(km$cluster == female_cluster)
n_M  <- sum(km$cluster == male_cluster)
S_F  <- cov(sex_mat[km$cluster == female_cluster, , drop = FALSE])
S_M  <- cov(sex_mat[km$cluster == male_cluster,   , drop = FALSE])
Sigma_pooled <- ((n_F - 1) * S_F + (n_M - 1) * S_M) / (n_F + n_M - 2)

# Ledoit-Wolf-style shrinkage to ensure invertibility:
# regularize toward diag(variances)
diag_target <- diag(diag(Sigma_pooled))
lambda <- 0.05  # mild shrinkage; pooled is already stable with N~1400
Sigma_reg <- (1 - lambda) * Sigma_pooled + lambda * diag_target

cat("Pooled covariance:\n"); print(round(Sigma_pooled, 4))
cat(sprintf("n_F=%d n_M=%d ; pooled cov shrunk with lambda=%.2f\n",
            n_F, n_M, lambda))

Sigma_inv <- solve(Sigma_reg)

# ---------- 6. Mahalanobis distance to assigned centroid ----------
mahal_d2 <- function(z, mu) {
  if (any(is.na(z))) return(NA_real_)
  d <- z - mu
  as.numeric(d %*% Sigma_inv %*% d)
}

z_mat <- as.matrix(sex_df[, .(xist_lcpm, y_logcpm_mean)])
colnames(z_mat) <- c("x","y")

d2_F <- apply(z_mat, 1, function(z) mahal_d2(z, mu_F))
d2_M <- apply(z_mat, 1, function(z) mahal_d2(z, mu_M))

d2_assigned <- ifelse(sex_df$assigned_sex == "F", d2_F, d2_M)
sex_df[, mahal_d2_assigned := d2_assigned]
sex_df[, mahal_d2_F := d2_F]
sex_df[, mahal_d2_M := d2_M]

# Logistic posterior in squared-distance gap
sex_df[, posterior_F := 1 / (1 + exp(d2_F - d2_M))]

# ---------- 7. Flag aneuploidy ----------
# Threshold from chi-sq 2 df at p<0.001:
TH_D2 <- 13.816

# Quadrant in z-score space (median-centered + MAD-scaled)
x_med <- median(sex_df$xist_lcpm,    na.rm = TRUE)
y_med <- median(sex_df$y_logcpm_mean, na.rm = TRUE)
x_mad <- mad(sex_df$xist_lcpm,    na.rm = TRUE) + 1e-6
y_mad <- mad(sex_df$y_logcpm_mean, na.rm = TRUE) + 1e-6
sex_df[, x_z := (xist_lcpm - x_med) / x_mad]
sex_df[, y_z := (y_logcpm_mean - y_med) / y_mad]

# Per-cluster centroid distances from medians (used to detect "between-cluster"
# samples)
# Quadrant decision: use cluster centroids + 1.5 SD shells
# - XIST_high & Y_high   -> XXY candidate (Klinefelter; ~0.1-0.2% prev)
# - XIST_low  & Y_high   -> typical M
# - XIST_high & Y_low    -> typical F (also catches XO carrier with reactivated XIST??)
# - XIST_low  & Y_low    -> XO_or_low_QC (Turner-like or library failure)
x_th_hi <- 0.5 * (mu_F["x"] + mu_M["x"])    # midpoint between clusters
y_th_hi <- 0.5 * (mu_F["y"] + mu_M["y"])
sex_df[, x_high := xist_lcpm >= x_th_hi]
sex_df[, y_high := y_logcpm_mean >= y_th_hi]

sex_df[, quadrant := fcase(
  is.na(x_high) | is.na(y_high), "missing",
  x_high == TRUE  & y_high == TRUE,  "XIST_high_Y_high",
  x_high == FALSE & y_high == TRUE,  "XIST_low_Y_high",
  x_high == TRUE  & y_high == FALSE, "XIST_high_Y_low",
  x_high == FALSE & y_high == FALSE, "XIST_low_Y_low",
  default = "ambiguous"
)]

sex_df[, flag_aneuploidy := !is.na(mahal_d2_assigned) & mahal_d2_assigned > TH_D2]

sex_df[, candidate_karyotype := fcase(
  !flag_aneuploidy, "typical",
  flag_aneuploidy & quadrant == "XIST_high_Y_high", "XXY_candidate",
  flag_aneuploidy & quadrant == "XIST_low_Y_low",   "XO_or_lowQC_candidate",
  flag_aneuploidy & quadrant == "XIST_high_Y_low" & assigned_sex == "M",
                                                     "XX_in_M_cluster_outlier",
  flag_aneuploidy & quadrant == "XIST_low_Y_high" & assigned_sex == "F",
                                                     "XY_in_F_cluster_outlier",
  flag_aneuploidy, "outlier_within_quadrant",
  default = "typical"
)]

# ---------- 8. Cross-check with pass_sex ----------
sex_df <- merge(sex_df, qc[, .(sample_id, pass_sex, pass_pca, pass_libsize,
                               total_counts, pass_technical, pass_all)],
                by = "sample_id", all.x = TRUE)
sex_df <- merge(sex_df, meta[, .(sample_id, dataset)],
                by = "sample_id", all.x = TRUE)

# Coincidence of aneuploidy flag with pass_sex==FALSE
overlap <- sex_df[, .N, by = .(flag_aneuploidy, pass_sex)]
cat("\n=== flag_aneuploidy x pass_sex ===\n")
print(overlap)

cat("\n=== candidate_karyotype counts ===\n")
print(sex_df[, .N, by = candidate_karyotype][order(-N)])

cat("\n=== Flagged samples (showing top 30) ===\n")
flagged <- sex_df[flag_aneuploidy == TRUE,
                  .(sample_id, dataset, xist_lcpm, y_logcpm_mean,
                    mahal_d2_assigned, candidate_karyotype, pass_sex,
                    pass_technical)]
print(flagged[order(-mahal_d2_assigned)][1:min(30, nrow(flagged))])

# ---------- 9. Write outputs ----------
out_path <- file.path(OUT, "aneuploidy_flags.csv")
fwrite(sex_df[, .(sample_id, dataset, xist_lcpm, y_logcpm_mean,
                  mahal_d2_assigned, mahal_d2_F, mahal_d2_M, posterior_F,
                  assigned_sex, quadrant, flag_aneuploidy, candidate_karyotype,
                  pass_sex, pass_technical, pass_all, total_counts)],
       out_path)
cat(sprintf("\nWritten: %s (%d rows)\n", out_path, nrow(sex_df)))

summary_dt <- sex_df[, .(
  n_samples = .N,
  n_flagged = sum(flag_aneuploidy, na.rm = TRUE),
  n_pass_sex_false = sum(!pass_sex, na.rm = TRUE),
  n_flag_and_pass_sex_false = sum(flag_aneuploidy & !pass_sex, na.rm = TRUE)
), by = candidate_karyotype]
fwrite(summary_dt, file.path(OUT, "aneuploidy_summary.tsv"), sep = "\t")
cat(sprintf("Written: %s\n", file.path(OUT, "aneuploidy_summary.tsv")))

# Final stats
cat(sprintf("\n=== Final ===\n"))
cat(sprintf("Total samples processed: %d\n", nrow(sex_df)))
cat(sprintf("Flagged (d^2 > %.3f): %d (%.2f%%)\n", TH_D2,
            sum(sex_df$flag_aneuploidy, na.rm = TRUE),
            100 * mean(sex_df$flag_aneuploidy, na.rm = TRUE)))
cat(sprintf("Of flagged, n with pass_sex==FALSE: %d\n",
            sum(sex_df$flag_aneuploidy & !sex_df$pass_sex, na.rm = TRUE)))

cat("\n=== T8 complete ===\n")
