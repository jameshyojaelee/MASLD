# 253_sc_donor_heterogeneity.R
# F3 sub-state (sc) — Steatohepatitis donor-level heterogeneity in hepatocyte
# subtype proportions.
#
# Hypothesis: if F3a/F3b sub-states exist (mixing-fraction interpretation),
# Steatohepatitis donors should differ in their hepatocyte subtype distribution
#   F3a-like donors  → high Disease-Neutral, lower Disease-Progressor
#   F3b-like donors  → high Disease-Progressor, lower Disease-Neutral
# Test:
#   (a) Donor-level subtype-proportion variance vs cell-level (donor effect ANOVA)
#   (b) Donor-level "Progressor fraction" — bimodal across donors?
#   (c) Cluster donors by subtype-proportion vector; do donor clusters emerge?
#
# Inputs:
#   Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv
#   Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/meta_subtype_mapping.csv
#
# Outputs (results/granular_staging/):
#   f3_substate_sc_donor_proportions.csv      — per-donor subtype proportions
#   f3_substate_sc_donor_diagnostics.csv      — bimodality + heterogeneity tests
#   f3_substate_sc_donor_clusters.csv         — donor-level k-means clustering
#   f3_substate_sc_summary.md                  — narrative

suppressPackageStartupMessages({
  library(matrixStats)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HEP_DIR <- file.path(PROJECT_ROOT, "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== F3 sub-state (sc) — Steatohepatitis donor heterogeneity ==\n")

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
hep_meta_path <- file.path(HEP_DIR, "hepatocyte_subtype_metadata.csv")
mapping_path <- file.path(HEP_DIR, "meta_subtype_mapping.csv")

cat("Reading hepatocyte metadata (large file)...\n")
hep <- read.csv(hep_meta_path, stringsAsFactors = FALSE)
mapping <- read.csv(mapping_path, stringsAsFactors = FALSE)

cat(sprintf("Total cells: %d, Total subtypes: %d\n",
            nrow(hep), length(unique(hep$hepatocyte_subtype))))
print(table(hep$disease_stage_coarse, useNA = "ifany"))

# Map subtype number → meta_subtype
hep$meta_subtype <- mapping$meta_subtype[match(hep$hepatocyte_subtype, mapping$subtype)]
print(table(hep$meta_subtype, useNA = "ifany"))

# Filter to Steatohepatitis donors
sh <- hep[hep$disease_stage_coarse == "Steatohepatitis", ]
cat(sprintf("\nSteatohepatitis cells: %d\n", nrow(sh)))
cat(sprintf("Donors (samples): %d\n", length(unique(sh$sample))))

# ---------------------------------------------------------------------------
# (1) Per-donor subtype proportions
# ---------------------------------------------------------------------------
donor_props <- as.data.frame.matrix(
  prop.table(table(sh$sample, sh$meta_subtype), margin = 1)
)
donor_props$sample <- rownames(donor_props)
donor_props$dataset <- sh$dataset[match(donor_props$sample, sh$sample)]
donor_props$n_cells <- as.integer(table(sh$sample))[match(donor_props$sample, names(table(sh$sample)))]
donor_props <- donor_props[donor_props$n_cells >= 50, ]

cat(sprintf("Donors with ≥50 hepatocytes: %d\n", nrow(donor_props)))
print(head(donor_props))

write.csv(donor_props, file.path(OUT_DIR, "f3_substate_sc_donor_proportions.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# (2) Bimodality of donor-level "Progressor fraction"
# ---------------------------------------------------------------------------
progressor_col <- "Disease-Progressor"
da_col <- "Disease-Associated"
neutral_col <- "Disease-Neutral"

if (!progressor_col %in% colnames(donor_props)) progressor_col <- grep("Progressor", colnames(donor_props), value = TRUE)[1]
if (!da_col %in% colnames(donor_props)) da_col <- grep("Associated", colnames(donor_props), value = TRUE)[1]
if (!neutral_col %in% colnames(donor_props)) neutral_col <- grep("Neutral", colnames(donor_props), value = TRUE)[1]

prog_frac <- donor_props[[progressor_col]]
da_frac <- donor_props[[da_col]]
neutral_frac <- donor_props[[neutral_col]]

cat(sprintf("\nProgressor fraction: median=%.3f IQR=(%.3f, %.3f)\n",
            median(prog_frac), quantile(prog_frac, 0.25), quantile(prog_frac, 0.75)))

# 2-Gaussian mixture vs 1-Gaussian on Progressor fraction
fit_gmm_2 <- function(x, max_iter = 200, tol = 1e-6) {
  n <- length(x)
  ord <- sort(x)
  mu1 <- mean(ord[1:floor(n/2)]); mu2 <- mean(ord[(floor(n/2)+1):n])
  sd1 <- max(sd(ord[1:floor(n/2)]), 1e-3); sd2 <- max(sd(ord[(floor(n/2)+1):n]), 1e-3)
  p1 <- 0.5; ll_old <- -Inf
  for (iter in 1:max_iter) {
    d1 <- p1 * dnorm(x, mu1, sd1); d2 <- (1 - p1) * dnorm(x, mu2, sd2)
    g1 <- d1 / (d1 + d2 + 1e-300)
    p1 <- mean(g1)
    mu1 <- sum(g1 * x) / sum(g1); mu2 <- sum((1 - g1) * x) / sum(1 - g1)
    sd1 <- max(sqrt(sum(g1 * (x - mu1)^2) / sum(g1)), 1e-3)
    sd2 <- max(sqrt(sum((1 - g1) * (x - mu2)^2) / sum(1 - g1)), 1e-3)
    ll <- sum(log(p1 * dnorm(x, mu1, sd1) + (1 - p1) * dnorm(x, mu2, sd2) + 1e-300))
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }
  list(mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, p1 = p1, ll = ll)
}
hartigan_dip <- function(x, n_mc = 999) {
  obs_dip <- {
    grid <- seq(min(x), max(x), length.out = 200)
    e <- sapply(grid, function(g) mean(x <= g))
    g_pdf <- pnorm(grid, mean(x), sd(x))
    max(abs(e - g_pdf))
  }
  null_dips <- replicate(n_mc, {
    y <- rnorm(length(x), mean(x), sd(x))
    grid <- seq(min(y), max(y), length.out = 200)
    e <- sapply(grid, function(g) mean(y <= g))
    g_pdf <- pnorm(grid, mean(y), sd(y))
    max(abs(e - g_pdf))
  })
  list(dip_obs = obs_dip,
       p_value = (sum(null_dips >= obs_dip) + 1) / (n_mc + 1))
}

# Bimodality on multiple axes
diag_rows <- list()
for (axis_name in c("progressor_fraction", "disease_associated_fraction",
                    "disease_neutral_fraction", "f3b_minus_f3a_proxy")) {
  if (axis_name == "progressor_fraction") x <- prog_frac
  if (axis_name == "disease_associated_fraction") x <- da_frac
  if (axis_name == "disease_neutral_fraction") x <- neutral_frac
  # F3b-like proxy = Progressor + DA; F3a-like proxy = Disease-Neutral
  if (axis_name == "f3b_minus_f3a_proxy") x <- (prog_frac + da_frac) - neutral_frac

  if (length(x) < 20) next
  n <- length(x)
  ll1 <- sum(dnorm(x, mean(x), sd(x), log = TRUE))
  bic1 <- -2 * ll1 + 2 * log(n)
  fit2 <- fit_gmm_2(x)
  bic2 <- -2 * fit2$ll + 5 * log(n)
  bic_diff <- bic1 - bic2
  sep <- abs(fit2$mu1 - fit2$mu2) / sqrt((fit2$sd1^2 + fit2$sd2^2) / 2)
  dip <- hartigan_dip(x, n_mc = 999)
  diag_rows[[axis_name]] <- data.frame(
    axis = axis_name, n_donors = n,
    median = median(x), sd = sd(x), range = max(x) - min(x),
    BIC_diff_2vs1 = bic_diff, GMM_sep_sigma = sep, GMM_p1 = fit2$p1,
    dipMC_p = dip$p_value,
    stringsAsFactors = FALSE
  )
}
diag_df <- do.call(rbind, diag_rows)
write.csv(diag_df, file.path(OUT_DIR, "f3_substate_sc_donor_diagnostics.csv"),
          row.names = FALSE)
print(diag_df)

# ---------------------------------------------------------------------------
# (3) Donor-level k-means / hierarchical clustering
# ---------------------------------------------------------------------------
ms_cols <- intersect(c("Healthy", "Neutral", "Disease-Neutral",
                        "Disease-Associated", "Disease-Progressor"),
                      colnames(donor_props))
prop_mat <- as.matrix(donor_props[, ms_cols])
rownames(prop_mat) <- donor_props$sample

# Ward hierarchical clustering
d <- dist(prop_mat)
hc <- hclust(d, method = "ward.D2")
# k = 2 cut
cl2 <- cutree(hc, k = 2)
cl3 <- cutree(hc, k = 3)
donor_props$cluster_k2 <- cl2[match(donor_props$sample, names(cl2))]
donor_props$cluster_k3 <- cl3[match(donor_props$sample, names(cl3))]

# Mean profile per cluster
cluster_summary <- aggregate(prop_mat, by = list(cluster = cl2), FUN = mean)
cat("\nDonor cluster (k=2) mean subtype proportions:\n")
print(cluster_summary)

# Cluster size balance check
cluster_sizes <- table(cl2)
cat(sprintf("Cluster sizes: %s\n", paste(cluster_sizes, collapse = " / ")))
balance_pass <- min(cluster_sizes) / sum(cluster_sizes) >= 0.15

# Cluster-stage Spearman if the cells have a continuous severity proxy
# Average disease_stage_numeric per donor
donor_stage <- aggregate(sh$disease_stage_numeric,
                          by = list(sample = sh$sample), FUN = function(x) mean(x, na.rm = TRUE))
colnames(donor_stage) <- c("sample", "mean_stage")
donor_props <- merge(donor_props, donor_stage, by = "sample", all.x = TRUE)

cl_means <- aggregate(donor_props$mean_stage,
                       by = list(cluster = donor_props$cluster_k2),
                       FUN = function(x) mean(x, na.rm = TRUE))
cat("Cluster mean disease_stage_numeric:\n"); print(cl_means)

write.csv(donor_props, file.path(OUT_DIR, "f3_substate_sc_donor_clusters.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
sink(file.path(OUT_DIR, "f3_substate_sc_summary.md"))
cat("# F3 sub-state (sc) — Steatohepatitis donor heterogeneity\n\n")
cat(sprintf("- Steatohepatitis cells: %d\n", nrow(sh)))
cat(sprintf("- Donors with ≥50 cells: %d\n\n", nrow(donor_props)))

cat("## Bimodality tests on donor-level proportions\n\n")
cat("| Axis | n donors | median | range | ΔBIC | sep σ | dip p |\n|---|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(diag_df))) {
  cat(sprintf("| %s | %d | %.3f | %.3f | %.1f | %.2f | %.3f |\n",
              diag_df$axis[i], diag_df$n_donors[i], diag_df$median[i],
              diag_df$range[i], diag_df$BIC_diff_2vs1[i],
              diag_df$GMM_sep_sigma[i], diag_df$dipMC_p[i]))
}
cat("\n## Donor-level Ward clustering (k=2)\n\n")
cat(sprintf("Cluster sizes: %s\n\n", paste(cluster_sizes, collapse = " / ")))
cat(sprintf("Balance (min cluster fraction): %.2f %s\n\n",
            min(cluster_sizes) / sum(cluster_sizes),
            ifelse(balance_pass, "(PASS — biologically meaningful)",
                   "(FAIL — outlier-driven, like bulk)")))
cat("Cluster mean subtype proportions:\n\n")
cat("| Cluster |", paste(ms_cols, collapse = " | "), "|\n")
cat("|---|", paste(rep(":---:", length(ms_cols)), collapse = " | "), "|\n")
for (i in seq_len(nrow(cluster_summary))) {
  vals <- sprintf("%.2f", as.numeric(cluster_summary[i, ms_cols]))
  cat(sprintf("| %d | %s |\n", cluster_summary$cluster[i], paste(vals, collapse = " | ")))
}
sink()

cat("\n== F3 sub-state (sc) complete. ==\n")
