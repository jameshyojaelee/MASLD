#!/usr/bin/env Rscript
# =============================================================================
# B3 NMF seed sensitivity (V5 verification)
# Re-runs NMF k=2 with 3 seeds, computes S2 fraction at F0-F2 vs F3-F4,
# and checks if doubling holds.
#
# Output: verification/controls/b3_nmf_seed_sensitivity.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(edgeR)
  library(limma)
  library(NMF)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
dge_path  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
meta_path <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
qc_path   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv")
model_meta_path <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier/modeling_metadata.csv")
nmf_orig_path <- file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv")

OUT_CSV <- file.path(BASE, "docs/manuscript/verification/controls/b3_nmf_seed_sensitivity.csv")
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

cat("=== B3 NMF seed sensitivity ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ---- Load ----
dge <- readRDS(dge_path)
meta <- fread(meta_path)
qc <- fread(qc_path)
nmf_orig <- fread(nmf_orig_path)  # columns: sample_id, subtype

# unified_metadata.csv already has fibrosis_stage — no need to join model_meta
meta <- meta %>%
  left_join(qc %>% select(sample_id, pass_technical), by = "sample_id")

disease_samples <- meta %>%
  filter(group_binary == "Disease" & pass_technical == TRUE) %>%
  pull(sample_id)
available <- intersect(disease_samples, colnames(dge))
cat("Disease QC-pass samples in DGE:", length(available), "\n")

dge_sub <- dge[, available]
meta_sub <- meta %>% filter(sample_id %in% available)

# ---- Normalize + batch-correct + select top 5000 ----
dge_sub <- calcNormFactors(dge_sub, method = "TMM")
v <- voom(dge_sub, design = NULL, plot = FALSE)
logcpm <- v$E
logcpm_corr <- removeBatchEffect(logcpm, batch = meta_sub$dataset)
gene_iqr <- apply(logcpm_corr, 1, IQR)
top_genes <- names(sort(gene_iqr, decreasing = TRUE))[1:5000]
mat <- logcpm_corr[top_genes, ]
mat_nn <- mat - apply(mat, 1, min)
mat_nn[mat_nn < 0] <- 0
cat("Matrix:", dim(mat_nn), "\n")

# ---- Identify S2 in original assignment for consistent labeling ----
# Use mean expression of top-var genes: subtype with more members in F3-F4 is "S2"
nmf_orig_stage <- nmf_orig %>%
  inner_join(meta_sub %>% select(sample_id, fibrosis_stage), by = "sample_id")

orig_stage_cross <- nmf_orig_stage %>%
  filter(!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4) %>%
  mutate(stage_group = ifelse(fibrosis_stage <= 2, "F0_F2", "F3_F4")) %>%
  count(nmf_subtype, stage_group)
cat("\nOriginal S1/S2 × stage_group:\n"); print(orig_stage_cross)

# Reference S2 fraction
orig_wide <- orig_stage_cross %>%
  tidyr::pivot_wider(names_from = stage_group, values_from = n, values_fill = 0L) %>%
  mutate(total = F0_F2 + F3_F4)
cat("\nReference (from nmf_assignments.csv):\n"); print(orig_wide)

# ---- Run NMF with 3 seeds ----
results <- list()
for (seed in c(42L, 2026L, 17L)) {
  cat(sprintf("\n--- NMF k=2 seed=%d ---\n", seed))
  set.seed(seed)
  # use single-seed integer so brunet is deterministic per seed
  res <- nmf(mat_nn, rank = 2, nrun = 1, method = "brunet",
             seed = seed, .options = list(verbose = FALSE))
  cluster <- predict(res)  # factor 1/2 per sample
  cluster_char <- as.character(cluster)

  # Map clusters to S1/S2 using concordance with original assignment
  # (S2 = cluster with higher F3-F4 fraction)
  sample_ids <- colnames(mat_nn)
  dt <- data.table(sample_id = sample_ids, cluster = cluster_char)
  dt <- merge(dt, as.data.table(meta_sub)[, .(sample_id, fibrosis_stage)],
              by = "sample_id")
  dt <- dt[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
  dt[, stage_group := ifelse(fibrosis_stage <= 2, "F0_F2", "F3_F4")]

  # Which cluster has higher F3-F4 share?
  cluster_stats <- dt[, .(n_F0F2 = sum(stage_group == "F0_F2"),
                          n_F3F4 = sum(stage_group == "F3_F4")), by = cluster]
  cluster_stats[, frac_F3F4 := n_F3F4 / (n_F0F2 + n_F3F4)]
  cat("Cluster stats:\n"); print(cluster_stats)

  s2_cluster <- cluster_stats[which.max(frac_F3F4), cluster]
  dt[, subtype := ifelse(cluster == s2_cluster, "S2", "S1")]

  # S2 fractions
  stage_tab <- dt[, .N, by = .(fibrosis_stage, subtype)]
  stage_wide <- dcast(stage_tab, fibrosis_stage ~ subtype, value.var = "N", fill = 0)
  stage_wide[, total := S1 + S2][, pct_S2 := 100 * S2 / total]
  cat("Per-stage S2 pct:\n"); print(stage_wide)

  s2_F0F2 <- dt[stage_group == "F0_F2", mean(subtype == "S2")] * 100
  s2_F3F4 <- dt[stage_group == "F3_F4", mean(subtype == "S2")] * 100
  ratio <- s2_F3F4 / s2_F0F2

  # Concordance with original S2
  orig_map <- setNames(nmf_orig$nmf_subtype, nmf_orig$sample_id)
  shared <- intersect(sample_ids, names(orig_map))
  conc <- mean(dt[sample_id %in% shared]$subtype ==
               orig_map[dt[sample_id %in% shared]$sample_id])

  results[[as.character(seed)]] <- data.table(
    seed = seed,
    n_samples = nrow(dt),
    s2_pct_F0_F2 = s2_F0F2,
    s2_pct_F3_F4 = s2_F3F4,
    ratio_F3F4_over_F0F2 = ratio,
    doubling_holds = (s2_F0F2 < 30 & s2_F3F4 > 40),
    concordance_with_orig_s2 = conc,
    cophenetic = tryCatch(cophcor(res), error = function(e) NA_real_)
  )
}

result_dt <- rbindlist(results)
cat("\n\n=== SUMMARY ===\n")
print(result_dt)

fwrite(result_dt, OUT_CSV)
cat("\nSaved:", OUT_CSV, "\n")

cat("\nEnd:", format(Sys.time()), "\n")
