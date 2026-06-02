# 253b_sc_donor_cluster_annotation.R
# F3 sub-state (sc) addendum — annotate the 35/26 donor cluster split with
# axis scores (lipid / inflammatory / ferroptosis) to distinguish:
#   - Severity interpretation: cluster 1 = more progressed (higher inflammatory)
#   - F3a/F3b interpretation: cluster 1 = F3b-like (higher ferroptosis), cluster 2
#     = F3a-like (higher lipid)

suppressPackageStartupMessages({})

PROJECT_ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")

clusters <- read.csv(file.path(OUT_DIR, "f3_substate_sc_donor_clusters.csv"),
                      stringsAsFactors = FALSE)
axes <- read.csv(file.path(OUT_DIR, "donor_axis_scores.csv"),
                  stringsAsFactors = FALSE)

merged <- merge(clusters, axes, by = "sample", all.x = TRUE)
cat(sprintf("Donors with cluster + axis: %d\n", nrow(merged)))

# Per-cluster mean axis values
agg <- aggregate(merged[, c("axis_inflammatory_avg", "axis_lipid_avg",
                            "axis_ferroptosis_avg")],
                  by = list(cluster_k2 = merged$cluster_k2),
                  FUN = function(x) mean(x, na.rm = TRUE))
cat("Cluster mean axis scores:\n")
print(agg)

# Tests: cluster 1 vs cluster 2
test_axis <- function(axis_name) {
  v1 <- merged[[axis_name]][merged$cluster_k2 == 1]
  v2 <- merged[[axis_name]][merged$cluster_k2 == 2]
  v1 <- v1[!is.na(v1)]; v2 <- v2[!is.na(v2)]
  if (length(v1) < 3 || length(v2) < 3) return(c(mean1 = NA, mean2 = NA, p = NA))
  wt <- wilcox.test(v1, v2)
  c(mean1 = mean(v1), mean2 = mean(v2),
    diff = mean(v1) - mean(v2),
    cohens_d = (mean(v1) - mean(v2)) / sqrt((var(v1) + var(v2)) / 2),
    p = wt$p.value)
}

results <- sapply(c("axis_inflammatory_avg", "axis_lipid_avg", "axis_ferroptosis_avg"),
                  test_axis)
res_df <- as.data.frame(t(results))
res_df$axis <- rownames(res_df)
res_df <- res_df[, c("axis", "mean1", "mean2", "diff", "cohens_d", "p")]
res_df$padj_bh <- p.adjust(res_df$p, method = "BH")
print(res_df)

write.csv(res_df, file.path(OUT_DIR, "f3_substate_sc_cluster_axis_annotation.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Interpret
# ---------------------------------------------------------------------------
sink(file.path(OUT_DIR, "f3_substate_sc_cluster_interpretation.md"))
cat("# Donor cluster interpretation: donor-level severity heterogeneity (orthogonal to F3a/F3b)\n\n")
cat("Cluster k=2 split: 35 (Progressor-dominant) / 26 (Associated/Neutral-dominant).\n\n")
cat("## Per-axis comparison\n\n")
cat("| Axis | Cluster 1 mean | Cluster 2 mean | Cohen's d | p (Wilcox) | BH padj |\n")
cat("|---|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(res_df))) {
  cat(sprintf("| %s | %.3f | %.3f | %.2f | %.2e | %.2e |\n",
              res_df$axis[i], res_df$mean1[i], res_df$mean2[i],
              res_df$cohens_d[i], res_df$p[i], res_df$padj_bh[i]))
}

cat("\n## Decision rule\n\n")
inf_d <- res_df$cohens_d[res_df$axis == "axis_inflammatory_avg"]
lip_d <- res_df$cohens_d[res_df$axis == "axis_lipid_avg"]
fer_d <- res_df$cohens_d[res_df$axis == "axis_ferroptosis_avg"]

cat(sprintf("- inflammatory: cluster1 - cluster2 = %.2f σ (positive = cluster1 higher)\n", inf_d))
cat(sprintf("- lipid:        cluster1 - cluster2 = %.2f σ\n", lip_d))
cat(sprintf("- ferroptosis:  cluster1 - cluster2 = %.2f σ\n\n", fer_d))

cat("**Interpretation guide:**\n\n")
cat("- Severity: Cluster 1 (Progressor-dominant) higher in BOTH inflammatory AND ferroptosis (more progressed)\n")
cat("- F3a/F3b: Cluster 2 has higher LIPID (F3a UPR/SREBP biology); Cluster 1 has higher ferroptosis but NOT inflammation (F3b ECM/senescence)\n")
cat("- Mixed: results don't cleanly fit either model\n\n")

verdict <- if (inf_d > 0.5 && fer_d > 0.5) {
  "SEVERITY interpretation likely (cluster 1 more progressed across both inflammation & ferroptosis)."
} else if (inf_d < 0.5 && lip_d < -0.5 && fer_d > 0.5) {
  "F3a/F3b interpretation supported (cluster 2 = F3a-like lipid-rich; cluster 1 = F3b-like ferroptosis-rich)."
} else {
  "MIXED — neither severity nor F3a/F3b interpretation fully fits. Caveat in discussion."
}
cat("**Verdict:** ", verdict, "\n", sep = "")
sink()

cat("Wrote f3_substate_sc_cluster_interpretation.md\n")
