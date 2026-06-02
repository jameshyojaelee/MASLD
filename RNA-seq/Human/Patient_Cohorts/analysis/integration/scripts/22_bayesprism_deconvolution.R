#!/usr/bin/env Rscript
# 22_bayesprism_deconvolution.R
# ---------------------------------------------------------------------------
# BayesPrism Deconvolution Concordance & Attribution Validation
#
# Purpose: Validate the 83% hepatocyte-intrinsic finding (from Script 25/MuSiC)
#          using BayesPrism as an independent second deconvolution method.
#
# Strategy:
#   1. Load per-cohort BayesPrism proportions (already computed)
#   2. Merge into unified matrix matching the 873 mega-analysis samples
#   3. Compare BayesPrism vs MuSiC fractions (target: Pearson r > 0.7)
#   4. Run deconvolution-adjusted dream with BayesPrism fractions
#   5. Classify genes as Hepatocyte_intrinsic / Composition_driven / Unmasked
#   6. Compare BayesPrism-based attribution with MuSiC-based attribution
#
# Inputs:
#   - Per-cohort BayesPrism proportions (Analysis/Deconvolution/results/)
#   - Per-cohort MuSiC proportions (same directory)
#   - merged_dge.rds, meta_matched.rds (from integration pipeline)
#   - dream_results.csv (unadjusted, from Script 05)
#
# Outputs:
#   - results/deconvolution/bayesprism/unified_bayesprism_proportions.csv
#   - results/deconvolution/bayesprism/method_concordance.csv
#   - results/deconvolution/bayesprism/bayesprism_attribution_scores.csv
#   - results/deconvolution/bayesprism/dream_results_bayesprism_adjusted.csv
#   - figures/bayesprism_concordance.pdf
#   - figures/bayesprism_attribution_comparison.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(yaml)
  library(edgeR)
})

# Force injection into lme4 namespace BEFORE loading variancePartition
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}

suppressPackageStartupMessages({
  library(variancePartition)
  library(BiocParallel)
  library(ggplot2)
  library(gridExtra)
})

cat("=== Script 22: BayesPrism Deconvolution Concordance & Attribution ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ============================================================
#  Paths
# ============================================================
BASE     <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT_DIR  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR     <- file.path(INT_DIR, "results/integration")
DECONV   <- file.path(BASE, "Analysis/Deconvolution/results")
OUT      <- file.path(INT_DIR, "results/deconvolution/bayesprism")
FIGDIR   <- file.path(BASE, "figures")
CAUSAL_DIR <- file.path(BASE, "RNA-seq/results/causal_inference")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
#  Parallel setup
# ============================================================
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
if (ncpus > 1) {
  register(SnowParam(ncpus))
} else {
  register(SerialParam())
}

# ============================================================
#  1. Load merged DGE and metadata (same as Script 25)
# ============================================================
cat("\n--- Step 1: Loading merged DGE object ---\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# Exclude datasets without healthy controls for mega-analysis
# Must match Script 05/25 exclusions for cross-method Jaccard to be valid
ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]
cat("Samples for mega-analysis (yaml mega cohorts:", paste(mega_cohorts, collapse = ", "), "):", ncol(dge_mega), "
")

# Metadata for inferred sex
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

# ============================================================
#  2. Load per-cohort BayesPrism proportions
# ============================================================
cat("\n--- Step 2: Loading per-cohort BayesPrism proportions ---\n")
datasets <- unique(dge_mega$samples$dataset)
cat("  Datasets:", paste(datasets, collapse = ", "), "\n")

bp_props <- lapply(datasets, function(ds) {
  prop_file <- file.path(DECONV, ds, paste0(ds, "_bayesprism_proportions.tsv"))
  if (file.exists(prop_file)) {
    df <- read.table(prop_file, header = TRUE, sep = "\t", check.names = FALSE)
    df$sample_id <- rownames(df)
    df$dataset <- ds
    df$method <- "BayesPrism"
    cat("  ", ds, ": loaded", nrow(df), "samples,",
        ncol(df) - 3, "cell types\n")
    return(as.data.table(df))
  } else {
    cat("  ", ds, ": BayesPrism file NOT FOUND\n")
    return(NULL)
  }
})
bp_dt <- rbindlist(bp_props, fill = TRUE)
cat("  Total BayesPrism samples:", nrow(bp_dt), "\n")

# ============================================================
#  3. Load per-cohort MuSiC proportions (for comparison)
# ============================================================
cat("\n--- Step 3: Loading per-cohort MuSiC proportions ---\n")
music_props <- lapply(datasets, function(ds) {
  prop_file <- file.path(DECONV, ds, paste0(ds, "_music_prop_weighted.tsv"))
  if (file.exists(prop_file)) {
    df <- read.table(prop_file, header = TRUE, sep = "\t", check.names = FALSE)
    df$sample_id <- rownames(df)
    df$dataset <- ds
    df$method <- "MuSiC"
    return(as.data.table(df))
  }
  return(NULL)
})
music_dt <- rbindlist(music_props, fill = TRUE)
cat("  Total MuSiC samples:", nrow(music_dt), "\n")

# ============================================================
#  4. Build unified proportions matching mega-analysis samples
# ============================================================
cat("\n--- Step 4: Building unified BayesPrism proportion matrix ---\n")
mega_samples <- colnames(dge_mega)

bp_hepa  <- bp_dt$Hepatocytes[match(mega_samples, bp_dt$sample_id)]
bp_macro <- bp_dt$Macrophages[match(mega_samples, bp_dt$sample_id)]
music_hepa  <- music_dt$Hepatocytes[match(mega_samples, music_dt$sample_id)]
music_macro <- music_dt$Macrophages[match(mega_samples, music_dt$sample_id)]

n_bp_matched <- sum(!is.na(bp_hepa))
n_music_matched <- sum(!is.na(music_hepa))
cat("  BayesPrism matched:", n_bp_matched, "/", length(mega_samples), "samples\n")
cat("  MuSiC matched:", n_music_matched, "/", length(mega_samples), "samples\n")

# Impute missing with median (same strategy as Script 25)
# Warn if imputation fraction is high (>20% of samples)
bp_pct_imputed <- 100 * sum(is.na(bp_hepa)) / length(mega_samples)
if (bp_pct_imputed > 20) {
  cat("  WARNING:", round(bp_pct_imputed, 1),
      "% of samples lack BayesPrism data — imputed with median.\n")
  cat("  Interpret BayesPrism-adjusted dream results with caution.\n")
}
if (any(is.na(bp_hepa)))  {
  cat("  Imputing", sum(is.na(bp_hepa)), "missing BayesPrism Hepatocyte fractions with median\n")
  bp_hepa[is.na(bp_hepa)]   <- median(bp_hepa, na.rm = TRUE)
}
if (any(is.na(bp_macro))) {
  cat("  Imputing", sum(is.na(bp_macro)), "missing BayesPrism Macrophage fractions with median\n")
  bp_macro[is.na(bp_macro)] <- median(bp_macro, na.rm = TRUE)
}
if (any(is.na(music_hepa)))  music_hepa[is.na(music_hepa)]   <- median(music_hepa, na.rm = TRUE)
if (any(is.na(music_macro))) music_macro[is.na(music_macro)] <- median(music_macro, na.rm = TRUE)

# Save unified proportions
unified_props <- data.table(
  sample_id       = mega_samples,
  dataset         = dge_mega$samples$dataset,
  bp_Hepatocytes  = bp_hepa,
  bp_Macrophages  = bp_macro,
  music_Hepatocytes = music_hepa,
  music_Macrophages = music_macro
)
fwrite(unified_props, file.path(OUT, "unified_bayesprism_proportions.csv"))
cat("  Saved: unified_bayesprism_proportions.csv\n")

# ============================================================
#  5. Method concordance: BayesPrism vs MuSiC
# ============================================================
cat("\n--- Step 5: Method concordance analysis ---\n")

# Only compare non-imputed samples
non_imputed <- !is.na(bp_dt$Hepatocytes[match(mega_samples, bp_dt$sample_id)]) &
               !is.na(music_dt$Hepatocytes[match(mega_samples, music_dt$sample_id)])

r_hepa  <- cor(bp_hepa[non_imputed], music_hepa[non_imputed], method = "pearson")
r_macro <- cor(bp_macro[non_imputed], music_macro[non_imputed], method = "pearson")
r_hepa_spearman  <- cor(bp_hepa[non_imputed], music_hepa[non_imputed], method = "spearman")
r_macro_spearman <- cor(bp_macro[non_imputed], music_macro[non_imputed], method = "spearman")

cat("  Hepatocyte fraction concordance:\n")
cat("    Pearson r =", round(r_hepa, 3), "\n")
cat("    Spearman rho =", round(r_hepa_spearman, 3), "\n")
cat("  Macrophage fraction concordance:\n")
cat("    Pearson r =", round(r_macro, 3), "\n")
cat("    Spearman rho =", round(r_macro_spearman, 3), "\n")

concordance_dt <- data.table(
  cell_type = c("Hepatocytes", "Macrophages"),
  pearson_r = c(r_hepa, r_macro),
  spearman_rho = c(r_hepa_spearman, r_macro_spearman),
  n_samples = c(sum(non_imputed), sum(non_imputed))
)
fwrite(concordance_dt, file.path(OUT, "method_concordance.csv"))
cat("  Saved: method_concordance.csv\n")

# ---- Concordance scatter plot ----
plot_dt <- data.table(
  bp_hepa = bp_hepa[non_imputed],
  music_hepa = music_hepa[non_imputed],
  bp_macro = bp_macro[non_imputed],
  music_macro = music_macro[non_imputed],
  dataset = dge_mega$samples$dataset[non_imputed]
)

pdf(file.path(FIGDIR, "bayesprism_concordance.pdf"), width = 12, height = 5)

p1 <- ggplot(plot_dt, aes(x = music_hepa, y = bp_hepa, color = dataset)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  annotate("text", x = 0.1, y = 0.95,
           label = paste0("r = ", round(r_hepa, 3)),
           size = 4, fontface = "bold") +
  labs(title = "Hepatocyte Fraction: BayesPrism vs MuSiC",
       x = "MuSiC Hepatocyte Fraction",
       y = "BayesPrism Hepatocyte Fraction") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

p2 <- ggplot(plot_dt, aes(x = music_macro, y = bp_macro, color = dataset)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  annotate("text", x = min(plot_dt$music_macro) + 0.02,
           y = max(plot_dt$bp_macro) * 0.95,
           label = paste0("r = ", round(r_macro, 3)),
           size = 4, fontface = "bold") +
  labs(title = "Macrophage Fraction: BayesPrism vs MuSiC",
       x = "MuSiC Macrophage Fraction",
       y = "BayesPrism Macrophage Fraction") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()
cat("  Saved: figures/bayesprism_concordance.pdf\n")

# ============================================================
#  6. Run BayesPrism-adjusted dream
# ============================================================
cat("\n--- Step 6: Running BayesPrism-adjusted dream ---\n")

info <- data.frame(
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(matched_sex),
  Hepatocytes  = bp_hepa,
  Macrophages  = bp_macro,
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge_mega)

cat("  Formula: ~ group_binary + inferred_sex + Hepatocytes + Macrophages + (1|dataset)\n")
form_adj <- ~ group_binary + inferred_sex + Hepatocytes + Macrophages + (1|dataset)

cat("  Running voomWithDreamWeights...\n")
v_adj <- suppressWarnings(voomWithDreamWeights(dge_mega, form_adj, info))

cat("  Running dream() with BayesPrism fractions...\n")
fit_adj <- suppressWarnings(dream(v_adj, form_adj, info))

res_adj <- topTable(fit_adj, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res_adj$gene <- rownames(res_adj)
dt_adj <- as.data.table(res_adj)
setnames(dt_adj, "adj.P.Val", "padj")

cat("  BayesPrism-adjusted DEGs (padj < 0.05, |logFC| > 0.5):",
    nrow(dt_adj[padj < 0.05 & abs(logFC) > 0.5]), "\n")

fwrite(dt_adj, file.path(OUT, "dream_results_bayesprism_adjusted.csv"))
cat("  Saved: dream_results_bayesprism_adjusted.csv\n")

# ============================================================
#  7. Attribution: Compare with unadjusted dream (Script 05)
# ============================================================
cat("\n--- Step 7: Attribution classification ---\n")
dt_unadj <- fread(file.path(RDIR, "dream_results.csv"))

common_genes <- intersect(dt_unadj$gene, dt_adj$gene)
cat("  Common genes:", length(common_genes), "\n")

dt_unadj_m <- dt_unadj[gene %in% common_genes]
dt_adj_m   <- dt_adj[gene %in% common_genes]
setkey(dt_unadj_m, gene)
setkey(dt_adj_m, gene)

# Classify using same logic as Script 25
PADJ_THR <- 0.05
LFC_THR  <- 0.5

attrib <- data.table(
  gene         = dt_unadj_m$gene,
  logFC_bp_adj = dt_adj_m[dt_unadj_m$gene, logFC],
  logFC_unadj  = dt_unadj_m$logFC,
  padj_bp_adj  = dt_adj_m[dt_unadj_m$gene, padj],
  padj_unadj   = dt_unadj_m$padj
)

sig_unadj <- attrib$padj_unadj < PADJ_THR & abs(attrib$logFC_unadj) > LFC_THR
sig_adj   <- attrib$padj_bp_adj < PADJ_THR & abs(attrib$logFC_bp_adj) > LFC_THR

attrib[, bp_class := "NS"]
attrib[sig_unadj & sig_adj,   bp_class := "Hepatocyte_intrinsic"]
attrib[sig_unadj & !sig_adj,  bp_class := "Composition_driven"]
attrib[!sig_unadj & sig_adj,  bp_class := "Unmasked"]

cat("  BayesPrism attribution results:\n")
print(attrib[, .N, by = bp_class][order(-N)])

fwrite(attrib, file.path(OUT, "bayesprism_attribution_scores.csv"))
cat("  Saved: bayesprism_attribution_scores.csv\n")

# ============================================================
#  8. Compare with MuSiC-based attribution (Script 25)
# ============================================================
cat("\n--- Step 8: Cross-method attribution concordance ---\n")
music_attrib_file <- file.path(CAUSAL_DIR, "deconv_attribution_scores.csv")

if (file.exists(music_attrib_file)) {
  music_attrib <- fread(music_attrib_file)
  cat("  MuSiC attribution loaded:", nrow(music_attrib), "genes\n")

  # Merge on gene — MuSiC attribution file uses "category" column
  music_class_col <- if ("category" %in% names(music_attrib)) "category" else
                     if ("class" %in% names(music_attrib)) "class" else
                     NA_character_

  if (!is.na(music_class_col)) {
    compare <- merge(
      attrib[, .(gene, bp_class)],
      music_attrib[, .(gene, music_class = get(music_class_col))],
      by = "gene"
    )

    cat("  Genes in both attribution sets:", nrow(compare), "\n")

    # Concordance table
    cat("\n  Cross-tabulation (BayesPrism vs MuSiC attribution):\n")
    cross_tab <- table(compare$bp_class, compare$music_class)
    print(cross_tab)

    # Overall agreement rate
    agree <- sum(compare$bp_class == compare$music_class)
    cat("\n  Overall agreement:", agree, "/", nrow(compare),
        "(", round(100 * agree / nrow(compare), 1), "%)\n")

    # Hepatocyte-intrinsic concordance
    both_hi <- sum(compare$bp_class == "Hepatocyte_intrinsic" &
                   compare$music_class == "Hepatocyte_intrinsic")
    either_hi <- sum(compare$bp_class == "Hepatocyte_intrinsic" |
                     compare$music_class == "Hepatocyte_intrinsic")
    jaccard_hi <- both_hi / either_hi
    cat("  Hepatocyte_intrinsic Jaccard:", round(jaccard_hi, 3), "\n")

    # Save comparison
    fwrite(compare, file.path(OUT, "attribution_cross_method_comparison.csv"))
    cat("  Saved: attribution_cross_method_comparison.csv\n")

    # ---- Comparison visualization ----
    pdf(file.path(FIGDIR, "bayesprism_attribution_comparison.pdf"),
        width = 10, height = 8)

    # Sankey/alluvial-style comparison (simplified as grouped bar)
    compare_summary <- as.data.table(compare)[, .N, by = .(bp_class, music_class)]
    p3 <- ggplot(compare_summary,
                 aes(x = music_class, y = N, fill = bp_class)) +
      geom_col(position = "dodge", width = 0.7) +
      scale_fill_brewer(palette = "Set2", name = "BayesPrism\nClassification") +
      labs(title = "Deconvolution Attribution: BayesPrism vs MuSiC",
           subtitle = paste0("Agreement: ", round(100 * agree / nrow(compare), 1),
                             "% | Hepatocyte-intrinsic Jaccard: ",
                             round(jaccard_hi, 3)),
           x = "MuSiC Classification",
           y = "Number of Genes") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            plot.title = element_text(face = "bold"))

    print(p3)
    dev.off()
    cat("  Saved: figures/bayesprism_attribution_comparison.pdf\n")
  } else {
    cat("  WARNING: Could not identify attribution class column in MuSiC results\n")
  }
} else {
  cat("  MuSiC attribution file not found at:", music_attrib_file, "\n")
  cat("  Skipping cross-method comparison. Run Script 25 first.\n")
}

# ============================================================
#  Summary
# ============================================================
cat("\n=== Script 22: BayesPrism Deconvolution Complete ===\n")
cat("  Method concordance (Hepatocytes): Pearson r =", round(r_hepa, 3), "\n")
cat("  Method concordance (Macrophages): Pearson r =", round(r_macro, 3), "\n")
bp_hi_count <- sum(attrib$bp_class == "Hepatocyte_intrinsic")
bp_cd_count <- sum(attrib$bp_class == "Composition_driven")
bp_um_count <- sum(attrib$bp_class == "Unmasked")
total_sig <- bp_hi_count + bp_cd_count + bp_um_count
cat("  BayesPrism attribution:\n")
cat("    Hepatocyte_intrinsic:", bp_hi_count,
    "(", round(100 * bp_hi_count / max(total_sig, 1), 1), "%)\n")
cat("    Composition_driven:", bp_cd_count,
    "(", round(100 * bp_cd_count / max(total_sig, 1), 1), "%)\n")
cat("    Unmasked:", bp_um_count,
    "(", round(100 * bp_um_count / max(total_sig, 1), 1), "%)\n")
cat("End time:", format(Sys.time()), "\n")
