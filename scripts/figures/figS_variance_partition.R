#!/usr/bin/env Rscript
# figS_variance_partition.R
# BMI sensitivity Prong 2: variance partition analysis
# Quantifies expression variance explained by disease, fibrosis, NAS, dataset, sex
# If NAS/steatosis (BMI proxy) explains little variance vs disease, BMI confounding is negligible

suppressPackageStartupMessages({
  library(variancePartition)
  library(edgeR)
  library(BiocParallel)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
})

set.seed(42)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

DGE_PATH   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
META_PATH  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier/modeling_metadata.csv")
DREAM_PATH <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv")

OUT_FIG  <- file.path(FIGS_SENS_DIR, "figS_variance_partition.pdf")
OUT_CSV  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/variance_partition_results.csv")
OUT_SUM  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/variance_partition_summary.csv")

dir.create(dirname(OUT_FIG), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGS_SENS_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "16"))
cat("Using", ncores, "cores\n")

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
cat("Loading DGE...\n")
dge <- readRDS(DGE_PATH)
cat("  Samples:", ncol(dge), " Genes:", nrow(dge), "\n")

meta <- read.csv(META_PATH, stringsAsFactors = FALSE)
dream <- read.csv(DREAM_PATH, stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# 2. Prepare metadata — merge into DGE$samples
# ---------------------------------------------------------------------------
idx <- match(colnames(dge), meta$sample_id)
stopifnot(all(!is.na(idx)))

dge$samples$fibrosis_stage  <- meta$fibrosis_stage[idx]
dge$samples$nas_score       <- meta$nas_score[idx]
dge$samples$steatosis_grade <- meta$steatosis_grade[idx]
# Recode steatosis -1 as NA (means unavailable)
dge$samples$steatosis_grade[dge$samples$steatosis_grade == -1] <- NA

# ---------------------------------------------------------------------------
# 3. Infer sex for samples with missing sex (k-means on XIST/DDX3Y)
# ---------------------------------------------------------------------------
cat("Inferring sex for samples with missing annotation...\n")
xist_id  <- grep("ENSG00000229807", rownames(dge), value = TRUE)
ddx3y_id <- grep("ENSG00000067048", rownames(dge), value = TRUE)
stopifnot(length(xist_id) == 1, length(ddx3y_id) == 1)

lcpm <- edgeR::cpm(dge, log = TRUE)
sex_mat <- cbind(XIST = lcpm[xist_id, ], DDX3Y = lcpm[ddx3y_id, ])
km <- kmeans(sex_mat, centers = 2, nstart = 25)

# Assign F to cluster with higher XIST
xist_means <- tapply(sex_mat[, "XIST"], km$cluster, mean)
female_cluster <- as.integer(names(which.max(xist_means)))
inferred_sex <- ifelse(km$cluster == female_cluster, "F", "M")

# Use annotated sex where available; inferred where missing
dge$samples$sex_covar <- ifelse(is.na(dge$samples$sex), inferred_sex,
                                as.character(dge$samples$sex))
cat("  Sex distribution (final):\n")
print(table(dge$samples$sex_covar))

# Validate: among samples with annotated sex, check concordance
annotated <- !is.na(dge$samples$sex)
concordance <- mean(dge$samples$sex_covar[annotated] == inferred_sex[annotated])
cat("  Sex inference concordance with annotation:", round(concordance, 4), "\n")

# ---------------------------------------------------------------------------
# 4. Filter and normalize
# ---------------------------------------------------------------------------
cat("Filtering genes...\n")
keep <- filterByExpr(dge, group = dge$samples$group_binary)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
cat("  Genes after filtering:", nrow(dge), "\n")

# Convert categorical variables to factors
dge$samples$group_binary <- factor(dge$samples$group_binary)
dge$samples$dataset      <- factor(dge$samples$dataset)
dge$samples$sex_covar    <- factor(dge$samples$sex_covar)

# ---------------------------------------------------------------------------
# 5. Model A1: All 1,444 samples (group + sex + dataset, no fibrosis)
# ---------------------------------------------------------------------------
cat("\n=== Model A1: Full cohort (N=", ncol(dge), ") ===\n")
form_a1 <- ~ (1|group_binary) + (1|dataset) + (1|sex_covar)

vobj_a1 <- voomWithDreamWeights(dge, form_a1, dge$samples)
param <- SnowParam(ncores, "SOCK", progressbar = TRUE)
cat("Fitting variance partition (A1)...\n")
t0 <- Sys.time()
varPart_a1 <- fitExtractVarPartModel(vobj_a1, form_a1, dge$samples, BPPARAM = param)
cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

# ---------------------------------------------------------------------------
# 6. Model A2: Fibrosis subset (~1,128 samples)
# ---------------------------------------------------------------------------
cat("\n=== Model A2: Fibrosis subset ===\n")
has_fib <- !is.na(dge$samples$fibrosis_stage)
cat("  Samples with fibrosis:", sum(has_fib), "\n")
dge_fib <- dge[, has_fib]
dge_fib$samples$fibrosis_stage <- factor(dge_fib$samples$fibrosis_stage)
dge_fib$samples <- droplevels(dge_fib$samples)

form_a2 <- ~ (1|group_binary) + (1|dataset) + (1|sex_covar) + (1|fibrosis_stage)
vobj_a2 <- voomWithDreamWeights(dge_fib, form_a2, dge_fib$samples)
cat("Fitting variance partition (A2)...\n")
t0 <- Sys.time()
varPart_a2 <- fitExtractVarPartModel(vobj_a2, form_a2, dge_fib$samples, BPPARAM = param)
cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

# ---------------------------------------------------------------------------
# 7. Model B: NAS subset (~660 samples)
# ---------------------------------------------------------------------------
cat("\n=== Model B: NAS subset ===\n")
has_nas <- !is.na(dge$samples$nas_score)
cat("  Samples with NAS:", sum(has_nas), "\n")
dge_nas <- dge[, has_nas]
dge_nas$samples$nas_score      <- factor(dge_nas$samples$nas_score)
dge_nas$samples$fibrosis_stage <- factor(dge_nas$samples$fibrosis_stage)
dge_nas$samples <- droplevels(dge_nas$samples)

# Check fibrosis availability in NAS subset
fib_in_nas <- sum(!is.na(dge_nas$samples$fibrosis_stage))
cat("  Fibrosis available in NAS subset:", fib_in_nas, "/", ncol(dge_nas), "\n")

if (fib_in_nas == ncol(dge_nas)) {
  form_b <- ~ (1|group_binary) + (1|dataset) + (1|sex_covar) + (1|fibrosis_stage) + (1|nas_score)
} else {
  form_b <- ~ (1|group_binary) + (1|dataset) + (1|sex_covar) + (1|nas_score)
}
vobj_b <- voomWithDreamWeights(dge_nas, form_b, dge_nas$samples)
cat("Fitting variance partition (B)...\n")
t0 <- Sys.time()
varPart_b <- fitExtractVarPartModel(vobj_b, form_b, dge_nas$samples, BPPARAM = param)
cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

# ---------------------------------------------------------------------------
# 8. Identify DEGs
# ---------------------------------------------------------------------------
if ("gene" %in% colnames(dream)) {
  gene_col <- "gene"
} else {
  gene_col <- colnames(dream)[1]
}
deg_genes <- dream[[gene_col]][!is.na(dream$adj.P.Val) & dream$adj.P.Val < 0.05 &
                                 !is.na(dream$logFC) & abs(dream$logFC) > 0.5]
cat("\nDEGs (padj<0.05, |logFC|>0.5):", length(deg_genes), "\n")

deg_in_vp   <- intersect(deg_genes, rownames(varPart_a2))
nondeg_in_vp <- setdiff(rownames(varPart_a2), deg_genes)
cat("DEGs in VP results:", length(deg_in_vp), "\n")
cat("Non-DEGs in VP results:", length(nondeg_in_vp), "\n")

# ---------------------------------------------------------------------------
# 9. Save per-gene results (Model A2 = most informative)
# ---------------------------------------------------------------------------
vp_df <- as.data.frame(varPart_a2)
vp_df$gene <- rownames(vp_df)
vp_df$is_deg <- vp_df$gene %in% deg_genes
write.csv(vp_df, OUT_CSV, row.names = FALSE)
cat("Saved per-gene results to:", OUT_CSV, "\n")

# ---------------------------------------------------------------------------
# 10. Summary table
# ---------------------------------------------------------------------------
make_summary <- function(vp_mat, model_name) {
  data.frame(
    model     = model_name,
    covariate = colnames(vp_mat),
    median_pct = apply(vp_mat, 2, median) * 100,
    mean_pct   = apply(vp_mat, 2, mean) * 100,
    q25_pct    = apply(vp_mat, 2, quantile, 0.25) * 100,
    q75_pct    = apply(vp_mat, 2, quantile, 0.75) * 100,
    stringsAsFactors = FALSE, row.names = NULL
  )
}

sum_a1       <- make_summary(varPart_a1, "A1_all_1444")
sum_a2       <- make_summary(varPart_a2, "A2_fibrosis_subset")
sum_a2_deg   <- make_summary(varPart_a2[deg_in_vp, , drop = FALSE], "A2_DEGs_only")
sum_a2_nondeg <- make_summary(varPart_a2[nondeg_in_vp, , drop = FALSE], "A2_nonDEGs")
sum_b        <- make_summary(varPart_b, "B_NAS_subset")

summary_all <- bind_rows(sum_a1, sum_a2, sum_a2_deg, sum_a2_nondeg, sum_b)
write.csv(summary_all, OUT_SUM, row.names = FALSE)
cat("Saved summary to:", OUT_SUM, "\n")

cat("\nModel A2 (fibrosis subset) median variance explained:\n")
print(sum_a2[, c("covariate", "median_pct")], row.names = FALSE)
cat("\nModel B (NAS subset) median variance explained:\n")
print(sum_b[, c("covariate", "median_pct")], row.names = FALSE)

# ---------------------------------------------------------------------------
# 11. Figures
# ---------------------------------------------------------------------------
cat("\nGenerating figures...\n")

covar_labels <- c(
  group_binary   = "Disease status",
  dataset        = "Dataset (batch)",
  sex_covar      = "Sex",
  fibrosis_stage = "Fibrosis stage",
  nas_score      = "NAS score",
  Residuals      = "Residual"
)

# --- Panel A: Violin plot (Model A2) ---
vp_long_a2 <- as.data.frame(varPart_a2) %>%
  mutate(gene = rownames(varPart_a2)) %>%
  pivot_longer(-gene, names_to = "covariate", values_to = "variance_fraction") %>%
  mutate(
    covariate = factor(covar_labels[covariate],
                       levels = c("Disease status", "Fibrosis stage",
                                  "Dataset (batch)", "Sex", "Residual"))
  )

covar_fills <- c(
  "Disease status"  = masld_colors$masld,
  "Fibrosis stage"  = masld_colors$fibrosis,
  "Dataset (batch)" = "#78909C",
  "Sex"             = "#7E57C2",
  "Residual"        = "#E0E0E0",
  "NAS score"       = "#FF8F00"
)

panel_a <- ggplot(vp_long_a2, aes(x = covariate, y = variance_fraction * 100,
                                   fill = covariate)) +
  geom_violin(scale = "width", alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = covar_fills) +
  labs(x = NULL, y = "Variance explained (%)",
       title = paste0("Expression variance decomposition (N=", sum(has_fib), ")")) +
  theme_masld() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

# --- Panel B: Stacked bar — DEGs vs non-DEGs ---
bar_data <- bind_rows(
  sum_a2_deg %>% mutate(gene_set = "DEGs"),
  sum_a2_nondeg %>% mutate(gene_set = "Non-DEGs")
) %>%
  mutate(
    covariate_label = covar_labels[covariate],
    covariate_label = factor(covariate_label,
                             levels = c("Residual", "Sex", "Dataset (batch)",
                                        "Fibrosis stage", "Disease status"))
  )

panel_b <- ggplot(bar_data, aes(x = gene_set, y = median_pct, fill = covariate_label)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = covar_fills, name = "Covariate") +
  labs(x = NULL, y = "Median variance explained (%)",
       title = "DEGs vs non-DEGs") +
  theme_masld() +
  theme(legend.position = "right")

# --- Panel C: Scatter — disease vs fibrosis variance ---
vp_scatter <- as.data.frame(varPart_a2) %>%
  mutate(gene = rownames(varPart_a2),
         is_deg = gene %in% deg_genes)

panel_c <- ggplot(vp_scatter, aes(x = group_binary * 100, y = fibrosis_stage * 100,
                                   color = is_deg)) +
  geom_point(size = 0.3, alpha = 0.3) +
  scale_color_manual(values = c("FALSE" = "#BDBDBD", "TRUE" = masld_colors$masld),
                     labels = c("Non-DEG", "DEG"), name = NULL) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Variance from disease status (%)",
       y = "Variance from fibrosis stage (%)",
       title = "Disease vs fibrosis variance") +
  theme_masld() +
  theme(legend.position = c(0.8, 0.9))

# --- Panel D: Top 20 genes by NAS variance (Model B) ---
vp_b_df <- as.data.frame(varPart_b) %>%
  mutate(gene = rownames(varPart_b)) %>%
  arrange(desc(nas_score))

# Map Ensembl IDs to gene symbols via dream results
gene_map <- NULL
if ("gene_name" %in% colnames(dream)) {
  gene_map <- dream[, c(gene_col, "gene_name")]
  colnames(gene_map) <- c("gene", "gene_name")
  gene_map <- gene_map[!duplicated(gene_map$gene), ]
}

top20 <- head(vp_b_df, 20)
if (!is.null(gene_map)) {
  top20 <- left_join(top20, gene_map, by = "gene")
  top20$display_name <- ifelse(is.na(top20$gene_name) | top20$gene_name == "",
                               sub("\\..*", "", top20$gene),
                               top20$gene_name)
} else {
  top20$display_name <- sub("\\..*", "", top20$gene)
}

# Pivot for stacked bar
var_cols <- intersect(colnames(varPart_b), colnames(top20))
top20_long <- top20 %>%
  select(display_name, all_of(var_cols)) %>%
  mutate(display_name = factor(display_name, levels = rev(display_name))) %>%
  pivot_longer(-display_name, names_to = "covariate", values_to = "variance_fraction") %>%
  mutate(
    covariate_label = covar_labels[covariate],
    covariate_label = factor(covariate_label,
                             levels = rev(c("Disease status", "Fibrosis stage", "NAS score",
                                            "Dataset (batch)", "Sex", "Residual")))
  )

panel_d <- ggplot(top20_long, aes(x = variance_fraction * 100, y = display_name,
                                   fill = covariate_label)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = covar_fills, name = "Covariate") +
  labs(x = "Variance explained (%)", y = NULL,
       title = "Top 20 genes by NAS variance") +
  theme_masld() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 7, face = "italic"))

# --- Compose 4-panel figure ---
fig <- plot_grid(panel_a, panel_b, panel_c, panel_d,
                 ncol = 2, labels = c("a", "b", "c", "d"),
                 label_size = 14, rel_widths = c(1, 1.2))

ggsave(OUT_FIG, fig, width = 14, height = 10, device = cairo_pdf)
cat("Saved figure to:", OUT_FIG, "\n")

# ---------------------------------------------------------------------------
# 12. Print key findings
# ---------------------------------------------------------------------------
cat("\n=== KEY FINDINGS ===\n")
cat("Model A1 (N=", ncol(dge), ", group + sex + dataset):\n")
for (i in seq_len(nrow(sum_a1))) {
  cat(sprintf("  %s: %.1f%% median\n", sum_a1$covariate[i], sum_a1$median_pct[i]))
}

cat("\nModel A2 (N=", sum(has_fib), ", group + fibrosis + sex + dataset):\n")
for (i in seq_len(nrow(sum_a2))) {
  cat(sprintf("  %s: %.1f%% median\n", sum_a2$covariate[i], sum_a2$median_pct[i]))
}

cat("\nModel B (N=", sum(has_nas), ", group + fibrosis + NAS + sex + dataset):\n")
for (i in seq_len(nrow(sum_b))) {
  cat(sprintf("  %s: %.1f%% median\n", sum_b$covariate[i], sum_b$median_pct[i]))
}

# Disease vs NAS comparison for DEGs
if ("nas_score" %in% colnames(varPart_b)) {
  deg_in_b <- intersect(deg_genes, rownames(varPart_b))
  if (length(deg_in_b) > 0) {
    cat("\nAmong DEGs in NAS model:\n")
    cat(sprintf("  Disease median: %.1f%%\n", median(varPart_b[deg_in_b, "group_binary"]) * 100))
    cat(sprintf("  NAS median: %.1f%%\n", median(varPart_b[deg_in_b, "nas_score"]) * 100))
    if ("fibrosis_stage" %in% colnames(varPart_b)) {
      cat(sprintf("  Fibrosis median: %.1f%%\n", median(varPart_b[deg_in_b, "fibrosis_stage"]) * 100))
    }

    # Genes where NAS > disease
    nas_dom <- sum(varPart_b[deg_in_b, "nas_score"] > varPart_b[deg_in_b, "group_binary"])
    cat(sprintf("  DEGs where NAS > disease: %d / %d (%.1f%%)\n",
                nas_dom, length(deg_in_b), nas_dom / length(deg_in_b) * 100))
  }
}

cat("\nDone.\n")
