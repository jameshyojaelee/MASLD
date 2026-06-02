#!/usr/bin/env Rscript
# Compositional Analysis: CLR + Linear Mixed Model
#
# Scientific question:
#   Do cell type proportions shift with MASLD?
#   Primary:   binary Healthy vs MASLD (all filtered samples)
#   Secondary: continuous NAS score + fibrosis stage (if >=30 paired samples)
#
# Key hypothesis: If hepatocytes are 83% intrinsic, their proportions
# should be roughly stable — a drop would instead suggest composition
# drives the bulk signal.
#
# Method: CLR (Aitchison-safe) + lme4::lmer with dataset as random effect
# Input:  results_gpu_v2/cell_type_proportions.csv
# Output: results_gpu_v2/compositional/

suppressPackageStartupMessages({
  library(lme4)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

sc_dir    <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2")
meta_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
out_dir   <- file.path(sc_dir, "compositional")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Load + filter proportions
# ---------------------------------------------------------------------------
message("Loading cell type proportions...")
props <- fread(file.path(sc_dir, "cell_type_proportions.csv"))
message("  Total samples: ", nrow(props))

# Filter: unsorted + nuclei only (comparable preparation methods)
keep_preps <- c("unsorted", "nuclei")
props <- props[preparation_method %in% keep_preps]
message("  After prep-method filter (", paste(keep_preps, collapse = "/"), "): ",
        nrow(props), " samples")

# Keep only Healthy vs MASLD
props <- props[condition_harmonized %in% c("Healthy", "MASLD")]
message("  After condition filter (Healthy/MASLD): ", nrow(props), " samples")
message("  Condition: ", paste(names(table(props$condition_harmonized)),
                               table(props$condition_harmonized), sep = "=",
                               collapse = ", "))
message("  Dataset:   ", paste(names(table(props$dataset)), table(props$dataset),
                               sep = "=", collapse = ", "))

# Identify cell type columns (all non-metadata columns)
meta_cols <- c("sample", "dataset", "condition", "preparation_method",
               "condition_harmonized")
ct_cols <- setdiff(colnames(props), meta_cols)
message("  Cell types (", length(ct_cols), "): ",
        paste(ct_cols, collapse = ", "))

# ---------------------------------------------------------------------------
# 2. Severity metadata (NAS / fibrosis stage) — exploratory
# ---------------------------------------------------------------------------
message("\nChecking for severity metadata (NAS/fibrosis)...")
nas_n <- 0
fib_n <- 0

if (file.exists(meta_path)) {
  bulk_meta <- fread(meta_path, select = c("sample_id", "dataset",
                                            "nas_score", "fibrosis_stage"))
  # The SC datasets are not in the bulk metadata
  sc_datasets <- unique(props$dataset)
  matched_meta <- bulk_meta[dataset %in% sc_datasets]
  if (nrow(matched_meta) > 0) {
    props <- merge(props, matched_meta, by.x = c("sample", "dataset"),
                   by.y = c("sample_id", "dataset"), all.x = TRUE)
    nas_n <- sum(!is.na(props$nas_score))
    fib_n <- sum(!is.na(props$fibrosis_stage))
  }
}

message("  Samples with NAS score: ", nas_n)
message("  Samples with fibrosis stage: ", fib_n)
message("  (Threshold for severity analysis: >=30 samples)")
run_nas_model <- nas_n >= 30
run_fib_model <- fib_n >= 30

if (!run_nas_model && !run_fib_model) {
  message("  SC datasets not in bulk metadata — severity analysis skipped.")
}

# ---------------------------------------------------------------------------
# 3. CLR transformation (Aitchison-safe)
# ---------------------------------------------------------------------------
message("\nApplying CLR transformation...")

# Extract proportion matrix
prop_mat <- as.matrix(props[, ..ct_cols])
K <- ncol(prop_mat)

# Pseudocount: 0.5/K (Aitchison recommendation)
pseudocount <- 0.5 / K
prop_mat_adj <- prop_mat + pseudocount

# CLR: log(x_i) - mean(log(x_j))  for j = 1..K
log_mat <- log(prop_mat_adj)
row_means <- rowMeans(log_mat)
clr_mat <- log_mat - row_means

colnames(clr_mat) <- ct_cols

# Save CLR values
clr_dt <- data.table(
  sample              = props$sample,
  dataset             = props$dataset,
  condition_harmonized = props$condition_harmonized,
  clr_mat
)
fwrite(clr_dt, file.path(out_dir, "clr_values.csv"))
message("  Saved: clr_values.csv (", nrow(clr_dt), " samples x ",
        length(ct_cols), " cell types)")

# ---------------------------------------------------------------------------
# 4. Binary LMM: clr_ct ~ condition_harmonized + (1|dataset)
# ---------------------------------------------------------------------------
message("\nFitting binary LMM (MASLD vs Healthy) per cell type...")

fit_lmm_binary <- function(ct_col, data) {
  df <- data.frame(
    clr_val   = data[[ct_col]],
    condition = factor(data$condition_harmonized,
                       levels = c("Healthy", "MASLD")),
    dataset   = factor(data$dataset)
  )
  n_datasets <- length(unique(df$dataset))

  tryCatch({
    if (n_datasets > 1) {
      fit <- lmer(clr_val ~ condition + (1 | dataset), data = df, REML = FALSE)
    } else {
      # Single dataset: no random effect needed
      fit <- lm(clr_val ~ condition, data = df)
    }

    # Extract coefficient for MASLD
    coef_tab <- tryCatch(coef(summary(fit)), error = function(e) NULL)
    if (is.null(coef_tab)) return(NULL)

    # Row for condition
    cond_row <- grep("conditionMASLD", rownames(coef_tab))
    if (length(cond_row) == 0) return(NULL)

    beta  <- coef_tab[cond_row, 1]
    se    <- coef_tab[cond_row, 2]
    t_val <- coef_tab[cond_row, 3]
    # p-value: lmer (class merMod) does NOT return p-values by default.
    # Use large-sample normal approximation (valid for n > 30 samples).
    # lm objects return a 4-column matrix with p in column 4.
    if (inherits(fit, "lm") && ncol(coef_tab) >= 4) {
      p_val <- coef_tab[cond_row, 4]
    } else {
      p_val <- 2 * pnorm(-abs(t_val))
    }
    if (is.na(p_val) || p_val < 0 || p_val > 1) {
      p_val <- 2 * pnorm(-abs(t_val))
    }

    data.frame(
      cell_type   = ct_col,
      beta        = beta,
      se          = se,
      t_stat      = t_val,
      p_value     = p_val,
      ci_lo       = beta - 1.96 * se,
      ci_hi       = beta + 1.96 * se,
      n_samples   = nrow(df),
      n_healthy   = sum(df$condition == "Healthy"),
      n_masld     = sum(df$condition == "MASLD"),
      n_datasets  = n_datasets,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("  Warning: LMM failed for ", ct_col, ": ", e$message)
    NULL
  })
}

binary_results <- rbindlist(lapply(ct_cols, fit_lmm_binary, data = clr_dt),
                             use.names = TRUE, fill = TRUE)

# BH correction across all cell types
binary_results[, padj := p.adjust(p_value, method = "BH")]
binary_results[, direction := ifelse(beta > 0, "Up in MASLD", "Down in MASLD")]
binary_results[, significant := padj < 0.05]

# Order by absolute effect size
binary_results[, abs_beta := abs(beta)]
setorder(binary_results, -abs_beta)
binary_results[, abs_beta := NULL]

fwrite(binary_results, file.path(out_dir, "lmm_binary_results.csv"))
message("  Saved: lmm_binary_results.csv")
message("  Significant (padj<0.05): ",
        sum(binary_results$significant, na.rm = TRUE), " cell types")

# Print summary
print(binary_results[, .(cell_type, beta, se, p_value, padj, significant,
                          direction)])

# ---------------------------------------------------------------------------
# 5. Severity models (NAS / fibrosis) — if enough samples
# ---------------------------------------------------------------------------
sev_results <- data.table()

fit_lmm_severity <- function(ct_col, data, sev_col) {
  df <- data.frame(
    clr_val  = data[[ct_col]],
    severity = data[[sev_col]],
    dataset  = factor(data$dataset)
  )
  df <- df[!is.na(df$severity), ]
  if (nrow(df) < 10) return(NULL)
  n_datasets <- length(unique(df$dataset))

  tryCatch({
    if (n_datasets > 1) {
      fit <- lmer(clr_val ~ severity + (1 | dataset), data = df, REML = FALSE)
    } else {
      fit <- lm(clr_val ~ severity, data = df)
    }
    coef_tab <- coef(summary(fit))
    sev_row  <- grep("severity", rownames(coef_tab))
    if (length(sev_row) == 0) return(NULL)

    beta  <- coef_tab[sev_row, 1]
    se    <- coef_tab[sev_row, 2]
    t_val <- coef_tab[sev_row, 3]
    if (inherits(fit, "lm") && ncol(coef_tab) >= 4) {
      p_val <- coef_tab[sev_row, 4]
    } else {
      p_val <- 2 * pnorm(-abs(t_val))
    }
    if (is.na(p_val) || p_val < 0 || p_val > 1) {
      p_val <- 2 * pnorm(-abs(t_val))
    }

    data.frame(
      cell_type    = ct_col,
      severity_var = sev_col,
      beta         = beta,
      se           = se,
      t_stat       = t_val,
      p_value      = p_val,
      ci_lo        = beta - 1.96 * se,
      ci_hi        = beta + 1.96 * se,
      n_samples    = nrow(df),
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
}

if (run_nas_model) {
  message("\nFitting NAS-score severity model...")
  nas_res <- rbindlist(lapply(ct_cols, fit_lmm_severity,
                               data = clr_dt, sev_col = "nas_score"),
                        use.names = TRUE, fill = TRUE)
  nas_res[, padj := p.adjust(p_value, method = "BH")]
  sev_results <- rbind(sev_results, nas_res)
}

if (run_fib_model) {
  message("Fitting fibrosis-stage severity model...")
  fib_res <- rbindlist(lapply(ct_cols, fit_lmm_severity,
                               data = clr_dt, sev_col = "fibrosis_stage"),
                        use.names = TRUE, fill = TRUE)
  fib_res[, padj := p.adjust(p_value, method = "BH")]
  sev_results <- rbind(sev_results, fib_res)
}

if (nrow(sev_results) == 0) {
  sev_results <- data.table(
    note = "Severity analysis skipped: SC datasets not in bulk metadata"
  )
}
fwrite(sev_results, file.path(out_dir, "lmm_severity_results.csv"))
message("  Saved: lmm_severity_results.csv")

# ---------------------------------------------------------------------------
# 6. Forest plot (binary model)
# ---------------------------------------------------------------------------
message("\nGenerating forest plot...")

# Prettify cell type names
binary_results[, ct_label := gsub("\\.", " ", cell_type)]
binary_results[, ct_label := gsub("Mono mono", "Mono+mono", ct_label)]

# Sort by beta for display
binary_results <- binary_results[order(beta)]
binary_results[, ct_label := factor(ct_label, levels = ct_label)]

fp <- ggplot(binary_results,
             aes(x = beta, y = ct_label,
                 color = significant, shape = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                orientation = "y", width = 0.3, linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#7F8C8D"),
                     labels = c("TRUE" = "padj < 0.05", "FALSE" = "NS"),
                     name = NULL) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "padj < 0.05", "FALSE" = "NS"),
                     name = NULL) +
  labs(
    title    = "Compositional Analysis: Cell Type Proportions in MASLD vs Healthy",
    subtitle = paste0("CLR-transformed proportions; lme4 LMM with dataset random effect\n",
                      "n=", nrow(clr_dt), " samples (unsorted/nuclei); ",
                      sum(binary_results$significant, na.rm = TRUE),
                      " cell types significant at BH-adj p<0.05"),
    x = "Beta (CLR units; positive = increased in MASLD)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major.y = element_line(color = "grey95"),
    legend.position   = "bottom",
    plot.title        = element_text(face = "bold", size = 11),
    plot.subtitle     = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(out_dir, "compositional_forest_plot.pdf"),
       fp, width = 7, height = 6)
message("  Saved: compositional_forest_plot.pdf")

# ---------------------------------------------------------------------------
# 7. Sample-level PCA of CLR proportions, colored by hepatocyte CLR score
# ---------------------------------------------------------------------------
message("\nGenerating PCA plot colored by hepatocyte CLR score...")

clr_only <- as.matrix(clr_dt[, ..ct_cols])

# PCA of CLR-transformed proportions
pca_res <- prcomp(clr_only, scale. = FALSE, center = TRUE)
pca_df  <- data.frame(
  PC1  = pca_res$x[, 1],
  PC2  = pca_res$x[, 2],
  condition = clr_dt$condition_harmonized,
  dataset   = clr_dt$dataset,
  hep_clr   = clr_only[, "Hepatocytes"]
)

pvar <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_plt <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = hep_clr, shape = condition), size = 2, alpha = 0.8) +
  scale_color_gradient2(low = "#2980B9", mid = "white", high = "#C0392B",
                        midpoint = 0, name = "Hepatocyte\nCLR score") +
  scale_shape_manual(values = c("Healthy" = 16, "MASLD" = 17), name = "Condition") +
  labs(
    title = "Sample PCA of CLR-transformed cell type proportions",
    subtitle = "Color = hepatocyte CLR score; shape = condition",
    x = paste0("PC1 (", pvar[1], "% var)"),
    y = paste0("PC2 (", pvar[2], "% var)")
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "compositional_umap_composition.pdf"),
       pca_plt, width = 6.5, height = 5)
message("  Saved: compositional_umap_composition.pdf (PCA of CLR proportions)")

# ---------------------------------------------------------------------------
# 8. Supplementary boxplot: hepatocyte + macrophage CLR by condition
# ---------------------------------------------------------------------------
key_cts <- c("Hepatocytes", "Macrophages", "Fibroblasts",
             "Endothelial.cells", "Cholangiocytes")
key_cts <- intersect(key_cts, ct_cols)

box_dt <- melt(clr_dt[, c("condition_harmonized", "dataset", key_cts),
                       with = FALSE],
               id.vars    = c("condition_harmonized", "dataset"),
               variable.name = "cell_type",
               value.name    = "clr_value")
box_dt[, cell_type := gsub("\\.", " ", cell_type)]
box_dt[, cell_type := gsub("Mono mono", "Mono+mono", cell_type)]

box_plt <- ggplot(box_dt,
                  aes(x = condition_harmonized, y = clr_value,
                      fill = condition_harmonized)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.7) +
  scale_fill_manual(values = c("Healthy" = "#3498DB", "MASLD" = "#E74C3C"),
                    guide = "none") +
  facet_wrap(~ cell_type, scales = "free_y", nrow = 1) +
  labs(
    title = "Key cell type CLR proportions: Healthy vs MASLD",
    x = NULL, y = "CLR score"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey95"),
    axis.text.x       = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(out_dir, "compositional_boxplot_key_celltypes.pdf"),
       box_plt, width = 10, height = 4)
message("  Saved: compositional_boxplot_key_celltypes.pdf")

message("\n=== Module 1 (CLR + LMM) complete ===")
message("Results in: ", out_dir)
message("Files:")
for (f in list.files(out_dir)) message("  ", f)
