#!/usr/bin/env Rscript
# Species Effect Size Normalization + Cross-Species Scatter
#
# Scientific question:
#   Do mouse hepatocytes show concordant disease programs to human hepatocytes?
#
# Method:
#   1. Human hepatocyte DE: sc pseudobulk (Hepatocytes_de.csv) — gene symbols
#   2. Mouse hepatocyte DE: bulk mouse dream (more powerful than sc n=2)
#   3. SD normalization: logFC_norm = logFC / sd(logFC) per study
#   4. Gene matching: human sc gene symbols → multi-evidence atlas
#      (which has human_symbol + ensembl_id + mouse_ortholog + mouse_meta_logFC)
#   5. Cross-species scatter: human logFC_norm vs mouse logFC_norm
#
# Note: The human sc pseudobulk DE uses gene symbols as identifiers
#   (same as CellRanger v3+ var_names from GENCODE). The multi-evidence
#   atlas joins by human_symbol, providing mouse_meta_logFC for matching.
#
# Key expected result: Spearman rho > 0.3 for hepatocyte concordance
#
# Input:
#   - results_gpu_v2/pseudobulk_de/Hepatocytes_de.csv  (human sc)
#   - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (cross-species)
#   - results_gpu_v2/mouse_sc/Hepatocytes_de.csv  (mouse sc, if available)
#
# Output: results_gpu_v2/mouse_sc/
#   - species_effect_norm_hepatocytes.csv
#   - cross_species_sc_scatter.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

sc_dir    <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2")
mouse_dir <- file.path(sc_dir, "mouse_sc")
out_dir   <- mouse_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Load human sc hepatocyte DE
# ---------------------------------------------------------------------------
human_sc_path <- file.path(sc_dir, "pseudobulk_de", "Hepatocytes_de.csv")
if (!file.exists(human_sc_path)) {
  stop("Human sc hepatocyte DE not found: ", human_sc_path)
}

human_sc <- fread(human_sc_path)
message("Human sc Hepatocytes DE: ", nrow(human_sc), " genes")

# Normalize column names
if ("logFC" %in% colnames(human_sc)) setnames(human_sc, "logFC", "human_sc_logFC")
if ("gene" %in% colnames(human_sc))  setnames(human_sc, "gene",  "human_symbol")
if ("padj" %in% colnames(human_sc))  setnames(human_sc, "padj",  "human_sc_padj")

# The gene column contains gene symbols (DDX11L2, ALB, ...) plus some ENSG IDs
# for genes without a proper symbol. Both are valid for matching to the atlas.
human_sc_logfc_sd <- sd(human_sc$human_sc_logFC, na.rm = TRUE)
human_sc[, human_sc_logFC_norm := human_sc_logFC / human_sc_logfc_sd]
message("  logFC SD: ", round(human_sc_logfc_sd, 3))

# ---------------------------------------------------------------------------
# 2. Load multi-evidence atlas (cross-species bridge)
#    Has: human_symbol, ensembl_id, mouse_ortholog, bulk_logFC,
#         mouse_meta_logFC, mouse_meta_padj
# ---------------------------------------------------------------------------
atlas_path <- file.path(BASE,
  "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
if (!file.exists(atlas_path)) {
  stop("Multi-evidence atlas not found: ", atlas_path)
}

atlas <- fread(atlas_path,
               select = c("human_symbol", "ensembl_id", "mouse_ortholog",
                          "bulk_logFC", "bulk_padj", "bulk_tstat",
                          "mouse_meta_logFC", "mouse_meta_padj",
                          "n_diets_sig", "is_conserved",
                          "primary_category"))
stopifnot(all(c("bulk_padj", "bulk_logFC") %in% names(atlas)))
message("Multi-evidence atlas: ", nrow(atlas), " genes")
message("  Columns: ", paste(colnames(atlas), collapse = ", "))

# Compute SD-normalized mouse bulk logFC
atlas_mouse_sd <- sd(atlas$mouse_meta_logFC, na.rm = TRUE)
atlas[, mouse_bulk_logFC_norm := mouse_meta_logFC / atlas_mouse_sd]

# Also compute SD-normalized human bulk logFC for comparison
atlas_human_sd <- sd(atlas$bulk_logFC, na.rm = TRUE)
atlas[, human_bulk_logFC_norm := bulk_logFC / atlas_human_sd]

message("  Mouse bulk logFC SD: ", round(atlas_mouse_sd, 3))
message("  Human bulk logFC SD: ", round(atlas_human_sd, 3))

# ---------------------------------------------------------------------------
# 3. Merge human sc DE → atlas (by human_symbol)
# ---------------------------------------------------------------------------
message("\nMerging human sc DE with multi-evidence atlas (by human_symbol)...")

merged <- merge(
  human_sc[, .(human_symbol, human_sc_logFC, human_sc_logFC_norm,
               human_sc_padj)],
  atlas[, .(human_symbol, ensembl_id, mouse_ortholog,
            bulk_logFC, human_bulk_logFC_norm,
            mouse_meta_logFC, mouse_bulk_logFC_norm,
            mouse_meta_padj, n_diets_sig, is_conserved,
            primary_category)],
  by = "human_symbol", all = FALSE
)
message("  Matched genes: ", nrow(merged))
message("  With mouse bulk logFC: ", sum(!is.na(merged$mouse_meta_logFC)))

# ---------------------------------------------------------------------------
# 4. Load mouse sc hepatocyte DE (optional/supplementary)
# ---------------------------------------------------------------------------
mouse_sc_path <- file.path(mouse_dir, "Hepatocytes_de.csv")
mouse_sc_available <- file.exists(mouse_sc_path)

if (mouse_sc_available) {
  mouse_sc_de <- fread(mouse_sc_path)
  message("\nMouse sc Hepatocytes DE: ", nrow(mouse_sc_de),
          " genes (n=2, exploratory)")
  if ("logFC" %in% colnames(mouse_sc_de)) {
    setnames(mouse_sc_de, "logFC", "mouse_sc_logFC")
  }
  gene_col_sc <- intersect(c("gene", "gene_id"), colnames(mouse_sc_de))[1]
  if (!is.na(gene_col_sc)) setnames(mouse_sc_de, gene_col_sc, "mouse_symbol")
  mouse_sc_sd <- sd(mouse_sc_de$mouse_sc_logFC, na.rm = TRUE)
  mouse_sc_de[, mouse_sc_logFC_norm := mouse_sc_logFC / mouse_sc_sd]
  # Mouse sc uses ENSMUSG IDs — match to atlas via mouse_ortholog column
  # (atlas$mouse_ortholog is the mouse gene symbol, not ENSMUSG)
  # TODO: if needed, join via ortholog map; for now this is supplementary
} else {
  message("\nMouse sc Hepatocytes DE not found — bulk comparison only.")
}

# ---------------------------------------------------------------------------
# 5. Save merged table
# ---------------------------------------------------------------------------
fwrite(merged, file.path(out_dir, "species_effect_norm_hepatocytes.csv"))
message("\nSaved: species_effect_norm_hepatocytes.csv (", nrow(merged), " genes)")

# ---------------------------------------------------------------------------
# 6. Cross-species scatter plots
# ---------------------------------------------------------------------------
message("\nGenerating cross-species scatter plots...")

plot_scatter <- function(df, x_col, y_col, x_label, y_label,
                         title, subtitle, out_path,
                         color_col = NULL) {
  df_plot <- df[!is.na(get(x_col)) & !is.na(get(y_col))]
  if (nrow(df_plot) < 10) {
    message("  Skipping ", basename(out_path),
            ": insufficient data (", nrow(df_plot), " genes)")
    return(invisible(NULL))
  }

  rho <- cor(df_plot[[x_col]], df_plot[[y_col]],
             method = "spearman", use = "complete.obs")
  r   <- cor(df_plot[[x_col]], df_plot[[y_col]],
             method = "pearson",  use = "complete.obs")
  n   <- nrow(df_plot)

  df_plot[, concordant := sign(get(x_col)) == sign(get(y_col))]
  prop_conc <- mean(df_plot$concordant, na.rm = TRUE)

  if (is.null(color_col)) {
    plt <- ggplot(df_plot, aes(x = get(x_col), y = get(y_col))) +
      geom_point(aes(color = concordant), size = 0.6, alpha = 0.5)
  } else {
    plt <- ggplot(df_plot, aes(x = get(x_col), y = get(y_col))) +
      geom_point(aes(color = get(color_col)), size = 0.7, alpha = 0.6)
  }

  plt <- plt +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_smooth(method = "lm", formula = y ~ x, color = "#2C3E50",
                se = TRUE, linewidth = 0.8) +
    scale_color_manual(
      values = if (is.null(color_col))
        c("TRUE" = "#27AE60", "FALSE" = "#E74C3C")
      else c("TRUE" = "#C0392B", "FALSE" = "#7F8C8D"),
      labels = if (is.null(color_col))
        c("TRUE" = "Concordant", "FALSE" = "Discordant")
      else c("TRUE" = "Conserved", "FALSE" = "Other"),
      name = NULL
    ) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = sprintf(
               "Spearman rho = %.3f\nPearson r = %.3f\nn = %d genes\n%.1f%% concordant",
               rho, r, n, 100 * prop_conc),
             size = 3.2, family = "mono") +
    labs(title = title, subtitle = subtitle,
         x = x_label, y = y_label) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")

  ggsave(out_path, plt, width = 6, height = 5.5)
  message("  Saved: ", basename(out_path),
          sprintf(" (rho=%.3f, r=%.3f, n=%d)", rho, r, n))
  invisible(rho)
}

# Primary: human sc vs mouse bulk
rho_main <- plot_scatter(
  df       = merged,
  x_col    = "human_sc_logFC_norm",
  y_col    = "mouse_bulk_logFC_norm",
  x_label  = "Human sc Hepatocyte logFC (SD-normalized)",
  y_label  = "Mouse bulk Hepatocyte dream logFC (SD-normalized)",
  title    = "Cross-species hepatocyte concordance",
  subtitle = "Human single-cell pseudobulk vs Mouse bulk dream pooled | Orthologous genes",
  out_path = file.path(out_dir, "cross_species_sc_scatter.pdf")
)

# Supplementary: human bulk vs mouse bulk (both from atlas)
plot_scatter(
  df       = merged[!is.na(bulk_logFC) & !is.na(mouse_meta_logFC)],
  x_col    = "human_bulk_logFC_norm",
  y_col    = "mouse_bulk_logFC_norm",
  x_label  = "Human bulk dream logFC (SD-normalized)",
  y_label  = "Mouse bulk dream logFC (SD-normalized)",
  title    = "Cross-species bulk concordance (reference)",
  subtitle = "Human bulk dream (10-cohort) vs Mouse bulk dream | Orthologous genes",
  out_path = file.path(out_dir, "cross_species_bulk_bulk_scatter.pdf")
)

# Conserved highlighted
if ("is_conserved" %in% colnames(merged)) {
  plot_scatter(
    df        = merged[!is.na(human_sc_logFC_norm) & !is.na(mouse_bulk_logFC_norm)],
    x_col     = "human_sc_logFC_norm",
    y_col     = "mouse_bulk_logFC_norm",
    x_label   = "Human sc Hepatocyte logFC (SD-normalized)",
    y_label   = "Mouse bulk Hepatocyte logFC (SD-normalized)",
    title     = "Conserved genes in cross-species scatter",
    subtitle  = "Red = Conserved (723 genes); cross-species conserved",
    out_path  = file.path(out_dir, "cross_species_sc_scatter_conserved.pdf"),
    color_col = "is_conserved"
  )
}

# ---------------------------------------------------------------------------
# 7. Effect size distribution comparison
# ---------------------------------------------------------------------------
plot_ef_dist <- function(merged, out_path) {
  df_long <- rbind(
    data.table(logFC_norm = merged[!is.na(human_sc_logFC_norm)]$human_sc_logFC_norm,
               source = "Human sc Hepatocytes"),
    data.table(logFC_norm = merged[!is.na(mouse_bulk_logFC_norm)]$mouse_bulk_logFC_norm,
               source = "Mouse bulk Hepatocytes")
  )
  if (nrow(df_long) == 0) return(invisible(NULL))

  hist_plt <- ggplot(df_long, aes(x = logFC_norm, fill = source)) +
    geom_histogram(binwidth = 0.05, alpha = 0.6, position = "identity") +
    scale_fill_manual(values = c("Human sc Hepatocytes" = "#3498DB",
                                 "Mouse bulk Hepatocytes" = "#E67E22"),
                      name = NULL) +
    labs(title = "SD-normalized effect size distributions",
         x = "logFC / SD(logFC)", y = "Count") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")

  ggsave(out_path, hist_plt, width = 6, height = 4)
  message("  Saved: ", basename(out_path))
}

plot_ef_dist(
  merged,
  out_path = file.path(out_dir, "species_effect_size_distribution.pdf")
)

# ---------------------------------------------------------------------------
# 8. Summary statistics per concordance tier
# ---------------------------------------------------------------------------
n_with_both <- sum(!is.na(merged$human_sc_logFC_norm) &
                   !is.na(merged$mouse_bulk_logFC_norm))
n_concordant <- sum(sign(merged$human_sc_logFC) == sign(merged$mouse_meta_logFC),
                    na.rm = TRUE)

message("\n=== Module 2 (Species Effect Size Normalization) complete ===")
message("Key stats:")
message("  Genes matched (human sc ∩ atlas): ", nrow(merged))
message("  Genes with both human sc + mouse bulk logFC: ", n_with_both)
message("  Concordant direction: ",
        round(100 * n_concordant / n_with_both, 1), "%")
if (!is.null(rho_main)) {
  message("  Spearman rho (human sc vs mouse bulk): ", round(rho_main, 3))
  if (rho_main > 0.3) {
    message("  -> Concordant hepatocyte disease programs detected")
  }
}
message("\nResults in: ", out_dir)
