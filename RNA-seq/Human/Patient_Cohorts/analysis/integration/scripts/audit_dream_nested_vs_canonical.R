#!/usr/bin/env Rscript
# audit_dream_nested_vs_canonical.R
# ---------------------------------------------------------------------------
# Agent B2 reviewer-defense: compare canonical dream_results.csv against
# dream_results_nested_re.csv and dream_results_ordinal.csv. Emit REPORT.md
# to RNA-seq/results/audit_sensitivity/dream_nested_re/.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

INPUT_ROOT  <- Sys.getenv("MASLD_INPUT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTPUT_ROOT <- Sys.getenv("MASLD_OUTPUT_ROOT", INPUT_ROOT)
IN_RDIR  <- file.path(INPUT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_RDIR <- file.path(OUTPUT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
# Canonical dream_results.csv is in IN_RDIR; nested + ordinal are in OUT_RDIR
RDIR <- IN_RDIR  # back-compat for any reference below
OUT_DIR <- file.path(OUTPUT_ROOT,
                     "RNA-seq/results/audit_sensitivity/dream_nested_re")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== B2 audit: dream nested-RE vs canonical ===\n")

canon_path <- file.path(IN_RDIR, "dream_results.csv")
nest_path  <- file.path(OUT_RDIR, "dream_results_nested_re.csv")
ord_path   <- file.path(OUT_RDIR, "dream_results_ordinal.csv")
meta_path  <- file.path(OUT_RDIR, "dream_results_nested_re_meta.json")

stopifnot(file.exists(canon_path), file.exists(nest_path))
canon <- fread(canon_path)
nest  <- fread(nest_path)
meta_info <- if (file.exists(meta_path)) fromJSON(meta_path) else list()

# Join on gene
m <- merge(canon[, .(gene, logFC_canon = logFC, t_canon = t,
                     padj_canon = padj)],
           nest[,  .(gene, logFC_nest  = logFC, t_nest  = t,
                     padj_nest  = padj, SE_nest = SE,
                     df_nest = df.total, isSingular)],
           by = "gene")
cat("Joined genes:", nrow(m), "\n")

# Per-gene t-stat ratio (nested / canonical)
m[, t_ratio := ifelse(is.finite(t_canon) & t_canon != 0,
                      t_nest / t_canon, NA_real_)]
med_ratio <- median(m$t_ratio, na.rm = TRUE)
q25       <- quantile(m$t_ratio, 0.25, na.rm = TRUE)
q75       <- quantile(m$t_ratio, 0.75, na.rm = TRUE)
n_deflated <- sum(m$t_ratio < 0.7, na.rm = TRUE)
n_inflated <- sum(m$t_ratio > 1.3, na.rm = TRUE)
n_have     <- sum(!is.na(m$t_ratio))

# Singular proportion
n_singular <- sum(m$isSingular, na.rm = TRUE)
n_sing_have <- sum(!is.na(m$isSingular))

# DEG overlap canonical vs nested (padj<0.1)
deg_canon <- canon[padj < 0.1, gene]
deg_nest  <- nest[padj < 0.1, gene]
inter <- length(intersect(deg_canon, deg_nest))
jacc  <- inter / length(union(deg_canon, deg_nest))

# Ordinal overlap (if exists)
ord_summary <- ""
if (file.exists(ord_path)) {
  ord <- fread(ord_path)
  trans_levels <- unique(ord$transition)
  ord_rows <- list()
  for (tr in trans_levels) {
    deg_t <- ord[transition == tr & padj < 0.1, gene]
    n_t <- length(deg_t)
    j_t <- length(intersect(deg_t, deg_canon)) /
           max(1L, length(union(deg_t, deg_canon)))
    ord_rows[[tr]] <- sprintf("| %s | %d | %.4f |", tr, n_t, j_t)
  }
  ord_summary <- paste0(
    "\n## Ordinal vs canonical binary DEG overlap\n\n",
    "| Transition | n DEG (padj<0.1) | Jaccard vs canonical binary |\n",
    "|---|---|---|\n",
    paste(unlist(ord_rows), collapse = "\n"), "\n")
}

# Flagged genes: ratio < 0.7
flagged <- m[t_ratio < 0.7][order(t_ratio)]
flagged_csv <- file.path(OUT_DIR, "flagged_deflated_genes.csv")
fwrite(flagged, flagged_csv)

# Top deflated for table
top_flagged <- head(flagged[, .(gene, logFC_canon, logFC_nest,
                                t_canon, t_nest, t_ratio,
                                padj_canon, padj_nest)], 20)

# Save merged comparison
fwrite(m, file.path(OUT_DIR, "merged_comparison.csv"))

# ---------------------------------------------------------------------------
# REPORT.md
# ---------------------------------------------------------------------------
report_path <- file.path(OUT_DIR, "REPORT.md")

lines <- c(
  "# Dream nested-RE sensitivity (Agent B2)",
  "",
  paste0("Generated: ", as.character(Sys.time())),
  "",
  "## Model metadata",
  "",
  if (length(meta_info) > 0)
    paste0("- Formula: `", meta_info$formula, "`")
  else "- Formula: (meta JSON missing)",
  if (length(meta_info) > 0)
    paste0("- use_nested: ", meta_info$use_nested) else "",
  if (length(meta_info) > 0)
    paste0("- batch_source: ", meta_info$batch_source) else "",
  if (length(meta_info) > 0)
    paste0("- n_samples: ", meta_info$n_samples) else "",
  "",
  "## t-statistic ratio (nested / canonical)",
  "",
  sprintf("- Genes compared: %d", n_have),
  sprintf("- Median ratio: %.4f", med_ratio),
  sprintf("- IQR: [%.4f, %.4f]", q25, q75),
  sprintf("- Genes with ratio < 0.7 (deflated, flagged): %d (%.2f%%)",
          n_deflated, 100 * n_deflated / n_have),
  sprintf("- Genes with ratio > 1.3 (inflated): %d (%.2f%%)",
          n_inflated, 100 * n_inflated / n_have),
  "",
  "## Singular fits",
  "",
  sprintf("- Singular: %d / %d (%.2f%%)", n_singular, n_sing_have,
          if (n_sing_have > 0) 100 * n_singular / n_sing_have else NA),
  "",
  "## DEG overlap (padj<0.1)",
  "",
  sprintf("- Canonical DEGs: %d", length(deg_canon)),
  sprintf("- Nested-RE DEGs: %d", length(deg_nest)),
  sprintf("- Intersection: %d", inter),
  sprintf("- Jaccard: %.4f", jacc),
  ord_summary,
  "",
  "## Top 20 most-deflated genes (lowest t_ratio)",
  "",
  "| gene | logFC_canon | logFC_nest | t_canon | t_nest | t_ratio | padj_canon | padj_nest |",
  "|---|---|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(top_flagged))) {
  r <- top_flagged[i]
  lines <- c(lines, sprintf("| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3g | %.3g |",
                            r$gene, r$logFC_canon, r$logFC_nest,
                            r$t_canon, r$t_nest, r$t_ratio,
                            r$padj_canon, r$padj_nest))
}
lines <- c(lines, "",
           sprintf("Full flagged list: `%s`", flagged_csv),
           sprintf("Merged comparison: `%s`", file.path(OUT_DIR, "merged_comparison.csv")),
           "")

writeLines(lines[nzchar(lines) | TRUE], report_path)
cat("Wrote:", report_path, "\n")
