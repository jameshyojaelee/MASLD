# 248_olink_substate_signature.R
# Phase 2.7 — Plasma signature detection (DOWNGRADED from classifier).
# Per power analysis (240): n_per_grp=15-25 + Bonferroni 1461 → power < 0.18 for d=1.
# Strategy: pre-specified panels (ECM, senescence, inflammation, UPR), single-protein
# tests within panel only (greatly reduced multiple-testing burden), report effect-size CIs.
#
# CONDITIONAL: Run only after Phase 1.7 GO decision.
#
# Inputs:
#   Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt
#   Analysis/Proteomics/results/gse276114_disease_metadata.csv (subject metadata)
#   results/granular_staging/f3_substate_pooled_labels.csv
#
# Outputs (results/granular_staging/):
#   olink_substate_panels.csv  — panel-level + per-protein results
#   olink_substate_summary.md  — pre-specified panel results

suppressPackageStartupMessages({
  library(matrixStats)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
OLINK_PATH <- file.path(PROJECT_ROOT, "Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt")
META_PATH <- file.path(PROJECT_ROOT, "Analysis/Proteomics/results/gse276114_disease_metadata.csv")

cat("== Phase 2.7 Olink plasma signature detection ==\n")

# ---------------------------------------------------------------------------
# Gate check
# ---------------------------------------------------------------------------
gates_path <- file.path(OUT_DIR, "phase1_gates.csv")
if (file.exists(gates_path)) {
  gates <- read.csv(gates_path, stringsAsFactors = FALSE)
  if (any(gates$gate == "OVERALL" & !gates$passes)) {
    cat("Phase 1.7 NO-GO. Phase 2.7 should not run unless overridden.\n")
    if (Sys.getenv("FORCE_PHASE2", "0") != "1") {
      cat("Set FORCE_PHASE2=1 to override.\n")
      quit(status = 0)
    }
  }
}

# ---------------------------------------------------------------------------
# Load Olink + metadata
# ---------------------------------------------------------------------------
if (!file.exists(OLINK_PATH)) {
  cat("Olink data not found; exiting.\n"); quit(status = 0)
}
olink <- read.delim(OLINK_PATH, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("Olink: %d proteins × %d subject columns\n", nrow(olink), ncol(olink) - 1))

# olink first column is "Assay" (protein name)
proteins <- olink[[1]]
expr_mat <- as.matrix(olink[, -1, drop = FALSE])
rownames(expr_mat) <- proteins

# Metadata: gse276114_disease_metadata.csv
if (!file.exists(META_PATH)) {
  cat("Olink metadata not found; exiting.\n"); quit(status = 0)
}
sub_meta <- read.csv(META_PATH, stringsAsFactors = FALSE)
cat(sprintf("Subject metadata: %d rows, columns %s\n",
            nrow(sub_meta), paste(colnames(sub_meta), collapse = ", ")))

# Map matrix_column_name to columns of expr_mat
col_in_meta <- intersect(colnames(expr_mat), sub_meta$matrix_column_name)
cat(sprintf("Subjects matched to metadata: %d / %d\n",
            length(col_in_meta), ncol(expr_mat)))

# ---------------------------------------------------------------------------
# Identify F3 subjects + map sub-state labels
# ---------------------------------------------------------------------------
# Olink subjects are different individuals from the bulk RNA-seq cohort —
# Yang et al. 2025 plasma data. They are NOT joinable per-subject to bulk
# F3 sub-state labels (different patients).
# Strategy: stratify Olink subjects by their own fibrosis/disease grouping
# and test whether F3-class subjects show the F3a/F3b signature axis.
# Then test whether the proportion / score of "F3b-like" plasma profiles
# correlates with disease severity.

if ("disease_group" %in% colnames(sub_meta)) {
  cat("\nDisease group distribution:\n"); print(table(sub_meta$disease_group, useNA = "ifany"))
}
if ("disease" %in% colnames(sub_meta)) {
  cat("Disease distribution:\n"); print(table(sub_meta$disease, useNA = "ifany"))
}

# Filter to MASLD F3-class subjects (where labels permit). The metadata uses
# coarse grouping (disease_group e.g., "MASLD F0-2", "MASLD F3", "MASLD F4", "HCC").
f3_mask <- if ("disease_group" %in% colnames(sub_meta)) {
  grepl("F3$|F3 ", sub_meta$disease_group, ignore.case = TRUE) |
  grepl("^MASLD F3", sub_meta$disease_group, ignore.case = TRUE)
} else rep(FALSE, nrow(sub_meta))

cat(sprintf("\nMASLD F3 plasma subjects: %d\n", sum(f3_mask)))

# ---------------------------------------------------------------------------
# Pre-specified panels (down from 1461 proteins to ~5-30 per panel)
# ---------------------------------------------------------------------------
panels <- list(
  ECM_collagen = c("COL1A1", "COL3A1", "COL5A1", "COL6A3", "FN1", "BGN", "LUM",
                   "DCN", "VCAN", "TIMP1", "TIMP2", "MMP2", "MMP9", "LOX",
                   "LOXL1", "LOXL2", "IGFBP7"),
  Senescence_SenMayo = c("CDKN1A", "CDKN2A", "SERPINE1", "CXCL8", "CCL2", "IL6",
                         "MMP1", "MMP3", "IGFBP3", "IGFBP4", "TNFRSF11B",
                         "ICAM1", "PLAUR"),
  Inflammation = c("IL6", "TNF", "CXCL10", "CCL2", "CXCL8", "CXCL12",
                   "CCR5", "CXCR2", "ICAM1", "VCAM1"),
  UPR_metabolic = c("HSPA5", "DDIT3", "ATF4", "FABP1", "PLIN2"),
  Hepatocyte_secretome = c("AFP", "HGF", "IGFBP1", "IGFBP2", "IGFBP3", "ALB"),
  GDF15_progression = c("GDF15", "PROC", "FBLN1", "MMP10", "GP2")
)

panel_results <- list()
per_protein_results <- list()

# F3 sub-state labels are not directly joinable to plasma subjects, so we
# stratify plasma subjects by disease severity within the available metadata
# and test panel-level associations.
if (sum(f3_mask) >= 6 && "disease_group" %in% colnames(sub_meta)) {
  # Define "F3" group vs "non-F3" disease group (F0-2 + F4 grouped together for power)
  group_f3 <- sub_meta[f3_mask, ]
  group_other_disease <- sub_meta[!f3_mask & sub_meta$disease %in% c("MASLD", "MASH"), ]
  # Match Olink columns
  cols_f3 <- intersect(colnames(expr_mat), group_f3$matrix_column_name)
  cols_other <- intersect(colnames(expr_mat), group_other_disease$matrix_column_name)
  cat(sprintf("F3 plasma columns: %d, non-F3 disease plasma columns: %d\n",
              length(cols_f3), length(cols_other)))

  if (length(cols_f3) >= 5 && length(cols_other) >= 5) {
    for (panel in names(panels)) {
      panel_genes <- panels[[panel]]
      panel_present <- intersect(panel_genes, rownames(expr_mat))
      if (length(panel_present) < 3) {
        cat(sprintf("Skip panel %s: only %d / %d proteins present\n",
                    panel, length(panel_present), length(panel_genes)))
        next
      }
      # Mean panel score per subject
      f3_scores <- colMeans(expr_mat[panel_present, cols_f3, drop = FALSE], na.rm = TRUE)
      other_scores <- colMeans(expr_mat[panel_present, cols_other, drop = FALSE], na.rm = TRUE)
      wt <- suppressWarnings(wilcox.test(f3_scores, other_scores))
      d <- (mean(f3_scores) - mean(other_scores)) /
           sqrt((var(f3_scores) * (length(f3_scores) - 1) +
                 var(other_scores) * (length(other_scores) - 1)) /
                (length(f3_scores) + length(other_scores) - 2))
      panel_results[[panel]] <- data.frame(
        panel = panel, n_proteins_in_panel = length(panel_present),
        n_f3 = length(cols_f3), n_other = length(cols_other),
        mean_f3 = mean(f3_scores), mean_other = mean(other_scores),
        cohens_d = d,
        wilcox_p = wt$p.value,
        stringsAsFactors = FALSE
      )

      # Per-protein within panel (Bonferroni n_proteins_in_panel)
      for (p in panel_present) {
        v_f3 <- as.numeric(expr_mat[p, cols_f3])
        v_other <- as.numeric(expr_mat[p, cols_other])
        if (sum(!is.na(v_f3)) < 3 || sum(!is.na(v_other)) < 3) next
        wt_p <- suppressWarnings(wilcox.test(v_f3, v_other))
        d_p <- (mean(v_f3, na.rm = TRUE) - mean(v_other, na.rm = TRUE)) /
               sqrt((var(v_f3, na.rm = TRUE) + var(v_other, na.rm = TRUE)) / 2)
        per_protein_results[[paste(panel, p, sep = "_")]] <- data.frame(
          panel = panel, protein = p,
          n_f3 = sum(!is.na(v_f3)), n_other = sum(!is.na(v_other)),
          mean_f3 = mean(v_f3, na.rm = TRUE),
          mean_other = mean(v_other, na.rm = TRUE),
          cohens_d = d_p,
          wilcox_p = wt_p$p.value,
          panel_bonferroni_alpha = 0.05 / length(panel_present),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

panel_df <- if (length(panel_results) > 0) do.call(rbind, panel_results) else NULL
prot_df <- if (length(per_protein_results) > 0) do.call(rbind, per_protein_results) else NULL

if (!is.null(panel_df)) {
  panel_df$padj_bh <- p.adjust(panel_df$wilcox_p, method = "BH")
}
if (!is.null(prot_df)) {
  # Within-panel BH adjustment
  prot_df$padj_within_panel <- ave(prot_df$wilcox_p, prot_df$panel,
                                    FUN = function(p) p.adjust(p, method = "BH"))
}

if (!is.null(panel_df)) {
  write.csv(panel_df, file.path(OUT_DIR, "olink_substate_panels.csv"), row.names = FALSE)
  cat("\nPanel-level results:\n"); print(panel_df)
}
if (!is.null(prot_df)) {
  write.csv(prot_df, file.path(OUT_DIR, "olink_substate_proteins.csv"), row.names = FALSE)
  cat(sprintf("\nWrote %d per-protein rows\n", nrow(prot_df)))
}

# ---------------------------------------------------------------------------
# Markdown summary
# ---------------------------------------------------------------------------
sink(file.path(OUT_DIR, "olink_substate_summary.md"))
cat("# Phase 2.7 — Olink plasma signature detection\n\n")
cat("**Pre-specified panel approach** (per power analysis: per-protein Bonferroni vs 1461 has power < 0.2 for d=1.0).\n\n")
cat("Subjects: F3 plasma (n=", if (!is.null(panel_df)) panel_df$n_f3[1] else 0,
    ") vs non-F3 MASLD/MASH (n=", if (!is.null(panel_df)) panel_df$n_other[1] else 0, ")\n\n", sep = "")
if (!is.null(panel_df)) {
  cat("## Panel-level Wilcoxon (BH-adjusted across panels)\n\n")
  cat("| Panel | n proteins | mean F3 | mean other | Cohen's d | Wilcox p | BH padj |\n")
  cat("|---|---:|---:|---:|---:|---:|---:|\n")
  for (i in seq_len(nrow(panel_df))) {
    cat(sprintf("| %s | %d | %.3f | %.3f | %.2f | %.3f | %.3f |\n",
                panel_df$panel[i], panel_df$n_proteins_in_panel[i],
                panel_df$mean_f3[i], panel_df$mean_other[i],
                panel_df$cohens_d[i], panel_df$wilcox_p[i], panel_df$padj_bh[i]))
  }
}
if (!is.null(prot_df)) {
  cat("\n## Top per-protein hits (within-panel BH < 0.05)\n\n")
  sig <- prot_df[prot_df$padj_within_panel < 0.05, ]
  if (nrow(sig) > 0) {
    cat("| Panel | Protein | mean F3 | mean other | d | wilcox p | within-panel padj |\n")
    cat("|---|---|---:|---:|---:|---:|---:|\n")
    sig <- sig[order(sig$padj_within_panel), ]
    for (i in seq_len(min(20, nrow(sig)))) {
      cat(sprintf("| %s | %s | %.3f | %.3f | %.2f | %.3f | %.3f |\n",
                  sig$panel[i], sig$protein[i], sig$mean_f3[i], sig$mean_other[i],
                  sig$cohens_d[i], sig$wilcox_p[i], sig$padj_within_panel[i]))
    }
  } else {
    cat("No per-protein hits at within-panel BH < 0.05.\n")
  }
}
cat("\n## Caveat\n\n")
cat("Yang et al. plasma cohort is independent of bulk RNA-seq cohorts; ",
    "F3 sub-state labels (from 241) cannot be transferred subject-by-subject. ",
    "This phase tests whether F3-staged plasma subjects show ECM / senescence / ",
    "inflammation / UPR-metabolic signal at the panel level vs other-stage MASLD ",
    "plasma subjects. It is a population-level, not subject-level, test.\n")
sink()

cat("\n== Phase 2.7 complete. ==\n")
