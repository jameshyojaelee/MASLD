# 244_transition_effect_sizes.R
# Phase 1.4 — Effect-size comparison: F1→F2 vs F3 sub-state-A → F3 sub-state-B.
# Determines narrative framing:
#   Effect ≥ 30% of F1→F2  → "hierarchical transitions" framing supported
#   Effect <  30% of F1→F2 → "F2 switch with F3 heterogeneity" framing only
#
# Inputs:
#   results/integration/merged_dge.rds, meta_matched.rds
#   results/granular_staging/f3_substate_pooled_labels.csv (from 241)
#
# Outputs (results/granular_staging/):
#   transition_effect_sizes.csv    — per-transition median |log2FC| of top-1000 DEGs
#   transition_pathway_counts.csv  — # significant pathways at each transition (fgsea)
#   transition_effect_sizes_summary.md  — narrative recommendation

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(matrixStats)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== Phase 1.4 Transition effect-size comparison ==\n")

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
dge <- readRDS(file.path(INT_DIR, "merged_dge.rds"))
meta <- readRDS(file.path(INT_DIR, "meta_matched.rds"))
meta <- meta[match(colnames(dge), meta$sample_id), ]
pooled_labels <- read.csv(file.path(OUT_DIR, "f3_substate_pooled_labels.csv"),
                          stringsAsFactors = FALSE)

# Map symbols using multi-evidence atlas
atlas_sym <- read.csv(file.path(PROJECT_ROOT, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
                      stringsAsFactors = FALSE)[, c("ensembl_id", "human_symbol")]
atlas_sym$ensembl_clean <- sub("\\..*", "", atlas_sym$ensembl_id)
ensg_to_sym <- function(ensg) {
  ensg_clean <- sub("\\..*", "", ensg)
  m <- atlas_sym$human_symbol[match(ensg_clean, atlas_sym$ensembl_clean)]
  ifelse(is.na(m) | m == "", ensg_clean, m)
}

# ---------------------------------------------------------------------------
# Helper: limma-voom contrast given two sample-index vectors
# ---------------------------------------------------------------------------
run_limma_contrast <- function(idx_a, idx_b, label_a = "A", label_b = "B") {
  use_idx <- c(idx_a, idx_b)
  group <- c(rep("A", length(idx_a)), rep("B", length(idx_b)))
  group <- factor(group, levels = c("A", "B"))
  cohort <- factor(meta$dataset[use_idx])

  dge_sub <- dge[, use_idx]
  keep <- filterByExpr(dge_sub, group = group, min.count = 5)
  dge_sub <- dge_sub[keep, , keep.lib.sizes = FALSE]
  dge_sub <- calcNormFactors(dge_sub)

  # Fit B vs A with cohort blocking if multi-cohort
  if (length(unique(cohort)) > 1) {
    design <- model.matrix(~ group + cohort)
  } else {
    design <- model.matrix(~ group)
  }
  v <- voom(dge_sub, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = "groupB", number = Inf, sort.by = "none")
  res$gene_id <- rownames(res)
  res$gene_symbol <- ensg_to_sym(res$gene_id)
  res$contrast <- paste0(label_b, "_vs_", label_a)
  res
}

summarise_effect <- function(res, top_n = 1000) {
  res_sorted <- res[order(res$P.Value), ]
  top <- head(res_sorted, top_n)
  tibble_like <- list(
    contrast = unique(res$contrast),
    n_total_genes = nrow(res),
    n_padj_05 = sum(!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05),
    n_padj_05_lfc05 = sum(!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05 &
                          abs(res$logFC) > 0.5),
    median_abs_logFC_top1000 = median(abs(top$logFC)),
    mean_abs_logFC_top1000 = mean(abs(top$logFC)),
    median_abs_logFC_padj05 = median(abs(res$logFC[!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05])),
    max_abs_logFC = max(abs(res$logFC), na.rm = TRUE)
  )
  as.data.frame(tibble_like, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Anchor 1: F1 vs F2 transition (the main "F2 switch")
# ---------------------------------------------------------------------------
cat("\n--- Computing F1 → F2 transition (the F2 switch) ---\n")
f1_idx <- which(meta$fibrosis_stage == 1)
f2_idx <- which(meta$fibrosis_stage == 2)
cat(sprintf("F1: %d  F2: %d\n", length(f1_idx), length(f2_idx)))

f1f2_res <- run_limma_contrast(f1_idx, f2_idx, "F1", "F2")
f1f2_summary <- summarise_effect(f1f2_res)
print(f1f2_summary)

# ---------------------------------------------------------------------------
# Anchor 2: F2 vs F3 transition (independent comparator)
# ---------------------------------------------------------------------------
cat("\n--- Computing F2 → F3 transition ---\n")
f3_idx <- which(meta$fibrosis_stage == 3)
f2f3_res <- run_limma_contrast(f2_idx, f3_idx, "F2", "F3")
f2f3_summary <- summarise_effect(f2f3_res)
print(f2f3_summary)

# ---------------------------------------------------------------------------
# Anchor 3: F3 sub-state A vs sub-state B (from 241 pooled labels)
# ---------------------------------------------------------------------------
substate_levels <- sort(unique(pooled_labels$cluster_pooled))
cat(sprintf("\n--- Computing F3 sub-state contrasts: %s ---\n",
            paste(substate_levels, collapse = " vs ")))

substate_summaries <- list()
substate_results <- list()

if (length(substate_levels) >= 2) {
  for (i in seq_along(substate_levels)[-length(substate_levels)]) {
    for (j in (i+1):length(substate_levels)) {
      ca <- substate_levels[i]; cb <- substate_levels[j]
      ids_a <- pooled_labels$sample_id[pooled_labels$cluster_pooled == ca]
      ids_b <- pooled_labels$sample_id[pooled_labels$cluster_pooled == cb]
      idx_a <- match(ids_a, meta$sample_id)
      idx_b <- match(ids_b, meta$sample_id)
      idx_a <- idx_a[!is.na(idx_a)]; idx_b <- idx_b[!is.na(idx_b)]
      if (length(idx_a) < 5 || length(idx_b) < 5) {
        cat(sprintf("Skip %s vs %s (n=%d/%d too small)\n", ca, cb, length(idx_a), length(idx_b)))
        next
      }
      cat(sprintf("%s (n=%d) vs %s (n=%d)\n", ca, length(idx_a), cb, length(idx_b)))
      res_ij <- run_limma_contrast(idx_a, idx_b, ca, cb)
      substate_results[[paste0(ca, "_vs_", cb)]] <- res_ij
      substate_summaries[[paste0(ca, "_vs_", cb)]] <- summarise_effect(res_ij)
    }
  }
}

# ---------------------------------------------------------------------------
# Compile effect-size table
# ---------------------------------------------------------------------------
all_summaries <- rbind(f1f2_summary, f2f3_summary)
if (length(substate_summaries) > 0) {
  all_summaries <- rbind(all_summaries, do.call(rbind, substate_summaries))
}
write.csv(all_summaries, file.path(OUT_DIR, "transition_effect_sizes.csv"), row.names = FALSE)
cat("\nEffect-size table:\n")
print(all_summaries)

# ---------------------------------------------------------------------------
# fgsea pathway count per transition
# ---------------------------------------------------------------------------
cat("\n--- fgsea pathway counts per transition ---\n")
msig <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) suppressWarnings(msigdbr(species = "Homo sapiens", category = "H"))
)
gene_col <- intersect(c("gene_symbol", "human_gene_symbol", "gs_symbol"), colnames(msig))[1]
set_col <- intersect(c("gs_name", "gs_id"), colnames(msig))[1]
hallmark_sets <- split(msig[[gene_col]], msig[[set_col]])
hallmark_sets <- lapply(hallmark_sets, unique)

run_fgsea <- function(res, label) {
  ranks <- res$logFC
  names(ranks) <- res$gene_symbol
  ranks <- ranks[!duplicated(names(ranks)) & !is.na(names(ranks)) & names(ranks) != ""]
  ranks <- sort(ranks, decreasing = TRUE)
  fgr <- suppressWarnings(fgsea(pathways = hallmark_sets, stats = ranks, minSize = 10, maxSize = 500))
  data.frame(
    contrast = label,
    n_pathways_total = nrow(fgr),
    n_pathways_padj_05 = sum(fgr$padj < 0.05, na.rm = TRUE),
    n_pathways_padj_01 = sum(fgr$padj < 0.01, na.rm = TRUE),
    median_abs_NES_top10 = median(abs(head(fgr[order(fgr$padj), ], 10)$NES), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

pwy_rows <- list(run_fgsea(f1f2_res, "F2_vs_F1"),
                 run_fgsea(f2f3_res, "F3_vs_F2"))
for (n in names(substate_results)) {
  pwy_rows[[n]] <- run_fgsea(substate_results[[n]], n)
}
pwy_df <- do.call(rbind, pwy_rows)
write.csv(pwy_df, file.path(OUT_DIR, "transition_pathway_counts.csv"), row.names = FALSE)
cat("Pathway counts:\n"); print(pwy_df)

# ---------------------------------------------------------------------------
# Decision: hierarchical vs single-switch framing
# ---------------------------------------------------------------------------
f1f2_eff <- f1f2_summary$median_abs_logFC_top1000
substate_eff <- if (length(substate_summaries) > 0) {
  max(sapply(substate_summaries, function(s) s$median_abs_logFC_top1000), na.rm = TRUE)
} else 0
ratio <- substate_eff / f1f2_eff
gate_pass <- ratio >= 0.30

cat("\n== Decision summary ==\n")
cat(sprintf("F1→F2 median |log2FC| (top 1000 DEGs): %.3f\n", f1f2_eff))
cat(sprintf("F3 sub-state max effect size:           %.3f\n", substate_eff))
cat(sprintf("Ratio (F3 sub / F1→F2):                 %.2f\n", ratio))
cat(sprintf("Gate (ratio ≥ 0.30):                    %s\n", ifelse(gate_pass, "PASS", "FAIL")))
cat(sprintf("Recommended narrative framing:          %s\n",
            ifelse(gate_pass, "HIERARCHICAL TRANSITIONS",
                   "F2 SWITCH WITH F3 HETEROGENEITY")))

# Markdown summary
sink(file.path(OUT_DIR, "transition_effect_sizes_summary.md"))
cat("# Phase 1.4 — Transition effect-size comparison\n\n")
cat("## Effect-size table (median |log2FC| of top 1000 DEGs)\n\n")
cat("| Contrast | n_padj<0.05 | n |LFC|>0.5 | median |LFC| top1000 | mean |LFC| top1000 |\n")
cat("|---|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(all_summaries))) {
  cat(sprintf("| %s | %d | %d | %.3f | %.3f |\n",
              all_summaries$contrast[i], all_summaries$n_padj_05[i],
              all_summaries$n_padj_05_lfc05[i],
              all_summaries$median_abs_logFC_top1000[i],
              all_summaries$mean_abs_logFC_top1000[i]))
}
cat("\n## fgsea Hallmark pathway count per transition\n\n")
cat("| Contrast | Total | padj<0.05 | padj<0.01 | median |NES| top10 |\n")
cat("|---|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(pwy_df))) {
  cat(sprintf("| %s | %d | %d | %d | %.2f |\n",
              pwy_df$contrast[i], pwy_df$n_pathways_total[i],
              pwy_df$n_pathways_padj_05[i], pwy_df$n_pathways_padj_01[i],
              pwy_df$median_abs_NES_top10[i]))
}
cat(sprintf("\n## Narrative gate\n\n"))
cat(sprintf("- F1→F2 effect size: %.3f\n", f1f2_eff))
cat(sprintf("- F3 sub-state max effect: %.3f\n", substate_eff))
cat(sprintf("- Ratio: %.2f (gate ≥ 0.30 → hierarchical framing)\n", ratio))
cat(sprintf("- Decision: **%s**\n",
            ifelse(gate_pass, "hierarchical transitions framing supported",
                   "single F2 switch with F3 heterogeneity framing only")))
sink()

cat(sprintf("\nWrote %s, %s, %s\n",
            file.path(OUT_DIR, "transition_effect_sizes.csv"),
            file.path(OUT_DIR, "transition_pathway_counts.csv"),
            file.path(OUT_DIR, "transition_effect_sizes_summary.md")))
cat("== Phase 1.4 complete. ==\n")
