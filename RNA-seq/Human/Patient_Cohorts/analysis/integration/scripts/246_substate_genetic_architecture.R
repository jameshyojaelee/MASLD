# 246_substate_genetic_architecture.R
# Phase 2.5 — Genetic architecture per transition.
# Test whether COLOC PP4>0.5 hits are enriched at each transition (F1→F2, F2→F3,
# F3 sub-state-A vs F3 sub-state-B). Uses canonical SuSiE-COLOC posterior column
# 'coloc_best_susie_pp4' from multi_evidence_atlas.csv (28-GWAS portfolio).
#
# CONDITIONAL: Run only after Phase 1.7 GO decision.
#
# Inputs:
#   results/multi_evidence/multi_evidence_atlas.csv (33943 × 218 cols)
#   results/granular_staging/transition_effect_sizes.csv (from 244)
#   results/integration/merged_dge.rds, meta_matched.rds
#   results/granular_staging/f3_substate_pooled_labels.csv
#
# Outputs (results/granular_staging/):
#   substate_genetic_enrichment.csv — Fisher OR per transition × COLOC class
#   substate_genetic_top_genes.csv  — gene-level calls per transition
#   substate_genetic_summary.md     — markdown summary

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
ATLAS <- file.path(PROJECT_ROOT, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

cat("== Phase 2.5 Genetic architecture per transition ==\n")

# ---------------------------------------------------------------------------
# Gate check
# ---------------------------------------------------------------------------
gates_path <- file.path(OUT_DIR, "phase1_gates.csv")
if (file.exists(gates_path)) {
  gates <- read.csv(gates_path, stringsAsFactors = FALSE)
  if (any(gates$gate == "OVERALL" & !gates$passes)) {
    cat("Phase 1.7 NO-GO. Phase 2.5 should not run unless overridden.\n")
    if (Sys.getenv("FORCE_PHASE2", "0") != "1") {
      cat("Set FORCE_PHASE2=1 to override.\n")
      quit(status = 0)
    }
  }
}

# ---------------------------------------------------------------------------
# Load atlas + DE/transition data
# ---------------------------------------------------------------------------
atlas <- read.csv(ATLAS, stringsAsFactors = FALSE)
cat(sprintf("Atlas: %d genes × %d cols\n", nrow(atlas), ncol(atlas)))

# Identify COLOC columns
coloc_cols <- intersect(c("coloc_best_susie_pp4", "coloc_pp4", "coloc_abf_best_pp4",
                          "broadaway_coloc_pp4", "ast_coloc_pp4", "ggt_coloc_pp4",
                          "pdff_coloc_pp4", "ukbb_alt_coloc_pp4",
                          "finngen_nafld_coloc_pp4", "finngen_nash_coloc_pp4",
                          "finngen_hcc_coloc_pp4", "bbj_alt_coloc_pp4"),
                        colnames(atlas))
cat("COLOC columns:", paste(coloc_cols, collapse = ", "), "\n")

# Canonical PP4>0.5 hit set (any coloc column ≥ 0.5)
coloc_mat <- as.matrix(atlas[, coloc_cols, drop = FALSE])
coloc_mat[is.na(coloc_mat)] <- 0
atlas$any_pp4_05 <- rowSums(coloc_mat >= 0.5) > 0
atlas$any_pp4_08 <- rowSums(coloc_mat >= 0.8) > 0
atlas$any_pp4_09 <- rowSums(coloc_mat >= 0.9) > 0
atlas$max_pp4 <- do.call(pmax, c(as.list(as.data.frame(coloc_mat)), list(na.rm = TRUE)))

cat(sprintf("PP4 ≥ 0.5: %d genes\n", sum(atlas$any_pp4_05)))
cat(sprintf("PP4 ≥ 0.8: %d genes\n", sum(atlas$any_pp4_08)))
cat(sprintf("PP4 ≥ 0.9: %d genes\n", sum(atlas$any_pp4_09)))

# ---------------------------------------------------------------------------
# Build transition DEG sets
# ---------------------------------------------------------------------------
# Re-derive transition DEG sets directly from contrasts to avoid dependence
# on intermediate files. This is essentially what 244 did, but more explicit.
atlas_sym <- atlas[, c("ensembl_id", "human_symbol")]
atlas_sym$ensembl_clean <- sub("\\..*", "", atlas_sym$ensembl_id)
ensg_to_sym <- function(ensg) {
  ensg_clean <- sub("\\..*", "", ensg)
  m <- atlas_sym$human_symbol[match(ensg_clean, atlas_sym$ensembl_clean)]
  ifelse(is.na(m) | m == "", ensg_clean, m)
}

dge <- readRDS(file.path(INT_DIR, "merged_dge.rds"))
meta <- readRDS(file.path(INT_DIR, "meta_matched.rds"))
meta <- meta[match(colnames(dge), meta$sample_id), ]
pooled_labels <- read.csv(file.path(OUT_DIR, "f3_substate_pooled_labels.csv"),
                          stringsAsFactors = FALSE)

transition_contrast <- function(idx_a, idx_b) {
  group <- factor(c(rep("A", length(idx_a)), rep("B", length(idx_b))), levels = c("A", "B"))
  cohort <- factor(meta$dataset[c(idx_a, idx_b)])
  dge_sub <- dge[, c(idx_a, idx_b)]
  keep <- filterByExpr(dge_sub, group = group, min.count = 5)
  dge_sub <- dge_sub[keep, , keep.lib.sizes = FALSE]
  dge_sub <- calcNormFactors(dge_sub)
  design <- if (length(unique(cohort)) > 1) model.matrix(~ group + cohort) else model.matrix(~ group)
  v <- voom(dge_sub, design)
  fit <- eBayes(lmFit(v, design))
  res <- topTable(fit, coef = "groupB", number = Inf, sort.by = "none")
  res$gene_id <- rownames(res)
  res$gene_symbol <- ensg_to_sym(res$gene_id)
  res
}

f1 <- which(meta$fibrosis_stage == 1)
f2 <- which(meta$fibrosis_stage == 2)
f3 <- which(meta$fibrosis_stage == 3)

transitions <- list()
cat("\nComputing F1→F2...\n")
transitions$F1F2 <- transition_contrast(f1, f2)
cat("Computing F2→F3...\n")
transitions$F2F3 <- transition_contrast(f2, f3)

substate_levels <- sort(unique(pooled_labels$cluster_pooled))
if (length(substate_levels) >= 2) {
  ca <- substate_levels[1]; cb <- substate_levels[2]
  ids_a <- pooled_labels$sample_id[pooled_labels$cluster_pooled == ca]
  ids_b <- pooled_labels$sample_id[pooled_labels$cluster_pooled == cb]
  idx_a <- match(ids_a, meta$sample_id); idx_a <- idx_a[!is.na(idx_a)]
  idx_b <- match(ids_b, meta$sample_id); idx_b <- idx_b[!is.na(idx_b)]
  cat(sprintf("Computing F3 sub-states: %s (n=%d) vs %s (n=%d)\n",
              ca, length(idx_a), cb, length(idx_b)))
  if (length(idx_a) >= 5 && length(idx_b) >= 5) {
    transitions$F3subAB <- transition_contrast(idx_a, idx_b)
  }
}

# Define transition gene sets (top-1000 by P-value, padj<0.05 fallback)
transition_genes <- lapply(transitions, function(res) {
  sig <- res[!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05, ]
  if (nrow(sig) < 100) {
    # Fallback to top 500 by p-value
    sig <- res[order(res$P.Value), ][1:500, ]
  }
  sig$gene_symbol[!is.na(sig$gene_symbol) & sig$gene_symbol != ""]
})
cat("\nTransition gene set sizes:\n")
print(sapply(transition_genes, length))

# ---------------------------------------------------------------------------
# Fisher enrichment per transition × PP4 class
# ---------------------------------------------------------------------------
all_atlas_genes <- atlas$human_symbol
n_atlas <- nrow(atlas)

run_fisher <- function(set_genes, hit_set, label) {
  set_in <- atlas$human_symbol %in% set_genes
  hit_in <- atlas$human_symbol %in% hit_set
  set_hits <- sum(set_in & hit_in)
  set_nonhits <- sum(set_in & !hit_in)
  nonset_hits <- sum(!set_in & hit_in)
  nonset_nonhits <- sum(!set_in & !hit_in)
  tab <- matrix(c(set_hits, set_nonhits, nonset_hits, nonset_nonhits), nrow = 2)
  ft <- tryCatch(fisher.test(tab, alternative = "greater"), error = function(e) NULL)
  data.frame(
    transition = label,
    set_size = sum(set_in), hit_size = sum(hit_in),
    set_hits = set_hits, set_nonhits = set_nonhits,
    nonset_hits = nonset_hits, nonset_nonhits = nonset_nonhits,
    OR = if (!is.null(ft)) as.numeric(ft$estimate) else NA_real_,
    OR_ci_low = if (!is.null(ft)) ft$conf.int[1] else NA_real_,
    OR_ci_high = if (!is.null(ft)) ft$conf.int[2] else NA_real_,
    p_value = if (!is.null(ft)) ft$p.value else NA_real_,
    stringsAsFactors = FALSE
  )
}

hit_classes <- list(
  "PP4>=0.5" = atlas$human_symbol[atlas$any_pp4_05],
  "PP4>=0.8" = atlas$human_symbol[atlas$any_pp4_08],
  "PP4>=0.9" = atlas$human_symbol[atlas$any_pp4_09]
)

enrich_rows <- list()
for (transition in names(transition_genes)) {
  for (cls in names(hit_classes)) {
    enrich_rows[[paste(transition, cls)]] <- run_fisher(
      transition_genes[[transition]],
      hit_classes[[cls]],
      paste0(transition, " : ", cls)
    )
  }
}
enrich <- do.call(rbind, enrich_rows)
enrich$padj <- p.adjust(enrich$p_value, method = "BH")
write.csv(enrich, file.path(OUT_DIR, "substate_genetic_enrichment.csv"), row.names = FALSE)
cat("\nGenetic enrichment per transition:\n")
print(enrich)

# ---------------------------------------------------------------------------
# Gene-level calls: which COLOC genes are most differential per transition?
# ---------------------------------------------------------------------------
gene_call_rows <- list()
hit_genes <- atlas$human_symbol[atlas$any_pp4_05]

for (transition in names(transitions)) {
  res <- transitions[[transition]]
  hit_in_res <- res[res$gene_symbol %in% hit_genes & !is.na(res$adj.P.Val) &
                    res$adj.P.Val < 0.1, ]
  if (nrow(hit_in_res) == 0) next
  hit_in_res <- hit_in_res[order(hit_in_res$adj.P.Val), ]
  top_hits <- head(hit_in_res, 30)
  top_hits$transition <- transition
  # Add max_pp4
  top_hits$max_pp4 <- atlas$max_pp4[match(top_hits$gene_symbol, atlas$human_symbol)]
  gene_call_rows[[transition]] <- top_hits[, c("transition", "gene_symbol", "logFC", "P.Value", "adj.P.Val", "max_pp4")]
}
gene_calls <- if (length(gene_call_rows) > 0) do.call(rbind, gene_call_rows) else NULL
if (!is.null(gene_calls)) {
  write.csv(gene_calls, file.path(OUT_DIR, "substate_genetic_top_genes.csv"), row.names = FALSE)
  cat(sprintf("\nWrote %d gene-level calls\n", nrow(gene_calls)))
}

# ---------------------------------------------------------------------------
# Markdown summary
# ---------------------------------------------------------------------------
sink(file.path(OUT_DIR, "substate_genetic_summary.md"))
cat("# Phase 2.5 — Genetic architecture per transition\n\n")
cat("## Fisher enrichment of COLOC hits\n\n")
cat("| Transition : Class | Set size | Hit overlap | OR | 95% CI | p | BH-padj |\n")
cat("|---|---:|---:|---:|---|---:|---:|\n")
for (i in seq_len(nrow(enrich))) {
  cat(sprintf("| %s | %d | %d | %.2f | (%.2f, %.2f) | %.2e | %.2e |\n",
              enrich$transition[i], enrich$set_size[i], enrich$set_hits[i],
              enrich$OR[i], enrich$OR_ci_low[i], enrich$OR_ci_high[i],
              enrich$p_value[i], enrich$padj[i]))
}
if (!is.null(gene_calls) && nrow(gene_calls) > 0) {
  cat("\n## Top COLOC-supported gene calls per transition (top 10 each)\n\n")
  for (transition in unique(gene_calls$transition)) {
    cat(sprintf("\n### %s\n\n", transition))
    sub <- head(gene_calls[gene_calls$transition == transition, ], 10)
    cat("| Gene | logFC | padj | max PP4 |\n|---|---:|---:|---:|\n")
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("| %s | %.2f | %.2e | %.2f |\n",
                  sub$gene_symbol[i], sub$logFC[i], sub$adj.P.Val[i], sub$max_pp4[i]))
    }
  }
}
sink()

cat("\n== Phase 2.5 complete. ==\n")
