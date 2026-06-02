#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# 57b_composite_score_sensitivity.R
#
# A7: Composite priority score sensitivity (ATAC improvement plan)
#
# Script 57 currently uses an arbitrary composite gene-level score:
#   gwas_atac_regulatory_score = pip_norm * (0.4 + 0.3 * da_norm + 0.3 * motif_norm)
#
# Here we evaluate 5 alternative scoring variants + 1 elastic-net data-driven
# weighting against COLOC PP.H4 > 0.5 as the supervised signal, and quantify
# top-50 list stability via pairwise Jaccard + Spearman.
#
# Outputs (RNA-seq/results/audit_sensitivity/atac/):
#   composite_score_sensitivity.csv     — per-gene 6 score columns + features
#   composite_score_top50_overlap.csv   — Jaccard/Spearman matrices, capture
#   elasticnet_weights.csv              — CV weights + 100-fold bootstrap CIs
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(glmnet)
})

set.seed(42)

# ── Paths ───────────────────────────────────────────────────────────────────
PROJ      <- Sys.getenv("MASLD_PROJECT_ROOT",
                        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATAC_DIR  <- file.path(PROJ, "GWAS/finemapping/results/gwas_atac")
COLOC_F   <- file.path(PROJ, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
OUT_DIR   <- file.path(PROJ, "RNA-seq/results/audit_sensitivity/atac")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat("A7: Composite priority score sensitivity\n")
cat("================================================================\n")
cat("OUT_DIR:", OUT_DIR, "\n\n")

# ── Load variant annotation + motif disruption ─────────────────────────────
ann <- fread(file.path(ATAC_DIR, "gwas_atac_variant_annotation.csv"))
cat("Loaded variant annotations:", nrow(ann), "variant×CT overlaps\n")

mot_file <- file.path(ATAC_DIR, "motif_disruption_scores.csv")
if (file.exists(mot_file)) {
  motif <- fread(mot_file)
  cat("Loaded motif disruption:", nrow(motif), "variant-motif pairs\n")
} else {
  stop("Missing motif_disruption_scores.csv")
}

# ── Build gene assignment (mirror Script 57 logic exactly) ──────────────────
ann_df <- as.data.frame(ann)
gene_variant_map <- ann_df %>%
  mutate(
    assigned_gene = case_when(
      !is.na(scenic_target_gene) & nzchar(scenic_target_gene) ~ scenic_target_gene,
      !is.na(linked_gene)        & nzchar(linked_gene)        ~ linked_gene,
      !is.na(nearest_gene)       & nzchar(nearest_gene) &
        !is.na(distance_to_tss)  & distance_to_tss <= 500000  ~ nearest_gene,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(assigned_gene), nzchar(assigned_gene))
cat("Gene-to-variant mapping:", nrow(gene_variant_map), "entries\n")

# ── Gene-level aggregation (raw features) ──────────────────────────────────
gene_pip <- gene_variant_map %>%
  group_by(assigned_gene) %>%
  summarise(max_pip = max(max_pip, na.rm = TRUE),
            n_variants = n_distinct(variant_id),
            .groups = "drop")

gene_da <- gene_variant_map %>%
  filter(!is.na(hep_da_logFC)) %>%
  group_by(assigned_gene) %>%
  summarise(max_da_abs = max(abs(hep_da_logFC), na.rm = TRUE), .groups = "drop")

# Motif → gene mapping (same logic as Script 57)
variant_gene <- ann_df %>%
  mutate(assigned_gene = case_when(
    !is.na(scenic_target_gene) & nzchar(scenic_target_gene) ~ scenic_target_gene,
    !is.na(linked_gene)        & nzchar(linked_gene)        ~ linked_gene,
    !is.na(nearest_gene)       & nzchar(nearest_gene) &
      !is.na(distance_to_tss)  & distance_to_tss <= 500000  ~ nearest_gene,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(assigned_gene)) %>%
  select(variant_id, assigned_gene) %>%
  distinct()

gene_motif <- motif %>%
  rename(variant_id = SNP_id) %>%
  inner_join(variant_gene, by = "variant_id") %>%
  group_by(assigned_gene) %>%
  summarise(max_motif_diff = max(abs(alleleDiff), na.rm = TRUE), .groups = "drop")

gene_tbl <- gene_pip %>%
  left_join(gene_da,    by = "assigned_gene") %>%
  left_join(gene_motif, by = "assigned_gene") %>%
  mutate(
    max_da_abs     = replace_na(max_da_abs, 0),
    max_motif_diff = replace_na(max_motif_diff, 0),
    pip_norm   = pmin(max_pip, 1),
    da_norm    = pmin(max_da_abs / 2, 1),
    motif_norm = pmin(max_motif_diff / 1, 1)
  )

cat("Gene-level table:", nrow(gene_tbl), "genes\n")

# ── Merge COLOC PP.H4 (for supervised target) ───────────────────────────────
coloc <- fread(COLOC_F, select = c("gene", "ensembl",
                                   "coloc_best_pp4", "coloc_best_susie_pp4"))
# join on symbol first, then fall back to ensembl
gene_tbl <- gene_tbl %>%
  left_join(coloc %>% select(gene, coloc_best_pp4, coloc_best_susie_pp4),
            by = c("assigned_gene" = "gene"))

# For rows that didn't match (ENSG IDs), try ensembl
need_ensg <- is.na(gene_tbl$coloc_best_pp4) & grepl("^ENSG", gene_tbl$assigned_gene)
if (any(need_ensg)) {
  gene_tbl_ens <- gene_tbl[need_ensg, ] %>%
    select(-coloc_best_pp4, -coloc_best_susie_pp4) %>%
    left_join(coloc %>% select(ensembl, coloc_best_pp4, coloc_best_susie_pp4),
              by = c("assigned_gene" = "ensembl"))
  gene_tbl[need_ensg, c("coloc_best_pp4", "coloc_best_susie_pp4")] <-
    gene_tbl_ens[, c("coloc_best_pp4", "coloc_best_susie_pp4")]
}

# Best-of-two posterior
gene_tbl <- gene_tbl %>%
  mutate(
    coloc_pp4 = pmax(coalesce(coloc_best_pp4, 0),
                     coalesce(coloc_best_susie_pp4, 0)),
    coloc_pos = coloc_pp4 > 0.5
  )

cat("Genes with COLOC PP4>0.5:", sum(gene_tbl$coloc_pos, na.rm = TRUE), "\n")

# ── 6 scoring variants ──────────────────────────────────────────────────────
gene_tbl <- gene_tbl %>%
  mutate(
    score_canonical   = pip_norm * (0.4 + 0.3 * da_norm + 0.3 * motif_norm),
    score_equal       = (pip_norm + da_norm + motif_norm) / 3,
    score_pip_only    = pip_norm,
    score_da_only     = da_norm,
    score_motif_only  = motif_norm
  )

# Elastic-net (data-driven) ─ COLOC pos as outcome
X <- as.matrix(gene_tbl[, c("pip_norm", "da_norm", "motif_norm")])
y <- as.integer(gene_tbl$coloc_pos)

cat("\nElastic-net training: n=", length(y),
    " pos=", sum(y), " (", round(100 * mean(y), 2), "%)\n", sep = "")

# Use grouped 10-fold CV; alpha=0.5
cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 0.5, nfolds = 10,
                    type.measure = "auc", standardize = TRUE)
best_lambda <- cv_fit$lambda.min
coefs <- as.numeric(coef(cv_fit, s = "lambda.min"))
names(coefs) <- c("(Intercept)", "pip_norm", "da_norm", "motif_norm")
cat("CV AUC at lambda.min:", round(max(cv_fit$cvm), 3), "\n")
cat("CV weights (lambda.min):\n"); print(round(coefs, 4))

# Predict probability on full data → use as enet score (in [0,1])
enet_score_full <- as.numeric(predict(cv_fit, newx = X,
                                      s = "lambda.min", type = "response"))
gene_tbl$score_enet <- enet_score_full

# ── 100-fold bootstrap CI on enet weights ──────────────────────────────────
cat("\nRunning 100 bootstrap replicates for enet weight CIs...\n")
B <- 100
boot_coefs <- matrix(NA_real_, nrow = B, ncol = 4)
colnames(boot_coefs) <- c("(Intercept)", "pip_norm", "da_norm", "motif_norm")
boot_auc <- numeric(B)
set.seed(2024)
for (b in seq_len(B)) {
  idx <- sample.int(nrow(X), replace = TRUE)
  tryCatch({
    fit_b <- cv.glmnet(X[idx, ], y[idx], family = "binomial",
                       alpha = 0.5, nfolds = 5, type.measure = "auc",
                       standardize = TRUE)
    boot_coefs[b, ] <- as.numeric(coef(fit_b, s = "lambda.min"))
    boot_auc[b]    <- max(fit_b$cvm, na.rm = TRUE)
  }, error = function(e) NULL)
  if (b %% 20 == 0) cat("  boot", b, "/", B, "\n")
}

enet_weights_out <- data.frame(
  feature   = colnames(boot_coefs),
  cv_weight = coefs,
  boot_mean = apply(boot_coefs, 2, mean, na.rm = TRUE),
  boot_lo   = apply(boot_coefs, 2, quantile, 0.025, na.rm = TRUE),
  boot_hi   = apply(boot_coefs, 2, quantile, 0.975, na.rm = TRUE),
  boot_nonzero_frac = apply(boot_coefs, 2,
                            function(z) mean(abs(z) > 1e-8, na.rm = TRUE))
)
cat("Bootstrap AUC: median=", round(median(boot_auc, na.rm = TRUE), 3),
    " [", round(quantile(boot_auc, 0.025, na.rm = TRUE), 3),
    "–", round(quantile(boot_auc, 0.975, na.rm = TRUE), 3), "]\n", sep = "")

fwrite(enet_weights_out, file.path(OUT_DIR, "elasticnet_weights.csv"))
cat("Wrote elasticnet_weights.csv\n")

# ── Write per-gene sensitivity table ────────────────────────────────────────
out_cols <- c("assigned_gene", "n_variants",
              "pip_norm", "da_norm", "motif_norm",
              "coloc_pp4", "coloc_pos",
              "score_canonical", "score_equal",
              "score_pip_only", "score_da_only", "score_motif_only",
              "score_enet")
fwrite(as.data.table(gene_tbl)[, ..out_cols],
       file.path(OUT_DIR, "composite_score_sensitivity.csv"))
cat("Wrote composite_score_sensitivity.csv (", nrow(gene_tbl), " genes)\n", sep = "")

# ── Top-50 comparison ───────────────────────────────────────────────────────
score_cols <- c("score_canonical", "score_equal",
                "score_pip_only", "score_da_only", "score_motif_only",
                "score_enet")

top_lists <- lapply(score_cols, function(sc) {
  ord <- order(gene_tbl[[sc]], decreasing = TRUE)
  gene_tbl$assigned_gene[ord][seq_len(min(50, nrow(gene_tbl)))]
})
names(top_lists) <- score_cols

# Pairwise Jaccard
nS <- length(score_cols)
jac <- matrix(NA_real_, nS, nS, dimnames = list(score_cols, score_cols))
for (i in seq_len(nS)) for (j in seq_len(nS)) {
  a <- top_lists[[i]]; b <- top_lists[[j]]
  jac[i, j] <- length(intersect(a, b)) / length(union(a, b))
}

# Pairwise Spearman (on score vectors, full gene set)
spr <- matrix(NA_real_, nS, nS, dimnames = list(score_cols, score_cols))
for (i in seq_len(nS)) for (j in seq_len(nS)) {
  spr[i, j] <- suppressWarnings(
    cor(gene_tbl[[score_cols[i]]], gene_tbl[[score_cols[j]]],
        method = "spearman", use = "complete.obs"))
}

# Drug-target capture
drug_targets <- c("THRB", "PPARA", "NR1H4", "PPARG")
drug_capture <- sapply(top_lists, function(g) sum(drug_targets %in% g))
drug_capture_genes <- sapply(top_lists,
                             function(g) paste(intersect(drug_targets, g),
                                               collapse = ";"))

# COLOC PP.H4 > 0.5 capture
coloc_pos_genes <- gene_tbl$assigned_gene[which(gene_tbl$coloc_pos)]
coloc_capture <- sapply(top_lists,
                        function(g) sum(g %in% coloc_pos_genes))

# Long-form overlap table
jac_long <- as.data.frame(as.table(jac))
colnames(jac_long) <- c("scoring_a", "scoring_b", "jaccard_top50")
spr_long <- as.data.frame(as.table(spr))
colnames(spr_long) <- c("scoring_a", "scoring_b", "spearman_full")
overlap_long <- merge(jac_long, spr_long,
                      by = c("scoring_a", "scoring_b"))

per_scoring <- data.frame(
  scoring               = score_cols,
  top50_n               = sapply(top_lists, length),
  drug_targets_captured = drug_capture,
  drug_targets_in_top50 = drug_capture_genes,
  coloc_pp4_05_in_top50 = coloc_capture
)

# Median pairwise Jaccard (off-diagonal)
off_diag <- jac[upper.tri(jac)]
med_jac  <- median(off_diag)
cat("\nMedian pairwise Jaccard (off-diagonal top-50):",
    round(med_jac, 3), "\n")
cat("Drug-target capture (THRB/PPARA/NR1H4/PPARG):\n")
print(per_scoring[, c("scoring", "drug_targets_captured",
                      "drug_targets_in_top50", "coloc_pp4_05_in_top50")])

# Combined CSV: overlap_long + per_scoring summary appended as comments
fwrite(overlap_long,
       file.path(OUT_DIR, "composite_score_top50_overlap.csv"))
fwrite(per_scoring,
       file.path(OUT_DIR, "composite_score_top50_capture.csv"))
cat("\nWrote composite_score_top50_overlap.csv + composite_score_top50_capture.csv\n")

cat("\n================================================================\n")
cat("A7 DONE\n")
cat("================================================================\n")
