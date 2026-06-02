# 261_held_out_replication.R
# Team M7 — Multi-step cascade hierarchy replication on a held-out cohort.
#
# Original numbers (from 250_two_transition_decomposition.R, all F0-F4 cohorts pooled):
#   F0->F1: 0.373   F1->F2: 0.245   F2->F3: 0.400   F3->F4: 0.466
#   (median |log2FC| of top 1000 DEGs by P-value; limma-voom + cohort blocking)
#
# Cohort situation: The original cascade used EVERY sample in meta_matched.rds
# with a non-NA fibrosis_stage. ALL F0-F4-staged cohorts (GSE130970, GSE213621,
# GSE162694, GSE174478, GSE193066, GSE240729, GSE135251) contributed to the
# original cascade. There is NO truly blinded held-out cohort.
#
# HONEST FALLBACK: Leave-one-cohort-out (LOCO) cross-validation.
#   For each eligible cohort C:
#     1. Compute cascade hierarchy on C alone (within-cohort transitions).
#     2. Spearman rho of C's cascade vs the literal original anchor.
#     3. COLOC OR per transition (top-500 DEGs vs PP4>=0.5 background) on C.
#   Aggregate across LOCO folds:
#     - Mean / median rho
#     - 2,000-bootstrap percentile 95% CI on mean rho
#     - F2->F3 OR >= F1->F2 OR consistency rate
#
# Pre-registered acceptance:
#   - Mean LOCO Spearman rho >= 0.6
#   - 95% CI excludes 0
#   - F2->F3 OR >= F1->F2 OR (rank consistent) on majority of folds
#
# Reference COLOC ORs (computed once on the full F0-F4 pool, the "original"):
#   - or_original_F0F1, or_original_F1F2, or_original_F2F3, or_original_F3F4
#
# Inputs:
#   results/integration/merged_dge.rds, meta_matched.rds
#   results/multi_evidence/multi_evidence_atlas.csv
#
# Outputs (results/granular_staging/):
#   held_out_replication.csv
#   held_out_replication.md

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(matrixStats)
})

set.seed(42)

# Force unbuffered stdout for live progress in SLURM logs
options(echo = FALSE)

cat0 <- function(...) { cat(..., sep = ""); flush.console() }
say  <- function(...) { cat(..., "\n", sep = ""); flush.console() }

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                           "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_DIR <- file.path(PROJECT_ROOT,
                     "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
ATLAS_PATH <- file.path(PROJECT_ROOT, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ORIGINAL_CASCADE <- c(F0F1 = 0.373, F1F2 = 0.245, F2F3 = 0.400, F3F4 = 0.466)
TRANSITIONS <- names(ORIGINAL_CASCADE)
N_BOOTSTRAP <- 2000

say("== 261 Held-out replication of multi-step cascade hierarchy ==")
say("Original cascade anchor:")
print(ORIGINAL_CASCADE); flush.console()

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
say("\nLoading DGE + meta...")
dge <- readRDS(file.path(INT_DIR, "merged_dge.rds"))
meta <- readRDS(file.path(INT_DIR, "meta_matched.rds"))
meta <- meta[match(colnames(dge), meta$sample_id), ]
say(sprintf("DGE: %d genes x %d samples", nrow(dge), ncol(dge)))

staged_idx <- which(!is.na(meta$fibrosis_stage))
meta_st <- meta[staged_idx, ]
dge_st <- dge[, staged_idx]
say(sprintf("Stage-annotated samples: %d", length(staged_idx)))
say("Per-cohort x stage table:")
print(table(meta_st$dataset, meta_st$fibrosis_stage)); flush.console()

# ---------------------------------------------------------------------------
# Eligibility
# ---------------------------------------------------------------------------
cohort_stage_counts <- table(meta_st$dataset, meta_st$fibrosis_stage)
cohorts_all <- rownames(cohort_stage_counts)
cohort_eligibility <- sapply(cohorts_all, function(c) {
  cnt <- cohort_stage_counts[c, ]
  has_pair <- FALSE
  for (i in 1:(length(cnt)-1)) {
    if (cnt[i] >= 5 && cnt[i+1] >= 5) { has_pair <- TRUE; break }
  }
  has_pair && sum(cnt) >= 10
})
candidate_cohorts <- cohorts_all[cohort_eligibility]
say(sprintf("\nLOCO candidate cohorts: %s", paste(candidate_cohorts, collapse = ", ")))

# ---------------------------------------------------------------------------
# Helper: limma-voom contrast on a subset
# ---------------------------------------------------------------------------
run_contrast <- function(idx_a, idx_b, meta_subset, dge_subset, require_n_each = 5) {
  if (length(idx_a) < require_n_each || length(idx_b) < require_n_each) return(NULL)
  use_idx <- c(idx_a, idx_b)
  group <- factor(c(rep("A", length(idx_a)), rep("B", length(idx_b))),
                  levels = c("A", "B"))
  cohort <- factor(meta_subset$dataset[use_idx])
  d <- dge_subset[, use_idx]
  keep <- filterByExpr(d, group = group, min.count = 5)
  d <- d[keep, , keep.lib.sizes = FALSE]
  if (nrow(d) < 200) return(NULL)
  d <- calcNormFactors(d)
  design <- if (length(unique(cohort)) > 1) model.matrix(~ group + cohort) else model.matrix(~ group)
  v <- voom(d, design)
  fit <- eBayes(lmFit(v, design))
  res <- topTable(fit, coef = "groupB", number = Inf, sort.by = "none")
  res$gene_id <- rownames(res)
  res
}

cascade_for_subset <- function(meta_sub, dge_sub) {
  out <- numeric(length(TRANSITIONS)); names(out) <- TRANSITIONS
  for (t in TRANSITIONS) {
    a <- as.integer(substr(t, 2, 2)); b <- as.integer(substr(t, 4, 4))
    idx_a <- which(meta_sub$fibrosis_stage == a)
    idx_b <- which(meta_sub$fibrosis_stage == b)
    res <- run_contrast(idx_a, idx_b, meta_sub, dge_sub)
    if (is.null(res)) {
      out[t] <- NA_real_
    } else {
      top <- head(res[order(res$P.Value), , drop = FALSE], 1000)
      out[t] <- median(abs(top$logFC), na.rm = TRUE)
    }
  }
  out
}

# ---------------------------------------------------------------------------
# COLOC ground truth
# ---------------------------------------------------------------------------
say("\nLoading multi-evidence atlas...")
atlas <- read.csv(ATLAS_PATH, stringsAsFactors = FALSE)
pp4_col <- if ("coloc_best_susie_pp4_polyfun" %in% colnames(atlas)) {
  "coloc_best_susie_pp4_polyfun"
} else if ("coloc_best_pp4_polyfun" %in% colnames(atlas)) {
  "coloc_best_pp4_polyfun"
} else "coloc_abf_best_pp4"
say(sprintf("PP4 column: %s", pp4_col))
atlas$ensembl_clean <- sub("\\..*", "", atlas$ensembl_id)
coloc_genes <- atlas$ensembl_clean[!is.na(atlas[[pp4_col]]) & atlas[[pp4_col]] >= 0.5]
say(sprintf("Atlas: %d genes; %d with PP4 >= 0.5", nrow(atlas), length(coloc_genes)))

odds_ratio_for_transition <- function(deg_res, top_n = 500) {
  if (is.null(deg_res)) return(c(or = NA_real_, p = NA_real_, n_top = NA_integer_,
                                  n_top_coloc = NA_integer_))
  top <- head(deg_res[order(deg_res$P.Value), , drop = FALSE], top_n)
  top_clean <- sub("\\..*", "", top$gene_id)
  background <- sub("\\..*", "", deg_res$gene_id)
  in_top <- background %in% top_clean
  in_coloc <- background %in% coloc_genes
  ct <- table(in_top, in_coloc)
  if (any(dim(ct) < 2)) return(c(or = NA_real_, p = NA_real_, n_top = length(top_clean),
                                  n_top_coloc = sum(top_clean %in% coloc_genes)))
  ft <- fisher.test(ct)
  c(or = unname(ft$estimate), p = ft$p.value,
    n_top = length(top_clean),
    n_top_coloc = sum(top_clean %in% coloc_genes))
}

# ---------------------------------------------------------------------------
# Reference (original): cascade + COLOC ORs on the full F0-F4 pool
# ---------------------------------------------------------------------------
say("\n--- Reference: full-pool (anchor) cascade + COLOC ORs ---")
ref_or <- list()
for (t in TRANSITIONS) {
  a <- as.integer(substr(t, 2, 2)); b <- as.integer(substr(t, 4, 4))
  idx_a <- which(meta_st$fibrosis_stage == a)
  idx_b <- which(meta_st$fibrosis_stage == b)
  say(sprintf("  Reference contrast %s: a=%d b=%d", t, length(idx_a), length(idx_b)))
  res <- run_contrast(idx_a, idx_b, meta_st, dge_st)
  ref_or[[t]] <- odds_ratio_for_transition(res)
}
say("Reference COLOC ORs:")
for (t in TRANSITIONS) {
  v <- ref_or[[t]]
  say(sprintf("  %s OR=%.2f  p=%.2e  n_top_coloc=%d", t,
              v["or"], v["p"], v["n_top_coloc"]))
}

# ---------------------------------------------------------------------------
# LOCO per-cohort
# ---------------------------------------------------------------------------
results_list <- list()
for (cohort in candidate_cohorts) {
  say(sprintf("\n--- LOCO held-out: %s ---", cohort))
  ho_keep <- which(meta_st$dataset == cohort)
  meta_ho <- meta_st[ho_keep, ]
  dge_ho <- dge_st[, ho_keep]
  say(sprintf("Held-out n=%d  stages present: %s", length(ho_keep),
              paste(sort(unique(meta_ho$fibrosis_stage)), collapse = ",")))

  # Held-out cascade
  ho_cascade <- tryCatch(cascade_for_subset(meta_ho, dge_ho),
                         error = function(e) { say("HO cascade error: ", conditionMessage(e)); rep(NA_real_, 4) })
  names(ho_cascade) <- TRANSITIONS
  say(sprintf("Held-out cascade: %s",
              paste(sprintf("%s=%s", names(ho_cascade),
                            ifelse(is.na(ho_cascade), "NA", sprintf("%.3f", ho_cascade))),
                    collapse = "  ")))

  # Spearman rho vs anchor
  valid <- !is.na(ho_cascade) & !is.na(ORIGINAL_CASCADE)
  rho_orig <- if (sum(valid) >= 2) {
    suppressWarnings(cor(ho_cascade[valid], ORIGINAL_CASCADE[valid], method = "spearman"))
  } else NA_real_
  say(sprintf("Spearman rho vs original anchor: %s",
              ifelse(is.na(rho_orig), "NA", sprintf("%.3f", rho_orig))))

  # Held-out COLOC ORs per transition
  trans_or <- list()
  for (t in TRANSITIONS) {
    a <- as.integer(substr(t, 2, 2)); b <- as.integer(substr(t, 4, 4))
    idx_a <- which(meta_ho$fibrosis_stage == a)
    idx_b <- which(meta_ho$fibrosis_stage == b)
    res <- run_contrast(idx_a, idx_b, meta_ho, dge_ho)
    trans_or[[t]] <- odds_ratio_for_transition(res)
  }
  or_f1f2 <- trans_or[["F1F2"]]["or"]
  or_f2f3 <- trans_or[["F2F3"]]["or"]
  rank_consistent <- !is.na(or_f1f2) && !is.na(or_f2f3) && or_f2f3 >= or_f1f2
  pass_rep <- !is.na(rho_orig) && rho_orig >= 0.6 && rank_consistent
  say(sprintf("F1F2 OR=%s  F2F3 OR=%s  rank_consistent=%s",
              ifelse(is.na(or_f1f2), "NA", sprintf("%.2f", or_f1f2)),
              ifelse(is.na(or_f2f3), "NA", sprintf("%.2f", or_f2f3)),
              rank_consistent))

  for (t in TRANSITIONS) {
    or_ho <- trans_or[[t]]
    or_ref <- ref_or[[t]]
    results_list[[paste(cohort, t)]] <- data.frame(
      held_out_cohort = cohort,
      transition = t,
      effect_size_held_out = unname(ho_cascade[t]),
      effect_size_original = unname(ORIGINAL_CASCADE[t]),
      or_held_out = unname(or_ho["or"]),
      or_held_out_p = unname(or_ho["p"]),
      or_original = unname(or_ref["or"]),
      or_original_p = unname(or_ref["p"]),
      n_top_coloc_held_out = unname(or_ho["n_top_coloc"]),
      n_top_coloc_original = unname(or_ref["n_top_coloc"]),
      rank_concordance_spearman = rho_orig,
      f2f3_ge_f1f2_held_out = rank_consistent,
      pass_replication = pass_rep,
      notes = "LOCO fallback; original cascade trained on union of all F0-F4 cohorts (no truly blinded cohort exists)",
      stringsAsFactors = FALSE
    )
  }
}

# ---------------------------------------------------------------------------
# Aggregate + bootstrap CI
# ---------------------------------------------------------------------------
out_df <- do.call(rbind, results_list)
write.csv(out_df, file.path(OUT_DIR, "held_out_replication.csv"), row.names = FALSE)
say(sprintf("\nWrote %s (%d rows)",
            file.path(OUT_DIR, "held_out_replication.csv"), nrow(out_df)))

per_cohort <- unique(out_df[, c("held_out_cohort", "rank_concordance_spearman", "f2f3_ge_f1f2_held_out")])
rho_vec <- per_cohort$rank_concordance_spearman
rho_vec <- rho_vec[!is.na(rho_vec)]
mean_rho <- if (length(rho_vec) > 0) mean(rho_vec) else NA_real_
median_rho <- if (length(rho_vec) > 0) median(rho_vec) else NA_real_

say(sprintf("\nPer-fold rho values:"))
for (i in seq_len(nrow(per_cohort))) {
  say(sprintf("  %-15s rho=%s  F2F3>=F1F2=%s",
              per_cohort$held_out_cohort[i],
              ifelse(is.na(per_cohort$rank_concordance_spearman[i]), "NA",
                     sprintf("%.3f", per_cohort$rank_concordance_spearman[i])),
              ifelse(is.na(per_cohort$f2f3_ge_f1f2_held_out[i]), "NA",
                     ifelse(per_cohort$f2f3_ge_f1f2_held_out[i], "YES", "no"))))
}
say(sprintf("Mean LOCO rho: %.3f  (n=%d folds)", mean_rho, length(rho_vec)))

if (length(rho_vec) >= 2) {
  boot_means <- replicate(N_BOOTSTRAP, mean(sample(rho_vec, length(rho_vec), replace = TRUE)))
  ci <- quantile(boot_means, c(0.025, 0.975))
} else {
  ci <- c(`2.5%` = NA_real_, `97.5%` = NA_real_)
}
ci_excludes_zero <- !is.na(ci[1]) && !is.na(ci[2]) && (ci[1] > 0 || ci[2] < 0)
say(sprintf("Bootstrap 95%% CI on mean rho: [%.3f, %.3f] (excludes 0: %s)",
            ci[1], ci[2], ci_excludes_zero))

n_consistent <- sum(per_cohort$f2f3_ge_f1f2_held_out, na.rm = TRUE)
n_total <- sum(!is.na(per_cohort$f2f3_ge_f1f2_held_out))
say(sprintf("F2->F3 OR >= F1->F2 OR consistent: %d/%d", n_consistent, n_total))

verdict_pass <- !is.na(mean_rho) && mean_rho >= 0.6 &&
                ci_excludes_zero &&
                n_consistent >= ceiling(n_total / 2)

# ---------------------------------------------------------------------------
# Markdown verdict
# ---------------------------------------------------------------------------
sink(file.path(OUT_DIR, "held_out_replication.md"))
cat("# Held-out replication of the multi-step cascade hierarchy (Team M7)\n\n")
cat("## Cohort situation\n\n")
cat("**No truly blinded held-out cohort exists.** The original cascade in `250_two_transition_decomposition.R` ")
cat("used every sample in `meta_matched.rds` with a non-NA `fibrosis_stage`. Inspection of the F0-F4 stage ")
cat("distribution per cohort confirms that all stage-bearing cohorts contributed to the original cascade.\n\n")
cat("**Honest fallback: leave-one-cohort-out (LOCO).** For each eligible cohort the cascade is recomputed ")
cat("from within-cohort transitions only and rank-correlated against the original anchor ")
cat("(0.373 / 0.245 / 0.400 / 0.466). COLOC ORs per transition use Fisher's exact test on top-500 DEGs vs ")
cat("PP4 >= 0.5 background; reference ORs are computed once on the full pooled F0-F4 cohort.\n\n")

cat("## Anchor (original 250 cascade)\n\n")
cat("| Transition | median |log2FC| top1000 |\n|---|---:|\n")
for (t in TRANSITIONS) cat(sprintf("| %s | %.3f |\n", t, ORIGINAL_CASCADE[t]))

cat("\n## LOCO per-cohort cascade and concordance\n\n")
cat("| Held-out | F0F1 | F1F2 | F2F3 | F3F4 | rho vs original | F2F3>=F1F2 OR? |\n")
cat("|---|---:|---:|---:|---:|---:|:---:|\n")
for (cohort in unique(out_df$held_out_cohort)) {
  rows <- out_df[out_df$held_out_cohort == cohort, ]
  vals <- setNames(rows$effect_size_held_out, rows$transition)
  fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
  cat(sprintf("| %s | %s | %s | %s | %s | %s | %s |\n",
              cohort, fmt(vals["F0F1"]), fmt(vals["F1F2"]),
              fmt(vals["F2F3"]), fmt(vals["F3F4"]),
              ifelse(is.na(rows$rank_concordance_spearman[1]), "NA",
                     sprintf("%.3f", rows$rank_concordance_spearman[1])),
              ifelse(is.na(rows$f2f3_ge_f1f2_held_out[1]), "NA",
                     ifelse(rows$f2f3_ge_f1f2_held_out[1], "YES", "no"))))
}

cat("\n## COLOC ORs (held-out vs full-pool reference)\n\n")
cat("| Held-out | F0F1 OR (HO/ref) | F1F2 OR (HO/ref) | F2F3 OR (HO/ref) | F3F4 OR (HO/ref) |\n")
cat("|---|---:|---:|---:|---:|\n")
for (cohort in unique(out_df$held_out_cohort)) {
  rows <- out_df[out_df$held_out_cohort == cohort, ]
  ho <- setNames(rows$or_held_out, rows$transition)
  ref <- setNames(rows$or_original, rows$transition)
  fmt <- function(x, y) sprintf("%s/%s",
                               ifelse(is.na(x), "NA", sprintf("%.2f", x)),
                               ifelse(is.na(y), "NA", sprintf("%.2f", y)))
  cat(sprintf("| %s | %s | %s | %s | %s |\n",
              cohort, fmt(ho["F0F1"], ref["F0F1"]),
              fmt(ho["F1F2"], ref["F1F2"]),
              fmt(ho["F2F3"], ref["F2F3"]),
              fmt(ho["F3F4"], ref["F3F4"])))
}

cat("\n## Aggregate concordance\n\n")
cat(sprintf("- LOCO folds with usable rho: **%d**\n", length(rho_vec)))
cat(sprintf("- Mean rho vs original anchor: **%.3f**\n", mean_rho))
cat(sprintf("- Median rho vs original anchor: **%.3f**\n", median_rho))
cat(sprintf("- Bootstrap (n=%d) 95%% CI on mean rho: **[%.3f, %.3f]**\n", N_BOOTSTRAP, ci[1], ci[2]))
cat(sprintf("- CI excludes 0: **%s**\n", ci_excludes_zero))
cat(sprintf("- F2->F3 OR >= F1->F2 OR consistent: **%d / %d** cohorts\n", n_consistent, n_total))

cat("\n## Pre-registered acceptance\n\n")
cat("| Criterion | Threshold | Observed | Pass? |\n|---|---|---|:---:|\n")
cat(sprintf("| Mean LOCO Spearman rho | >= 0.6 | %s | %s |\n",
            ifelse(is.na(mean_rho), "NA", sprintf("%.3f", mean_rho)),
            ifelse(!is.na(mean_rho) && mean_rho >= 0.6, "YES", "no")))
cat(sprintf("| 95%% CI excludes 0 | yes | [%.3f, %.3f] | %s |\n",
            ci[1], ci[2], ifelse(ci_excludes_zero, "YES", "no")))
cat(sprintf("| F2F3 OR >= F1F2 OR | majority | %d/%d | %s |\n",
            n_consistent, n_total,
            ifelse(n_consistent >= ceiling(n_total / 2), "YES", "no")))

cat("\n## Verdict\n\n")
if (verdict_pass) {
  cat("**Replicates** under LOCO fallback. The cascade hierarchy F0->F1 / F1->F2 / F2->F3 / F3->F4 is ")
  cat("recoverable from independent within-cohort transitions, and the F2->F3 vs F1->F2 OR ordering is ")
  cat("rank-consistent across the LOCO folds.\n\n")
} else {
  cat("**Does NOT cleanly replicate** under LOCO fallback. The original cascade may reflect partly ")
  cat("between-cohort composition rather than a within-cohort biological gradient. Disclose the LOCO ")
  cat("result as a sensitivity check in the methods text.\n\n")
}

cat("## Methods notes\n\n")
cat("- Effect size: median |log2FC| of top 1000 DEGs ranked by P-value (limma-voom + cohort blocking when applicable). Identical contrast logic to `250_two_transition_decomposition.R`.\n")
cat(sprintf("- COLOC ground truth: `%s` >= 0.5 from `multi_evidence_atlas.csv`.\n", pp4_col))
cat("- COLOC OR: Fisher's exact test, top-500 DEGs vs background of all expressed genes per contrast, COLOC PP4 >= 0.5 as positive label.\n")
cat("- Reference COLOC ORs: computed once on the full F0-F4 pool (the closest analog to the original anchor's COLOC enrichment).\n")
cat("- Bootstrap: 2,000 percentile resamples of LOCO fold rho values with replacement.\n")
cat("- Eligibility: cohort must have >=2 adjacent stages with >=5 samples each AND >=10 staged samples total.\n")
sink()

say(sprintf("\nWrote %s", file.path(OUT_DIR, "held_out_replication.md")))
say("\n== 261 complete. ==")
