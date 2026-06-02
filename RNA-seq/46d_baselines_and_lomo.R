#!/usr/bin/env Rscript
##############################################################################
# 46d_baselines_and_lomo.R
#
# Wave 1 Agents C + D consolidated (run in-process, no SLURM).
#
# Computes:
#   (C) baseline competitors vs. 46d convergence score on Govaere/NIDDK/OT panels
#         1. Nearest gene to GWAS lead SNP (approximate via top PP.H4 gene per
#            GWAS locus from per-GWAS coloc columns — used as a PP4-free proxy)
#         2. Top |dream_logFC| among padj<0.05
#         3. Top coloc_best_pp4 (max SuSiE PP4 across 28 GWAS)
#         4. 46d convergence_score (reference score)
#       → AUROC + bootstrap 95% CI + PR-AUC + Wilcoxon p per panel per score
#
#   (D) LOMO-M AUROC (S1..S8) + expression-decile-matched permutation null on 46d
#       → identifies which modalities drive the headline AUROC; produces the
#         conservative LOMO-S1 (bulk-DE-stripped) number that replaces the
#         full-stack 0.93 in the revised abstract per adversarial review.
#
# Outputs:
#   RNA-seq/results/multi_evidence/convergence_evidence_baseline_benchmark.csv
#   RNA-seq/results/multi_evidence/convergence_evidence_lomo_validation.csv
#   RNA-seq/results/multi_evidence/convergence_evidence_permutation_null_matched.csv
#   RNA-seq/results/multi_evidence/convergence_evidence_validation_summary.md
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
})
# PR-AUC: manual trapezoidal integration so we don't depend on PRROC
pr_auc_manual <- function(scores, labels) {
  labels <- as.logical(labels)
  if (sum(labels) < 5 || sum(!labels) < 5) return(NA_real_)
  ok <- !is.na(scores); scores <- scores[ok]; labels <- labels[ok]
  ord <- order(scores, decreasing = TRUE)
  s <- scores[ord]; l <- labels[ord]
  tp <- cumsum(l); fp <- cumsum(!l)
  precision <- tp / (tp + fp)
  recall    <- tp / sum(l)
  # Prepend (recall=0, precision=1)
  recall <- c(0, recall); precision <- c(1, precision)
  sum(diff(recall) * (precision[-1] + precision[-length(precision)]) / 2)
}

set.seed(42)
BOOT_REPS <- 1000L
PERM_REPS <- 1000L
PANELS_DIR <- "data/published_gene_panels"
ME         <- "RNA-seq/results/multi_evidence"
ATLAS_F    <- file.path(ME, "multi_evidence_atlas.csv")
BAYEV_F    <- file.path(ME, "convergence_evidence.csv")

cat("=== Loading inputs ===\n")
bayev <- fread(BAYEV_F)
cat(sprintf("  convergence_evidence: %d rows\n", nrow(bayev)))

atlas <- fread(ATLAS_F,
  select = c("human_symbol","dream_logFC","dream_padj","coloc_susie_best_pp4",
             "coloc_abf_best_pp4"))
setnames(atlas, "human_symbol", "symbol")
cat(sprintf("  atlas columns loaded: %d rows\n", nrow(atlas)))

# Merge panel-relevant scores onto bayev on human_symbol.
# Note: 2026-04-23 overhaul replaced log_BF_S2a + log_BF_S2b with a single
# log_BF_S2_intact (Okamoto 2023 INTACT consolidation).
dt <- merge(bayev[, .(symbol = human_symbol, convergence_score, excluded_from_ranking,
                       log_BF_S1, log_BF_S2_intact, log_BF_S3,
                       log_BF_S4, log_BF_S5, log_BF_S6, log_BF_S7, log_BF_S8,
                       evidence_unsigned, evidence_signed_magnitude,
                       evidence_signed_sum, concordance_state)],
            atlas, by = "symbol", all.x = TRUE)
dt <- dt[excluded_from_ranking == FALSE]
cat(sprintf("  non-excluded merged frame: %d rows\n", nrow(dt)))

# ---- Panel loading (shares logic with 46d's load_panel) -----------------------
load_panel <- function(path, name) {
  if (!file.exists(path)) { cat("  SKIP", name, ": missing\n"); return(NULL) }
  raw <- readLines(path)
  hdr <- which(startsWith(raw, "gene_symbol"))[1]
  if (is.na(hdr)) { cat("  SKIP", name, ": no header\n"); return(NULL) }
  d <- fread(path, skip = hdr - 1, sep = "\t", header = TRUE)
  g <- unique(d$gene_symbol); g <- g[!is.na(g) & nchar(g) > 0]
  cat(sprintf("  %s panel: %d unique genes (%d in atlas)\n",
              name, length(g), sum(g %in% dt$symbol)))
  g
}

panels <- list(
  Govaere     = load_panel(file.path(PANELS_DIR, "govaere_2020_panel.tsv"), "Govaere"),
  NIDDK       = load_panel(file.path(PANELS_DIR, "niddk_pipeline_2024.tsv"), "NIDDK"),
  OpenTargets = load_panel(file.path(PANELS_DIR, "opentargets_masld_2025.tsv"), "OpenTargets")
)
panels <- panels[!sapply(panels, is.null)]

# ---- Baseline scores ----------------------------------------------------------
# 1. Nearest-gene proxy: use max across individual per-GWAS PP4 columns from
#    the full atlas, which effectively ranks genes by their strongest single-
#    GWAS colocalization — a PP4-based proxy for "which gene did the finemap
#    credible set prefer." We separately record top_pp4 via atlas column.
#
# 2. Top |dream_logFC| with padj<0.05 filter: sig-DEG-rank score
#
# 3. Top coloc_best_pp4: genetic-only baseline
dt[, score_top_logFC := ifelse(!is.na(dream_padj) & dream_padj < 0.05,
                                abs(dream_logFC), 0)]
dt[is.na(score_top_logFC), score_top_logFC := 0]
dt[, score_top_pp4 := pmax(coloc_susie_best_pp4, coloc_abf_best_pp4, na.rm = TRUE)]
dt[is.na(score_top_pp4) | is.infinite(score_top_pp4), score_top_pp4 := 0]
dt[, score_nearest_gene := score_top_pp4]  # same proxy (documented limitation)

# ---- LOMO scores (set one modality log_BF to 0, re-rank) ----------------------
# Approximate LOMO: set log_BF_M = 0, recompute evidence_total as
#   E_unsigned + |E_signed_sum − sign_contribution_M|
# For simplicity, recompute as rowSums of all other pmax(log_BF, 0) for
# unsigned modalities and |signed sum minus contribution|. We approximate by
# recomputing evidence_total = sum(pmax(log_BF_i, 0)) across kept modalities.
# This is a coarse approximation but captures the magnitude effect.
compute_lomo_scores <- function(dt) {
  bf_cols <- c("log_BF_S1","log_BF_S2_intact","log_BF_S3","log_BF_S4",
               "log_BF_S5","log_BF_S6","log_BF_S7","log_BF_S8")
  lomo <- list()
  base <- rowSums(sapply(bf_cols, function(c) pmax(dt[[c]], 0, na.rm = TRUE)))
  lomo[["full"]] <- base
  for (m in c("S1","S2_intact","S3","S4","S5","S6","S7","S8")) {
    cols_keep <- setdiff(bf_cols, paste0("log_BF_", m))
    ev <- rowSums(sapply(cols_keep, function(c) pmax(dt[[c]], 0, na.rm = TRUE)))
    lomo[[paste0("LOMO_", m)]] <- ev
  }
  lomo
}
lomo_scores <- compute_lomo_scores(dt)

# ---- AUROC + bootstrap CI + PR-AUC + Wilcoxon p -------------------------------
wilcox_auroc <- function(scores, labels) {
  scores <- scores[!is.na(scores)]; labels <- labels[!is.na(scores) & seq_along(labels) %in% seq_along(scores)]
  labels <- as.logical(labels)
  if (sum(labels) < 5 || sum(!labels) < 5) return(list(auroc=NA, p=NA))
  r <- rank(scores)
  nP <- sum(labels); nN <- sum(!labels)
  auroc <- (sum(r[labels]) - nP*(nP+1)/2) / (nP*nN)
  p <- tryCatch(
    wilcox.test(scores[labels], scores[!labels], alternative="greater")$p.value,
    error = function(e) NA_real_
  )
  list(auroc = auroc, p = p)
}

pr_auc <- pr_auc_manual

bootstrap_auroc <- function(scores, labels, reps = BOOT_REPS) {
  labels <- as.logical(labels)
  if (sum(labels) < 5 || sum(!labels) < 5) return(c(NA, NA))
  n <- length(scores)
  aurocs <- numeric(reps)
  for (i in seq_len(reps)) {
    idx <- sample.int(n, n, replace = TRUE)
    s <- scores[idx]; l <- labels[idx]
    if (sum(l) < 2 || sum(!l) < 2) { aurocs[i] <- NA; next }
    r <- rank(s); nP <- sum(l); nN <- sum(!l)
    aurocs[i] <- (sum(r[l]) - nP*(nP+1)/2) / (nP*nN)
  }
  quantile(aurocs, c(0.025, 0.975), na.rm = TRUE)
}

# ---- Benchmark table per panel per score -------------------------------------
cat("\n=== (C) Baseline benchmark ===\n")
baseline_rows <- list()
for (pname in names(panels)) {
  panel_genes <- panels[[pname]]
  labels <- dt$symbol %in% panel_genes
  n_panel <- sum(labels)
  cat(sprintf("  %s: n_in_atlas=%d\n", pname, n_panel))
  for (sc_name in c("convergence_score", "score_top_logFC", "score_top_pp4", "score_nearest_gene")) {
    scores <- dt[[sc_name]]
    valid <- !is.na(scores)
    res <- wilcox_auroc(scores[valid], labels[valid])
    ci <- bootstrap_auroc(scores[valid], labels[valid])
    pr <- pr_auc(scores[valid], labels[valid])
    baseline_rows[[paste(pname, sc_name, sep="|")]] <- data.table(
      panel = pname, score = sc_name, n_panel = n_panel,
      auroc = round(res$auroc, 4),
      auroc_ci_low = round(ci[1], 4), auroc_ci_high = round(ci[2], 4),
      pr_auc = round(pr, 4),
      wilcoxon_p = formatC(res$p, format="e", digits=2)
    )
  }
}
baseline_dt <- rbindlist(baseline_rows)
fwrite(baseline_dt, file.path(ME, "convergence_evidence_baseline_benchmark.csv"))
print(baseline_dt)

# ---- LOMO per panel -----------------------------------------------------------
cat("\n=== (D) LOMO validation ===\n")
lomo_rows <- list()
for (pname in names(panels)) {
  panel_genes <- panels[[pname]]
  labels <- dt$symbol %in% panel_genes
  for (lm in names(lomo_scores)) {
    scores <- lomo_scores[[lm]]
    valid <- !is.na(scores)
    res <- wilcox_auroc(scores[valid], labels[valid])
    ci <- bootstrap_auroc(scores[valid], labels[valid])
    pr <- pr_auc(scores[valid], labels[valid])
    lomo_rows[[paste(pname, lm, sep="|")]] <- data.table(
      panel = pname, lomo_variant = lm,
      auroc = round(res$auroc, 4),
      auroc_ci_low = round(ci[1], 4), auroc_ci_high = round(ci[2], 4),
      pr_auc = round(pr, 4),
      wilcoxon_p = formatC(res$p, format="e", digits=2)
    )
  }
}
lomo_dt <- rbindlist(lomo_rows)
fwrite(lomo_dt, file.path(ME, "convergence_evidence_lomo_validation.csv"))
print(lomo_dt)

# ---- Expression-decile-matched permutation null -------------------------------
cat("\n=== (D) Expression-matched permutation null ===\n")
# Decile-bin genes by |dream_logFC|
dt[, abs_lfc := abs(dream_logFC)]
dt[is.na(abs_lfc), abs_lfc := 0]
dt[, lfc_decile := cut(abs_lfc, quantile(abs_lfc, probs = seq(0, 1, 0.1), na.rm=TRUE),
                       include.lowest = TRUE, labels = FALSE)]
dt[is.na(lfc_decile), lfc_decile := 1L]

perm_rows <- list()
for (pname in names(panels)) {
  panel_genes <- panels[[pname]]
  panel_idx <- which(dt$symbol %in% panel_genes)
  n_panel <- length(panel_idx)
  if (n_panel < 5) next
  # Observed
  obs_res <- wilcox_auroc(dt$convergence_score, dt$symbol %in% panel_genes)
  # Null: sample n_panel genes matched on LFC decile, 1000 times
  null_aurocs <- numeric(PERM_REPS)
  decile_tab <- dt[panel_idx, .N, by = lfc_decile]
  for (i in seq_len(PERM_REPS)) {
    sampled_idx <- unlist(lapply(seq_len(nrow(decile_tab)), function(r) {
      d <- decile_tab$lfc_decile[r]; n <- decile_tab$N[r]
      pool <- which(dt$lfc_decile == d)
      sample(pool, n, replace = FALSE)
    }))
    null_lbl <- seq_len(nrow(dt)) %in% sampled_idx
    res <- wilcox_auroc(dt$convergence_score, null_lbl)
    null_aurocs[i] <- res$auroc
  }
  perm_rows[[pname]] <- data.table(
    panel = pname,
    n_panel = n_panel,
    observed_auroc = round(obs_res$auroc, 4),
    null_mean = round(mean(null_aurocs, na.rm=TRUE), 4),
    null_sd = round(sd(null_aurocs, na.rm=TRUE), 4),
    null_95pct = round(quantile(null_aurocs, 0.95, na.rm=TRUE), 4),
    observed_percentile = round(100 * mean(null_aurocs <= obs_res$auroc, na.rm=TRUE), 2),
    perm_p = round(mean(null_aurocs >= obs_res$auroc, na.rm=TRUE), 4)
  )
  cat(sprintf("  %s: obs=%.3f null_mean=%.3f null_95pct=%.3f perm_p=%.4f\n",
              pname, obs_res$auroc, mean(null_aurocs, na.rm=TRUE),
              quantile(null_aurocs, 0.95, na.rm=TRUE),
              mean(null_aurocs >= obs_res$auroc, na.rm=TRUE)))
}
perm_dt <- rbindlist(perm_rows)
fwrite(perm_dt, file.path(ME, "convergence_evidence_permutation_null_matched.csv"))
print(perm_dt)

# ---- Summary markdown ---------------------------------------------------------
cat("\n=== Writing summary markdown ===\n")
summary_lines <- c(
  "# 46d validation summary — baselines + LOMO + expression-matched null",
  "",
  "**Generated by `RNA-seq/46d_baselines_and_lomo.R` on the existing `convergence_evidence.csv` outputs (pre-script-overhaul Wave 1).**",
  "",
  "## Headline numbers for the revised abstract",
  "",
  sprintf("- **46d full-stack Govaere AUROC:** %s (current abstract number, adversarial review consensus says this is in-domain leakage)",
          baseline_dt[panel=="Govaere" & score=="convergence_score", auroc]),
  sprintf("- **46d LOMO-S1 (bulk-DE-stripped) Govaere AUROC:** %s ← **this is the revised headline number**",
          lomo_dt[panel=="Govaere" & lomo_variant=="LOMO_S1", auroc]),
  sprintf("- **Top |dream_logFC| baseline Govaere AUROC:** %s (Patel 2025 concern: how much does 46d add over a trivial |logFC| rank?)",
          baseline_dt[panel=="Govaere" & score=="score_top_logFC", auroc]),
  sprintf("- **Top COLOC PP4 baseline Govaere AUROC:** %s",
          baseline_dt[panel=="Govaere" & score=="score_top_pp4", auroc]),
  sprintf("- **Expression-matched permutation null Govaere p:** %s (1,000 LFC-decile-matched random 25-gene panels)",
          if (nrow(perm_dt[panel=="Govaere"]) > 0) perm_dt[panel=="Govaere", perm_p] else NA),
  "",
  "## Full baseline benchmark (per panel × score)",
  "",
  paste0("See `convergence_evidence_baseline_benchmark.csv`"),
  "",
  "## LOMO per-modality AUROC",
  "",
  paste0("See `convergence_evidence_lomo_validation.csv`"),
  "",
  "## Expression-matched permutation null (1,000 reps)",
  "",
  paste0("See `convergence_evidence_permutation_null_matched.csv`"),
  "",
  "## Implications for revised paper claims",
  "",
  "(Draft — fill in with actual numbers after this script runs.)",
  "",
  "- The revised abstract should lead with LOMO-S1 AUROC, not the full-stack 0.93.",
  "- Baseline comparison shows whether 46d genuinely adds value over trivial |logFC| or top-PP4 ranks.",
  "- Permutation-null percentile establishes whether the score is base-rate-inflated at the panel's small n.",
  "",
  "## Known limitations",
  "",
  "- Nearest-gene baseline is approximated via `score_top_pp4` (max COLOC PP.H4 across 28 GWAS) — true genomic-distance nearest-gene would require GENCODE + SuSiE lead-SNP positions; punted for Wave 1.",
  "- LOMO recomputes evidence_total as sum of kept pmax(log_BF, 0); does NOT re-run the full 4-state concordance classifier. This is an approximation; fine for AUROC perturbation magnitude but not exact.",
  "- Sample-disjoint Govaere rerun (dropping GSE135251) deferred to Wave 2 — would require full pipeline re-run."
)
writeLines(summary_lines, file.path(ME, "convergence_evidence_validation_summary.md"))

cat("\n=== Done. Outputs:\n")
cat("  ", file.path(ME, "convergence_evidence_baseline_benchmark.csv"), "\n")
cat("  ", file.path(ME, "convergence_evidence_lomo_validation.csv"), "\n")
cat("  ", file.path(ME, "convergence_evidence_permutation_null_matched.csv"), "\n")
cat("  ", file.path(ME, "convergence_evidence_validation_summary.md"), "\n")
