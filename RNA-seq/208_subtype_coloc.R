#!/usr/bin/env Rscript
# 208_subtype_coloc.R — k-Program Genetic Architecture
#
# REFACTORED 2026-04-20: replaces binary S1/S2 COLOC analysis with per-program
# enrichment + multinomial logistic LOCO-CV across k biological programs from
# Script 44 (Metabolic / Inflammatory / Fibrotic / Sex-* / HCC-like / ...).
#
# Analysis:
#   1. Per-program Fisher enrichment: COLOC PP.H4>0.5/0.8/0.9 genes overlap
#      each program's upregulated marker set.
#   2. Per-program causal score: PP4 * |program_logFC| (signed).
#   3. LOCO-CV MULTINOMIAL logistic (nnet::multinom) predicting
#      dominant_program from COLOC+TWAS features → one-vs-rest AUROC per
#      program.
#   4. Sex-program decomposition: per-program Fisher(program-up ∩ female-biased).
#
# Outputs (to RNA-seq/results/stratified_causal/):
#   program_coloc_enrichment.csv  (long: program x gene_set x coloc_threshold)
#   program_causal_scores.csv     (per-gene x per-program)
#   program_multinom_auroc.csv    (one-vs-rest AUROC per program)
#   program_multinom_predictions.csv
#   program_sex_decomposition.csv
#   subtype_coloc_enrichment.csv  (DEPRECATED ALIAS — points to fibrotic program rows)
#   subtype_causal_scores.csv     (DEPRECATED ALIAS)
#   s2_genetic_risk.csv           (DEPRECATED ALIAS — Fibrotic program AUROC)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(nnet)
  library(pROC)
  library(fgsea)
  library(edgeR)
  library(limma)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("=== 208_subtype_coloc.R: k-Program Genetic Architecture ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ============================================================
# 1. Load inputs
# ============================================================
cat("=== Loading inputs ===\n")
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("ensembl_id","human_symbol"))
atlas <- atlas[human_symbol != "" & !is.na(human_symbol)]
atlas[, ens_base := sub("\\..*", "", ensembl_id)]

coloc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]
# T0.4 swap (review B3#5): prefer SuSiE PP4 over legacy ABF — same canonical pattern as
# 207/210/213/217b. Without this, subtype COLOC enrichment was ABF-based (317 genes change
# their PP4>0.5 classification) and non-comparable to the sibling stratified scripts.
if ("coloc_best_susie_pp4" %in% names(coloc)) {
  coloc[, coloc_best_abf_pp4 := coloc_best_pp4]
  coloc[, coloc_best_pp4 := fifelse(!is.na(coloc_best_susie_pp4),
                                    coloc_best_susie_pp4, coloc_best_pp4)]
  cat("  T0.4: coloc_best_pp4 now sourced from SuSiE with ABF fallback\n")
}
coloc_05 <- coloc[coloc_best_pp4 > 0.5, unique(gene)]
coloc_08 <- coloc[coloc_best_pp4 > 0.8, unique(gene)]
coloc_09 <- coloc[coloc_best_pp4 > 0.9, unique(gene)]
cat(sprintf("  COLOC genes: total=%d, PP4>0.5=%d, PP4>0.8=%d, PP4>0.9=%d\n",
            nrow(coloc), length(coloc_05), length(coloc_08), length(coloc_09)))

nmf <- fread(file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv"))
stopifnot(all(c("sample_id","dominant_program_code","dominant_program") %in% names(nmf)))
cat(sprintf("  NMF assignments: %d samples\n", nrow(nmf)))
cat("  Programs:\n"); print(nmf[, .N, by = .(dominant_program_code, dominant_program)])

prog_labels <- fread(file.path(BASE, "RNA-seq/results/subtypes/program_labels.csv"))
K <- nrow(prog_labels)
program_codes <- prog_labels$program_code
program_bios  <- prog_labels$biological_label
cat(sprintf("  k = %d programs\n", K))

# Identify Fibrogenic program (for deprecated-alias output)
fib_code <- prog_labels$program_code[prog_labels$label_category == "Fibrogenic"]
stopifnot(length(fib_code) == 1)
if (length(fib_code) == 0) {
  cat("  NOTE: no Fibrogenic-labeled program — legacy S2 aliases will be empty.\n")
  fib_code <- ""
} else if (length(fib_code) > 1) {
  fib_code <- fib_code[1]
  cat(sprintf("  WARNING: multiple Fibrogenic programs, using %s for legacy aliases\n", fib_code))
}

markers <- fread(file.path(BASE, "RNA-seq/results/subtypes/subtype_markers.csv"))
stopifnot(all(c("gene","program_code","logFC","padj","direction") %in% names(markers)))
markers[, ens_base := sub("\\..*", "", gene)]
markers <- merge(markers, atlas[, .(ens_base, human_symbol)], by = "ens_base", all.x = TRUE)
cat(sprintf("  Marker DE rows: %d (genes x programs)\n", nrow(markers)))

twas <- fread(file.path(BASE, "RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv"))
twas_best <- twas[, .(twas_z = zscore[which.min(pvalue)],
                      twas_p = min(pvalue)), by = .(gene_name)]
cat(sprintf("  TWAS genes: %d\n", nrow(twas_best)))

sex_file <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
sex_degs <- fread(sex_file)
sex_degs[, ens_base := sub("\\..*", "", gene)]
sex_degs <- merge(sex_degs, atlas[, .(ens_base, human_symbol)], by = "ens_base", all.x = TRUE)
sex_classes <- unique(sex_degs$sex_class)
use_v2 <- "Female_biased" %in% sex_classes
female_lbl <- if (use_v2) "Female_biased" else "Female_specific"
male_lbl   <- if (use_v2) "Male_biased"   else "Male_specific"
female_specific <- sex_degs[sex_class == female_lbl & !is.na(human_symbol), unique(human_symbol)]
male_specific   <- sex_degs[sex_class == male_lbl   & !is.na(human_symbol), unique(human_symbol)]
cat(sprintf("  Female-biased=%d, Male-biased=%d\n",
            length(female_specific), length(male_specific)))

# ============================================================
# 2. Per-program up-marker gene sets
# ============================================================
cat("\n=== Building per-program up-marker gene sets ===\n")
program_up_sets <- list()
program_up_strict_sets <- list()  # padj<0.05 & logFC>0.3
for (pc in program_codes) {
  up_all <- markers[program_code == pc & direction == "up" & padj < 0.05 &
                    !is.na(human_symbol), unique(human_symbol)]
  up_strict <- markers[program_code == pc & direction == "up" & padj < 0.05 &
                       logFC > 0.3 & !is.na(human_symbol), unique(human_symbol)]
  program_up_sets[[pc]] <- up_all
  program_up_strict_sets[[pc]] <- up_strict
  cat(sprintf("  %s (%s): up padj<0.05=%d, up+LFC>0.3=%d\n",
              pc, prog_labels$biological_label[prog_labels$program_code == pc],
              length(up_all), length(up_strict)))
}

# ============================================================
# 3. Per-program Fisher enrichment x COLOC thresholds
# ============================================================
cat("\n=== Per-program Fisher enrichment vs COLOC ===\n")
universe <- unique(coloc$gene); universe <- universe[universe != "" & !is.na(universe)]
cat(sprintf("  COLOC universe: %d genes\n", length(universe)))

run_fisher <- function(coloc_set, test_set, univ, prog, thr) {
  coloc_u <- intersect(coloc_set, univ)
  test_u  <- intersect(test_set, univ)
  a <- length(intersect(coloc_u, test_u))
  b <- length(setdiff(coloc_u, test_u))
  c_ <- length(setdiff(test_u, coloc_u))
  d <- length(univ) - a - b - c_
  ft <- fisher.test(matrix(c(a,b,c_,d), nrow=2))
  data.table(
    program = prog, coloc_threshold = thr,
    n_universe = length(univ), n_coloc = length(coloc_u),
    n_program_up = length(test_u), n_overlap = a,
    expected_overlap = round(length(coloc_u)*length(test_u)/length(univ), 2),
    odds_ratio = as.numeric(ft$estimate),
    ci_lower = ft$conf.int[1], ci_upper = ft$conf.int[2],
    pvalue = ft$p.value,
    direction = ifelse(ft$estimate > 1, "enriched", "depleted")
  )
}

enrich_rows <- list()
for (pc in program_codes) {
  for (thr in c(0.5, 0.8, 0.9)) {
    cs <- if (thr == 0.5) coloc_05 else if (thr == 0.8) coloc_08 else coloc_09
    enrich_rows[[paste(pc, thr, sep="_")]] <- run_fisher(
      cs, program_up_sets[[pc]], universe, pc, thr
    )
  }
}
enrich_dt <- rbindlist(enrich_rows)
enrich_dt[, padj := p.adjust(pvalue, method = "BH")]
enrich_dt[, biological_label := prog_labels$biological_label[match(program, prog_labels$program_code)]]
setcolorder(enrich_dt, c("program","biological_label","coloc_threshold"))
fwrite(enrich_dt, file.path(outdir, "program_coloc_enrichment.csv"))
cat(sprintf("  program_coloc_enrichment.csv: %d rows\n", nrow(enrich_dt)))
print(enrich_dt[coloc_threshold == 0.5, .(program, biological_label,
                                           n_overlap, odds_ratio = round(odds_ratio,2),
                                           pvalue = signif(pvalue,3),
                                           padj = signif(padj,3))], row.names=FALSE)

# DEPRECATED ALIAS: subtype_coloc_enrichment.csv = Fibrotic program rows
if (fib_code != "") {
  fib_rows <- enrich_dt[program == fib_code]
  fib_rows[, gene_set := "S2_up_markers"]
  fwrite(fib_rows, file.path(outdir, "subtype_coloc_enrichment.csv"))
} else {
  fwrite(enrich_dt[0], file.path(outdir, "subtype_coloc_enrichment.csv"))
}

# ============================================================
# 4. Per-program causal scores (PP4 x |logFC|)
# ============================================================
cat("\n=== Per-program causal scores ===\n")
causal_scores <- copy(markers[, .(gene, ens_base, human_symbol,
                                   program_code, program, logFC, padj, direction)])
causal_scores <- merge(
  causal_scores,
  coloc[, .(human_symbol = gene, coloc_best_pp4, coloc_best_gwas)],
  by = "human_symbol", all.x = TRUE
)
causal_scores[is.na(coloc_best_pp4), coloc_best_pp4 := 0]
causal_scores[, program_causal_score := coloc_best_pp4 * abs(logFC)]
causal_scores[, program_causal_signed := coloc_best_pp4 * logFC]  # + = up in prog

causal_scores <- merge(causal_scores, twas_best[, .(human_symbol = gene_name, twas_z, twas_p)],
                       by = "human_symbol", all.x = TRUE)
causal_scores <- merge(causal_scores,
                       sex_degs[!is.na(human_symbol), .(human_symbol, sex_class,
                                                        logFC_F, logFC_M)],
                       by = "human_symbol", all.x = TRUE)
setorder(causal_scores, program_code, -program_causal_score)
causal_scores[, rank_in_program := seq_len(.N), by = program_code]
fwrite(causal_scores, file.path(outdir, "program_causal_scores.csv"))
cat(sprintf("  program_causal_scores.csv: %d rows\n", nrow(causal_scores)))

# Top 5 per program preview
cat("  Top-5 causal per program:\n")
for (pc in program_codes) {
  top5 <- causal_scores[program_code == pc & rank_in_program <= 5,
                        .(human_symbol, pp4 = round(coloc_best_pp4, 3),
                          logFC = round(logFC, 2),
                          score = round(program_causal_score, 3))]
  cat(sprintf("  %s (%s):\n", pc, prog_labels$biological_label[prog_labels$program_code == pc]))
  print(top5, row.names = FALSE)
}

# DEPRECATED ALIAS: subtype_causal_scores.csv = Fibrotic program rows
if (fib_code != "") {
  fib_causal <- causal_scores[program_code == fib_code]
  setnames(fib_causal, c("program_causal_score","program_causal_signed"),
           c("s2_causal_score","s2_causal_signed"))
  fwrite(fib_causal, file.path(outdir, "subtype_causal_scores.csv"))
} else {
  fwrite(causal_scores[0], file.path(outdir, "subtype_causal_scores.csv"))
}

# ============================================================
# 5. Multinomial LOCO-CV
# ============================================================
cat("\n=== Multinomial LOCO-CV: predict dominant_program from COLOC+TWAS ===\n")
dge <- readRDS(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds"))
meta <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"))

nmf_meta <- merge(nmf, meta[, .(sample_id, dataset)], by = "sample_id")
nmf_samples <- intersect(nmf_meta$sample_id, colnames(dge))
nmf_meta <- nmf_meta[sample_id %in% nmf_samples]
dge_sub <- dge[, nmf_samples]
dge_sub <- calcNormFactors(dge_sub, method = "TMM")
v <- voom(dge_sub, design = NULL, plot = FALSE)
logcpm <- v$E
# NOTE (2026-04-22, T2.14): removeBatchEffect moved inside the LOCO-CV loop
# (see run_loco_multinom below) so batch correction is re-estimated per fold,
# excluding the held-out dataset. The prior pre-correction used all-dataset
# batch estimates, which leaked test-set information into training folds.
base_rows <- sub("\\..*", "", rownames(logcpm))
ens2sym <- setNames(atlas$human_symbol, atlas$ens_base)
mapped <- ens2sym[base_rows]
keep <- !is.na(mapped) & !duplicated(mapped)
expr_sym <- logcpm[keep, ]
rownames(expr_sym) <- mapped[keep]
cat(sprintf("  Expression matrix (symbols, UNCORRECTED; per-fold batch correction): %d genes x %d samples\n",
            nrow(expr_sym), ncol(expr_sym)))

coloc_feat <- intersect(coloc_05, rownames(expr_sym))
twas_sig   <- twas_best[twas_p < 0.05, gene_name]
coloc_twas <- union(coloc_feat, intersect(twas_sig, rownames(expr_sym)))
cat(sprintf("  COLOC features: %d; +TWAS: %d\n", length(coloc_feat), length(coloc_twas)))

y <- nmf_meta[match(colnames(expr_sym), sample_id), dominant_program_code]
names(y) <- colnames(expr_sym)
y_fac <- factor(y, levels = program_codes)
datasets <- nmf_meta[match(colnames(expr_sym), sample_id), dataset]
unique_ds <- unique(datasets)

run_loco_multinom <- function(feat_genes, label) {
  if (length(feat_genes) < 5) return(NULL)
  X_raw <- t(expr_sym[feat_genes, , drop = FALSE])
  # Output per-sample class probability matrix
  prob_mat <- matrix(NA_real_, nrow = length(y), ncol = length(program_codes),
                     dimnames = list(names(y), program_codes))
  for (ds in unique_ds) {
    test_idx  <- which(datasets == ds)
    train_idx <- which(datasets != ds)
    y_tr <- y_fac[train_idx]
    if (length(test_idx) < 5 || length(unique(y_tr)) < 2) next
    tr_tab <- table(y_tr)
    if (any(tr_tab < 3)) next

    # Per-fold batch correction (T2.14 fix): fit removeBatchEffect on
    # training samples only, then apply the learned batch coefficients to
    # the held-out fold. We re-run on the feature-subset matrix (transposed
    # back to genes-x-samples) so batch effects are re-estimated per fold.
    train_mat <- t(X_raw[train_idx, , drop = FALSE])  # features x samples
    test_mat  <- t(X_raw[test_idx,  , drop = FALSE])
    train_batch <- datasets[train_idx]
    # removeBatchEffect on training fold (batch has >=2 levels for LOCO since
    # train excludes exactly one dataset)
    train_corr <- tryCatch(
      limma::removeBatchEffect(train_mat, batch = train_batch),
      error = function(e) train_mat
    )
    # Test fold is a single held-out dataset: re-centre by subtracting the
    # held-out-fold gene means (matches what removeBatchEffect would produce
    # if the test batch had been included as a factor level).
    test_corr <- test_mat - rowMeans(test_mat, na.rm = TRUE) +
                   rowMeans(train_corr, na.rm = TRUE)
    Xtr <- t(train_corr)
    Xte <- t(test_corr)
    # Drop constant features in training fold
    var_tr <- apply(Xtr, 2, var)
    keep_f <- which(var_tr > 1e-10)
    if (length(keep_f) < 5) next
    Xtr2 <- Xtr[, keep_f, drop = FALSE]
    Xte2 <- Xte[, keep_f, drop = FALSE]
    # Scale using training stats
    mu <- colMeans(Xtr2); sdv <- apply(Xtr2, 2, sd)
    Xtr2 <- scale(Xtr2, center = mu, scale = sdv)
    Xte2 <- scale(Xte2, center = mu, scale = sdv)
    fit <- tryCatch(
      nnet::multinom(y_tr ~ ., data = data.frame(y_tr = y_tr, Xtr2),
                     trace = FALSE, MaxNWts = 1e5, maxit = 100),
      error = function(e) { cat(sprintf("    multinom failed fold %s: %s\n", ds, e$message)); NULL })
    if (is.null(fit)) next
    preds <- tryCatch(predict(fit, newdata = data.frame(Xte2), type = "probs"),
                      error = function(e) NULL)
    if (is.null(preds)) next
    if (is.null(dim(preds))) preds <- matrix(preds, nrow = 1)
    # nnet::multinom returns one column per class in sorted levels order
    cls <- fit$lev
    for (ci in seq_along(cls)) {
      prob_mat[test_idx, cls[ci]] <- preds[, ci]
    }
  }
  # Per-program one-vs-rest AUROC
  auroc_rows <- list()
  for (pc in program_codes) {
    y_bin <- as.integer(y == pc)
    probs <- prob_mat[, pc]
    ok <- !is.na(probs)
    if (sum(ok) < 20 || length(unique(y_bin[ok])) < 2) {
      auroc_rows[[pc]] <- data.table(
        model = label, program = pc,
        biological_label = prog_labels$biological_label[prog_labels$program_code == pc],
        auroc = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        n_positives = sum(y_bin[ok] == 1), n_negatives = sum(y_bin[ok] == 0),
        n_features = length(feat_genes)
      )
      next
    }
    roc_ <- suppressMessages(roc(y_bin[ok], probs[ok], quiet = TRUE))
    ci_  <- as.numeric(ci.auc(roc_, conf.level = 0.95))
    auroc_rows[[pc]] <- data.table(
      model = label, program = pc,
      biological_label = prog_labels$biological_label[prog_labels$program_code == pc],
      auroc = as.numeric(auc(roc_)),
      ci_lower = ci_[1], ci_upper = ci_[3],
      n_positives = sum(y_bin[ok] == 1),
      n_negatives = sum(y_bin[ok] == 0),
      n_features = length(feat_genes)
    )
  }
  list(auroc = rbindlist(auroc_rows), prob_mat = prob_mat)
}

result_coloc <- run_loco_multinom(coloc_feat, "COLOC_only")
result_coltw <- run_loco_multinom(coloc_twas, "COLOC_plus_TWAS")
set.seed(42)
random_genes <- sample(rownames(expr_sym), min(length(coloc_feat), nrow(expr_sym)))
result_rand  <- run_loco_multinom(random_genes, "random_baseline")

auroc_out <- rbindlist(list(result_coloc$auroc, result_coltw$auroc, result_rand$auroc), fill = TRUE)
fwrite(auroc_out, file.path(outdir, "program_multinom_auroc.csv"))
cat("\n  Multinomial LOCO-CV AUROC summary:\n")
print(auroc_out[, .(model, program, biological_label,
                    auroc = round(auroc,3),
                    ci = sprintf("[%.2f-%.2f]", ci_lower, ci_upper),
                    n_pos = n_positives, n_neg = n_negatives)],
      row.names = FALSE)

# Save per-sample predictions (from COLOC+TWAS)
if (!is.null(result_coltw$prob_mat)) {
  pm <- as.data.table(result_coltw$prob_mat, keep.rownames = "sample_id")
  pm[, model := "COLOC_plus_TWAS"]
  pm[, y_true := y[match(pm$sample_id, names(y))]]
  fwrite(pm, file.path(outdir, "program_multinom_predictions.csv"))
}

# DEPRECATED ALIAS: s2_genetic_risk.csv = Fibrotic-program AUROC rows
if (fib_code != "") {
  fib_auroc <- auroc_out[program == fib_code, .(model, auroc,
                                                auroc_ci_lower = ci_lower,
                                                auroc_ci_upper = ci_upper,
                                                n_features, n_samples = n_positives + n_negatives)]
  fwrite(fib_auroc, file.path(outdir, "s2_genetic_risk.csv"))
} else {
  fwrite(data.table(), file.path(outdir, "s2_genetic_risk.csv"))
}

# ============================================================
# 6. Per-program sex decomposition
# ============================================================
cat("\n=== Per-program sex decomposition ===\n")
sex_rows <- list()
for (pc in program_codes) {
  up_set <- program_up_sets[[pc]]
  up_coloc_set <- intersect(up_set, coloc_05)
  n_up_female <- length(intersect(up_set, female_specific))
  n_up_male   <- length(intersect(up_set, male_specific))
  n_up_coloc_female <- length(intersect(up_coloc_set, female_specific))
  n_up_coloc_male   <- length(intersect(up_coloc_set, male_specific))
  # Fisher: is up-coloc set enriched for female-biased DEGs?
  if (length(up_coloc_set) >= 3) {
    a <- n_up_coloc_female
    b <- length(up_coloc_set) - a
    c_ <- length(intersect(coloc_05, female_specific)) - a
    d <- length(coloc_05) - a - b - c_
    c_ <- max(c_, 0); d <- max(d, 0)
    ft <- suppressWarnings(fisher.test(matrix(c(a,b,c_,d), nrow=2)))
    or_f <- as.numeric(ft$estimate); p_f <- ft$p.value
  } else { or_f <- NA_real_; p_f <- NA_real_ }
  sex_rows[[pc]] <- data.table(
    program = pc,
    biological_label = prog_labels$biological_label[prog_labels$program_code == pc],
    n_up = length(up_set), n_up_coloc = length(up_coloc_set),
    n_up_female = n_up_female, n_up_male = n_up_male,
    n_up_coloc_female = n_up_coloc_female, n_up_coloc_male = n_up_coloc_male,
    female_or = or_f, female_p = p_f
  )
}
sex_decomp <- rbindlist(sex_rows)
sex_decomp[, female_padj := p.adjust(female_p, "BH")]
fwrite(sex_decomp, file.path(outdir, "program_sex_decomposition.csv"))
print(sex_decomp, row.names = FALSE)

# DEPRECATED ALIAS: s2_sex_decomposition.csv
if (fib_code != "") {
  fwrite(sex_decomp[program == fib_code], file.path(outdir, "s2_sex_decomposition.csv"))
} else {
  fwrite(data.table(), file.path(outdir, "s2_sex_decomposition.csv"))
}

# ============================================================
# 7. Summary
# ============================================================
cat("\n=== Final summary ===\n")
sig_enrich <- enrich_dt[padj < 0.05 & direction == "enriched"]
cat(sprintf("  Significant program x COLOC enrichments (padj<0.05): %d\n", nrow(sig_enrich)))
if (nrow(sig_enrich)) {
  print(sig_enrich[, .(program, biological_label, coloc_threshold, odds_ratio = round(odds_ratio,2),
                       padj = signif(padj,3))], row.names = FALSE)
}
best_auroc <- auroc_out[model == "COLOC_plus_TWAS"][which.max(auroc)]
if (nrow(best_auroc)) {
  cat(sprintf("  Best program predictability (COLOC+TWAS): %s (%s) AUROC=%.3f\n",
              best_auroc$program, best_auroc$biological_label, best_auroc$auroc))
}
cat(sprintf("\n=== 208_subtype_coloc.R complete (k=%d programs) ===\n", K))
cat("End:", format(Sys.time()), "\n")
