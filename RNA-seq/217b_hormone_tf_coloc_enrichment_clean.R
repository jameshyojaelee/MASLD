#!/usr/bin/env Rscript
# 217b_hormone_tf_coloc_enrichment_clean.R -- non-circular Fisher + permutation null
#
# Agent: B5b (C5 critique fix)
#
# CRITICAL CHANGES vs 217b:
#   1. Universe restricted to genes simultaneously tested in BOTH
#      sex_deg_classification.csv AND gene_level_coloc.csv,
#      AND protein-coding or lncRNA biotype (per gencode_v49_gene_metadata.tsv.gz).
#      Previous universe (~34K) over-counted untested / pseudogene genes.
#   2. Permutation null: B=1000 random target-label shuffles per TF.
#      Empirical p-value = fraction of permutations with OR >= observed.
#   3. Comparison logging vs B5 outputs (delta_OR per TF).
#
# Inputs:
#   - RNA-seq/results/stratified_causal/hormone_tf_targetsets_long_clean.csv
#   - RNA-seq/Human/.../results/integration/sex_deg_classification.csv
#   - GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#   - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (ID mapping)
#   - data/gencode_v49_gene_metadata.tsv.gz (biotype filter)
#
# Outputs:
#   - RNA-seq/results/stratified_causal/hormone_tf_sex_coloc_clean.csv (primary)
#   - RNA-seq/results/stratified_causal/hormone_tf_loglin_clean.csv
#   - RNA-seq/results/stratified_causal/hormone_tf_triple_genes_clean.csv
#   - RNA-seq/results/stratified_causal/hormone_tf_negative_controls.csv
#   - RNA-seq/results/stratified_causal/B5_vs_B5b_comparison.csv

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(42)
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "1000"))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

HORMONE_TF_PANEL <- c("AR", "ESR1", "ESR2", "FOXA1", "FOXA2",
                      "STAT5A", "STAT5B", "BCL6", "CUX2", "HNF4A")
NEG_CTRL_TFS <- c("CTCF", "MYC", "RFX5")

cat("=== 217b_clean: NON-CIRCULAR Fisher + permutation null ===\n")
cat("  N permutations:", N_PERM, "\n")

# ===========================================================================
# 1. Load inputs
# ===========================================================================
tf_long <- fread(file.path(outdir, "hormone_tf_targetsets_long_clean.csv"))
cat("  TF target edges (clean):", nrow(tf_long), "\n")

atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("ensembl_id", "human_symbol"))

# Map TF target symbols to Ensembl
tf_long[, is_ens := grepl("^ENSG", target_gene)]
tf_sym <- merge(tf_long[!(is_ens)], atlas, by.x = "target_gene",
                by.y = "human_symbol", all.x = TRUE)
tf_ens <- tf_long[(is_ens)]
tf_ens[, ensembl_id := target_gene]
tf_map <- rbind(tf_sym, tf_ens, fill = TRUE)
tf_map <- tf_map[!is.na(ensembl_id) & ensembl_id != ""]

# Sex DEGs (v3 mashr Bayesian preferred; v2 fallback — v3 provides v2-compat `sex_class`)
sex_v3_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/sex_deg_classification_v3.csv")
sex_v2_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
sex_file <- if (file.exists(sex_v3_path)) sex_v3_path else sex_v2_path
cat("  Source:", basename(dirname(sex_file)), "/", basename(sex_file), "\n")
sex_degs <- fread(sex_file)
sex_degs[, ensembl_id := sub("\\..*", "", gene)]
cat("  Sex DEGs:", nrow(sex_degs), "rows\n")

# COLOC
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]
if ("coloc_best_susie_pp4" %in% names(coloc)) {
  coloc[, coloc_best_abf_pp4 := coloc_best_pp4]
  coloc[, coloc_best_pp4 := fifelse(!is.na(coloc_best_susie_pp4),
                                    coloc_best_susie_pp4, coloc_best_pp4)]
}
coloc <- merge(coloc, atlas, by.x = "gene", by.y = "human_symbol", all.x = TRUE)
coloc[is.na(ensembl_id), ensembl_id := ensembl]
coloc_best <- coloc[, .(coloc_best_pp4 = max(coloc_best_pp4, na.rm = TRUE),
                         coloc_best_gwas = coloc_best_gwas[which.max(coloc_best_pp4)]),
                    by = ensembl_id]
cat("  COLOC gene-rows:", nrow(coloc_best), "\n")

# Biotype metadata
biotype_path <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
biotype <- fread(biotype_path,
                 select = c("ensembl_base", "gene_biotype"))
setnames(biotype, "ensembl_base", "ensembl_id")
biotype <- unique(biotype)
keep_biotypes <- c("protein_coding", "lncRNA")
biotype_keep <- biotype[gene_biotype %in% keep_biotypes, ensembl_id]
cat("  protein_coding+lncRNA ensembl_ids:", length(biotype_keep), "\n")

# ===========================================================================
# 2. RESTRICTED universe: tested in BOTH sex_degs AND coloc, AND
#    protein-coding or lncRNA biotype
# ===========================================================================
sex_tested <- unique(sex_degs$ensembl_id)
coloc_tested <- unique(coloc_best$ensembl_id)
both_tested <- intersect(sex_tested, coloc_tested)
universe_ids <- intersect(both_tested, biotype_keep)
cat("  Universe (sex AND coloc AND PC/lncRNA):", length(universe_ids), "genes\n")
cat("    (compare to B5's 27,638-gene tested universe; STAR -s 2, 2026-05-28 P0-J)\n")

universe <- merge(sex_degs[ensembl_id %in% universe_ids,
                            .(ensembl_id, sex_class, logFC_M, logFC_F,
                              padj_M, padj_F, interaction_padj)],
                  coloc_best, by = "ensembl_id", all.x = TRUE)
universe[is.na(coloc_best_pp4) | !is.finite(coloc_best_pp4),
         coloc_best_pp4 := 0]
cat("  Universe rows:", nrow(universe), "\n")

universe[, sex_biased := sex_class %in% c("Female_biased", "Male_biased", "Divergent",
                                           "Female_specific", "Male_specific")]
universe[, female_biased := sex_class %in% c("Female_biased", "Female_specific")]
universe[, male_biased := sex_class %in% c("Male_biased", "Male_specific")]

# ===========================================================================
# 3. Per-TF Fisher + permutation null
# ===========================================================================
ALL_TF <- sort(unique(tf_long$tf))
cat("\n--- Per-TF Fisher + permutation ---\n")

fisher_results <- list()
loglin_results <- list()
triple_gene_list <- list()
neg_ctrl_results <- list()
PP4_THRESHOLDS <- c(0.5, 0.8)

permute_OR <- function(is_target_logical, is_test_logical, B = N_PERM) {
  # Permute the target labels B times, compute OR each time
  n_uni <- length(is_target_logical)
  n_tgt <- sum(is_target_logical)
  observed <- {
    tbl <- table(factor(is_target_logical, levels = c(FALSE, TRUE)),
                 factor(is_test_logical, levels = c(FALSE, TRUE)))
    # OR with Haldane correction
    a <- tbl[2,2] + 0.5; b <- tbl[2,1] + 0.5
    c <- tbl[1,2] + 0.5; d <- tbl[1,1] + 0.5
    (a / b) / (c / d)
  }
  perm_or <- numeric(B)
  for (i in seq_len(B)) {
    perm_tgt <- logical(n_uni)
    perm_tgt[sample.int(n_uni, n_tgt)] <- TRUE
    tbl <- table(factor(perm_tgt, levels = c(FALSE, TRUE)),
                 factor(is_test_logical, levels = c(FALSE, TRUE)))
    a <- tbl[2,2] + 0.5; b <- tbl[2,1] + 0.5
    c <- tbl[1,2] + 0.5; d <- tbl[1,1] + 0.5
    perm_or[i] <- (a / b) / (c / d)
  }
  list(observed = observed, perm = perm_or,
       emp_p = mean(perm_or >= observed),
       emp_p_two_sided = mean(abs(log(perm_or)) >= abs(log(observed))))
}

for (tf_i in ALL_TF) {
  target_ens <- unique(tf_map[tf == tf_i, ensembl_id])
  is_target <- universe$ensembl_id %in% target_ens
  n_tgt <- sum(is_target)
  is_neg <- tf_i %in% NEG_CTRL_TFS
  if (n_tgt < 5) {
    cat(sprintf("  %-10s  n_targets=%d  SKIPPED (too few)\n", tf_i, n_tgt))
    next
  }

  for (pp4_thr in PP4_THRESHOLDS) {
    is_coloc <- universe$coloc_best_pp4 >= pp4_thr
    n_coloc <- sum(is_coloc)

    f_sex <- fisher.test(table(is_target, universe$sex_biased))
    f_fem <- if (sum(universe$female_biased) >= 5)
      fisher.test(table(is_target, universe$female_biased)) else NULL
    f_mal <- if (sum(universe$male_biased) >= 5)
      fisher.test(table(is_target, universe$male_biased)) else NULL
    f_col <- if (n_coloc >= 5)
      fisher.test(table(is_target, is_coloc)) else NULL
    triple <- is_target & universe$sex_biased & is_coloc
    f_trp <- if (sum(triple) >= 1 && n_coloc >= 5)
      tryCatch(fisher.test(table(is_target, universe$sex_biased & is_coloc)),
               error = function(e) NULL) else NULL

    # Permutation null for the key OR (sex AND coloc)
    if (n_coloc >= 5 && pp4_thr == 0.5) {
      perm_triple <- permute_OR(is_target, universe$sex_biased & is_coloc, B = N_PERM)
      emp_p_triple <- perm_triple$emp_p
      emp_p_triple_two <- perm_triple$emp_p_two_sided
      perm_sex <- permute_OR(is_target, universe$sex_biased, B = N_PERM)
      emp_p_sex <- perm_sex$emp_p_two_sided
      perm_col <- permute_OR(is_target, is_coloc, B = N_PERM)
      emp_p_coloc <- perm_col$emp_p_two_sided
    } else {
      emp_p_triple <- NA_real_; emp_p_triple_two <- NA_real_
      emp_p_sex <- NA_real_; emp_p_coloc <- NA_real_
    }

    fisher_results[[paste(tf_i, pp4_thr, sep = "_")]] <- data.table(
      tf = tf_i,
      is_negative_control = is_neg,
      pp4_threshold = pp4_thr,
      n_targets = n_tgt,
      n_sex_biased_targets = sum(is_target & universe$sex_biased),
      n_female_biased_targets = sum(is_target & universe$female_biased),
      n_male_biased_targets = sum(is_target & universe$male_biased),
      n_coloc_targets = sum(is_target & is_coloc),
      n_triple = sum(triple),
      OR_sex = if (!is.null(f_sex)) f_sex$estimate else NA_real_,
      P_sex = if (!is.null(f_sex)) f_sex$p.value else NA_real_,
      emp_P_sex = emp_p_sex,
      OR_female = if (!is.null(f_fem)) f_fem$estimate else NA_real_,
      P_female = if (!is.null(f_fem)) f_fem$p.value else NA_real_,
      OR_male = if (!is.null(f_mal)) f_mal$estimate else NA_real_,
      P_male = if (!is.null(f_mal)) f_mal$p.value else NA_real_,
      OR_coloc = if (!is.null(f_col)) f_col$estimate else NA_real_,
      P_coloc = if (!is.null(f_col)) f_col$p.value else NA_real_,
      emp_P_coloc = emp_p_coloc,
      OR_triple = if (!is.null(f_trp)) f_trp$estimate else NA_real_,
      P_triple = if (!is.null(f_trp)) f_trp$p.value else NA_real_,
      emp_P_triple = emp_p_triple,
      emp_P_triple_two_sided = emp_p_triple_two
    )

    if (pp4_thr == 0.5) {
      trip_genes <- universe[triple,
                              .(tf = tf_i, ensembl_id, sex_class,
                                logFC_M, logFC_F, interaction_padj,
                                coloc_best_pp4, coloc_best_gwas)]
      if (nrow(trip_genes) > 0) {
        triple_gene_list[[tf_i]] <- trip_genes
      }
    }

    if (pp4_thr == 0.5) {
      tbl3 <- table(target = is_target,
                    sex_biased = universe$sex_biased,
                    coloc = is_coloc)
      if (all(dim(tbl3) == c(2, 2, 2)) && all(tbl3 > 0)) {
        fit_homog <- tryCatch(loglin(tbl3, list(c(1,2), c(1,3), c(2,3)),
                                     fit = FALSE, print = FALSE),
                              error = function(e) NULL)
        if (!is.null(fit_homog)) {
          lr <- fit_homog$lrt
          df_h <- fit_homog$df
          p_lr <- pchisq(lr, df_h, lower.tail = FALSE)
          loglin_results[[tf_i]] <- data.table(
            tf = tf_i,
            is_negative_control = is_neg,
            lrt_3way = lr,
            df = df_h,
            p_3way = p_lr,
            n_target = n_tgt,
            n_sex_biased = sum(universe$sex_biased),
            n_coloc = n_coloc,
            n_triple = sum(triple))
        }
      }
    }
  }

  fr <- fisher_results[[paste(tf_i, 0.5, sep = "_")]]
  if (!is.null(fr)) {
    cat(sprintf("  %-10s%s  n_tgt=%-4d  sex_OR=%.2f(emp_p=%.3f)  coloc_OR=%.2f(emp_p=%.3f)  triple=%d  emp_p_triple=%.3f\n",
                tf_i, ifelse(is_neg, " [NEG]", "     "),
                fr$n_targets, fr$OR_sex, fr$emp_P_sex,
                fr$OR_coloc, fr$emp_P_coloc,
                fr$n_triple, fr$emp_P_triple))
  }
}

# ===========================================================================
# 4. Combine + BH correct (split hormone vs negative control)
# ===========================================================================
fisher_dt <- rbindlist(fisher_results, use.names = TRUE, fill = TRUE)
for (col in c("P_sex", "P_female", "P_male", "P_coloc", "P_triple")) {
  padj_col <- sub("^P_", "padj_", col)
  fisher_dt[, (padj_col) := p.adjust(get(col), method = "BH"),
            by = pp4_threshold]
}
setcolorder(fisher_dt,
            c("tf", "is_negative_control", "pp4_threshold", "n_targets",
              "n_sex_biased_targets", "n_female_biased_targets",
              "n_male_biased_targets", "n_coloc_targets", "n_triple",
              "OR_sex", "P_sex", "padj_sex", "emp_P_sex",
              "OR_female", "P_female", "padj_female",
              "OR_male", "P_male", "padj_male",
              "OR_coloc", "P_coloc", "padj_coloc", "emp_P_coloc",
              "OR_triple", "P_triple", "padj_triple",
              "emp_P_triple", "emp_P_triple_two_sided"))
fisher_dt <- fisher_dt[order(is_negative_control, pp4_threshold, padj_triple,
                              P_triple, na.last = TRUE)]
fwrite(fisher_dt, file.path(outdir, "hormone_tf_sex_coloc_clean.csv"))
cat("\n  Wrote hormone_tf_sex_coloc_clean.csv (", nrow(fisher_dt), "rows)\n")

# Negative-control summary
neg_ctrl_dt <- fisher_dt[is_negative_control == TRUE]
fwrite(neg_ctrl_dt, file.path(outdir, "hormone_tf_negative_controls.csv"))
cat("  Wrote hormone_tf_negative_controls.csv (", nrow(neg_ctrl_dt), "rows)\n")

loglin_dt <- rbindlist(loglin_results, use.names = TRUE, fill = TRUE)
if (nrow(loglin_dt) > 0) {
  loglin_dt[, padj_3way := p.adjust(p_3way, method = "BH")]
  loglin_dt <- loglin_dt[order(padj_3way, p_3way)]
}
fwrite(loglin_dt, file.path(outdir, "hormone_tf_loglin_clean.csv"))
cat("  Wrote hormone_tf_loglin_clean.csv (", nrow(loglin_dt), "rows)\n")

triple_dt <- rbindlist(triple_gene_list, use.names = TRUE, fill = TRUE)
if (nrow(triple_dt) > 0) {
  triple_dt <- merge(triple_dt, atlas, by = "ensembl_id", all.x = TRUE)
  triple_dt <- triple_dt[order(tf, -coloc_best_pp4)]
}
fwrite(triple_dt, file.path(outdir, "hormone_tf_triple_genes_clean.csv"))
cat("  Wrote hormone_tf_triple_genes_clean.csv (", nrow(triple_dt), "rows)\n")

# ===========================================================================
# 5. B5 vs B5b comparison table -- the headline diagnostic
# ===========================================================================
b5_path <- file.path(outdir, "hormone_tf_sex_coloc.csv")
if (file.exists(b5_path)) {
  b5_dt <- fread(b5_path)
  b5_05 <- b5_dt[pp4_threshold == 0.5,
                  .(tf, OR_sex_B5 = OR_sex, OR_coloc_B5 = OR_coloc,
                    OR_triple_B5 = OR_triple,
                    n_triple_B5 = n_triple)]
  b5b_05 <- fisher_dt[pp4_threshold == 0.5,
                       .(tf, is_negative_control,
                         OR_sex_B5b = OR_sex, OR_coloc_B5b = OR_coloc,
                         OR_triple_B5b = OR_triple,
                         n_triple_B5b = n_triple,
                         emp_P_triple_B5b = emp_P_triple)]
  cmp <- merge(b5b_05, b5_05, by = "tf", all = TRUE)
  cmp[, delta_OR_sex := OR_sex_B5 - OR_sex_B5b]
  cmp[, delta_OR_coloc := OR_coloc_B5 - OR_coloc_B5b]
  cmp[, delta_OR_triple := OR_triple_B5 - OR_triple_B5b]
  # circularity impact score: relative inflation of triple OR by disease regulons
  cmp[, circularity_impact_score := ifelse(
    !is.na(OR_triple_B5) & !is.na(OR_triple_B5b) & OR_triple_B5b > 0,
    (OR_triple_B5 - OR_triple_B5b) / OR_triple_B5b, NA_real_)]
  setcolorder(cmp,
              c("tf", "is_negative_control",
                "OR_sex_B5", "OR_sex_B5b", "delta_OR_sex",
                "OR_coloc_B5", "OR_coloc_B5b", "delta_OR_coloc",
                "OR_triple_B5", "OR_triple_B5b", "delta_OR_triple",
                "n_triple_B5", "n_triple_B5b",
                "emp_P_triple_B5b",
                "circularity_impact_score"))
  cmp <- cmp[order(is_negative_control, -abs(circularity_impact_score),
                    na.last = TRUE)]
  fwrite(cmp, file.path(outdir, "B5_vs_B5b_comparison.csv"))
  cat("  Wrote B5_vs_B5b_comparison.csv (", nrow(cmp), "rows)\n")
  cat("\nCirculary impact summary (top by |delta_OR_triple|):\n")
  print(cmp[!is.na(delta_OR_triple), .(tf, is_negative_control,
                                         OR_triple_B5, OR_triple_B5b,
                                         delta_OR_triple,
                                         circularity_impact_score)])
} else {
  cat("  B5 hormone_tf_sex_coloc.csv not found -- comparison skipped\n")
}

cat("\n=== 217b_clean complete ===\n")
