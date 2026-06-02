#!/usr/bin/env Rscript
# 217b_hormone_tf_coloc_enrichment.R — Fisher + log-linear: TF target x sex-class x COLOC
#
# Agent: B5 (Sex-Hormone TF x COLOC) -- HEADLINE A7
#
# For each sex-hormone TF, ask:
#   (a) Are TF targets enriched among sex-biased DEGs?
#   (b) Are TF targets enriched among high-PP4 COLOC genes?
#   (c) Triple-intersection: TF target AND sex-biased AND COLOC -- log-linear
#   (d) Stratified by sex direction (Female_biased vs Male_biased)
#
# Inputs:
#   - RNA-seq/results/stratified_causal/hormone_tf_targetsets_long.csv (from 217a)
#   - RNA-seq/Human/.../results/integration/sex_deg_classification.csv
#   - GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#   - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (ID mapping)
#
# Outputs:
#   - RNA-seq/results/stratified_causal/hormone_tf_sex_coloc.csv (primary)
#   - RNA-seq/results/stratified_causal/hormone_tf_loglin.csv
#   - RNA-seq/results/stratified_causal/hormone_tf_triple_genes.csv

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("=== 217b: Hormone TF x sex x COLOC enrichment ===\n")

# ===========================================================================
# 1. Load inputs
# ===========================================================================
tf_long <- fread(file.path(outdir, "hormone_tf_targetsets_long.csv"))
cat("  TF target edges:", nrow(tf_long), "\n")

atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("ensembl_id", "human_symbol"))

# Map TF target symbols to Ensembl
# Some targets may already be Ensembl IDs (in SCENIC); detect and pass-through
tf_long[, is_ens := grepl("^ENSG", target_gene)]
tf_sym <- merge(tf_long[!(is_ens)], atlas, by.x = "target_gene",
                by.y = "human_symbol", all.x = TRUE)
tf_ens <- tf_long[(is_ens)]
tf_ens[, ensembl_id := target_gene]
tf_map <- rbind(tf_sym, tf_ens, fill = TRUE)
tf_map <- tf_map[!is.na(ensembl_id) & ensembl_id != ""]

# Sex DEGs
sex_file <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
sex_degs <- fread(sex_file)
sex_degs[, ensembl_id := sub("\\..*", "", gene)]
cat("  Sex DEGs:", nrow(sex_degs), "rows\n")
cat("  Sex class table:\n"); print(sex_degs[, .N, by = sex_class][order(-N)])

# COLOC
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]
# SuSiE preferred (T0.4 convention)
if ("coloc_best_susie_pp4" %in% names(coloc)) {
  coloc[, coloc_best_abf_pp4 := coloc_best_pp4]
  coloc[, coloc_best_pp4 := fifelse(!is.na(coloc_best_susie_pp4),
                                    coloc_best_susie_pp4, coloc_best_pp4)]
}
# Map coloc symbols -> ensembl via atlas
coloc <- merge(coloc, atlas, by.x = "gene", by.y = "human_symbol", all.x = TRUE)
coloc[is.na(ensembl_id), ensembl_id := ensembl]
coloc_best <- coloc[, .(coloc_best_pp4 = max(coloc_best_pp4, na.rm = TRUE),
                         coloc_best_gwas = coloc_best_gwas[which.max(coloc_best_pp4)]),
                    by = ensembl_id]
cat("  COLOC genes:", nrow(coloc_best), "\n")

# ===========================================================================
# 2. Build universe = sex_degs (full transcriptome tested) joined with COLOC
# ===========================================================================
universe <- merge(sex_degs[, .(ensembl_id, sex_class, logFC_M, logFC_F,
                                padj_M, padj_F, interaction_padj)],
                  coloc_best, by = "ensembl_id", all.x = TRUE)
universe[is.na(coloc_best_pp4) | !is.finite(coloc_best_pp4),
         coloc_best_pp4 := 0]
cat("  Universe:", nrow(universe), "genes\n")

# Sex class flags
universe[, sex_biased := sex_class %in% c("Female_biased", "Male_biased", "Divergent",
                                           "Female_specific", "Male_specific")]
universe[, female_biased := sex_class %in% c("Female_biased", "Female_specific")]
universe[, male_biased := sex_class %in% c("Male_biased", "Male_specific")]

# ===========================================================================
# 3. Per-TF Fisher tests
# ===========================================================================
TF_PANEL <- sort(unique(tf_long$tf))
cat("\n--- Per-TF Fisher tests ---\n")

fisher_results <- list()
loglin_results <- list()
triple_gene_list <- list()

PP4_THRESHOLDS <- c(0.5, 0.8)

for (tf_i in TF_PANEL) {
  # Note: avoid `tf` as loop var since column is also named `tf`
  target_ens <- unique(tf_map[tf == tf_i, ensembl_id])
  tf <- tf_i  # legacy variable preserved for downstream string assembly
  is_target <- universe$ensembl_id %in% target_ens
  n_tgt <- sum(is_target)
  if (n_tgt < 5) {
    cat(sprintf("  %-10s  n_targets=%d  SKIPPED (too few)\n", tf, n_tgt))
    next
  }

  for (pp4_thr in PP4_THRESHOLDS) {
    is_coloc <- universe$coloc_best_pp4 >= pp4_thr
    n_coloc <- sum(is_coloc)

    # 2x2: target x sex_biased
    f_sex <- fisher.test(table(is_target, universe$sex_biased))
    # 2x2: target x female_biased
    f_fem <- if (sum(universe$female_biased) >= 5)
      fisher.test(table(is_target, universe$female_biased)) else NULL
    # 2x2: target x male_biased
    f_mal <- if (sum(universe$male_biased) >= 5)
      fisher.test(table(is_target, universe$male_biased)) else NULL
    # 2x2: target x coloc
    f_col <- if (n_coloc >= 5)
      fisher.test(table(is_target, is_coloc)) else NULL
    # 2x2: target x (sex_biased AND coloc)  -- the triple
    triple <- is_target & universe$sex_biased & is_coloc
    f_trp <- if (sum(triple) >= 1 && n_coloc >= 5)
      tryCatch(fisher.test(table(is_target, universe$sex_biased & is_coloc)),
               error = function(e) NULL) else NULL

    fisher_results[[paste(tf, pp4_thr, sep = "_")]] <- data.table(
      tf = tf,
      pp4_threshold = pp4_thr,
      n_targets = n_tgt,
      n_sex_biased_targets = sum(is_target & universe$sex_biased),
      n_female_biased_targets = sum(is_target & universe$female_biased),
      n_male_biased_targets = sum(is_target & universe$male_biased),
      n_coloc_targets = sum(is_target & is_coloc),
      n_triple = sum(triple),
      OR_sex = if (!is.null(f_sex)) f_sex$estimate else NA_real_,
      P_sex = if (!is.null(f_sex)) f_sex$p.value else NA_real_,
      OR_female = if (!is.null(f_fem)) f_fem$estimate else NA_real_,
      P_female = if (!is.null(f_fem)) f_fem$p.value else NA_real_,
      OR_male = if (!is.null(f_mal)) f_mal$estimate else NA_real_,
      P_male = if (!is.null(f_mal)) f_mal$p.value else NA_real_,
      OR_coloc = if (!is.null(f_col)) f_col$estimate else NA_real_,
      P_coloc = if (!is.null(f_col)) f_col$p.value else NA_real_,
      OR_triple = if (!is.null(f_trp)) f_trp$estimate else NA_real_,
      P_triple = if (!is.null(f_trp)) f_trp$p.value else NA_real_
    )

    # Triple-intersection gene list (PP4 >= 0.5)
    if (pp4_thr == 0.5) {
      trip_genes <- universe[triple,
                              .(tf = tf, ensembl_id, sex_class,
                                logFC_M, logFC_F, interaction_padj,
                                coloc_best_pp4, coloc_best_gwas)]
      if (nrow(trip_genes) > 0) {
        triple_gene_list[[tf]] <- trip_genes
      }
    }

    # ===== 3-way log-linear (target x sex_biased x coloc) =====
    if (pp4_thr == 0.5) {
      tbl3 <- table(target = is_target,
                    sex_biased = universe$sex_biased,
                    coloc = is_coloc)
      if (all(dim(tbl3) == c(2, 2, 2)) && all(tbl3 > 0)) {
        # Saturated 3-way model fit; report higher-order interaction LRT
        fit_full <- tryCatch(loglin(tbl3, list(c(1,2,3)),
                                    fit = FALSE, print = FALSE),
                             error = function(e) NULL)
        fit_homog <- tryCatch(loglin(tbl3, list(c(1,2), c(1,3), c(2,3)),
                                     fit = FALSE, print = FALSE),
                              error = function(e) NULL)
        if (!is.null(fit_homog)) {
          lr <- fit_homog$lrt
          df_h <- fit_homog$df
          p_lr <- pchisq(lr, df_h, lower.tail = FALSE)
          loglin_results[[tf]] <- data.table(
            tf = tf,
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

  fr <- fisher_results[[paste(tf, 0.5, sep = "_")]]
  if (!is.null(fr)) {
    cat(sprintf("  %-10s  n_tgt=%-4d  sex_OR=%.2f p=%.2e  coloc_OR=%.2f p=%.2e  triple=%d\n",
                tf, fr$n_targets, fr$OR_sex, fr$P_sex,
                fr$OR_coloc, fr$P_coloc, fr$n_triple))
  }
}

# ===========================================================================
# 4. Combine + BH correct
# ===========================================================================
fisher_dt <- rbindlist(fisher_results, use.names = TRUE, fill = TRUE)
for (col in c("P_sex", "P_female", "P_male", "P_coloc", "P_triple")) {
  padj_col <- sub("^P_", "padj_", col)
  fisher_dt[, (padj_col) := p.adjust(get(col), method = "BH"),
            by = pp4_threshold]
}
setcolorder(fisher_dt,
            c("tf", "pp4_threshold", "n_targets",
              "n_sex_biased_targets", "n_female_biased_targets",
              "n_male_biased_targets", "n_coloc_targets", "n_triple",
              "OR_sex", "P_sex", "padj_sex",
              "OR_female", "P_female", "padj_female",
              "OR_male", "P_male", "padj_male",
              "OR_coloc", "P_coloc", "padj_coloc",
              "OR_triple", "P_triple", "padj_triple"))
fisher_dt <- fisher_dt[order(pp4_threshold, padj_triple, P_triple, na.last = TRUE)]
fwrite(fisher_dt, file.path(outdir, "hormone_tf_sex_coloc.csv"))
cat("\n  Wrote hormone_tf_sex_coloc.csv (", nrow(fisher_dt), "rows)\n")

loglin_dt <- rbindlist(loglin_results, use.names = TRUE, fill = TRUE)
if (nrow(loglin_dt) > 0) {
  loglin_dt[, padj_3way := p.adjust(p_3way, method = "BH")]
  loglin_dt <- loglin_dt[order(padj_3way, p_3way)]
}
fwrite(loglin_dt, file.path(outdir, "hormone_tf_loglin.csv"))
cat("  Wrote hormone_tf_loglin.csv (", nrow(loglin_dt), "rows)\n")

triple_dt <- rbindlist(triple_gene_list, use.names = TRUE, fill = TRUE)
# Annotate with human_symbol
triple_dt <- merge(triple_dt, atlas, by = "ensembl_id", all.x = TRUE)
triple_dt <- triple_dt[order(tf, -coloc_best_pp4)]
fwrite(triple_dt, file.path(outdir, "hormone_tf_triple_genes.csv"))
cat("  Wrote hormone_tf_triple_genes.csv (", nrow(triple_dt), "rows)\n")

cat("\n=== 217b complete ===\n")
