#!/usr/bin/env Rscript
# 218d_sex_progression_coloc.R — A9 — Sex × progression × COLOC
#
# For each fibrosis-stage transition (F0→F1, F1→F2, F2→F3, F3→F4):
#   2×2 fisher: sex_class (Female_biased | Male_biased) × coloc (yes | no)
# Tests whether female-biased COLOC concentrates in early vs late transitions
# (or vice versa).
#
# Outputs:
#   RNA-seq/results/stratified_causal/sex_progression_coloc.csv

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n", sep = " ")
log_msg("218d sex × progression × COLOC starting")

# ─── 1. Inputs ─────────────────────────────────────────────────────────────────
log_msg("Loading sex_deg classification")
# v3 mashr Bayesian preferred; v2 fallback. v3 CSV provides v2-compatible `sex_class` alias.
sex_v3_path <- file.path(BASE,
   "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/sex_deg_classification_v3.csv")
sex_v2_path <- file.path(BASE,
   "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
sex_deg_file <- if (file.exists(sex_v3_path)) sex_v3_path else sex_v2_path
log_msg(sprintf("sex_deg source: %s", basename(dirname(sex_deg_file))))
sex_deg <- fread(sex_deg_file)
sex_deg[, ensembl_base := sub("\\..*", "", gene)]
gencode <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
sex_deg <- merge(sex_deg, gencode[, .(ensembl_base, gene_name)],
                 by = "ensembl_base", all.x = TRUE)

log_msg("Loading progression_coloc_by_transition")
prog <- fread(file.path(BASE,
   "RNA-seq/results/stratified_causal/progression_coloc_by_transition.csv"))

# prog has columns: symbol, in_F0_F1, in_F1_F2, in_F2_F3, in_F3_F4 (logical),
#                   is_coloc, coloc_best_pp4, progression_class
log_msg(sprintf("Progression table: %d genes", nrow(prog)))

# Merge sex class onto progression genes
prog_sx <- merge(prog,
                 sex_deg[, .(symbol = gene_name, sex_class)],
                 by = "symbol", all.x = TRUE)
log_msg(sprintf("Genes with sex_class: %d / %d",
                sum(!is.na(prog_sx$sex_class)), nrow(prog_sx)))

# ─── 2. Per-transition 2×2 fisher ──────────────────────────────────────────────
transitions <- c("F0_F1" = "in_F0_F1",
                 "F1_F2" = "in_F1_F2",
                 "F2_F3" = "in_F2_F3",
                 "F3_F4" = "in_F3_F4")

results <- list()

# Background for non-coloc enrichment: all sex_deg genes with classification
bg_all <- sex_deg[!is.na(gene_name) & nchar(gene_name) > 0]

for (tname in names(transitions)) {
  col <- transitions[[tname]]
  trans_genes <- prog_sx[get(col) == TRUE]
  if (nrow(trans_genes) == 0) next

  # Test 1: per sex_class — is COLOC enriched within sex_class in this transition?
  for (sxc in c("Female_biased", "Male_biased", "Divergent",
                "Concordant", "Not_significant")) {
    in_class <- trans_genes[sex_class == sxc]
    if (nrow(in_class) < 5) next
    n_class    <- nrow(in_class)
    n_coloc    <- sum(in_class$is_coloc, na.rm = TRUE)
    n_other    <- sum(trans_genes$sex_class != sxc | is.na(trans_genes$sex_class),
                      na.rm = TRUE)
    n_other_coloc <- sum((trans_genes$sex_class != sxc | is.na(trans_genes$sex_class)) &
                         trans_genes$is_coloc, na.rm = TRUE)
    if (n_class < 5 || n_other < 5) next
    m <- matrix(c(n_coloc, n_class - n_coloc,
                  n_other_coloc, n_other - n_other_coloc), nrow = 2)
    if (any(rowSums(m) == 0) || any(colSums(m) == 0)) next
    ft <- fisher.test(m)
    results[[length(results) + 1L]] <- data.table(
      transition       = tname,
      sex_class        = sxc,
      n_class          = n_class,
      n_coloc_in_class = n_coloc,
      n_other          = n_other,
      n_coloc_in_other = n_other_coloc,
      odds_ratio       = unname(ft$estimate),
      ci_lower         = ft$conf.int[1],
      ci_upper         = ft$conf.int[2],
      p_value          = ft$p.value
    )
  }
}

res_dt <- rbindlist(results, fill = TRUE)
res_dt[, padj := p.adjust(p_value, "BH")]
log_msg(sprintf("Per-transition × sex_class fisher: %d tests", nrow(res_dt)))

# ─── 3. Stage-direction analysis ───────────────────────────────────────────────
# For each sex_class, count COLOC genes per transition and test
# distribution skew (early F0-F2 vs late F2-F4)
log_msg("Stage-direction skew (early vs late)")

stage_results <- list()
for (sxc in c("Female_biased", "Male_biased")) {
  in_class <- prog_sx[sex_class == sxc & is_coloc == TRUE]
  if (nrow(in_class) == 0) next
  early <- sum(in_class$in_F0_F1 | in_class$in_F1_F2, na.rm = TRUE)
  late  <- sum(in_class$in_F2_F3 | in_class$in_F3_F4, na.rm = TRUE)
  # Same totals for background (genes regardless of sex_class)
  bg_in <- prog_sx[is_coloc == TRUE]
  bg_early <- sum(bg_in$in_F0_F1 | bg_in$in_F1_F2, na.rm = TRUE)
  bg_late  <- sum(bg_in$in_F2_F3 | bg_in$in_F3_F4, na.rm = TRUE)
  m <- matrix(c(early, late,
                bg_early - early, bg_late - late), nrow = 2,
              dimnames = list(c("sex_class","background"), c("early","late")))
  if (any(m < 0) || any(rowSums(m) == 0) || any(colSums(m) == 0)) next
  ft <- fisher.test(m)
  stage_results[[sxc]] <- data.table(
    sex_class = sxc,
    n_early_class = early,
    n_late_class  = late,
    n_early_bg    = bg_early,
    n_late_bg     = bg_late,
    odds_ratio    = unname(ft$estimate),
    ci_lower      = ft$conf.int[1],
    ci_upper      = ft$conf.int[2],
    p_value       = ft$p.value,
    interpretation = ifelse(ft$estimate > 1,
                            paste0(sxc, " COLOC skews EARLY"),
                            paste0(sxc, " COLOC skews LATE"))
  )
}
stage_dt <- rbindlist(stage_results, fill = TRUE)
if (nrow(stage_dt) > 0) {
  stage_dt[, padj := p.adjust(p_value, "BH")]
}

log_msg("Stage skew results:")
print(stage_dt)

# ─── 4. Write outputs ─────────────────────────────────────────────────────────-
out_path <- file.path(OUT_DIR, "sex_progression_coloc.csv")
fwrite(res_dt, out_path)
log_msg(sprintf("Wrote %s", out_path))

stage_path <- file.path(OUT_DIR, "sex_progression_coloc_stage_skew.csv")
fwrite(stage_dt, stage_path)
log_msg(sprintf("Wrote %s", stage_path))

# Print summary
log_msg("Top per-transition signals (padj<0.1):")
top <- res_dt[padj < 0.1 | p_value < 0.01]
if (nrow(top) > 0) print(top[order(p_value),
                              .(transition, sex_class, n_class,
                                n_coloc_in_class, odds_ratio, p_value, padj)])

log_msg("218d complete")
