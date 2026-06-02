#!/usr/bin/env Rscript
# 14.5d_T11_loo_cohort.R
# ---------------------------------------------------------------------------
# A7 T11: Per-cohort LOO forest plot of top-100 sex-dimorphic genes
#
# For each of 5 mega-cohorts, refit the sex-interaction dream model
# EXCLUDING that cohort. For the top-100 canonical sex-dimorphic genes,
# report interaction logFC +/- SE per LOO fit.
#
# BONUS (PART1/PART2 consistency, cold-critique driven):
# Also refit sex_dimorphic flag on the matched-dataset subset only
# (datasets with controls in BOTH sexes) and compare to canonical
# (which uses unmatched 5-cohort set for PART1 flagging).
#
# Outputs:
#   audit_sensitivity/sex_per_cohort_lfc/forest_input.csv
#   audit_sensitivity/sex_per_cohort_lfc/loo_summary.csv
#   audit_sensitivity/sex_part1_part2_consistency/matched_only_sex_dimorphic.csv
#   audit_sensitivity/sex_part1_part2_consistency/comparison_summary.csv
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas); library(lme4); library(data.table); library(edgeR); library(yaml)
})

ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}
suppressPackageStartupMessages({ library(variancePartition); library(BiocParallel) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUT  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/sex_per_cohort_lfc")
OUT2 <- file.path(BASE, "RNA-seq/results/audit_sensitivity/sex_part1_part2_consistency")
dir.create(OUT,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT2, showWarnings = FALSE, recursive = TRUE)

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("CPUs:", ncpus, "\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# --- Load ---
cat("Loading data...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("Mega cohorts:", paste(mega_cohorts, collapse = ", "), "\n")

keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]
qc_report <- fread(file.path(INT, "qc/sample_qc_report.csv"))
sex_pass_ids <- qc_report[pass_sex == TRUE, sample_id]
sex_fail <- !colnames(dge_mega) %in% sex_pass_ids
if (any(sex_fail)) dge_mega <- dge_mega[, !sex_fail]

meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

info <- data.frame(
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge_mega)
na_sex <- is.na(info$inferred_sex)
if (any(na_sex)) { dge_mega <- dge_mega[, !na_sex]; info <- info[!na_sex, , drop = FALSE] }
cat("Final samples:", nrow(info), "\n")

# --- Top-100 canonical sex-dimorphic ---
canon <- fread(file.path(RDIR, "sex_interaction_dream.csv"))
top100_canon <- canon[!is.na(padj) & padj < 0.05][order(-abs(logFC))][1:min(100, .N), gene]
cat("Top-100 canonical sex-dimorphic selected (length:", length(top100_canon), ")\n")

# --- Helper: fit interaction model on a subset and return LFC + SE per gene ---
fit_interaction <- function(dge_in, info_in, BPPARAM) {
  keep_g <- filterByExpr(dge_in, group = info_in$group_binary)
  d <- dge_in[keep_g, , keep.lib.sizes = FALSE]
  d <- calcNormFactors(d, method = "TMM")
  info_in$dataset <- droplevels(info_in$dataset)
  form <- ~ group_binary * inferred_sex + (1 | dataset)
  v <- suppressWarnings(voomWithDreamWeights(d, form, info_in, BPPARAM = BPPARAM))
  fit <- suppressWarnings(dream(v, form, info_in, BPPARAM = BPPARAM))
  cf  <- colnames(fit$coefficients)
  int_coef <- grep("group_binary.*inferred_sex|inferred_sex.*group_binary", cf, value = TRUE)[1]
  # logFC + SE: dream returns t-stat => SE = logFC / t for each gene
  tt <- topTable(fit, coef = int_coef, number = Inf, sort.by = "none")
  dt <- as.data.table(tt)
  dt[, gene := rownames(tt)]
  setnames(dt, "adj.P.Val", "padj", skip_absent = TRUE)
  dt[, SE := ifelse(t == 0, NA_real_, logFC / t)]
  dt
}

# ============================================================================
# PART A: LOO per cohort (top-100 forest)
# ============================================================================
cat("\n===== PART A: Per-cohort LOO interaction fits =====\n")
forest <- list()
loo_summary <- list()
for (excl in mega_cohorts) {
  cat("\n--- LOO excluding:", excl, "---\n")
  keep_loo <- info$dataset != excl
  info_loo <- info[keep_loo, , drop = FALSE]
  dge_loo  <- dge_mega[, keep_loo]
  cat("  Samples remaining:", nrow(info_loo), "\n")
  cat("  group x sex:\n"); print(table(info_loo$group_binary, info_loo$inferred_sex))

  res_loo <- tryCatch(fit_interaction(dge_loo, info_loo, param),
                      error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
  if (is.null(res_loo)) next

  sub <- res_loo[gene %in% top100_canon, .(gene, logFC, SE, t, padj)]
  sub[, cohort_excluded := excl]
  forest[[excl]] <- sub

  # Summary stats vs canonical
  canon_sub <- canon[gene %in% top100_canon, .(gene, logFC_canon = logFC, padj_canon = padj)]
  mrg <- merge(sub, canon_sub, by = "gene")
  rho <- cor(mrg$logFC, mrg$logFC_canon, method = "spearman", use = "pairwise.complete.obs")
  sign_consistent <- sum(sign(mrg$logFC) == sign(mrg$logFC_canon), na.rm = TRUE)
  loo_summary[[excl]] <- data.table(
    cohort_excluded = excl,
    n_samples_remaining = nrow(info_loo),
    n_top100_recovered = nrow(sub),
    spearman_rho_to_canon = rho,
    n_sign_consistent = sign_consistent,
    n_sex_dimorphic_padj05 = sum(res_loo$padj < 0.05, na.rm = TRUE)
  )
  cat("  rho-to-canonical:", round(rho, 3),
      " sign-consistent:", sign_consistent, "/", nrow(mrg),
      " n_dim(padj<0.05):", loo_summary[[excl]]$n_sex_dimorphic_padj05, "\n")
}

forest_dt <- rbindlist(forest, use.names = TRUE, fill = TRUE)
fwrite(forest_dt, file.path(OUT, "forest_input.csv"))
cat("\nSaved forest_input.csv (rows:", nrow(forest_dt), ")\n")

loo_dt <- rbindlist(loo_summary, use.names = TRUE, fill = TRUE)
fwrite(loo_dt, file.path(OUT, "loo_summary.csv"))
cat("Saved loo_summary.csv\n"); print(loo_dt)

# ============================================================================
# PART B: PART1/PART2 consistency — matched-only sex_dimorphic flag
# ============================================================================
cat("\n===== PART B: Matched-dataset PART1 (cold-critique) =====\n")
sex_levels <- levels(info$inferred_sex)
male_label   <- sex_levels[grepl("^M", sex_levels)][1]
female_label <- sex_levels[grepl("^F", sex_levels)][1]
ds_m_ctrl <- unique(as.character(info$dataset[info$group_binary == "Control" & info$inferred_sex == male_label]))
ds_f_ctrl <- unique(as.character(info$dataset[info$group_binary == "Control" & info$inferred_sex == female_label]))
matched_ds <- intersect(ds_m_ctrl, ds_f_ctrl)
cat("Matched datasets (controls in both sexes):", paste(matched_ds, collapse = ", "), "\n")

keep_matched <- info$dataset %in% matched_ds
info_matched <- info[keep_matched, , drop = FALSE]
dge_matched  <- dge_mega[, keep_matched]
cat("Samples after matching:", nrow(info_matched), "\n")
cat("group x sex (matched):\n"); print(table(info_matched$group_binary, info_matched$inferred_sex))

res_matched <- tryCatch(fit_interaction(dge_matched, info_matched, param),
                        error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
if (!is.null(res_matched)) {
  fwrite(res_matched, file.path(OUT2, "matched_only_sex_dimorphic.csv"))
  cat("Saved matched_only_sex_dimorphic.csv\n")

  # Comparison summary
  canon2 <- canon[, .(gene, logFC_canon = logFC, padj_canon = padj)]
  mrg2   <- merge(res_matched[, .(gene, logFC_matched = logFC, padj_matched = padj)],
                  canon2, by = "gene")
  rho_all  <- cor(mrg2$logFC_matched, mrg2$logFC_canon, method = "spearman", use = "pairwise.complete.obs")
  n_canon  <- sum(canon$padj < 0.05, na.rm = TRUE)
  n_match  <- sum(res_matched$padj < 0.05, na.rm = TRUE)
  overlap  <- length(intersect(
    canon[!is.na(padj) & padj < 0.05, gene],
    res_matched[!is.na(padj) & padj < 0.05, gene]
  ))
  comp_sum <- data.table(
    metric = c("spearman_rho_all_genes", "n_sex_dimorphic_canonical_padj05",
               "n_sex_dimorphic_matched_padj05", "overlap_padj05",
               "jaccard_padj05"),
    value  = c(rho_all, n_canon, n_match, overlap,
               overlap / max(length(union(
                 canon[!is.na(padj) & padj < 0.05, gene],
                 res_matched[!is.na(padj) & padj < 0.05, gene])), 1L))
  )
  fwrite(comp_sum, file.path(OUT2, "comparison_summary.csv"))
  cat("\nMatched vs canonical:\n"); print(comp_sum)
} else {
  cat("Matched fit failed.\n")
}

cat("\nDone:", as.character(Sys.time()), "\n")
