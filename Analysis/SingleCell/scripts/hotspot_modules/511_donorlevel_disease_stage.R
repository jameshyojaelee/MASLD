#!/usr/bin/env Rscript
# ============================================================================
# Donor-level recompute of the Hotspot module disease-stage statistic.
#
# FIXES a run-vs-donor pseudoreplication bug. Script 505 fits
#   lmer(score ~ disease_stage_ordinal + (1|dataset))
# on RUN-level per-sample module scores (donor_scores_all.tsv is keyed on the
# atlas `sample` = a sequencing run / sort-fraction library, NOT a biological
# donor). Four datasets carry multiple runs per donor
#   GSE244832 117 runs -> 18 donors, GSE202379 67 -> 46, GSE185477 21 -> 3,
#   GSE136103 20 GSMs -> 10 donors.
# In the actual disease-stage fit (after the protocol-contamination filter) this
# is 211 runs standing in for only 73 true donors -> n inflated ~2.9x, SE biased
# down, significance inflated. `(1|dataset)` does NOT absorb within-donor
# correlation (two runs of one donor are more correlated than two donors of one
# dataset), so the bias is real.
#
# THE FIX: collapse per-run module scores to ONE score per true biological donor
# BEFORE the regression, keeping (1|dataset). Module DEFINITIONS (Hotspot
# co-expression membership) are a co-expression structure unaffected by
# pseudoreplication and are NOT recomputed.
#   CANONICAL  = unweighted mean of a donor's per-run module scores.
#   SENSITIVITY= pooled-cell (cell-count-weighted) donor mean from
#                cell_scores_all.parquet (via 511_pooled_cell_donor_scores.py) =
#                the exact donor-level statistic. Reported to confirm robustness.
#
# Outputs (NEW paths; canonical run-level files untouched):
#   donor_collapse/phenotype_correlations_donor.tsv           (unweighted)
#   donor_collapse/phenotype_correlations_donor_weighted.tsv  (pooled-cell, if avail)
#   donor_collapse/all_modules_donor.tsv                      (unweighted, drives Fig3G)
#   donor_collapse/disease_stage_before_after.tsv             (per-CT sig counts)
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(lme4); library(lmerTest)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
RES <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
OUT <- file.path(RES, "donor_collapse"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CTS_FIG <- c("hepatocytes", "fibroblasts", "macrophages", "cholangiocytes", "tcells")  # Fig3G cell types

srr_to_donor <- build_srr_to_donor_map(BASE)

# ---------------------------------------------------------------------------
# (1) donor-level metadata (one row per biological donor)
# ---------------------------------------------------------------------------
donor_meta <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"))
donor_meta <- unique(donor_meta, by = "sample")
if (!"exclude_stage_analysis" %in% names(donor_meta)) donor_meta[, exclude_stage_analysis := FALSE]
if (!"F_stage_inferred" %in% names(donor_meta)) donor_meta[, F_stage_inferred := NA_real_]
if (!"F_stage_source" %in% names(donor_meta)) donor_meta[, F_stage_source := NA_character_]
if (!"nas_score" %in% names(donor_meta)) donor_meta[, nas_score := NA_real_]
donor_meta[, disease_stage_ordinal := fcase(
  disease_stage_coarse == "Healthy", 0, disease_stage_coarse == "Steatosis", 1,
  disease_stage_coarse == "Steatohepatitis", 2, disease_stage_coarse == "Cirrhosis", 3,
  default = NA_real_)]
donor_meta[, donor := ifelse(sample %in% names(srr_to_donor), srr_to_donor[sample], sample)]

# stage / dataset / exclude are invariant within a donor (verified: 0 disagreements);
# F_stage_inferred and nas_score are averaged across the donor's runs.
mode_chr <- function(v) { t <- table(v); if (length(t) == 0L) NA_character_ else names(sort(t, decreasing = TRUE))[1] }
meta_donor <- donor_meta[, .(
  disease_stage_coarse   = disease_stage_coarse[1],
  disease_stage_ordinal  = disease_stage_ordinal[1],
  dataset                = dataset[1],
  exclude_stage_analysis = as.logical(any(exclude_stage_analysis %in% TRUE)),
  F_stage_inferred = if (all(is.na(F_stage_inferred))) NA_real_ else mean(F_stage_inferred, na.rm = TRUE),
  F_stage_source   = mode_chr(F_stage_source),
  nas_score        = if (all(is.na(nas_score))) NA_real_ else mean(nas_score, na.rm = TRUE)
), by = donor]

# ---------------------------------------------------------------------------
# (2) fit engine (matches 505 fit_one: >=20 non-NA on the axis, >=2 datasets;
#     lmerTest::lmer(score ~ axis + (1|dataset), REML=FALSE))
# ---------------------------------------------------------------------------
fit_one <- function(d, axis_col) {
  v <- d[[axis_col]]
  if (sum(!is.na(v)) < 20 || uniqueN(d$dataset) < 2)
    return(data.table(beta = NA_real_, SE = NA_real_, t = NA_real_, p = NA_real_, n = sum(!is.na(v))))
  fm <- as.formula(sprintf("score ~ %s + (1|dataset)", axis_col))
  m <- tryCatch(lmerTest::lmer(fm, data = d, REML = FALSE), error = function(e) NULL)
  if (is.null(m)) return(data.table(beta = NA_real_, SE = NA_real_, t = NA_real_, p = NA_real_, n = nrow(d)))
  cf <- tryCatch(summary(m)$coefficients, error = function(e) NULL)
  if (is.null(cf) || !(axis_col %in% rownames(cf)))
    return(data.table(beta = NA_real_, SE = NA_real_, t = NA_real_, p = NA_real_, n = nrow(d)))
  row <- cf[axis_col, ]
  p_val <- if ("Pr(>|t|)" %in% colnames(cf)) row["Pr(>|t|)"] else NA_real_
  data.table(beta = unname(row["Estimate"]), SE = unname(row["Std. Error"]),
             t = unname(row["t value"]), p = unname(p_val), n = nrow(d))
}

AXES <- list(
  disease_stage           = list(col = "disease_stage_ordinal", filt = function(d) d),
  F_stage                 = list(col = "F_stage_inferred",      filt = function(d) d),
  F_stage_documented_only = list(col = "F_stage_inferred",      filt = function(d) d[F_stage_source == "documented"]),
  NAS                     = list(col = "nas_score",             filt = function(d) d)
)

compute_pheno <- function(scores_donor, label) {
  # scores_donor: cell_type, module, donor, score
  joined <- merge(scores_donor, meta_donor, by = "donor")
  joined <- joined[!is.na(score)]
  n_pre <- uniqueN(joined$donor)
  joined <- joined[is.na(exclude_stage_analysis) | exclude_stage_analysis == FALSE]
  n_post <- uniqueN(joined$donor)
  cat(sprintf("[%s] protocol-contamination filter: %d donors retained (was %d)\n", label, n_post, n_pre))
  res <- rbindlist(lapply(names(AXES), function(ax) {
    cfg <- AXES[[ax]]; d_ax <- cfg$filt(joined)
    per <- d_ax[, fit_one(.SD, cfg$col), by = .(cell_type, module)]
    per[, axis := ax][, q := p.adjust(p, method = "BH")][]
  }))
  res
}

# ---------------------------------------------------------------------------
# (3) CANONICAL unweighted donor collapse
# ---------------------------------------------------------------------------
scores <- fread(file.path(RES, "donor_scores_all.tsv"))          # sample, module, score, cell_type
scores[, donor := ifelse(sample %in% names(srr_to_donor), srr_to_donor[sample], sample)]
scores_donor <- scores[, .(score = mean(score, na.rm = TRUE)), by = .(cell_type, module, donor)]
cat(sprintf("[unweighted] run-level rows %d -> donor-level rows %d\n", nrow(scores), nrow(scores_donor)))
pheno_uw <- compute_pheno(scores_donor, "unweighted")
fwrite(pheno_uw, file.path(OUT, "phenotype_correlations_donor.tsv"), sep = "\t")

# ---------------------------------------------------------------------------
# (4) SENSITIVITY pooled-cell (weighted) donor collapse, if available
# ---------------------------------------------------------------------------
wfile <- file.path(OUT, "donor_scores_all_weighted.tsv")
pheno_w <- NULL
if (file.exists(wfile)) {
  sw <- fread(wfile)                                             # sample(=donor), module, score, cell_type
  setnames(sw, "sample", "donor")
  pheno_w <- compute_pheno(sw[, .(cell_type, module, donor, score)], "weighted")
  fwrite(pheno_w, file.path(OUT, "phenotype_correlations_donor_weighted.tsv"), sep = "\t")
} else {
  cat("[weighted] donor_scores_all_weighted.tsv not found - skipping pooled-cell sensitivity arm\n")
}

# ---------------------------------------------------------------------------
# (5) build all_modules_donor.tsv (508-style merge with donor-level pheno)
# ---------------------------------------------------------------------------
am <- fread(file.path(RES, "all_modules.tsv"))
pheno_cols <- c("disease_stage_beta", "F_stage_beta", "F_stage_documented_only_beta", "NAS_beta",
                "disease_stage_q", "F_stage_q", "F_stage_documented_only_q", "NAS_q",
                "progression_module", "activity_module", "bulk_replicated")
am_base <- am[, setdiff(names(am), pheno_cols), with = FALSE]     # keep novelty/loo/bulk_module_score/names etc

pw <- copy(pheno_uw)
pw[, module := as.integer(sub("^.*__", "", module))]             # strip "<celltype>__" namespace (matches 508)
pheno_wide <- dcast(pw[axis %in% c("disease_stage", "F_stage", "F_stage_documented_only", "NAS")],
                    cell_type + module ~ axis, value.var = c("beta", "q"))
setnames(pheno_wide,
  c("beta_disease_stage", "beta_F_stage", "beta_F_stage_documented_only", "beta_NAS",
    "q_disease_stage", "q_F_stage", "q_F_stage_documented_only", "q_NAS"),
  c("disease_stage_beta", "F_stage_beta", "F_stage_documented_only_beta", "NAS_beta",
    "disease_stage_q", "F_stage_q", "F_stage_documented_only_q", "NAS_q"))
pheno_wide[, progression_module := (!is.na(disease_stage_q) & disease_stage_q < 0.05 &
                                    !is.na(F_stage_q) & F_stage_q < 0.05 &
                                    sign(disease_stage_beta) == sign(F_stage_beta))]
pheno_wide[, activity_module := (!is.na(disease_stage_q) & disease_stage_q < 0.05 &
                                 (!is.na(F_stage_q) & F_stage_q >= 0.10))]

all_modules_donor <- merge(am_base, pheno_wide, by = c("cell_type", "module"), all.x = TRUE)
all_modules_donor[, bulk_replicated := !is.na(bulk_module_score) &
                    sign(bulk_module_score) == sign(disease_stage_beta)]
setcolorder(all_modules_donor, intersect(names(am), names(all_modules_donor)))
fwrite(all_modules_donor, file.path(OUT, "all_modules_donor.tsv"), sep = "\t")
cat(sprintf("Wrote all_modules_donor.tsv (%d rows)\n", nrow(all_modules_donor)))

# ---------------------------------------------------------------------------
# (6) BEFORE (run-level) vs AFTER (donor-level) disease-significant module counts
# ---------------------------------------------------------------------------
CTS_ALL <- c("global", "hepatocytes", "macrophages", "fibroblasts",
             "endothelial_cells", "cholangiocytes", "tcells")
tot <- function(dt) dt[!is.na(disease_stage_q) & disease_stage_q < 0.05, .N]
sig_by_ct <- function(dt) {
  x <- dt[!is.na(disease_stage_q) & disease_stage_q < 0.05, .(n_sig = .N), by = cell_type]
  merge(data.table(cell_type = CTS_ALL), x, by = "cell_type", all.x = TRUE)[is.na(n_sig), n_sig := 0][]
}
before <- sig_by_ct(am);               setnames(before, "n_sig", "n_sig_run_level")
after  <- sig_by_ct(all_modules_donor); setnames(after,  "n_sig", "n_sig_donor_level")
ntested <- merge(
  am[!is.na(disease_stage_q), .(n_tested = .N), by = cell_type],
  data.table(cell_type = CTS_ALL), by = "cell_type", all.y = TRUE)
cmp <- Reduce(function(a, b) merge(a, b, by = "cell_type", all = TRUE), list(ntested, before, after))
cmp[is.na(n_tested), n_tested := 0][is.na(n_sig_run_level), n_sig_run_level := 0][is.na(n_sig_donor_level), n_sig_donor_level := 0]
cmp[, cell_type := factor(cell_type, levels = CTS_ALL)]; setorder(cmp, cell_type)
if (!is.null(pheno_w)) {
  pw2 <- copy(pheno_w); pw2 <- pw2[axis == "disease_stage"]
  pw2[, module := as.integer(sub("^.*__", "", module))]
  wct <- pw2[!is.na(q) & q < 0.05, .(n_sig_donor_weighted = .N), by = cell_type]
  cmp <- merge(cmp, wct, by = "cell_type", all.x = TRUE)[is.na(n_sig_donor_weighted), n_sig_donor_weighted := 0]
}
fwrite(cmp, file.path(OUT, "disease_stage_before_after.tsv"), sep = "\t")

cat("\n=================  disease-significant (disease_stage_q < 0.05) modules  =================\n")
print(cmp)
figCTs <- cmp[cell_type %in% CTS_FIG]
cat(sprintf("\nFig3G 5 cell types: BEFORE (run-level) = %d of %d ; AFTER (donor, unweighted) = %d of %d",
            sum(figCTs$n_sig_run_level), sum(figCTs$n_tested),
            sum(figCTs$n_sig_donor_level), sum(figCTs$n_tested)))
if (!is.null(pheno_w)) cat(sprintf(" ; AFTER (donor, pooled-cell) = %d", sum(figCTs$n_sig_donor_weighted)))
cat("\n")
