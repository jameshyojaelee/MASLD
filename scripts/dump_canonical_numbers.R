#!/usr/bin/env Rscript
# dump_canonical_numbers.R
# Single source of truth for the GWAS/eQTL/COLOC/atlas headline numbers.
# Recomputes everything from the on-disk canonical files so docs (NUMBERS.md,
# progress.md, paper_outline.md, MEMORY.md) can be updated from a reproducible dump
# rather than stale hand-typed values. Written 2026-05-30 (GWAS review fix campaign).
#
# Usage:  Rscript scripts/dump_canonical_numbers.R   [> docs/review/canonical_numbers_<date>.txt]
# Read-only. Safe to run any time; reflects whatever is currently on disk.

suppressMessages({ library(data.table) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

hr <- function() cat(strrep("-", 72), "\n")
section <- function(s) { cat("\n"); hr(); cat("== ", s, "\n", sep = ""); hr() }
val <- function(label, x) cat(sprintf("  %-46s %s\n", label, x))

# count distinct genes over a threshold in a numeric column, NA-safe
n_over <- function(dt, col, thr) {
  if (!col %in% names(dt)) return(NA_integer_)
  v <- suppressWarnings(as.numeric(dt[[col]]))
  sum(v > thr, na.rm = TRUE)
}

cat("CANONICAL NUMBERS DUMP  —  recomputed from on-disk files\n")
cat("Project:", BASE, "\n")
# NOTE: timestamp intentionally omitted (Sys.time avoided for reproducibility); add date when filing.

# ---------------------------------------------------------------------------
# 1. SuSiE-COLOC / ABF-COLOC  (gene_level_coloc.csv)
# ---------------------------------------------------------------------------
section("1. COLOC (gene_level_coloc.csv — canonical atlas source)")
clf <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
if (file.exists(clf)) {
  cl <- fread(clf)
  val("file", clf)
  val("eGenes tested (rows)", nrow(cl))
  if ("coloc_best_gwas" %in% names(cl))
    val("distinct GWAS contributing (coloc_best_gwas)", uniqueN(cl$coloc_best_gwas[cl$coloc_best_gwas != ""]))
  cat("  -- SuSiE (coloc_best_susie_pp4) --\n")
  val("  PP.H4.susie > 0.5", n_over(cl, "coloc_best_susie_pp4", 0.5))
  val("  PP.H4.susie > 0.8", n_over(cl, "coloc_best_susie_pp4", 0.8))
  val("  PP.H4.susie > 0.9", n_over(cl, "coloc_best_susie_pp4", 0.9))
  if ("coloc_best_susie_pp4" %in% names(cl))
    val("  genes with non-NA SuSiE PP4 (converged)", sum(!is.na(cl$coloc_best_susie_pp4)))
  cat("  -- ABF (coloc_best_pp4) --\n")
  val("  PP.H4.abf > 0.5", n_over(cl, "coloc_best_pp4", 0.5))
  val("  PP.H4.abf > 0.8", n_over(cl, "coloc_best_pp4", 0.8))
  val("  PP.H4.abf > 0.9", n_over(cl, "coloc_best_pp4", 0.9))
} else {
  val("MISSING", clf)
}

# ---------------------------------------------------------------------------
# 2. Multi-evidence atlas dimensions + causal-column presence
# ---------------------------------------------------------------------------
section("2. Multi-evidence atlas (multi_evidence_atlas.csv)")
af <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
if (file.exists(af)) {
  ahead <- names(fread(af, nrows = 0L))
  nrow_a <- nrow(fread(af, select = 1L))  # row count via first column only (fast, robust)
  val("rows (genes)", nrow_a)
  val("columns", length(ahead))
  val("ctwas_pip present? (should be FALSE — dropped)", "ctwas_pip" %in% ahead)
  val("coloc_best_pp4 + coloc_best_susie_pp4 both present?",
      all(c("coloc_best_pp4", "coloc_best_susie_pp4") %in% ahead) ||
      all(c("coloc_susie_best_pp4") %in% ahead))
  val("hyprcoloc_method present? (method-gating)", "hyprcoloc_method" %in% ahead)
} else {
  val("MISSING (atlas not yet rebuilt?)", af)
}

# ---------------------------------------------------------------------------
# 3. Convergence ranking + Tier-1 (convergence_evidence.csv)
# ---------------------------------------------------------------------------
section("3. Convergence ranking (convergence_evidence.csv)")
cvf <- file.path(BASE, "RNA-seq/results/multi_evidence/convergence_evidence.csv")
if (file.exists(cvf)) {
  cv <- fread(cvf)
  val("rows", nrow(cv))
  tier_col <- intersect(c("tier", "convergence_tier"), names(cv))[1]
  if (!is.na(tier_col)) { cat("  tier breakdown:\n"); print(cv[, .N, by = tier_col]) }
  if ("coloc_best_gwas" %in% names(cv))
    { cat("  coloc_best_gwas (genetic-channel provenance) breakdown:\n"); print(cv[, .N, by = coloc_best_gwas][order(-N)][1:min(10,.N)]) }
} else {
  val("MISSING (46d not yet re-run?)", cvf)
}

# ---------------------------------------------------------------------------
# 4. INTACT (intact_scores.csv) — multi vs ct_single
# ---------------------------------------------------------------------------
section("4. INTACT (intact_scores.csv) — RETIRED 2026-06-19, no longer feeds 46d")
inf <- file.path(BASE, "RNA-seq/results/gwas_rna_integration/intact_scores.csv")
if (file.exists(inf)) {
  ic <- fread(inf)
  val("rows", nrow(ic))
  if ("intact_score_source" %in% names(ic)) { cat("  source:\n"); print(ic[, .N, by = intact_score_source]) }
  val("multi-product intact_score_bulk > 0.5 (Tier-1 eligible)", n_over(ic, "intact_score_bulk", 0.5))
  val("single-product intact_score_ct > 0.5 (supplementary)", n_over(ic, "intact_score_ct", 0.5))
} else {
  val("MISSING", inf)
}

# ---------------------------------------------------------------------------
# 5. Positive-control recovery (sanity check)
# ---------------------------------------------------------------------------
section("5. MASLD positive-control gene recovery (COLOC)")
pc <- c("PNPLA3", "TM6SF2", "HSD17B13", "GCKR", "MARC1", "MTARC1", "MBOAT7", "TRIB1", "GPAM", "APOE")
if (file.exists(clf)) {
  cl <- fread(clf)
  gcol <- intersect(c("gene", "human_symbol"), names(cl))[1]
  cat(sprintf("  %-10s %-12s %-12s\n", "gene", "best_susie", "best_abf"))
  for (g in pc) {
    r <- cl[get(gcol) == g]
    if (nrow(r) == 0) { cat(sprintf("  %-10s %-12s %-12s\n", g, "ABSENT", "ABSENT")); next }
    su <- if ("coloc_best_susie_pp4" %in% names(r)) max(r$coloc_best_susie_pp4, na.rm = TRUE) else NA
    ab <- if ("coloc_best_pp4" %in% names(r)) max(r$coloc_best_pp4, na.rm = TRUE) else NA
    cat(sprintf("  %-10s %-12s %-12s\n", g,
                ifelse(is.finite(su), round(su, 3), "NA"),
                ifelse(is.finite(ab), round(ab, 3), "NA")))
  }
  cat("  (PNPLA3 is the strongest NAFLD locus but is an exonic/missense effect — it is\n")
  cat("   commonly NOT an eQTL-COLOC hit; absence here is expected, not a pipeline failure.)\n")
} else {
  val("cannot check — gene_level_coloc.csv missing", clf)
}

# ---------------------------------------------------------------------------
# 6. Cross-ancestry validated (if present)
# ---------------------------------------------------------------------------
section("6. Cross-ancestry validated targets")
caf <- file.path(BASE, "RNA-seq/results/causal_inference/cross_ancestry_validated_targets.csv")
if (!file.exists(caf)) caf <- file.path(BASE, "RNA-seq/results/stratified_causal/cross_ancestry_validated_targets.csv")
if (file.exists(caf)) { ca <- fread(caf); val("file", caf); val("validated genes (rows)", nrow(ca)) } else {
  val("not found (adjust path)", "cross_ancestry_validated_targets.csv")
}

cat("\n")
hr()
cat("DONE. Update NUMBERS.md / progress.md / paper_outline.md / MEMORY.md from the above.\n")
hr()
