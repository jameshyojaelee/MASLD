#!/usr/bin/env Rscript
# 04c_disease_master_regulators.R
# =============================================================================
# Build a DEFENSIBLE set of disease-associated hepatic master regulators by
# CROSS-MODALITY convergence, replacing the FDR-gated single-cohort regulon set.
#
# WHY THIS SCRIPT EXISTS
# ----------------------
# The donor-level (n=18) SCENIC+ regulon-activity DE test honestly returns ZERO
# regulons at FDR < 0.05 -> `scenic_plus/disease_regulons.csv` is empty (header
# only). The historical "12/24 disrupted TFs" figure came from a CELL-LEVEL test
# (pseudoreplication: cells from the same donor treated as independent), which
# inflates significance and must NOT be restored.
#
# The PI decision is to SEPARATE the claims:
#   * Define disease-associated hepatic master regulators by WELL-POWERED,
#     cross-modality disease evidence (bulk MASLD DE, n=846; or genetic COLOC).
#   * Report the single-cohort regulon-activity DE HONESTLY as underpowered
#     (n=18, n.s. at FDR) -- carried only as a supporting annotation, NOT a gate.
#
# DEFINITION (non-circular w.r.t. the motif-disruption layer in Script 56)
# ------------------------------------------------------------------------
#   A "disease master regulator" is a hepatocyte SCENIC+ regulon TF (i.e. a TF
#   that forms a hepatocyte GRN module) that is ALSO:
#       (bulk MASLD DEG)  OR  (genetic COLOC hit).
#   Neither evidence layer touches the GWAS-variant motif-disruption test that
#   Script 56 runs, so cross-referencing the two is not circular.
#
# Externally validated positives expected in the output (sanity check):
#   HNF4A, RORA, THRB (THRB = resmetirom target; RORA/THRB COLOC PP4 ~ 1).
#
# ENV: rnaseq  |  Run: micromamba run -n rnaseq Rscript 04c_disease_master_regulators.R
# =============================================================================

suppressMessages({
  library(data.table)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATAC_DIR <- file.path(BASE_DIR, "Analysis/ATAC/Human_Multiome")
SCENIC_DIR <- file.path(ATAC_DIR, "scenic_plus")

HEP_REGULON_FILE <- file.path(SCENIC_DIR, "hepatocyte_regulons.csv")
DREAM_ASHR_FILE  <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
DREAM_FILE       <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
COLOC_FILE       <- file.path(BASE_DIR,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
OUT_FILE         <- file.path(SCENIC_DIR, "disease_master_regulators.csv")

# Well-powered evidence thresholds
LFSR_CUT  <- 0.05   # ashr local false sign rate (or padj fallback)
PADJ_CUT  <- 0.05
COLOC_CUT <- 0.5    # canonical PP.H4 threshold (SuSiE or ABF)

cat("============================================================\n")
cat("04c_disease_master_regulators.R\n")
cat("Cross-modality disease master regulators (bulk DE OR COLOC)\n")
cat("============================================================\n\n")

# ── 1. Candidate universe: hepatocyte SCENIC+ regulon TFs ────────────────────
stopifnot(file.exists(HEP_REGULON_FILE))
hep <- fread(HEP_REGULON_FILE)
cat("Loaded", nrow(hep), "hepatocyte regulon rows from", basename(HEP_REGULON_FILE), "\n")

# Collapse to one row per TF. regulon_activity_diff / activity_padj are constant
# within a TF (donor-level test is per-regulon), so take the first.
reg <- hep[, .(
  regulon_activity_diff = regulon_activity_diff[1],
  activity_padj         = activity_padj[1],
  n_target_genes        = if ("n_target_genes" %in% names(hep)) n_target_genes[1] else NA_integer_,
  n_donors              = if ("n_donors" %in% names(hep)) n_donors[1] else NA_integer_
), by = .(tf_name)]
reg[, tf_upper := toupper(tf_name)]
cat("Candidate hepatocyte regulon TFs (universe):", nrow(reg), "\n\n")

# ── 2a. Bulk MASLD differential expression (n=846, well-powered) ─────────────
if (file.exists(DREAM_ASHR_FILE)) {
  dream <- fread(DREAM_ASHR_FILE)
  dream_source <- "dream_results_ashr.csv"
  dream[, dream_logFC := if ("shrunk_logFC" %in% names(dream)) shrunk_logFC else logFC]
  dream[, dream_sig_stat := if ("lfsr" %in% names(dream)) lfsr else padj]
  dream_stat_name <- if ("lfsr" %in% names(dream)) "lfsr" else "padj"
  dream[, dream_symbol := if ("symbol" %in% names(dream)) symbol else gene]
} else if (file.exists(DREAM_FILE)) {
  dream <- fread(DREAM_FILE)
  dream_source <- "dream_results.csv"
  dream[, dream_logFC := logFC]
  dream[, dream_sig_stat := padj]
  dream_stat_name <- "padj"
  # dream_results.csv keys on gene; symbol may live in `gene` already
  dream[, dream_symbol := if ("symbol" %in% names(dream)) symbol else gene]
} else {
  stop("No dream results file found (ashr or plain).")
}
cat("Bulk DE source:", dream_source, "(significance stat:", dream_stat_name, ")\n")

dream[, dream_symbol_upper := toupper(dream_symbol)]
# One row per symbol: keep most significant
dream_tf <- dream[order(dream_sig_stat),
                  .SD[1], by = dream_symbol_upper,
                  .SDcols = c("dream_logFC", "dream_sig_stat")]
setnames(dream_tf, "dream_sig_stat", "dream_lfsr")  # column name kept generic

reg <- merge(reg, dream_tf, by.x = "tf_upper", by.y = "dream_symbol_upper", all.x = TRUE)
reg[, is_bulk_deg := !is.na(dream_lfsr) & dream_lfsr < ifelse(dream_stat_name == "lfsr", LFSR_CUT, PADJ_CUT)]
cat("  bulk DEG TFs (", dream_stat_name, " < ",
    ifelse(dream_stat_name == "lfsr", LFSR_CUT, PADJ_CUT), "): ",
    sum(reg$is_bulk_deg, na.rm = TRUE), "\n", sep = "")

# ── 2b. Genetic causal colocalization (well-powered) ─────────────────────────
if (file.exists(COLOC_FILE)) {
  coloc <- fread(COLOC_FILE)
  coloc[, gene_upper := toupper(gene)]
  has_susie <- "coloc_best_susie_pp4" %in% names(coloc)
  has_abf   <- "coloc_best_pp4" %in% names(coloc)
  coloc_keep <- c("gene_upper",
                  if (has_susie) "coloc_best_susie_pp4",
                  if (has_abf)   "coloc_best_pp4")
  coloc_tf <- unique(coloc[, ..coloc_keep], by = "gene_upper")
  reg <- merge(reg, coloc_tf, by.x = "tf_upper", by.y = "gene_upper", all.x = TRUE)
} else {
  cat("WARNING: coloc file not found — is_coloc set FALSE for all\n")
  reg[, coloc_best_susie_pp4 := NA_real_]
  reg[, coloc_best_pp4 := NA_real_]
}
if (!"coloc_best_susie_pp4" %in% names(reg)) reg[, coloc_best_susie_pp4 := NA_real_]
if (!"coloc_best_pp4"       %in% names(reg)) reg[, coloc_best_pp4 := NA_real_]

reg[, coloc_pp4 := pmax(coloc_best_susie_pp4, coloc_best_pp4, na.rm = TRUE)]
reg[is.infinite(coloc_pp4), coloc_pp4 := NA_real_]
reg[, is_coloc := (!is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 > COLOC_CUT) |
                  (!is.na(coloc_best_pp4)       & coloc_best_pp4       > COLOC_CUT)]
reg[is.na(is_coloc), is_coloc := FALSE]
cat("  COLOC TFs (PP.H4 SuSiE or ABF > ", COLOC_CUT, "): ",
    sum(reg$is_coloc, na.rm = TRUE), "\n\n", sep = "")

# ── 3. Define the master-regulator set: bulk DEG OR COLOC ────────────────────
reg[is.na(is_bulk_deg), is_bulk_deg := FALSE]
reg[, n_evidence := as.integer(is_bulk_deg) + as.integer(is_coloc)]
reg[, is_master_regulator := is_bulk_deg | is_coloc]

# Honest single-cohort annotation (NOT a gate)
reg[, sc_regulon_dir := sign(regulon_activity_diff)]
reg[, sc_underpowered := TRUE]   # n=18 donor-level test; n.s. at FDR

master <- reg[is_master_regulator == TRUE]
setorder(master, -n_evidence, dream_lfsr, -coloc_pp4)

DEFINITION <- paste0(
  "cross-modality: hepatocyte SCENIC+ regulon TF that is a bulk MASLD DEG ",
  "or COLOC hit; single-cohort regulon-activity DE underpowered (n=18)")

out <- master[, .(
  tf_name,
  regulon_activity_diff,
  activity_padj,
  is_bulk_deg,
  dream_logFC,
  dream_lfsr,
  is_coloc,
  coloc_pp4,
  n_evidence,
  sc_regulon_dir,
  sc_underpowered,
  definition = DEFINITION
)]

fwrite(out, OUT_FILE)

cat("============================================================\n")
cat("Disease master regulators (cross-modality definition)\n")
cat("============================================================\n")
cat("Candidate hepatocyte regulon TFs : ", nrow(reg), "\n")
cat("  bulk DEG                       : ", sum(reg$is_bulk_deg), "\n")
cat("  COLOC                          : ", sum(reg$is_coloc), "\n")
cat("  master regulators (DEG OR COLOC): ", nrow(out), "\n")
cat("Wrote: ", OUT_FILE, "\n\n")

cat("Full master-regulator TF set (n =", nrow(out), "):\n")
print(out[, .(tf_name, is_bulk_deg, dream_logFC = round(dream_logFC, 3),
              dream_lfsr = signif(dream_lfsr, 3), is_coloc,
              coloc_pp4 = round(coloc_pp4, 3), n_evidence, sc_regulon_dir)],
      nrows = nrow(out))

cat("\n--- Sanity check: externally validated positives ---\n")
for (tf in c("HNF4A", "RORA", "THRB")) {
  r <- out[toupper(tf_name) == tf]
  if (nrow(r) == 0) {
    cat(sprintf("  %-6s: ABSENT from master set (UNEXPECTED)\n", tf))
  } else {
    cat(sprintf("  %-6s: PRESENT | is_bulk_deg=%s dream_logFC=%.3f dream_lfsr=%.2e | is_coloc=%s coloc_pp4=%.3f | n_evidence=%d\n",
                tf, r$is_bulk_deg, r$dream_logFC, r$dream_lfsr, r$is_coloc, r$coloc_pp4, r$n_evidence))
  }
}
cat("\nDone.\n")
