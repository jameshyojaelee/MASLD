#!/usr/bin/env Rscript
# 31_spatial_immune_exclusion.R
#
# Analysis I2 (v2) — Spatial immune exclusion zones.
#
# Strategy:
#   1. Load cell2location PER-SPOT cell-type abundances (~6,546 spots).
#   2. Define immune abundance = T cells + NK + cDC1s + cDC2s + B cells.
#   3. Define stromal/fibrotic abundance = Fibroblasts + Macrophages + Endothelial.
#   4. Classify spots: immune-excluded = low immune, high stromal.
#   5. Test: are repulsive ligands (IDO1, TGFB1, VEGFA, CXCL12) elevated in
#      excluded vs non-excluded spots? (donor-aware, BH over found ligands)
#
# 2026-06-01 rigor pass:
#   * F120 — input is now the TRUE per-spot abundance matrix (ontrac per-spot
#     composition, 6,546 spots) instead of the donor-mean
#     spatial_cell_type_proportions.csv (5 rows). Spots are classified within
#     donors; donor-level summaries treat the ~5 donors as the unit.
#   * F121 — the repulsive-ligand elevation test (docstring step 5) is now
#     IMPLEMENTED on log1p-CPM spot expression (extracted by
#     31a_extract_repulsive_ligand_expr.py via spatial_stats.ensure_lognorm),
#     donor-aware, BH-corrected over the ligands actually present.
#   * F122 — pct_immune_excluded uses a branch-agnostic `excluded_label`, so it
#     is no longer identically 0 in the small-n branch.
#   * F123 — cell-type co-occurrence is computed on the per-spot matrix (per
#     donor, then aggregated) with a documented donor-level caveat, NOT on n=5
#     donor means.
#
# Outputs: Analysis/Spatial/results/immune_exclusion/
# Env: rnaseq
# Prereq (for step 5): spatial-env extractor
#   micromamba run -n spatial python 31a_extract_repulsive_ligand_expr.py

suppressPackageStartupMessages({library(data.table)})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SCRIPTS <- file.path(BASE, "Analysis/Spatial/scripts")
# F120: per-spot composition + metadata (one row per spot, aligned 1:1 to the
# deconvolved h5ad), NOT the donor-mean spatial_cell_type_proportions.csv.
COMP_CSV <- file.path(BASE, "Analysis/Spatial/results/ontrac/input/ontrac_cell_type_composition.csv")
META_CSV <- file.path(BASE, "Analysis/Spatial/results/ontrac/input/ontrac_metadata.csv")
LIG_CSV  <- file.path(BASE, "Analysis/Spatial/results/immune_exclusion/repulsive_ligand_spot_expr.csv")
OUTDIR <- file.path(BASE, "Analysis/Spatial/results/immune_exclusion")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[1] Loading PER-SPOT cell-type abundances...")
prop <- fread(COMP_CSV)                 # first col = spot barcode (unnamed -> V1)
setnames(prop, 1, "spot")
meta <- fread(META_CSV)                 # spot, Cell_ID, Sample, Cell_Type, x, y, condition
setnames(meta, 1, "spot")
meta <- meta[, .(spot, sample_id = Sample, condition)]
prop <- merge(prop, meta, by = "spot")
message(sprintf("  %d spots across %d donors: %s", nrow(prop), uniqueN(prop$sample_id),
                paste(sort(unique(prop$sample_id)), collapse = ", ")))

immune_cts <- c("T cells","Circulating NK/NKT","Resident NK","cDC1s","cDC2s","B cells","Plasma cells")
stromal_cts <- c("Fibroblasts","Macrophages","Endothelial cells")

immune_cols <- intersect(immune_cts, names(prop))
stromal_cols <- intersect(stromal_cts, names(prop))
stopifnot(length(immune_cols) > 0, length(stromal_cols) > 0)

message(sprintf("  Immune columns: %s", paste(immune_cols, collapse = ", ")))
message(sprintf("  Stromal columns: %s", paste(stromal_cols, collapse = ", ")))

prop[, immune_score := rowSums(.SD, na.rm = TRUE), .SDcols = immune_cols]
prop[, stromal_score := rowSums(.SD, na.rm = TRUE), .SDcols = stromal_cols]
prop[, hep_score := Hepatocytes]

# Classify SPOTS (per-spot, F120). Quantile split over all spots; with ~6.5k
# spots the quantile branch is the one taken.
n_rows <- nrow(prop)
if (n_rows < 20) {
  excluded_label <- "low_immune_high_stromal"     # F122: define per branch
  prop[, spot_class := fifelse(immune_score <= median(immune_score, na.rm = TRUE) &
                                stromal_score >= median(stromal_score, na.rm = TRUE),
                                excluded_label,
                                fifelse(immune_score > median(immune_score, na.rm = TRUE) &
                                        stromal_score < median(stromal_score, na.rm = TRUE),
                                        "high_immune_low_stromal", "intermediate"))]
} else {
  excluded_label <- "immune_excluded_fibrotic"     # F122: define per branch
  immune_q25 <- quantile(prop$immune_score, 0.25, na.rm = TRUE)
  stromal_q75 <- quantile(prop$stromal_score, 0.75, na.rm = TRUE)
  prop[, spot_class := fifelse(immune_score <= immune_q25 & stromal_score >= stromal_q75,
                                excluded_label,
                                fifelse(immune_score > immune_q25 & stromal_score < stromal_q75,
                                        "non_excluded", "intermediate"))]
}

message(sprintf("[2] Spot classes: %s",
                paste(sprintf("%s=%d", names(table(prop$spot_class)), table(prop$spot_class)),
                      collapse = ", ")))

message("[3] Per-donor immune-exclusion fraction (donors are the unit, n~5)...")
# F122: branch-agnostic excluded_label so pct_immune_excluded is not always 0.
sample_summary <- prop[, .(n_spots = .N,
                            pct_immune_excluded = mean(spot_class == excluded_label),
                            mean_immune = mean(immune_score, na.rm = TRUE),
                            mean_stromal = mean(stromal_score, na.rm = TRUE),
                            mean_hep = mean(hep_score, na.rm = TRUE)),
                        by = sample_id]
fwrite(sample_summary, file.path(OUTDIR, "per_sample_exclusion_summary.csv"))

# -------------------------------------------------------------------------- #
# [4] Cell-type spatial co-occurrence (F123): per-spot, per-donor, aggregated.
# -------------------------------------------------------------------------- #
# Compute Spearman co-occurrence WITHIN each donor over its thousands of spots
# (so within-donor spatial structure drives it, not n=5 donor means), then
# average rho across donors. CAVEAT: only ~5 donors, and spots within a slice
# are spatially autocorrelated, so the across-donor mean is descriptive — no
# across-donor significance is claimed.
ct_all <- intersect(c(immune_cts, stromal_cts, "Hepatocytes","Cholangiocytes"), names(prop))
donors <- sort(unique(prop$sample_id))
cor_list <- lapply(donors, function(d) {
  sub <- prop[sample_id == d, ..ct_all]
  # need variance in each column within the donor
  cor(sub, use = "pairwise.complete.obs", method = "spearman")
})
# Average per-donor correlation matrices (NaN-skipping).
cor_arr <- simplify2array(cor_list)
cor_mat <- apply(cor_arr, c(1, 2), function(v) mean(v, na.rm = TRUE))
dimnames(cor_mat) <- list(ct_all, ct_all)
# Output keeps the same 3-column shape the figure consumes; `rho` is now the
# across-donor MEAN of per-donor spot-level Spearman (F123), not an n=5 value.
fwrite(as.data.table(as.table(cor_mat)),
       file.path(OUTDIR, "celltype_cooccurrence_spearman.csv"))

# Per cell-type: mean abundance in each class (per-spot means).
abund_by_class <- prop[, lapply(.SD, mean, na.rm = TRUE),
                        .SDcols = ct_all, by = spot_class]
fwrite(abund_by_class, file.path(OUTDIR, "celltype_mean_by_spot_class.csv"))

# -------------------------------------------------------------------------- #
# [5] Repulsive-ligand elevation test (F121): excluded vs non-excluded spots,
#     donor-aware, on log1p-CPM spot expression.
# -------------------------------------------------------------------------- #
message("[5] Repulsive-ligand elevation (excluded vs non-excluded spots)...")
lig_res <- NULL
if (!file.exists(LIG_CSV)) {
  # Generate the per-spot lognorm ligand matrix on demand (spatial env helper).
  message("  repulsive_ligand_spot_expr.csv missing - generating via 31a extractor...")
  rc <- tryCatch(
    system2("micromamba",
            c("run", "-n", "spatial", "python",
              file.path(SCRIPTS, "31a_extract_repulsive_ligand_expr.py")),
            stdout = TRUE, stderr = TRUE),
    error = function(e) { message("  extractor failed: ", conditionMessage(e)); NULL })
}
if (file.exists(LIG_CSV)) {
  lig <- fread(LIG_CSV)                 # spot, sample_id, condition, <ligands...>
  ligand_cols <- setdiff(names(lig), c("spot", "sample_id", "condition"))
  message(sprintf("  Ligands available: %s", paste(ligand_cols, collapse = ", ")))
  # Attach spot class; only excluded vs non-excluded spots are compared.
  lig <- merge(lig, prop[, .(spot, spot_class)], by = "spot")
  lig <- lig[spot_class %in% c(excluded_label, "non_excluded")]
  lig[, is_excluded := spot_class == excluded_label]

  # Donor-aware: per donor, mean lognorm expr in excluded vs non-excluded spots;
  # paired Wilcoxon across donors (each donor contributes both classes), BH over
  # the ligands present. With ~5 donors this is the honest test (descriptive).
  rows <- lapply(ligand_cols, function(g) {
    pd <- lig[, .(mu_excl = mean(get(g)[is_excluded], na.rm = TRUE),
                  mu_nonexcl = mean(get(g)[!is_excluded], na.rm = TRUE),
                  n_excl = sum(is_excluded), n_nonexcl = sum(!is_excluded)),
              by = sample_id]
    pd <- pd[is.finite(mu_excl) & is.finite(mu_nonexcl)]
    n_donors <- nrow(pd)
    # effect = mean over donors of (excluded - non-excluded); >0 means elevated
    # in immune-excluded spots (the hypothesized direction).
    delta <- mean(pd$mu_excl - pd$mu_nonexcl)
    pval <- NA_real_
    if (n_donors >= 3 && length(unique(pd$mu_excl - pd$mu_nonexcl)) > 1) {
      pval <- tryCatch(
        suppressWarnings(wilcox.test(pd$mu_excl, pd$mu_nonexcl, paired = TRUE)$p.value),
        error = function(e) NA_real_)
    }
    data.table(ligand = g, n_donors = n_donors,
               mean_excluded = mean(pd$mu_excl), mean_non_excluded = mean(pd$mu_nonexcl),
               delta_excluded_minus_non = delta, pval = pval)
  })
  lig_res <- rbindlist(rows)
  lig_res[, padj_bh := p.adjust(pval, method = "BH")]
  fwrite(lig_res, file.path(OUTDIR, "repulsive_ligand_elevation.csv"))
  message(sprintf("  Tested %d ligand(s), donor-paired Wilcoxon (n~%d donors).",
                  nrow(lig_res), if (nrow(lig_res)) max(lig_res$n_donors) else 0L))
} else {
  message("  SKIPPED step 5: repulsive_ligand_spot_expr.csv unavailable. ",
          "Run 31a_extract_repulsive_ligand_expr.py in the spatial env first. ",
          "Do NOT cite a repulsive-ligand result until this CSV exists.")
}

summary_lines <- c(
  sprintf("Per-spot abundance source: %s", COMP_CSV),
  sprintf("Spots total: %d  |  Donors: %d", nrow(prop), uniqueN(prop$sample_id)),
  sprintf("Excluded class label: %s", excluded_label),
  sprintf("Median immune score: %.3f", median(prop$immune_score, na.rm = TRUE)),
  sprintf("Median stromal score: %.3f", median(prop$stromal_score, na.rm = TRUE)),
  "",
  "=== Spot class distribution ===",
  capture.output(print(table(prop$spot_class))),
  "",
  "=== Per-donor immune exclusion % (donor is the unit, n~5) ===",
  capture.output(print(head(sample_summary[order(-pct_immune_excluded)], 15L))),
  "",
  "=== Cell-type mean abundance by spot class ===",
  capture.output(print(abund_by_class)),
  "",
  "CAVEAT (F123): co-occurrence rho is the across-donor MEAN of per-donor",
  "spot-level Spearman; only ~5 donors and spots within a slice are spatially",
  "autocorrelated, so it is descriptive — no across-donor p-value is claimed.",
  "=== Spearman co-occurrence (top pairs, donor-mean of per-spot rho) ===",
  capture.output({
    cm <- as.data.table(as.table(cor_mat))
    # as.data.table(as.table(...)) yields columns V1, V2, N
    setnames(cm, names(cm), c("ct1","ct2","rho"))
    cm <- cm[ct1 != ct2]
    print(head(cm[order(-abs(rho))], 20L))
  })
)
if (!is.null(lig_res)) {
  summary_lines <- c(summary_lines, "",
    "=== Repulsive-ligand elevation (excluded vs non-excluded spots, donor-paired) ===",
    "delta_excluded_minus_non > 0 = elevated in immune-excluded spots; BH over ligands.",
    capture.output(print(lig_res)))
}
writeLines(summary_lines, file.path(OUTDIR, "immune_exclusion_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
