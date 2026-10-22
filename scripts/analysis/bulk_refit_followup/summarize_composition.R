#!/usr/bin/env Rscript
# Compare the InstaPrism composition rebuilt on the 774 corrected-count units with the
# accepted composition table built on the 844-sample version (kallisto counts).
# Environment: COMP_ROOT (the corrected_774_<TS> root from run_composition_rebuild.sbatch),
#              COMP_SUMMARY_ROOT (must not exist).
#
# Three comparisons, because two things change at once (count source and sample set):
#   old_844  accepted_composition.tsv.gz rows of the five cohorts (844 samples)
#   old_774  the same table restricted to the 774 units (set effect only)
#   new_774  rebuilt on corrected fragment counts (count-source effect, paired by sample)
# Acceptance diagnostics reuse the thresholds frozen in
# scripts/manuscript/program_observability_map/config.R. The estimator-acceptance step itself
# (nine cohorts, 1,257 samples) is not rerun.

suppressPackageStartupMessages(library(data.table))
set.seed(20260923)
fail <- function(...) stop(..., call. = FALSE)
COMP_ROOT <- Sys.getenv("COMP_ROOT", "")
OUT <- Sys.getenv("COMP_SUMMARY_ROOT", "")
if (!nzchar(COMP_ROOT) || !nzchar(OUT)) fail("COMP_ROOT and COMP_SUMMARY_ROOT are required")
if (file.exists(OUT)) fail("Refusing to overwrite: ", OUT)
dir.create(OUT, recursive = TRUE)

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
V6 <- file.path(BASE, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10")
OLD <- file.path(V6, "workstreams/BULK-PROGRAM-MAP-v9/precontract/composition_acceptance/accepted_composition.tsv.gz")
OLD_META <- file.path(V6, "inputs/BG001-DECISION/arms/F_five/results/integration/meta_matched.rds")
NEW_META <- file.path(BASE, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds")
COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
ZERO_MAX <- 0.05; FLOOR <- 1e-4; HEP_RANGE <- c(0.60, 0.95)
EXPECTED_SIGN <- c(Hepatocytes = -1, Fibroblasts = 1, Cholangiocytes = 1)

read_props <- function(path) {
  v <- fread(path, check.names = FALSE)
  setnames(v, 1L, "sample_id")
  v[, sample_id := as.character(sample_id)]
  v
}

# ---- control: unchanged kallisto input must reproduce the stored GSE126848 proportions
stored <- read_props(file.path(BASE, "Analysis/Deconvolution/results/GSE126848/GSE126848_bayesprism_proportions.tsv"))
ctrl <- read_props(file.path(COMP_ROOT, "control_kallisto_input/GSE126848/GSE126848_bayesprism_proportions.tsv"))
lin_all <- setdiff(names(stored), "sample_id")
if (!setequal(names(ctrl), names(stored))) fail("Control lineage columns differ from the stored table")
ctrl <- ctrl[match(stored$sample_id, sample_id)]
control <- data.table(check = "unchanged_kallisto_GSE126848_vs_stored_proportions",
                      n_samples = nrow(stored),
                      max_abs_diff = max(abs(as.matrix(ctrl[, ..lin_all]) - as.matrix(stored[, ..lin_all]))),
                      tolerance = 1e-6)
control[, pass := !is.na(max_abs_diff) & max_abs_diff <= tolerance]
fwrite(control, file.path(OUT, "reproduction_control.tsv"), sep = "\t")
print(control)

# ---- composition tables
old <- fread(OLD)[dataset %in% COHORTS]
lineages <- setdiff(names(old), c("sample_id", "dataset"))
if (nrow(old) != 844L) fail("Accepted composition does not hold 844 samples for the five cohorts")
new <- rbindlist(lapply(COHORTS, function(ds) {
  v <- read_props(file.path(COMP_ROOT, "results", ds, paste0(ds, "_bayesprism_proportions.tsv")))
  if (!setequal(setdiff(names(v), "sample_id"), lineages)) fail("Lineage columns differ for ", ds)
  v[, dataset := ds]
  v[, c("sample_id", "dataset", lineages), with = FALSE]
}))
meta_new <- as.data.table(readRDS(NEW_META))
meta_old <- as.data.table(readRDS(OLD_META))
if (nrow(new) != 774L || !setequal(new$sample_id, meta_new$sample_id)) fail("Rebuilt composition is not the 774 units")
if (max(abs(rowSums(new[, ..lineages]) - 1)) > 1e-6) fail("Rebuilt fractions do not sum to 1")
fwrite(new, file.path(OUT, "composition_new_774.tsv.gz"), sep = "\t")
old774 <- old[sample_id %in% new$sample_id]
paired <- merge(old774, new, by = c("sample_id", "dataset"), suffixes = c(".old", ".new"))
if (nrow(paired) != 774L) fail("Old and new composition do not pair on 774 samples")

# ---- shift by cohort and lineage
shift <- rbindlist(lapply(COHORTS, function(ds) rbindlist(lapply(lineages, function(l) {
  p <- paired[dataset == ds]
  o <- p[[paste0(l, ".old")]]; n <- p[[paste0(l, ".new")]]
  data.table(cohort = ds, lineage = l, n_774 = nrow(p),
             mean_old_844set = old[dataset == ds, mean(get(l))],
             mean_old_774 = mean(o), mean_new_774 = mean(n), delta_mean_new_minus_old = mean(n) - mean(o),
             median_old_774 = median(o), median_new_774 = median(n),
             paired_pearson = if (sd(o) > 0 && sd(n) > 0) cor(o, n) else NA_real_,
             paired_spearman = if (sd(o) > 0 && sd(n) > 0) cor(o, n, method = "spearman") else NA_real_,
             mean_abs_paired_diff = mean(abs(n - o)))
}))))
fwrite(shift, file.path(OUT, "composition_shift_by_cohort_lineage.tsv"), sep = "\t")

# ---- hepatocyte and non-parenchymal means by group (unadjusted participant means, as in Figure S3)
long <- function(comp, meta, version) {
  m <- merge(comp[, .(sample_id, dataset, hepatocyte = Hepatocytes)], meta[, .(sample_id, group_binary, fibrosis_stage)], by = "sample_id")
  m[, non_parenchymal := 1 - hepatocyte]
  rbindlist(list(
    m[, .(grouping = "group_binary", group = group_binary, sample_id, dataset, hepatocyte, non_parenchymal)],
    m[!is.na(fibrosis_stage), .(grouping = "fibrosis_stage", group = paste0("F", fibrosis_stage), sample_id, dataset, hepatocyte, non_parenchymal)]
  ))[, version := version]
}
g <- rbindlist(list(
  long(old, meta_old, "old_844_own_labels"),
  long(old774, meta_new, "old_774"),
  long(new, meta_new, "new_774")))
by_group <- g[, .(n = .N,
                  mean_hepatocyte = mean(hepatocyte), se_hepatocyte = sd(hepatocyte) / sqrt(.N),
                  mean_non_parenchymal = mean(non_parenchymal)), by = .(version, dataset, grouping, group)]
by_group <- dcast(by_group, dataset + grouping + group ~ version,
                  value.var = c("n", "mean_hepatocyte", "mean_non_parenchymal"))
setorder(by_group, dataset, grouping, group)
fwrite(by_group, file.path(OUT, "composition_hepatocyte_npc_by_group.tsv"), sep = "\t")

# ---- acceptance diagnostics on the five cohorts (thresholds from config.R)
diag_one <- function(comp, meta, version) {
  m <- merge(comp, meta[, .(sample_id, fibrosis_stage)], by = "sample_id")
  rbindlist(lapply(COHORTS, function(ds) {
    d <- m[dataset == ds]
    M <- as.matrix(d[, ..lineages])
    p05 <- apply(M, 2L, quantile, probs = 0.05, names = FALSE)
    rho <- vapply(names(EXPECTED_SIGN), function(l) {
      k <- is.finite(d$fibrosis_stage)
      if (sum(k) > 10L && sd(d[[l]][k]) > 0) cor(d[[l]][k], d$fibrosis_stage[k], method = "spearman") else NA_real_
    }, numeric(1))
    data.table(version = version, cohort = ds, n = nrow(d),
               max_zero_fraction = max(colMeans(M == 0)),
               median_hepatocyte = median(d$Hepatocytes),
               n_lineages_p05_below_1e4 = sum(p05 < FLOOR),
               spearman_hepatocytes_stage = rho[["Hepatocytes"]],
               spearman_fibroblasts_stage = rho[["Fibroblasts"]],
               spearman_cholangiocytes_stage = rho[["Cholangiocytes"]],
               direction_ok = sum(sign(rho) == EXPECTED_SIGN[names(rho)], na.rm = TRUE))
  }))
}
diag <- rbindlist(list(diag_one(old, meta_old, "old_844"),
                       diag_one(old774, meta_new, "old_774"),
                       diag_one(new, meta_new, "new_774")))
diag[, passes_zero := max_zero_fraction <= ZERO_MAX]
diag[, passes_hepatocyte_range := median_hepatocyte >= HEP_RANGE[[1]] & median_hepatocyte <= HEP_RANGE[[2]]]
fwrite(diag, file.path(OUT, "acceptance_diagnostics_five_cohorts.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
options(width = 230, digits = 4)
cat("\n== mean fraction, five main lineages (old_774 -> new_774) ==\n")
print(shift[lineage %in% c("Hepatocytes", "Cholangiocytes", "Fibroblasts", "Macrophages", "Endothelial cells", "T cells"),
            .(cohort, lineage, n_774, mean_old_844set, mean_old_774, mean_new_774, paired_spearman)])
cat("\n== hepatocyte mean by group_binary ==\n")
print(by_group[grouping == "group_binary"])
cat("\n== acceptance diagnostics ==\n"); print(diag)
cat("COMPOSITION_SUMMARY_COMPLETE", OUT, "\n")
