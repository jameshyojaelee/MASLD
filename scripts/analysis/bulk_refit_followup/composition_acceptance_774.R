#!/usr/bin/env Rscript
# Bulk composition acceptance and Figure S3 composition summaries on the 774-unit roster.
#
# Estimator: BayesPrism (InstaPrism implementation), a Bayesian deconvolution of bulk RNA-seq
# counts against a single-cell RNA-seq reference. It is unchanged from the accepted table:
# Analysis/Deconvolution/scripts/12_run_bayesprism.R with reference_human_scalesc_sce.rds.
# This script runs no deconvolution. It reads
#   five recounted, one-record-per-person cohorts: the rebuilt proportions (job 22135553)
#   four other cohorts (GSE167523, GSE174478, GSE193066, GSE240729): the accepted table, unchanged
# and applies the acceptance criteria frozen in scripts/manuscript/program_observability_map/
# config.R to three bases:
#   old_1257         accepted table as stored (844 + 413 records); must reproduce the stored audit
#   old_values_774   accepted values, five cohorts restricted to the 774 units (roster effect only)
#   new_774          rebuilt values on the 774 units (roster and count-source effect)
# The four unchanged cohorts are identical in all three bases.
# Only BayesPrism is evaluated. MuSiC and Rectangle were not rerun on the corrected counts, so
# the three-estimator choice is not repeated here.
#
# Environment: COMPOSITION_CANDIDATE_ROOT (must not exist).

suppressPackageStartupMessages(library(data.table))
options(digits = 17)
SEED <- 20260923L
set.seed(SEED)  # no random step is used; set for the record
fail <- function(...) stop(..., call. = FALSE)
out <- Sys.getenv("COMPOSITION_CANDIDATE_ROOT", "")
if (!nzchar(out)) fail("COMPOSITION_CANDIDATE_ROOT is required")
if (file.exists(out)) fail("Refusing to overwrite: ", out)
for (d in c("composition", "acceptance", "comparison")) dir.create(file.path(out, d), recursive = TRUE)

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
V6 <- file.path(BASE, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10")
ACC <- file.path(V6, "workstreams/BULK-PROGRAM-MAP-v9/precontract/composition_acceptance")
OLD_COMP <- file.path(ACC, "accepted_composition.tsv.gz")
OLD_META <- file.path(V6, "inputs/BG001-DECISION/arms/F_legacy/results/integration/meta_matched.rds")
NEW_ROOT <- file.path(BASE, "Analysis/Deconvolution/corrected_774_20261009T141331Z")
NEW_META <- file.path(BASE, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds")
STAGE_NEW <- file.path(BASE, "figures/candidates/bulk-release-corrected-controls-20261008T162308Z")
STAGE_OLD <- file.path(BASE, "figures/candidates/bulk-release-corrected-controls-20260925T175215Z")
SCREEN <- file.path(BASE, "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_primary_sample_ids.tsv")

FIVE <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
FOUR <- c("GSE167523", "GSE174478", "GSE193066", "GSE240729")
# Frozen in config.R
DISCOVERY <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478")
HOLDOUT <- "GSE193066"
REQUIRED <- c(DISCOVERY, HOLDOUT)
ZERO_MAX <- 0.05; FLOOR <- 1e-4; HEP_RANGE <- c(0.60, 0.95); MIN_DIRECTION <- 3L
COVARIATE_LINEAGES <- c("Cholangiocytes", "Fibroblasts", "Macrophages", "T cells")
EXPECTED_SIGN <- c(Hepatocytes = -1, Fibroblasts = 1, Cholangiocytes = 1)

sha256_file <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE)[[1L]], "[[:space:]]+")[[1L]][[1L]]
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")

# ------------------------------------------------------------------ inputs
old <- fread(OLD_COMP, check.names = FALSE)
lineages <- setdiff(names(old), c("sample_id", "dataset"))
if (nrow(old) != 1257L || anyDuplicated(old$sample_id) || length(lineages) != 16L) fail("Accepted table is not 1,257 x 16")
if (!setequal(unique(old$dataset), c(FIVE, FOUR))) fail("Accepted table cohorts differ from the nine expected")

new_five <- rbindlist(lapply(FIVE, function(ds) {
  path <- file.path(NEW_ROOT, "results", ds, paste0(ds, "_bayesprism_proportions.tsv"))
  v <- read.delim(path, row.names = 1L, check.names = FALSE)
  if (!setequal(colnames(v), lineages)) fail("Lineage columns differ for ", ds)
  data.table(sample_id = rownames(v), dataset = ds, as.data.table(v[, lineages, drop = FALSE]))
}))
meta_new <- as.data.table(readRDS(NEW_META))
meta_old <- as.data.table(readRDS(OLD_META))
if (nrow(new_five) != 774L || anyDuplicated(new_five$sample_id) || !setequal(new_five$sample_id, meta_new$sample_id) ||
    anyDuplicated(meta_new$analysis_unit_id)) fail("Rebuilt proportions are not the 774 one-per-person units")
m_new <- as.matrix(new_five[, ..lineages])
if (anyNA(m_new) || any(m_new < 0 | m_new > 1) || max(abs(rowSums(m_new) - 1)) > 1e-6) fail("Rebuilt fractions are not a composition")
if (!all(new_five$dataset == meta_new$dataset[match(new_five$sample_id, meta_new$sample_id)])) fail("Cohort labels disagree with the roster")
if (!all(new_five$sample_id %in% old$sample_id)) fail("A 774 unit is absent from the accepted table")

unchanged <- old[dataset %in% FOUR]
if (nrow(unchanged) != 413L) fail("The four unchanged cohorts do not hold 413 records")
tables <- list(
  old_1257 = old,
  old_values_774 = rbind(old[sample_id %in% new_five$sample_id], unchanged),
  new_774 = rbind(new_five, unchanged)
)
if (nrow(tables$old_values_774) != 1187L || nrow(tables$new_774) != 1187L) fail("Hybrid tables are not 774 + 413 records")

# Stage labels: the 774 roster for the five corrected cohorts in the two 774 bases, the
# accepted-run metadata elsewhere.
phenotype <- function(version) {
  p_old <- meta_old[, .(sample_id, fibrosis_stage = as.numeric(fibrosis_stage))]
  if (version == "old_1257") return(p_old)
  rbind(meta_new[, .(sample_id, fibrosis_stage = as.numeric(fibrosis_stage))],
        p_old[sample_id %in% unchanged$sample_id])
}

hybrid <- copy(tables$new_774)
fwrite(hybrid, file.path(out, "composition", "accepted_composition_774_plus_unchanged4.tsv.gz"), sep = "\t", quote = FALSE)
basis <- hybrid[, .(sample_id, dataset,
  basis = fifelse(dataset %in% FIVE, "rebuilt_corrected_fragment_counts_one_record_per_person",
                  "unchanged_accepted_table_kallisto_counts_roster_not_identity_resolved"))]
write_tsv(basis, file.path(out, "composition", "composition_basis.tsv"))
write_tsv(basis[, .(n_records = .N), by = .(dataset, basis)][order(dataset)], file.path(out, "composition", "composition_basis_by_cohort.tsv"))

# ------------------------------------------------------------------ acceptance criteria
audit <- function(version) {
  comp <- tables[[version]]
  ph <- phenotype(version)
  M <- as.matrix(comp[, ..lineages])
  lineage <- data.table(version = version, celltype = lineages, zero_fraction = colMeans(M == 0),
                        median_fraction = apply(M, 2L, median), max_fraction = apply(M, 2L, max))
  dynamic <- rbindlist(lapply(lineages, function(l) comp[dataset %in% REQUIRED, .(
    version = version, celltype = l, n = .N, p05 = quantile(get(l), 0.05, names = FALSE), p50 = median(get(l)),
    fraction_below_floor = mean(get(l) < FLOOR)), by = dataset]))
  # Same statistic as the accepted audit: per-cohort Spearman of fraction against Kleiner stage,
  # more than 10 staged records. GSE213621 is reported but is not a discovery cohort in the rule.
  direction <- rbindlist(lapply(names(EXPECTED_SIGN), function(l) rbindlist(lapply(c(DISCOVERY, "GSE213621"), function(cohort) {
    d <- merge(comp[dataset == cohort, c("sample_id", l), with = FALSE], ph, by = "sample_id")
    d <- d[is.finite(fibrosis_stage)]
    ok <- nrow(d) > 10L && sd(d[[l]]) > 0 && sd(d$fibrosis_stage) > 0
    test <- if (ok) suppressWarnings(cor.test(d[[l]], d$fibrosis_stage, method = "spearman", exact = FALSE)) else NULL
    rho <- if (ok) unname(test$estimate) else NA_real_
    data.table(version = version, celltype = l, cohort = cohort, in_acceptance_rule = cohort %in% DISCOVERY,
               n = nrow(d), spearman_fibrosis = rho, p_nominal = if (ok) test$p.value else NA_real_,
               expected_sign = EXPECTED_SIGN[[l]], agrees = isTRUE(is.finite(rho) && sign(rho) == EXPECTED_SIGN[[l]]))
  }))))
  # p values are a companion to the frozen sign rule. BH family: every lineage-by-cohort
  # correlation computed for this basis (3 lineages x 5 cohorts).
  direction[, p_bh := p.adjust(p_nominal, method = "BH")]
  direction[, bh_family_size := sum(!is.na(p_nominal))]
  eligibility <- dynamic[, .(n_cohorts = .N, min_p05 = min(p05), passes_dynamic_range = all(p05 >= FLOOR)),
                         by = .(version, celltype)]
  elig <- eligibility[celltype %in% COVARIATE_LINEAGES & passes_dynamic_range == TRUE, celltype]
  dir_rule <- direction[in_acceptance_rule == TRUE, .(ok = sum(agrees) >= MIN_DIRECTION), by = celltype]
  hep <- lineage[celltype == "Hepatocytes", median_fraction]
  criteria <- data.table(
    version = version, method = "bayesprism", n_records = nrow(comp),
    min_required_cohort_coverage = 1,  # every record of each basis has a complete row (checked above)
    zero_inflation_max = max(lineage$zero_fraction), passes_zero_inflation = max(lineage$zero_fraction) <= ZERO_MAX,
    median_hepatocyte = hep, passes_plausibility = hep >= HEP_RANGE[[1L]] && hep <= HEP_RANGE[[2L]],
    n_direction_lineages_ok = sum(dir_rule$ok), passes_direction = all(dir_rule$ok),
    n_eligible_covariate_lineages = length(elig), eligible_covariate_lineages = paste(elig, collapse = ";"))
  criteria[, eligible := passes_zero_inflation & passes_plausibility & passes_direction]
  cohort <- comp[, .(version = version, n_records = .N, median_hepatocyte = median(Hepatocytes),
                     mean_hepatocyte = mean(Hepatocytes), max_zero_fraction = max(colMeans(as.matrix(.SD) == 0))),
                 by = dataset, .SDcols = lineages]
  list(lineage = lineage, dynamic = dynamic, direction = direction, eligibility = eligibility,
       criteria = criteria, cohort = cohort)
}
audits <- lapply(names(tables), audit); names(audits) <- names(tables)
bind <- function(name) rbindlist(lapply(audits, `[[`, name))
write_tsv(bind("lineage"), file.path(out, "acceptance", "lineage_audit.tsv"))
write_tsv(bind("dynamic"), file.path(out, "acceptance", "dynamic_range_audit.tsv"))
write_tsv(bind("direction"), file.path(out, "acceptance", "stage_direction_audit.tsv"))
write_tsv(bind("eligibility"), file.path(out, "acceptance", "lineage_eligibility.tsv"))
write_tsv(bind("criteria"), file.path(out, "acceptance", "method_criteria.tsv"))
write_tsv(bind("cohort"), file.path(out, "acceptance", "cohort_summary.tsv"))

# The old_1257 basis must reproduce the stored acceptance audit. This checks this script, not the data.
stored <- list(
  lineage = fread(file.path(ACC, "lineage_audit.tsv"))[method == "bayesprism"],
  direction = fread(file.path(ACC, "stage_direction_audit.tsv"))[method == "bayesprism"],
  eligibility = fread(file.path(ACC, "lineage_eligibility.tsv"))[method == "bayesprism"],
  criteria = fread(file.path(ACC, "method_criteria.tsv"))[method == "bayesprism"])
mine <- audits$old_1257
j_lin <- merge(stored$lineage, mine$lineage, by = "celltype", suffixes = c(".stored", ".mine"))
j_dir <- merge(stored$direction, mine$direction, by = c("celltype", "cohort"), suffixes = c(".stored", ".mine"))
j_el <- merge(stored$eligibility, mine$eligibility, by = "celltype", suffixes = c(".stored", ".mine"))
reproduction <- data.table(
  check = c("lineage_median_fraction", "lineage_zero_fraction", "stage_direction_spearman", "stage_direction_n",
            "dynamic_range_min_p05", "median_hepatocyte", "eligible_covariate_lineages"),
  n_compared = c(nrow(j_lin), nrow(j_lin), nrow(j_dir), nrow(j_dir), nrow(j_el), 1L, 1L),
  max_abs_diff = c(max(abs(j_lin$median_fraction.stored - j_lin$median_fraction.mine)),
                   max(abs(j_lin$zero_fraction.stored - j_lin$zero_fraction.mine)),
                   max(abs(j_dir$spearman_fibrosis.stored - j_dir$spearman_fibrosis.mine)),
                   max(abs(j_dir$n.stored - j_dir$n.mine)),
                   max(abs(j_el$min_p05.stored - j_el$min_p05.mine)),
                   abs(stored$criteria$median_hepatocyte - mine$criteria$median_hepatocyte),
                   as.numeric(stored$criteria$eligible_covariate_lineages != mine$criteria$eligible_covariate_lineages)))
reproduction[, pass := n_compared > 0L & max_abs_diff <= 1e-12]
write_tsv(reproduction, file.path(out, "acceptance", "reproduction_of_stored_audit.tsv"))
if (nrow(j_lin) != 16L || nrow(j_dir) != 12L || nrow(j_el) != 16L || !all(reproduction$pass)) {
  fail("old_1257 does not reproduce the stored acceptance audit; see acceptance/reproduction_of_stored_audit.tsv")
}

# Identity-screen sensitivity for the four unchanged cohorts. A genotype screen exists for them
# (all_screened_primary_sample_ids.tsv) but no recount. GSE193066 has no qualified unit.
screen_ids <- fread(SCREEN, header = TRUE)[[1L]]
screened <- unchanged[, .(n_accepted_table = .N, n_passing_identity_screen = sum(sample_id %in% screen_ids)), by = dataset]
write_tsv(screened, file.path(out, "acceptance", "unchanged_cohort_identity_screen_counts.tsv"))
sens <- rbind(new_five, unchanged[sample_id %in% screen_ids])
sens_dyn <- rbindlist(lapply(COVARIATE_LINEAGES, function(l) sens[dataset %in% REQUIRED, .(
  celltype = l, n = .N, p05 = quantile(get(l), 0.05, names = FALSE)), by = dataset]))
sens_dyn[, basis := "new_774_plus_screened_units_of_unchanged_cohorts"]
write_tsv(sens_dyn, file.path(out, "acceptance", "dynamic_range_identity_screen_sensitivity.tsv"))

# ------------------------------------------------------------------ Figure S3 composition summaries
# Same estimand as progression_cascade.R: unadjusted mean and SE of the hepatocyte fraction and of
# 1 - hepatocyte, one record per analysis unit, by fibrosis stage and by NAS group.
read_manifests <- function(root) {
  ext <- file.path(root, "analysis", "stage_extensions")
  fib <- fread(file.path(ext, "stage_sample_manifest.tsv"))
  if (!"analysis_unit_id" %in% names(fib)) fib[, analysis_unit_id := paste(dataset, sample_id, sep = "::")]
  fib <- unique(fib[, .(sample_id, analysis_unit_id, dataset, group = paste0("F", fibrosis_stage))], by = "analysis_unit_id")
  nas <- fread(file.path(ext, "stage_extension_sample_manifest.tsv"))[axis == "NAS"]
  if (!"analysis_unit_id" %in% names(nas)) nas[, analysis_unit_id := paste(dataset, sample_id, sep = "::")]
  nas <- unique(nas[, .(sample_id, analysis_unit_id, dataset, group)], by = "analysis_unit_id")
  rbind(fib[, axis := "fibrosis"], nas[, axis := "NAS"])
}
summarize <- function(comp, samples, version) {
  d <- merge(comp[, .(sample_id, dataset, hepatocyte = Hepatocytes)], samples, by = c("sample_id", "dataset"))
  if (nrow(d) != nrow(samples)) fail("Composition join is incomplete for ", version)
  d[, non_parenchymal := 1 - hepatocyte]
  z <- melt(d, id.vars = c("analysis_unit_id", "dataset", "axis", "group"), measure.vars = c("hepatocyte", "non_parenchymal"),
            variable.name = "compartment", value.name = "fraction")
  list(pooled = z[, .(version = version, mean = mean(fraction), SE = sd(fraction) / sqrt(.N), n_participants = .N),
                  by = .(axis, group, compartment)],
       by_cohort = z[, .(version = version, mean = mean(fraction), SE = sd(fraction) / sqrt(.N), n_participants = .N),
                     by = .(axis, dataset, group, compartment)],
       units = d)
}
man_new <- read_manifests(STAGE_NEW)
man_old <- read_manifests(STAGE_OLD)
s <- list(old_844 = summarize(old, man_old, "old_844_records_accepted_values"),
          old_values_774 = summarize(old, man_new, "774_units_accepted_values"),
          new_774 = summarize(new_five, man_new, "774_units_rebuilt_values"))
pooled <- rbindlist(lapply(s, `[[`, "pooled"))
write_tsv(pooled, file.path(out, "comparison", "s3_composition_group_means_long.tsv"))
write_tsv(rbindlist(lapply(s, `[[`, "by_cohort")), file.path(out, "comparison", "s3_composition_group_means_by_cohort_long.tsv"))
wide <- dcast(pooled[compartment == "hepatocyte"], axis + group ~ version, value.var = c("mean", "SE", "n_participants"))
write_tsv(wide, file.path(out, "comparison", "s3_hepatocyte_fraction_old_vs_new.tsv"))

# The two existing Figure S3 composition source tables must be reproduced by this code.
check_table <- function(path, mine, label) {
  ref <- fread(path)[compartment == "hepatocyte", .(axis, group = as.character(group), mean_ref = mean, n_ref = n_participants)]
  j <- merge(ref, mine[compartment == "hepatocyte", .(axis, group = as.character(group), mean, n_participants)], by = c("axis", "group"))
  data.table(check = label, n_groups = nrow(j), max_abs_mean_diff = max(abs(j$mean_ref - j$mean)),
             n_mismatch = sum(j$n_ref != j$n_participants))
}
s3_checks <- rbind(
  check_table(file.path(STAGE_OLD, "source_tables", "figs3_composition_summary.tsv"), s$old_844$pooled,
              "844-record Figure S3 composition table reproduced"),
  check_table(file.path(STAGE_NEW, "source_tables", "figs3_composition_summary.tsv"), s$old_values_774$pooled,
              "existing 774-unit Figure S3 composition table is the accepted (old) values on 774 units"))
s3_checks[, pass := n_groups == 9L & max_abs_mean_diff <= 1e-12 & n_mismatch == 0L]
write_tsv(s3_checks, file.path(out, "comparison", "s3_source_table_reproduction.tsv"))

# Stage trend on the pooled S3 roster, as a descriptive companion: Spearman of hepatocyte
# fraction against the ordered group, one record per unit. Rank-based, so identical for the
# hepatocyte fraction, its logit and the hepatocyte / non-parenchymal log-ratio.
trend <- rbindlist(lapply(names(s), function(v) s[[v]]$units[, {
  ord <- if (axis[[1L]] == "fibrosis") match(group, paste0("F", 0:4)) else match(group, c("NAS0", "NAS1-2", "NAS3-4", "NAS5-8"))
  t <- suppressWarnings(cor.test(hepatocyte, ord, method = "spearman", exact = FALSE))
  .(version = v, n = .N, n_cohorts = uniqueN(dataset), spearman_hepatocyte_vs_ordered_group = unname(t$estimate), p_nominal = t$p.value)
}, by = axis]))
trend[, p_bh := p.adjust(p_nominal, method = "BH"), by = version]  # family: the two axes of one basis
write_tsv(trend, file.path(out, "comparison", "s3_hepatocyte_stage_trend.tsv"))

# Paired change of the same 774 libraries between accepted and rebuilt values.
paired <- merge(old[, c("sample_id", "dataset", lineages), with = FALSE], new_five, by = c("sample_id", "dataset"), suffixes = c(".old", ".new"))
if (nrow(paired) != 774L) fail("Paired table is not 774")
shift <- rbindlist(lapply(c(FIVE, "all_five"), function(ds) rbindlist(lapply(c("Hepatocytes", COVARIATE_LINEAGES, "Endothelial cells"), function(l) {
  p <- if (ds == "all_five") paired else paired[dataset == ds]
  o <- p[[paste0(l, ".old")]]; n <- p[[paste0(l, ".new")]]
  data.table(cohort = ds, lineage = l, n = nrow(p), mean_old = mean(o), mean_new = mean(n), mean_paired_diff = mean(n - o),
             median_old = median(o), median_new = median(n), paired_spearman = cor(o, n, method = "spearman"))
}))))
write_tsv(shift, file.path(out, "comparison", "paired_shift_774_units.tsv"))

inputs <- c(accepted_composition_1257 = OLD_COMP, accepted_run_metadata = OLD_META, roster_774_metadata = NEW_META,
            identity_screen_ids = SCREEN,
            stage_manifest_774 = file.path(STAGE_NEW, "analysis/stage_extensions/stage_sample_manifest.tsv"),
            nas_manifest_774 = file.path(STAGE_NEW, "analysis/stage_extensions/stage_extension_sample_manifest.tsv"),
            stage_manifest_844 = file.path(STAGE_OLD, "analysis/stage_extensions/stage_sample_manifest.tsv"),
            setNames(file.path(NEW_ROOT, "results", FIVE, paste0(FIVE, "_bayesprism_proportions.tsv")), paste0("rebuilt_", FIVE)))
write_tsv(data.table(input = names(inputs), path = inputs, sha256 = vapply(inputs, sha256_file, character(1))),
          file.path(out, "input_manifest.tsv"))
writeLines(c(paste("seed", SEED, "(set; the script has no random step)"), capture.output(sessionInfo())), file.path(out, "sessionInfo.txt"))

options(width = 220, digits = 5)
cat("\n== criteria ==\n"); print(bind("criteria"))
cat("\n== stage direction ==\n"); print(bind("direction"))
cat("\n== covariate lineage dynamic range ==\n"); print(bind("eligibility")[celltype %in% COVARIATE_LINEAGES])
cat("\n== Figure S3 hepatocyte fraction ==\n"); print(wide)
cat("\n== S3 source table checks ==\n"); print(s3_checks)
cat("\n== trend ==\n"); print(trend)
cat("COMPOSITION_ACCEPTANCE_774_COMPLETE", out, "\n")
