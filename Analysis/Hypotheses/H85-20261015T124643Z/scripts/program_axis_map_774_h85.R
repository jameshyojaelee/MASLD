#!/usr/bin/env Rscript

# Figure 4E program-axis models on the adopted 774-unit bulk roster.
#
# Reuses, unchanged, the scoring, CLR-component, HC3 and meta-analysis functions
# of scripts/manuscript/program_observability_map/analysis_lib.R and the frozen
# constants in config.R. No program is rediscovered, reweighted or filtered by
# outcome. Four arms are fitted with the same code path:
#
#   old469_reproduction      stored BULK-PROGRAM-MAP-v9 inputs, four cohorts
#                            (positive control; must reproduce the stored map)
#   old_three_cohorts        stored inputs without GSE174478 (isolates the loss
#                            of the fourth cohort from the change of roster)
#   new774_three_cohorts     PRIMARY: adopted 774-unit DGE and rebuilt
#                            composition; cohorts carrying both fibrosis and NAS
#   new774_plus_GSE174478    primary arm plus the unchanged GSE174478 records
#                            (outside the 774 units; old four-cohort design)
#
# H85 copy (Analysis/Hypotheses/H85-20261015T124643Z/PRESPEC_H85.md, "Compute").
# Changes from the stored producer: output root from H85_OUT; seed 20260923;
# input sha256 checks (S1); fit_arm also returns its roster and effect tables;
# the stored-map assertion feeds stop rule S2 instead of aborting; after the
# reproduction block, a second build_composition_scores call with the earlier
# fractions on the 774-unit metadata (pcs_earlier), the arms
# new774_plus_GSE174478_earlierfrac, new774_three_cohorts_earlierfrac (D2) and
# D1, PC-1a, diagnostics 2 to 9, NC-1, D3, and at the end P1/P2, PC-2, PC-3,
# D4, D5 and RESULT_H85.md. H85_DEBUG=1 replaces every earlier-set input of the
# outcome with within-cohort permuted rows (null data) for debugging.

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
})

repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
Sys.setenv(MASLD_PROJECT_ROOT = repo)
lib_dir <- file.path(repo, "scripts/manuscript/program_observability_map")
source(file.path(lib_dir, "config.R"))
source(file.path(lib_dir, "analysis_lib.R"))

out_root <- Sys.getenv("H85_OUT")
if (!nzchar(out_root)) fail("H85_OUT is required")
analysis_root <- file.path(out_root, "analysis")
if (dir.exists(analysis_root)) fail("Refusing to overwrite: ", analysis_root)
dir.create(analysis_root, recursive = TRUE)

# No step below draws random numbers (OLS, HC3, prcomp, REML). The seed is the
# frozen config seed, set so that any later stochastic addition is reproducible.
h85_seed <- 20260923L
set.seed(h85_seed)

# ---- H85 run controls, result writer, stop rules ---------------------------
h85_dir <- file.path(repo, "Analysis/Hypotheses/H85-20261015T124643Z")
h85_prespec_sha256 <- "5cb77b8edcbdfc7b7007d107ed7fd96d0c69041e4646a7a491c9cc56672ef3ef"
h85_debug <- identical(Sys.getenv("H85_DEBUG"), "1")
nc1_draws <- as.integer(Sys.getenv("H85_NC1_B", "100"))
if (!h85_debug && !identical(nc1_draws, 100L)) fail("The final run uses exactly 100 NC-1 draws")
result_path <- file.path(out_root, if (h85_debug) "RESULT_H85_DEBUG.md" else "RESULT_H85.md")
diagnostics_path <- file.path(out_root, "diagnostics.tsv")

fmt_value <- function(v) {
  if (is.numeric(v)) v <- signif(v, 8)
  paste(as.character(v), collapse = "; ")
}
diag_rows <- list()
add_diag <- function(diagnostic, item, value, expected = NA, pass = NA) {
  row <- data.table(diagnostic = as.character(diagnostic), item = item,
                    value = fmt_value(value), expected = fmt_value(expected),
                    pass = as.logical(pass))
  diag_rows[[length(diag_rows) + 1L]] <<- row
  message(sprintf("[diagnostic %s] %s: %s (expected %s; pass %s)",
                  diagnostic, item, row$value, row$expected, row$pass))
  write_tsv(rbindlist(diag_rows), diagnostics_path)
  invisible(row)
}
result_header <- function(verdict) {
  c("# RESULT H85", "",
    if (h85_debug) c("DEBUG RUN. Every earlier-set input of the outcome was replaced by within-cohort permuted rows (null data). Not a result.", "") else NULL,
    paste0("Verdict: ", verdict), "",
    paste0("Prespecification: ", file.path(h85_dir, "PRESPEC_H85.md"), " (sha256 ", h85_prespec_sha256, ")"),
    paste0("Run directory: ", out_root),
    paste0("Written: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("Slurm job: ", Sys.getenv("SLURM_JOB_ID", "none")), "")
}
stop_rule <- function(rule, reason) {
  verdict <- if (rule == "S5") "UNTESTABLE (S5)" else paste0("STOPPED (", rule, ")")
  writeLines(c(result_header(verdict),
               paste0("Stop rule ", rule, ": ", reason), "",
               "No number from this run is used. Finding R05-10 stays open with this reason.", "",
               paste0("Diagnostics written so far: ", diagnostics_path)), result_path)
  message("H85_STOP_RULE ", rule, ": ", reason)
  quit(save = "no", status = 0L)
}
options(error = function() {
  msg <- geterrmessage()
  try(writeLines(c(result_header("FAILED (unexpected error; no number from this run is used)"),
                   "```", msg, "```"), result_path), silent = TRUE)
  quit(save = "no", status = 1L)
})
message("H85 run: debug=", h85_debug, ", NC-1 draws=", nc1_draws, ", seed=", h85_seed)

# ---- H85 diagnostic 1 and stop rule S1: input hashes ------------------------
expected_inputs <- fread(file.path(h85_dir, "scripts", "inputs_expected_sha256.tsv"))
prespec_copy <- file.path(h85_dir, "PRESPEC_H85.md")
prespec_original <- file.path(repo, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/PRESPEC_H85.md")
prespec_text <- paste(readLines(prespec_copy, warn = FALSE), collapse = "\n")
hash_check <- expected_inputs[, .(
  role, path, expected_sha256 = sha256,
  observed_sha256 = vapply(file.path(repo, path), function(p)
    if (file.exists(p)) sha256_file(p) else NA_character_, character(1)),
  listed_in_prespec = vapply(sha256, function(h) grepl(h, prespec_text, fixed = TRUE), logical(1))
)]
hash_check[, matches := !is.na(observed_sha256) & observed_sha256 == expected_sha256 & listed_in_prespec]
write_tsv(hash_check, file.path(out_root, "input_hash_check.tsv"))
prespec_hashes <- c(copy = sha256_file(prespec_copy), original = sha256_file(prespec_original))
prespec_ok <- all(prespec_hashes == h85_prespec_sha256)
add_diag("1", "prespec sha256 (copy; original)", unname(prespec_hashes), h85_prespec_sha256, prespec_ok)
add_diag("1", "inputs whose sha256 matches the prespec table", sum(hash_check$matches), 26L,
         nrow(hash_check) == 26L && all(hash_check$matches))
producer_copy <- file.path(h85_dir, "scripts", "program_axis_map_774_h85.R")
producer_original <- file.path(repo, "figures/candidates/fig4e-program-axis-774-20261012T010830Z",
                               "producers/program_axis_map_774.R")
diff_path <- file.path(out_root, "producer_copy.diff")
suppressWarnings(system2("diff", c("-u", producer_original, producer_copy), stdout = diff_path))
diff_lines <- readLines(diff_path, warn = FALSE)
add_diag("1", "producer copy diff: original lines removed / lines added",
         paste(sum(grepl("^-", diff_lines) & !grepl("^---", diff_lines)),
               sum(grepl("^\\+", diff_lines) & !grepl("^\\+\\+\\+", diff_lines)), sep = " / "),
         "changes limited to the Compute list (reviewed in DEVIATIONS.md)", NA)
if (!prespec_ok || nrow(hash_check) != 26L || !all(hash_check$matches)) {
  stop_rule("S1", paste("sha256 mismatch:", paste(c(
    if (!prespec_ok) "PRESPEC_H85.md" else NULL,
    hash_check[matches == FALSE, role]), collapse = ", ")))
}

new_root <- file.path(
  repo, "RNA-seq/results/manuscript_release/candidates",
  "bulk-corrected-controls-20261008T162308Z/primary_source_controls"
)
new_dge_path <- file.path(new_root, "merged_dge.rds")
new_meta_path <- file.path(new_root, "meta_matched.rds")
new_composition_path <- file.path(
  repo, "figures/candidates/composition-774-20261011T151616Z/composition",
  "accepted_composition_774_plus_unchanged4.tsv.gz"
)
unit_manifest_path <- file.path(
  repo, "figures/candidates/bulk-release-corrected-controls-20261008T162308Z",
  "analysis/stage_extensions/five_cohort_sample_manifest.tsv"
)
stored_map_path <- file.path(discovery_root, "program_axis_map.tsv")
stored_meta_effects_path <- file.path(discovery_root, "meta_effects.tsv")

stromal_ecm <- "hotspot_hepatocytes_f05c535ae5bbc0b9"
ductular <- "hotspot_hepatocytes_48f39dd4d817a10e"
three_cohorts <- c("GSE130970", "GSE135251", "GSE162694")
extra_cohort <- "GSE174478"

registry <- fread(program_registry)
membership <- fread(program_membership)
assert_true(nrow(registry) == 117L, "Frozen program registry must contain 117 programs")
assert_true(all(c(stromal_ecm, ductular) %in% registry$program_uid),
            "Named programs missing from the registry")
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))

score_substrate <- function(dge, meta, composition, label) {
  message("[", label, "] samples ", ncol(dge), ", genes ", nrow(dge))
  assert_true(setequal(colnames(dge), meta$sample_id), "DGE and metadata rosters differ")
  meta <- meta[match(colnames(dge), sample_id)]
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
  rm(logcpm)
  scores <- score_programs(symbol_matrix, meta, membership, registry, program_weight_coverage)
  composition <- composition[match(meta$sample_id, sample_id)]
  assert_true(identical(composition$sample_id, meta$sample_id),
              "Composition does not cover every record")
  lineages <- setdiff(names(composition)[vapply(composition, is.numeric, logical(1))],
                      c("sample_id", "dataset"))
  assert_true(length(lineages) == 16L, "Composition must carry 16 lineages")
  fractions <- as.matrix(composition[, ..lineages])
  zero_audit <- data.table(
    substrate = label, n_records = nrow(fractions), n_lineages = ncol(fractions),
    n_exact_zero = sum(fractions == 0), n_below_pseudocount = sum(fractions < composition_pseudocount),
    min_fraction = min(fractions), pseudocount = composition_pseudocount,
    max_abs_rowsum_minus_1 = max(abs(rowSums(fractions) - 1))
  )
  comp <- build_composition_scores(
    composition, meta, lineages, composition_pseudocount,
    eligible_lineages = character(0), n_components = composition_n_clr_components
  )
  list(meta = meta, primary = scores$primary, coverage = scores$coverage,
       cohort_coverage = scores$cohort_coverage, pcs = comp$pcs,
       pc_variance = comp$pc_variance, zero_audit = zero_audit)
}

fit_arm <- function(arm, S, meta, pcs, coverage, cohorts, unit_col) {
  message("[", arm, "] fitting")
  common <- c("sample_id", "dataset", "sex_final", "fibrosis_stage", "nas_score", unit_col)
  discovery <- meta[dataset %in% cohorts & !is.na(fibrosis_stage) & !is.na(nas_score), ..common]
  assert_true(all(!is.na(discovery$sex_final)), "Missing sex_final in the model roster")
  assert_true(uniqueN(discovery[[unit_col]]) == nrow(discovery), "More than one record per unit")
  composition_meta <- merge(discovery, pcs, by = "sample_id", all.x = TRUE, sort = FALSE)
  setkey(composition_meta, NULL)
  testable <- coverage[testable == TRUE, feature_id]
  S <- S[testable, , drop = FALSE]
  linear <- fit_feature_models(S, discovery, cohorts, "primary")
  adjusted <- fit_feature_models(S, composition_meta, cohorts, "composition_pc")
  cohort_effects <- rbindlist(list(linear, adjusted), use.names = TRUE, fill = TRUE)
  meta_effects <- rbindlist(list(meta_analyze(linear), meta_analyze(adjusted)),
                            use.names = TRUE, fill = TRUE)
  # Same multiplicity rule as run_discovery.R: one BH family per model, pooling
  # the fibrosis and NAS terms of every eligible program.
  meta_effects[, `:=`(
    multiplicity_family = paste("hotspot_program::weighted_mean_z", model_kind, sep = "::"),
    family_size = .N,
    q_value = p.adjust(p_value_meta, method = "BH", n = .N)
  ), by = model_kind]
  value_cols <- c("beta_meta", "se_meta", "ci_lower", "ci_upper", "p_value_meta",
                  "q_value", "tau2", "n_cohorts", "n_same_direction")
  linear_wide <- dcast(meta_effects[model_kind == "primary"], feature_id ~ axis,
                       value.var = value_cols)
  comp_wide <- dcast(meta_effects[model_kind == "composition_pc"], feature_id ~ axis,
                     value.var = value_cols)
  setnames(comp_wide, setdiff(names(comp_wide), "feature_id"),
           paste0("composition_pc_", setdiff(names(comp_wide), "feature_id")))
  axis_map <- Reduce(function(x, y) merge(x, y, by = "feature_id", all.x = TRUE),
                     list(registry[, .(feature_id = program_uid, cell_type, module, module_name)],
                          coverage[, .(feature_id, weight_coverage, testable)],
                          linear_wide, comp_wide))
  roster <- discovery[, .(n_units = .N, n_female = sum(sex_final == "F"),
                          n_male = sum(sex_final == "M"),
                          fibrosis_levels = paste(sort(unique(fibrosis_stage)), collapse = ","),
                          nas_levels = paste(sort(unique(nas_score)), collapse = ",")),
                      by = dataset]
  arm_dir <- file.path(analysis_root, arm)
  dir.create(arm_dir)
  write_tsv(cohort_effects, file.path(arm_dir, "cohort_effects.tsv"))
  write_tsv(meta_effects, file.path(arm_dir, "meta_effects.tsv"))
  write_tsv(axis_map, file.path(arm_dir, "program_axis_map.tsv"))
  write_tsv(roster, file.path(arm_dir, "model_roster_by_cohort.tsv"))
  write_tsv(discovery, file.path(arm_dir, "model_roster.tsv"))
  list(axis_map = axis_map, roster = roster, n = nrow(discovery), discovery = discovery,
       cohort_effects = cohort_effects, meta_effects = meta_effects)
}

# ---- old substrate -------------------------------------------------------
old_meta <- as.data.table(readRDS(nonholdout_meta))
old_meta[, unit := sample_id]
old <- score_substrate(readRDS(nonholdout_dge), old_meta,
                       fread(nonholdout_composition, check.names = FALSE), "old_1097_records")
gc()

# ---- new substrate -------------------------------------------------------
new_meta <- as.data.table(readRDS(new_meta_path))
manifest <- fread(unit_manifest_path)
assert_true(nrow(new_meta) == 774L && uniqueN(new_meta$analysis_unit_id) == 774L,
            "The adopted roster must be 774 units, one record each")
assert_true(setequal(manifest$sample_id, new_meta$sample_id), "Unit manifest roster differs")
mm <- manifest[match(new_meta$sample_id, sample_id)]
assert_true(identical(mm$fibrosis_stage, new_meta$fibrosis_stage) &&
              identical(mm$nas_score, new_meta$nas_score),
            "Manifest stage or NAS differs from the model metadata")
new_meta[, unit := analysis_unit_id]
new <- score_substrate(readRDS(new_dge_path), new_meta,
                       fread(new_composition_path, check.names = FALSE), "new_774_units")
gc()

write_tsv(rbind(old$zero_audit, new$zero_audit), file.path(analysis_root, "clr_zero_audit.tsv"))
write_tsv(rbind(cbind(substrate = "old_1097_records", old$pc_variance),
                cbind(substrate = "new_774_units", new$pc_variance)),
          file.path(analysis_root, "composition_pc_variance.tsv"))
write_tsv(merge(old$coverage[, .(feature_id, old_weight_coverage = weight_coverage,
                                 old_testable = testable)],
                new$coverage[, .(feature_id, new_weight_coverage = weight_coverage,
                                 new_testable = testable)], by = "feature_id"),
          file.path(analysis_root, "program_weight_coverage_old_vs_new.tsv"))

# ---- arms ----------------------------------------------------------------
arms <- list()
arms$old469_reproduction <- fit_arm("old469_reproduction", old$primary, old$meta, old$pcs,
                                    old$coverage, discovery_cohorts, "unit")
arms$old_three_cohorts <- fit_arm("old_three_cohorts", old$primary, old$meta, old$pcs,
                                  old$coverage, three_cohorts, "unit")
arms$new774_three_cohorts <- fit_arm("new774_three_cohorts", new$primary, new$meta, new$pcs,
                                     new$coverage, three_cohorts, "unit")

# Primary scores for the three 774-unit cohorts plus the unchanged GSE174478
# records scored in their own (old) expression matrix. Eligibility is the
# 774-unit weight coverage; a program not scored in GSE174478 is fitted on the
# remaining three cohorts.
extra_ids <- old$meta[dataset == extra_cohort, sample_id]
combined_scores <- cbind(
  new$primary,
  old$primary[rownames(new$primary), extra_ids, drop = FALSE]
)
keep <- c("sample_id", "dataset", "sex_final", "fibrosis_stage", "nas_score", "unit")
combined_meta <- rbind(new$meta[, ..keep], old$meta[dataset == extra_cohort, ..keep])
combined_pcs <- rbind(new$pcs, old$pcs[sample_id %in% extra_ids])
assert_true(identical(colnames(combined_scores), combined_meta$sample_id),
            "Combined score and metadata order mismatch")
arms$new774_plus_GSE174478 <- fit_arm(
  "new774_plus_GSE174478", combined_scores, combined_meta, combined_pcs,
  new$coverage, c(three_cohorts, extra_cohort), "unit"
)

# ---- reproduction of the stored map --------------------------------------
stored <- fread(stored_map_path)
repro <- merge(stored, arms$old469_reproduction$axis_map, by = "feature_id",
               suffixes = c("_stored", "_rerun"))
check_cols <- c("beta_meta_fibrosis", "beta_meta_nas", "q_value_fibrosis", "q_value_nas",
                "composition_pc_beta_meta_fibrosis", "composition_pc_beta_meta_nas",
                "composition_pc_q_value_fibrosis", "composition_pc_q_value_nas")
reproduction <- rbindlist(lapply(check_cols, function(col) {
  a <- repro[[paste0(col, "_stored")]]; b <- repro[[paste0(col, "_rerun")]]
  data.table(quantity = col, n_compared = sum(is.finite(a) & is.finite(b)),
             n_na_mismatch = sum(is.finite(a) != is.finite(b)),
             max_abs_difference = max(abs(a - b), na.rm = TRUE))
}))
reproduction[, reproduced := n_na_mismatch == 0L & max_abs_difference < 1e-8]
write_tsv(reproduction, file.path(analysis_root, "stored_map_reproduction.tsv"))
pc1b_assertion <- tryCatch({
  assert_true(all(reproduction$reproduced) &&
                identical(sort(stored[testable == TRUE, feature_id]),
                          sort(arms$old469_reproduction$axis_map[testable == TRUE, feature_id])),
              "The stored 469-record map was not reproduced")
  TRUE
}, error = function(e) FALSE)

# ===========================================================================
# H85 block: diagnostics 2 to 9, stop rules S2 to S5, earlier-set arms, NC-1, D3
# ===========================================================================
stored_fig4e_analysis <- file.path(repo, "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis")
pc_cols <- paste0("composition_pc", seq_len(composition_n_clr_components))
four_cohorts <- c(three_cohorts, extra_cohort)
dataset_of <- c(setNames(new$meta$dataset, new$meta$sample_id),
                setNames(old$meta[dataset == extra_cohort, dataset], extra_ids))

# ---- diagnostic 2 and stop rule S3: roster --------------------------------
run_roster <- arms$new774_plus_GSE174478$discovery
stored_roster <- fread(file.path(stored_fig4e_analysis, "new774_plus_GSE174478/model_roster.tsv"))
gate_roster <- readLines(file.path(repo, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z",
                                   "round4-20261015T112148Z/gate/R05-10/roster457.txt"))
expected_counts <- c(GSE135251 = 207L, GSE162694 = 86L, GSE174478 = 93L, GSE130970 = 71L)
run_counts <- run_roster[, .N, by = dataset]
run_counts_vec <- setNames(run_counts$N, run_counts$dataset)[names(expected_counts)]
roster_ok <- nrow(run_roster) == 457L && uniqueN(run_roster$sample_id) == 457L &&
  setequal(run_roster$sample_id, stored_roster$sample_id) &&
  setequal(run_roster$sample_id, gate_roster) &&
  identical(as.integer(run_counts_vec), unname(expected_counts))
add_diag("2", "model roster records", nrow(run_roster), 457L, nrow(run_roster) == 457L)
add_diag("2", "per-cohort records GSE135251 / GSE162694 / GSE174478 / GSE130970",
         as.integer(run_counts_vec), unname(expected_counts),
         identical(as.integer(run_counts_vec), unname(expected_counts)))
add_diag("2", "sample IDs identical to the stored model_roster.tsv and the gate roster",
         setequal(run_roster$sample_id, stored_roster$sample_id) && setequal(run_roster$sample_id, gate_roster),
         TRUE, setequal(run_roster$sample_id, stored_roster$sample_id) && setequal(run_roster$sample_id, gate_roster))

# ---- diagnostic 3 and stop rule S3: earlier-fraction coverage ---------------
earlier_composition <- fread(nonholdout_composition, check.names = FALSE)
rebuilt_composition <- fread(new_composition_path, check.names = FALSE)
lineage_cols <- function(x) setdiff(names(x)[vapply(x, is.numeric, logical(1))], c("sample_id", "dataset"))
lineages <- lineage_cols(rebuilt_composition)
lineages_ok <- identical(lineages, lineage_cols(earlier_composition)) && length(lineages) == 16L
add_diag("3", "16 lineages, identical names and order in both fraction files", length(lineages), 16L, lineages_ok)
pca_records <- new$meta[dataset %in% three_cohorts, sample_id]
model_records <- run_roster$sample_id
cov_pca <- sum(pca_records %in% earlier_composition$sample_id)
cov_model <- sum(model_records %in% earlier_composition$sample_id)
cov_774 <- sum(new$meta$sample_id %in% earlier_composition$sample_id)
add_diag("3", "PCA records of the three cohorts (count; covered by the earlier file)",
         c(length(pca_records), cov_pca), c(418L, 418L), length(pca_records) == 418L && cov_pca == 418L)
add_diag("3", "model records covered by the earlier file", cov_model, 457L, cov_model == 457L)
add_diag("3", "774-unit records covered by the earlier file (needed by the CLR call)", cov_774, 774L, cov_774 == 774L)
frac_matrix <- function(comp, ids) as.matrix(comp[match(ids, sample_id), ..lineages])
for (set_name in c("earlier", "rebuilt")) {
  comp <- if (set_name == "earlier") earlier_composition else rebuilt_composition
  for (record_set in c("pca_418", "model_457")) {
    ids <- if (record_set == "pca_418") pca_records else model_records
    fr <- frac_matrix(comp, ids)
    add_diag("3", paste0(set_name, " fractions, ", record_set, ": exact zeros / values below 1e-6"),
             c(sum(fr == 0, na.rm = TRUE), sum(fr < 1e-6, na.rm = TRUE)),
             if (record_set == "model_457") c(0L, if (set_name == "earlier") 1336L else 1404L) else c(0L, NA),
             sum(fr == 0, na.rm = TRUE) == 0L)
  }
}
if (!roster_ok || !lineages_ok || cov_pca != 418L || cov_model != 457L || cov_774 != 774L) {
  stop_rule("S3", "model roster or earlier-fraction coverage differs from the prespecification (diagnostics 2 and 3)")
}

# ---- diagnostic 4 and stop rule S2: PC-1b and PC-1a --------------------------
old469_map <- arms$old469_reproduction$axis_map
pc1b_counts <- c(sum(old469_map$composition_pc_q_value_fibrosis < 0.05, na.rm = TRUE),
                 sum(old469_map$composition_pc_q_value_nas < 0.05, na.rm = TRUE))
pc1b_ok <- isTRUE(pc1b_assertion) && all(pc1b_counts == 0L)
add_diag("4", "PC-1b stored 469-record map, max abs difference over 8 columns",
         max(reproduction$max_abs_difference), "< 1e-8, no NA mismatch, same eligible set", pc1b_assertion)
add_diag("4", "PC-1b old469 adjusted calls fibrosis / NAS", pc1b_counts, c(0L, 0L), all(pc1b_counts == 0L))

stored_457_map <- fread(file.path(stored_fig4e_analysis, "new774_plus_GSE174478/program_axis_map.tsv"))
run_457_map <- arms$new774_plus_GSE174478$axis_map
pc1a_base <- as.vector(outer(c("beta_meta", "se_meta", "q_value"), c("fibrosis", "nas"), paste, sep = "_"))
pc1a_cols <- c(pc1a_base, paste0("composition_pc_", pc1a_base))
pc1a_merge <- merge(stored_457_map[testable == TRUE], run_457_map[testable == TRUE],
                    by = "feature_id", suffixes = c("_stored", "_rerun"))
pc1a <- rbindlist(lapply(pc1a_cols, function(col) {
  a <- pc1a_merge[[paste0(col, "_stored")]]; b <- pc1a_merge[[paste0(col, "_rerun")]]
  data.table(quantity = col, n_compared = sum(is.finite(a) & is.finite(b)),
             n_na_mismatch = sum(is.finite(a) != is.finite(b)),
             max_abs_difference = max(abs(a - b), na.rm = TRUE))
}))
pc1a[, reproduced := n_compared == 112L & n_na_mismatch == 0L & max_abs_difference < 1e-8]
write_tsv(pc1a, file.path(out_root, "pc1a_stored_457_reproduction.tsv"))
pc1a_same_set <- identical(sort(stored_457_map[testable == TRUE, feature_id]),
                           sort(run_457_map[testable == TRUE, feature_id])) &&
  sum(run_457_map$testable) == 112L
count_q05 <- function(map, prefix, axis) sum(map[[paste0(prefix, "q_value_", axis)]] < 0.05, na.rm = TRUE)
pc1a_counts <- c(count_q05(run_457_map, "", "fibrosis"), count_q05(run_457_map, "", "nas"),
                 count_q05(run_457_map, "composition_pc_", "fibrosis"),
                 count_q05(run_457_map, "composition_pc_", "nas"))
pc1a_ok <- all(pc1a$reproduced) && pc1a_same_set && identical(pc1a_counts, c(16L, 1L, 3L, 0L))
add_diag("4", "PC-1a stored 457 arm, max abs difference over 12 columns (112 programs)",
         max(pc1a$max_abs_difference), "< 1e-8, no NA mismatch, same 112 eligible",
         all(pc1a$reproduced) && pc1a_same_set)
add_diag("4", "PC-1a calls unadjusted fibrosis / NAS, adjusted fibrosis / NAS", pc1a_counts,
         c(16L, 1L, 3L, 0L), identical(pc1a_counts, c(16L, 1L, 3L, 0L)))
if (!pc1a_ok || !pc1b_ok) {
  stop_rule("S2", paste0("PC-1a ", if (pc1a_ok) "passed" else "failed", "; PC-1b ",
                         if (pc1b_ok) "passed" else "failed"))
}

# ---- earlier-set CLR components on the 774-unit metadata ---------------------
if (h85_debug) {
  # Null data for debugging: the earlier fraction rows are shuffled among records
  # within each dataset, so no record keeps its own earlier fractions.
  set.seed(1L)
  earlier_composition[, sample_id := sample_id[sample.int(.N)], by = dataset]
}
comp_earlier <- build_composition_scores(
  earlier_composition, new$meta, lineages, composition_pseudocount,
  eligible_lineages = character(0), n_components = composition_n_clr_components
)
pcs_earlier <- comp_earlier$pcs
pcs_earlier_three <- as.matrix(pcs_earlier[sample_id %in% pca_records, ..pc_cols])
if (any(!is.finite(pcs_earlier_three)) || nrow(pcs_earlier_three) != 418L) {
  stop_rule("S5", "earlier-set CLR components are NA for at least one PCA record of the three cohorts")
}
write_tsv(rbind(cbind(set = "earlier", comp_earlier$pc_variance),
                cbind(set = "rebuilt", new$pc_variance)),
          file.path(out_root, "pc_variance_by_set.tsv"))

# ---- earlier-set arms -------------------------------------------------------
combined_pcs_earlier <- rbind(pcs_earlier, old$pcs[sample_id %in% extra_ids])
arms$new774_plus_GSE174478_earlierfrac <- fit_arm(
  "new774_plus_GSE174478_earlierfrac", combined_scores, combined_meta, combined_pcs_earlier,
  new$coverage, four_cohorts, "unit"
)
arms$new774_three_cohorts_earlierfrac <- fit_arm(
  "new774_three_cohorts_earlierfrac", new$primary, new$meta, pcs_earlier,
  new$coverage, three_cohorts, "unit"
)
d1_pcs <- copy(old$pcs)
if (h85_debug) {
  d1_dataset <- old$meta$dataset[match(d1_pcs$sample_id, old$meta$sample_id)]
  for (cohort in three_cohorts) {
    rows <- which(d1_dataset == cohort)
    shuffled <- rows[sample.int(length(rows))]
    for (col in pc_cols) set(d1_pcs, i = rows, j = col, value = d1_pcs[[col]][shuffled])
  }
}
arms$D1_new774_plus_GSE174478_earlier_pcs_as_computed <- fit_arm(
  "D1_new774_plus_GSE174478_earlier_pcs_as_computed", combined_scores, combined_meta, d1_pcs,
  new$coverage, four_cohorts, "unit"
)

# ---- diagnostic 5 and stop rule S4: invariance -------------------------------
arm_r <- arms$new774_plus_GSE174478
arm_e <- arms$new774_plus_GSE174478_earlierfrac
max_diff_table <- function(x, y, keys, cols) {
  m <- merge(x, y, by = keys, suffixes = c("_r", "_e"))
  diffs <- vapply(cols, function(col) {
    a <- m[[paste0(col, "_r")]]; b <- m[[paste0(col, "_e")]]
    if (any(is.finite(a) != is.finite(b))) return(Inf)
    d <- abs(a - b)[is.finite(a)]
    if (length(d)) max(d) else 0
  }, numeric(1))
  list(n_x = nrow(x), n_y = nrow(y), n_merged = nrow(m), max_diff = max(diffs))
}
inv_unadj_cohort <- max_diff_table(arm_r$cohort_effects[model_kind == "primary"],
                                   arm_e$cohort_effects[model_kind == "primary"],
                                   c("feature_id", "cohort", "axis"),
                                   c("beta", "se_hc3", "statistic_hc3", "p_value"))
inv_unadj_meta <- max_diff_table(arm_r$meta_effects[model_kind == "primary"],
                                 arm_e$meta_effects[model_kind == "primary"],
                                 c("feature_id", "axis"),
                                 c("beta_meta", "se_meta", "p_value_meta", "q_value", "tau2"))
inv_174478 <- max_diff_table(arm_r$cohort_effects[model_kind == "composition_pc" & cohort == extra_cohort],
                             arm_e$cohort_effects[model_kind == "composition_pc" & cohort == extra_cohort],
                             c("feature_id", "cohort", "axis"),
                             c("beta", "se_hc3", "statistic_hc3", "p_value"))
inv_ok <- function(x) x$n_x == x$n_y && x$n_x == x$n_merged && x$n_x > 0L && x$max_diff < 1e-12
add_diag("5", "unadjusted cohort effects, earlier vs same-run rebuilt arm: rows / max abs difference",
         c(inv_unadj_cohort$n_merged, inv_unadj_cohort$max_diff), "896 rows; < 1e-12", inv_ok(inv_unadj_cohort))
add_diag("5", "unadjusted meta effects, earlier vs same-run rebuilt arm: rows / max abs difference",
         c(inv_unadj_meta$n_merged, inv_unadj_meta$max_diff), "224 rows; < 1e-12", inv_ok(inv_unadj_meta))
add_diag("5", "GSE174478 adjusted cohort effects, earlier vs same-run rebuilt arm: rows / max abs difference",
         c(inv_174478$n_merged, inv_174478$max_diff), "224 rows; < 1e-12", inv_ok(inv_174478))
if (!inv_ok(inv_unadj_cohort) || !inv_ok(inv_unadj_meta) || !inv_ok(inv_174478)) {
  stop_rule("S4", "unadjusted results or GSE174478 adjusted results differ between the earlier and rebuilt arms")
}

# ---- diagnostic 6: PCs under each set ----------------------------------------
pcs_rebuilt_all <- combined_pcs
variance_for <- function(pc_variance, ch) {
  v <- pc_variance[pc_variance$cohort == ch, ]
  v$proportion_variance[match(pc_cols, v$component)]
}
pc_compare <- rbindlist(lapply(four_cohorts, function(ch) {
  ids <- pcs_rebuilt_all$sample_id[dataset_of[pcs_rebuilt_all$sample_id] == ch]
  A <- as.matrix(combined_pcs_earlier[match(ids, sample_id), ..pc_cols])
  B <- as.matrix(pcs_rebuilt_all[match(ids, sample_id), ..pc_cols])
  cc <- stats::cancor(A, B)$cor
  r <- vapply(seq_along(pc_cols), function(j) stats::cor(A[, j], B[, j]), numeric(1))
  var_e <- variance_for(if (ch == extra_cohort) old$pc_variance else comp_earlier$pc_variance, ch)
  var_r <- variance_for(if (ch == extra_cohort) old$pc_variance else new$pc_variance, ch)
  data.table(cohort = ch, n_pca_records = length(ids), component = pc_cols,
             variance_explained_earlier = var_e, variance_explained_rebuilt = var_r,
             canonical_correlation = cc, pearson_r_raw = r, pearson_r_sign_aligned = abs(r),
             sign_flipped = r < 0)
}))
write_tsv(pc_compare, file.path(out_root, "diagnostic6_pc_comparison.tsv"))
for (ch in four_cohorts) {
  pc <- pc_compare[pc_compare$cohort == ch]
  add_diag("6", paste0(ch, ": PCA records; variance PC1-3 earlier | rebuilt; canonical r; sign-aligned r PC1-3"),
           c(pc$n_pca_records[1], round(pc$variance_explained_earlier, 4), round(pc$variance_explained_rebuilt, 4),
             round(pc$canonical_correlation, 4), round(pc$pearson_r_sign_aligned, 4)), NA, NA)
}

# ---- diagnostic 7: axis variance absorbed by each set's PCs -----------------
r2_rows <- rbindlist(lapply(four_cohorts, function(ch) {
  md <- run_roster[run_roster$dataset == ch]
  rbindlist(lapply(c("earlier", "rebuilt"), function(set_name) {
    pcs <- if (set_name == "earlier") combined_pcs_earlier else pcs_rebuilt_all
    d <- merge(md, pcs, by = "sample_id")
    data.table(cohort = ch, set = set_name, n_model_records = nrow(d),
               r2_fibrosis = summary(stats::lm(fibrosis_stage ~ composition_pc1 + composition_pc2 + composition_pc3, data = d))$r.squared,
               r2_nas = summary(stats::lm(nas_score ~ composition_pc1 + composition_pc2 + composition_pc3, data = d))$r.squared)
  }))
}))
write_tsv(r2_rows, file.path(out_root, "diagnostic7_pc_r2_on_axes.tsv"))
for (i in seq_len(nrow(r2_rows))) {
  add_diag("7", paste0(r2_rows$cohort[i], " ", r2_rows$set[i], ": n; R2 fibrosis; R2 NAS"),
           c(r2_rows$n_model_records[i], round(r2_rows$r2_fibrosis[i], 4), round(r2_rows$r2_nas[i], 4)), NA, NA)
}

# ---- diagnostic 8 and stop rule S5: eligible set, family, estimability -------
h85_arm_names <- c("new774_plus_GSE174478", "new774_plus_GSE174478_earlierfrac",
                   "new774_three_cohorts", "new774_three_cohorts_earlierfrac",
                   "D1_new774_plus_GSE174478_earlier_pcs_as_computed")
family_ok_all <- TRUE
for (arm in h85_arm_names) {
  a <- arms[[arm]]
  fam <- a$meta_effects[, .(family_size = unique(family_size)), by = model_kind]
  ok <- sum(a$axis_map$testable) == 112L && nrow(fam) == 2L && all(fam$family_size == 224L)
  family_ok_all <- family_ok_all && ok
  add_diag("8", paste0(arm, ": eligible programs; BH family size primary / composition_pc"),
           c(sum(a$axis_map$testable), fam[model_kind == "primary", family_size],
             fam[model_kind == "composition_pc", family_size]), c(112L, 224L, 224L), ok)
}
earlier_estimable <- arm_e$meta_effects[model_kind == "composition_pc" & estimable == TRUE, .N, by = axis]
earlier_estimable_vec <- setNames(earlier_estimable$N, earlier_estimable$axis)[c("fibrosis", "nas")]
add_diag("8", "earlier-set adjusted estimable meta-analyses fibrosis / NAS",
         as.integer(earlier_estimable_vec), c(112L, 112L), all(earlier_estimable_vec == 112L, na.rm = FALSE))
if (!family_ok_all || anyNA(earlier_estimable_vec) || any(earlier_estimable_vec < 112L) ||
    sum(new$coverage$testable) != 112L) {
  stop_rule("S5", "eligible set is not 112, a family is not 224, or the earlier-set adjusted model has fewer than 112 estimable meta-analyses on an axis")
}

# ---- NC-1 (diagnostic 9): shuffled earlier-set components --------------------
nc1_base <- copy(combined_pcs_earlier)
nc1_dataset <- unname(dataset_of[nc1_base$sample_id])
nc1_order <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478")
nc1_discovery <- arm_e$discovery
nc1_S <- combined_scores[new$coverage[testable == TRUE, feature_id], , drop = FALSE]
nc1_unadjusted_fib <- arm_e$meta_effects[model_kind == "primary" & axis == "fibrosis",
                                         .(feature_id, beta_unadjusted = beta_meta)]
nc1_indices <- vector("list", nc1_draws)
nc1_rows <- vector("list", nc1_draws)
set.seed(h85_seed)
for (b in seq_len(nc1_draws)) {
  perm <- copy(nc1_base)
  idx_b <- list()
  for (cohort in nc1_order) {
    rows <- which(nc1_dataset == cohort)
    idx <- sample.int(length(rows))
    for (col in pc_cols) set(perm, i = rows, j = col, value = nc1_base[[col]][rows[idx]])
    idx_b[[cohort]] <- idx
  }
  nc1_indices[[b]] <- idx_b
  cm <- merge(nc1_discovery, perm, by = "sample_id", all.x = TRUE, sort = FALSE)
  setkey(cm, NULL)
  adj <- fit_feature_models(nc1_S, cm, four_cohorts, "composition_pc")
  me <- meta_analyze(adj)
  me[, q_value := p.adjust(p_value_meta, method = "BH", n = .N)]
  fib <- merge(me[axis == "fibrosis", .(feature_id, beta_meta, q_value)], nc1_unadjusted_fib, by = "feature_id")
  ok <- is.finite(fib$beta_meta) & is.finite(fib$beta_unadjusted)
  nc1_rows[[b]] <- data.table(
    draw = b, family_size = nrow(me),
    n_estimable_fibrosis = sum(me[axis == "fibrosis", estimable]),
    n_estimable_nas = sum(me[axis == "nas", estimable]),
    n_q05_fibrosis = sum(me[axis == "fibrosis", q_value] < 0.05, na.rm = TRUE),
    n_q05_nas = sum(me[axis == "nas", q_value] < 0.05, na.rm = TRUE),
    pearson_r_fibrosis_vs_unadjusted = stats::cor(fib$beta_meta[ok], fib$beta_unadjusted[ok])
  )
  if (b %% 10L == 0L || b == nc1_draws) message("[NC-1] draw ", b, " of ", nc1_draws)
}
nc1 <- rbindlist(nc1_rows)
write_tsv(nc1, file.path(out_root, "nc1_draws.tsv"))
saveRDS(list(seed = h85_seed, rng_kind = RNGkind(), cohort_order = nc1_order,
             sample_order = nc1_base$sample_id, dataset = nc1_dataset, indices = nc1_indices),
        file.path(out_root, "nc1_permutation_indices.rds"))
nc1_median_count <- stats::median(nc1$n_q05_fibrosis)
nc1_median_r <- stats::median(nc1$pearson_r_fibrosis_vs_unadjusted)
nc1_failed <- !(isTRUE(nc1_median_count >= 10) && isTRUE(nc1_median_r >= 0.95))
add_diag("9", "NC-1 draws; all families 224", c(nrow(nc1), all(nc1$family_size == 224L)),
         c(nc1_draws, TRUE), nrow(nc1) == nc1_draws && all(nc1$family_size == 224L))
add_diag("9", "NC-1 median adjusted fibrosis calls (range)",
         c(nc1_median_count, min(nc1$n_q05_fibrosis), max(nc1$n_q05_fibrosis)),
         "12 to 18 predicted; fail if < 10", nc1_median_count >= 10)
add_diag("9", "NC-1 median NAS calls (range)",
         c(stats::median(nc1$n_q05_nas), min(nc1$n_q05_nas), max(nc1$n_q05_nas)), NA, NA)
add_diag("9", "NC-1 median Pearson r shuffled-adjusted vs unadjusted fibrosis betas (range)",
         c(nc1_median_r, min(nc1$pearson_r_fibrosis_vs_unadjusted), max(nc1$pearson_r_fibrosis_vs_unadjusted)),
         ">= 0.95", nc1_median_r >= 0.95)
add_diag("9", "NC-1 draws with 3 or fewer adjusted fibrosis calls", sum(nc1$n_q05_fibrosis <= 3L), NA, NA)
add_diag("9", "NC-1 outcome", if (nc1_failed) "FAILED" else "PASSED", "PASSED", !nc1_failed)
message("H85 diagnostics 1 to 9 written: ", diagnostics_path)

# ---- D3: test variants for both models under both fraction sets --------------
meta_analyze_test <- function(cohort_effects, test, min_cohorts = 3L) {
  groups <- unique(cohort_effects[, .(feature_id, axis, model_kind)])
  rbindlist(lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i]
    d <- cohort_effects[feature_id == g$feature_id & axis == g$axis & model_kind == g$model_kind &
                          estimable == TRUE & is.finite(beta) & is.finite(se_hc3) & se_hc3 > 0]
    if (nrow(d) < min_cohorts) {
      return(data.table(feature_id = g$feature_id, axis = g$axis, model_kind = g$model_kind,
                        estimable = FALSE, failure_reason = "fewer_than_three_cohorts"))
    }
    fit <- tryCatch(metafor::rma(yi = d$beta, sei = d$se_hc3, method = "REML", test = test),
                    error = function(e) e)
    if (inherits(fit, "error")) {
      return(data.table(feature_id = g$feature_id, axis = g$axis, model_kind = g$model_kind,
                        estimable = FALSE, failure_reason = conditionMessage(fit)))
    }
    data.table(feature_id = g$feature_id, axis = g$axis, model_kind = g$model_kind,
               estimable = TRUE, failure_reason = NA_character_, beta_meta = as.numeric(fit$b),
               se_meta = fit$se, p_value_meta = fit$pval)
  }), use.names = TRUE, fill = TRUE)
}
d3_rows <- list()
d3_tables <- list()
for (set_name in c("rebuilt", "earlier")) {
  arm <- if (set_name == "rebuilt") arm_r else arm_e
  for (model in c("primary", "composition_pc")) {
    for (test in c("knha", "adhoc", "z")) {
      me <- meta_analyze_test(arm$cohort_effects[model_kind == model], test)
      if (!"p_value_meta" %in% names(me)) me[, p_value_meta := NA_real_]
      me[, q_value := p.adjust(p_value_meta, method = "BH", n = .N)]
      me[, `:=`(set = set_name, test = test)]
      d3_tables[[length(d3_tables) + 1L]] <- me
      d3_rows[[length(d3_rows) + 1L]] <- data.table(
        set = set_name, model = model, test = test, family_size = nrow(me),
        n_estimable = sum(me$estimable),
        n_q05_fibrosis = sum(me[axis == "fibrosis", q_value] < 0.05, na.rm = TRUE),
        n_q05_nas = sum(me[axis == "nas", q_value] < 0.05, na.rm = TRUE),
        min_q_fibrosis = suppressWarnings(min(me[axis == "fibrosis", q_value], na.rm = TRUE)),
        failure_reasons = paste(unique(na.omit(me$failure_reason)), collapse = "; "))
    }
  }
}
d3 <- rbindlist(d3_rows)
write_tsv(d3, file.path(out_root, "D3_test_variant_counts.tsv"))
write_tsv(rbindlist(d3_tables, use.names = TRUE, fill = TRUE), file.path(out_root, "D3_test_variant_effects.tsv"))
# knha in D3 must equal the producer's meta_analyze, a check on meta_analyze_test.
d3_knha_check <- c(
  d3[set == "rebuilt" & model == "primary" & test == "knha", c(n_q05_fibrosis, n_q05_nas)],
  d3[set == "rebuilt" & model == "composition_pc" & test == "knha", c(n_q05_fibrosis, n_q05_nas)])

# ---- summaries -----------------------------------------------------------
count_calls <- function(map, prefix, axis) {
  q <- map[[paste0(prefix, "q_value_", axis)]]
  beta <- map[[paste0(prefix, "beta_meta_", axis)]]
  list(n_estimable = sum(is.finite(q)), n_q05 = sum(q < 0.05, na.rm = TRUE),
       n_q05_positive = sum(q < 0.05 & beta > 0, na.rm = TRUE),
       n_q05_negative = sum(q < 0.05 & beta < 0, na.rm = TRUE),
       n_nominal_p05 = sum(map[[paste0(prefix, "p_value_meta_", axis)]] < 0.05, na.rm = TRUE),
       min_q = min(q, na.rm = TRUE),
       median_abs_beta = median(abs(beta), na.rm = TRUE),
       median_se = median(map[[paste0(prefix, "se_meta_", axis)]], na.rm = TRUE))
}
summary_rows <- rbindlist(lapply(names(arms), function(arm) {
  map <- arms[[arm]]$axis_map
  rbindlist(lapply(c("", "composition_pc_"), function(prefix) {
    rbindlist(lapply(c("fibrosis", "nas"), function(axis) {
      c(list(arm = arm, n_units = arms[[arm]]$n, n_cohorts = nrow(arms[[arm]]$roster),
             model = if (prefix == "") "unadjusted" else "composition_adjusted",
             axis = axis, n_registry = nrow(map), n_eligible = sum(map$testable),
             bh_family_size = 2L * sum(map$testable)),
        count_calls(map, prefix, axis))
    }))
  }))
}))
write_tsv(summary_rows, file.path(analysis_root, "arm_summary.tsv"))
write_tsv(rbindlist(lapply(names(arms), function(arm) cbind(arm = arm, arms[[arm]]$roster))),
          file.path(analysis_root, "arm_rosters.tsv"))

compare <- function(reference, target) {
  a <- arms[[reference]]$axis_map; b <- arms[[target]]$axis_map
  m <- merge(a, b, by = "feature_id", suffixes = c("_ref", "_new"))
  rbindlist(lapply(c("", "composition_pc_"), function(prefix) {
    rbindlist(lapply(c("fibrosis", "nas"), function(axis) {
      x <- m[[paste0(prefix, "beta_meta_", axis, "_ref")]]
      y <- m[[paste0(prefix, "beta_meta_", axis, "_new")]]
      qx <- m[[paste0(prefix, "q_value_", axis, "_ref")]]
      qy <- m[[paste0(prefix, "q_value_", axis, "_new")]]
      ok <- is.finite(x) & is.finite(y)
      data.table(
        reference = reference, target = target,
        model = if (prefix == "") "unadjusted" else "composition_adjusted", axis = axis,
        n_programs_both = sum(ok),
        n_same_sign = sum(sign(x[ok]) == sign(y[ok])),
        sign_agreement = mean(sign(x[ok]) == sign(y[ok])),
        pearson_r = stats::cor(x[ok], y[ok]),
        spearman_rho = stats::cor(x[ok], y[ok], method = "spearman"),
        slope_new_on_ref = unname(stats::coef(stats::lm(y[ok] ~ x[ok]))[2]),
        n_q05_ref = sum(qx[ok] < 0.05), n_q05_new = sum(qy[ok] < 0.05),
        n_q05_both = sum(qx[ok] < 0.05 & qy[ok] < 0.05),
        n_q05_ref_same_sign_in_new = sum(qx[ok] < 0.05 & sign(x[ok]) == sign(y[ok]))
      )
    }))
  }))
}
comparisons <- rbindlist(list(
  compare("old469_reproduction", "new774_three_cohorts"),
  compare("old469_reproduction", "new774_plus_GSE174478"),
  compare("old469_reproduction", "old_three_cohorts"),
  compare("old_three_cohorts", "new774_three_cohorts")
))
write_tsv(comparisons, file.path(analysis_root, "old_vs_new_agreement.tsv"))

named <- rbindlist(lapply(names(arms), function(arm) {
  me <- fread(file.path(analysis_root, arm, "meta_effects.tsv"))
  me[feature_id %in% c(stromal_ecm, ductular),
     .(arm = arm, program = fifelse(feature_id == stromal_ecm, "Stromal ECM", "Ductular injury"),
       feature_id, model = fifelse(model_kind == "primary", "unadjusted", "composition_adjusted"),
       axis, beta_meta, ci_lower, ci_upper, p_value_meta, q_value, family_size,
       n_cohorts, n_same_direction, tau2, i2)]
}))
setorder(named, program, axis, model, arm)
write_tsv(named, file.path(analysis_root, "stromal_ecm_ductular_by_arm.tsv"))
named_cohort <- rbindlist(lapply(names(arms), function(arm) {
  ce <- fread(file.path(analysis_root, arm, "cohort_effects.tsv"))
  ce[feature_id %in% c(stromal_ecm, ductular),
     .(arm = arm, program = fifelse(feature_id == stromal_ecm, "Stromal ECM", "Ductular injury"),
       cohort, model_kind, axis, beta, se_hc3, ci_lower, ci_upper, p_value, n)]
}))
write_tsv(named_cohort, file.path(analysis_root, "stromal_ecm_ductular_cohort_effects.tsv"))

per_program <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), lapply(
  names(arms), function(arm) {
    cols <- c("feature_id", grep("^(composition_pc_)?(beta_meta|ci_lower|ci_upper|p_value_meta|q_value|n_cohorts)_",
                                 names(arms[[arm]]$axis_map), value = TRUE), "testable")
    d <- arms[[arm]]$axis_map[, ..cols]
    setnames(d, setdiff(names(d), "feature_id"), paste0(arm, "::", setdiff(names(d), "feature_id")))
    d
  }))
per_program <- merge(registry[, .(feature_id = program_uid, cell_type, module, module_name)],
                     per_program, by = "feature_id")
write_tsv(per_program, file.path(analysis_root, "per_program_all_arms.tsv"))

inputs <- c(old_dge = nonholdout_dge, old_meta = nonholdout_meta,
            old_composition = nonholdout_composition, new_dge = new_dge_path,
            new_meta = new_meta_path, new_composition = new_composition_path,
            unit_manifest = unit_manifest_path, program_registry = program_registry,
            program_membership = program_membership, gene_annotation = gene_annotation,
            stored_program_axis_map = stored_map_path,
            analysis_lib = file.path(lib_dir, "analysis_lib.R"),
            config = file.path(lib_dir, "config.R"))
write_tsv(data.table(input = names(inputs), path = inputs,
                     sha256 = vapply(inputs, sha256_file, character(1))),
          file.path(analysis_root, "input_manifest.tsv"))
writeLines(c(paste("seed:", h85_seed, "(set at start and again immediately before the NC-1 loop; RNGkind",
                   paste(RNGkind(), collapse = " / "), ")"),
             capture.output(sessionInfo())), file.path(analysis_root, "sessionInfo.txt"))
print(summary_rows)
print(comparisons)
print(named)
print(reproduction)
cat("FIG4E_774_ANALYSIS_COMPLETE\n")

# ===========================================================================
# H85 primary output P1, decision, P2, PC-2, PC-3, D4, D5, RESULT_H85.md
# ===========================================================================
map_r <- arm_r$axis_map[testable == TRUE]
map_e <- arm_e$axis_map[testable == TRUE]
call_set <- function(map, axis) {
  q <- map[[paste0("composition_pc_q_value_", axis)]]
  b <- map[[paste0("composition_pc_beta_meta_", axis)]]
  sel <- which(q < 0.05)
  out <- setNames(as.integer(sign(b[sel])), map$feature_id[sel])
  out[order(names(out))]
}
stored_fib_calls <- c(hotspot_cholangiocytes_576c469934529fdd = 1L,
                      hotspot_hepatocytes_48f39dd4d817a10e = 1L,
                      hotspot_hepatocytes_676d360db6c1efb1 = 1L)
calls <- list(
  rebuilt_fibrosis = call_set(map_r, "fibrosis"), rebuilt_nas = call_set(map_r, "nas"),
  earlier_fibrosis = call_set(map_e, "fibrosis"), earlier_nas = call_set(map_e, "nas"))
rebuilt_matches_stored <- identical(calls$rebuilt_fibrosis, stored_fib_calls[order(names(stored_fib_calls))]) &&
  length(calls$rebuilt_nas) == 0L
same_fibrosis <- setequal(names(calls$earlier_fibrosis), names(stored_fib_calls)) &&
  all(calls$earlier_fibrosis[names(stored_fib_calls)] == stored_fib_calls)
stable <- same_fibrosis && length(calls$earlier_nas) == 0L
write_tsv(rbindlist(lapply(names(calls), function(k) {
  if (!length(calls[[k]])) return(data.table(set_axis = k, feature_id = NA_character_, sign = NA_integer_))
  data.table(set_axis = k, feature_id = names(calls[[k]]), sign = unname(calls[[k]]))
})), file.path(out_root, "P1_adjusted_call_sets.tsv"))

# P2 (co-reported)
signed_ids <- function(x) if (length(x)) paste(names(x), x, sep = ":") else character(0)
p2 <- rbindlist(lapply(c("fibrosis", "nas"), function(axis) {
  m <- merge(map_r[, .(feature_id, br = get(paste0("composition_pc_beta_meta_", axis)),
                       sr = get(paste0("composition_pc_se_meta_", axis)),
                       qr = get(paste0("composition_pc_q_value_", axis)))],
             map_e[, .(feature_id, be = get(paste0("composition_pc_beta_meta_", axis)),
                       qe = get(paste0("composition_pc_q_value_", axis)))], by = "feature_id")
  ok <- is.finite(m$br) & is.finite(m$be)
  d <- abs(m$be - m$br)[ok]
  A <- signed_ids(calls[[paste0("rebuilt_", axis)]]); B <- signed_ids(calls[[paste0("earlier_", axis)]])
  data.table(axis = axis, n_programs = sum(ok),
             pearson_r = stats::cor(m$br[ok], m$be[ok]),
             spearman_rho = stats::cor(m$br[ok], m$be[ok], method = "spearman"),
             median_abs_beta_difference = stats::median(d), max_abs_beta_difference = max(d),
             median_abs_difference_in_rebuilt_se = stats::median(d / m$sr[ok]),
             max_abs_difference_in_rebuilt_se = max(d / m$sr[ok]),
             n_called_rebuilt_only = sum(m$qr[ok] < 0.05 & !(m$qe[ok] < 0.05)),
             n_called_earlier_only = sum(m$qe[ok] < 0.05 & !(m$qr[ok] < 0.05)),
             n_calls_rebuilt = length(A), n_calls_earlier = length(B),
             jaccard_signed_call_sets = if (length(union(A, B))) length(intersect(A, B)) / length(union(A, B)) else NA_real_)
}))
write_tsv(p2, file.path(out_root, "P2_between_set_agreement.tsv"))

# PC-2 and PC-3
duct_e <- map_e[feature_id == ductular]
duct_r <- map_r[feature_id == ductular]
pc2_pass <- isTRUE(duct_e$composition_pc_beta_meta_fibrosis > 0 && duct_e$composition_pc_ci_lower_fibrosis > 0)
ecm_e <- map_e[feature_id == stromal_ecm]
ecm_r <- map_r[feature_id == stromal_ecm]
pc3_pass <- isTRUE(ecm_e$composition_pc_beta_meta_fibrosis > ecm_r$composition_pc_beta_meta_fibrosis)

# D4 focal table and D5 per-program table
focal <- c(names(stored_fib_calls), stromal_ecm)
focal_label <- c(hotspot_cholangiocytes_576c469934529fdd = "Rho-GTPase / adhesion (cholangiocytes)",
                 hotspot_hepatocytes_48f39dd4d817a10e = "Ductular injury (BICC1)",
                 hotspot_hepatocytes_676d360db6c1efb1 = "Hepatocyte stress (SLC38A1)",
                 hotspot_hepatocytes_f05c535ae5bbc0b9 = "Stromal ECM (IGFBP7)")
d4 <- rbindlist(lapply(c("rebuilt", "earlier"), function(set_name) {
  arm <- if (set_name == "rebuilt") arm_r else arm_e
  me <- arm$meta_effects[feature_id %in% focal]
  ce <- dcast(arm$cohort_effects[feature_id %in% focal], feature_id + model_kind + axis ~ cohort,
              value.var = "beta")
  setnames(ce, intersect(four_cohorts, names(ce)), paste0("beta_", intersect(four_cohorts, names(ce))))
  out <- merge(me[, .(feature_id, model_kind, axis, beta_meta, se_meta, ci_lower, ci_upper, p_value_meta,
                      q_value, tau2, i2, q_heterogeneity, n_cohorts, n_same_direction)],
               ce, by = c("feature_id", "model_kind", "axis"))
  cbind(set = set_name, program = focal_label[out$feature_id], out)
}))
setorder(d4, program, axis, model_kind, set)
write_tsv(d4, file.path(out_root, "D4_focal_programs.tsv"))
d5 <- rbindlist(lapply(c("fibrosis", "nas"), function(axis) {
  m <- merge(map_r, map_e, by = "feature_id", suffixes = c("_r", "_e"))
  data.table(feature_id = m$feature_id, cell_type = m$cell_type_r, module_name = m$module_name_r, axis = axis,
             unadjusted_beta = m[[paste0("beta_meta_", axis, "_r")]],
             unadjusted_se = m[[paste0("se_meta_", axis, "_r")]],
             unadjusted_q = m[[paste0("q_value_", axis, "_r")]],
             adjusted_beta_rebuilt = m[[paste0("composition_pc_beta_meta_", axis, "_r")]],
             adjusted_se_rebuilt = m[[paste0("composition_pc_se_meta_", axis, "_r")]],
             adjusted_q_rebuilt = m[[paste0("composition_pc_q_value_", axis, "_r")]],
             adjusted_beta_earlier = m[[paste0("composition_pc_beta_meta_", axis, "_e")]],
             adjusted_se_earlier = m[[paste0("composition_pc_se_meta_", axis, "_e")]],
             adjusted_q_earlier = m[[paste0("composition_pc_q_value_", axis, "_e")]])
}))
d5[, `:=`(beta_difference_earlier_minus_rebuilt = adjusted_beta_earlier - adjusted_beta_rebuilt,
          difference_in_rebuilt_se = (adjusted_beta_earlier - adjusted_beta_rebuilt) / adjusted_se_rebuilt)]
write_tsv(d5, file.path(out_root, "D5_per_program_both_sets.tsv"))

# Decision and reading, computed by the prespecified rule
n_e_fib <- length(calls$earlier_fibrosis); n_e_nas <- length(calls$earlier_nas)
decision <- if (stable) "STABLE" else "DEPENDS ON THE FRACTION SET"
verdict <- if (nc1_failed) {
  paste0("NC-1 FAILED: no composition-adjusted count is proposed (call-set decision would have been ", decision, ")")
} else decision
legend_old <- "Sixteen fibrosis and one NAS association meet BH q <0.05; after composition adjustment, three fibrosis and zero NAS associations remain."
legend_new <- if (nc1_failed) {
  "Sixteen fibrosis and one NAS association meet BH q <0.05. The composition-adjusted count is not reported because shuffled composition components also change it."
} else if (stable) {
  "Sixteen fibrosis and one NAS association meet BH q <0.05; after composition adjustment, three fibrosis and zero NAS associations remain with either deconvolution input."
} else {
  sprintf("Sixteen fibrosis and one NAS association meet BH q <0.05; after adjustment for three cell-type composition components, the count depends on the deconvolution input: three fibrosis and zero NAS associations with the rebuilt fractions, %s fibrosis and %s NAS with the earlier fractions; no adjusted program is named as a finding.",
          n_e_fib, n_e_nas)
}
p2_fib <- p2[axis == "fibrosis"]; p2_nas <- p2[axis == "nas"]
results_line <- sprintf(
  "  Figure 4E four-cohort arm (457 records, BH over 224 tests): composition-adjusted calls are 3 fibrosis and 0 NAS with the rebuilt fractions and %s fibrosis and %s NAS with the earlier fractions; Pearson r of adjusted betas between sets %.3f (fibrosis) and %.3f (NAS) (H85).",
  n_e_fib, n_e_nas, p2_fib$pearson_r, p2_nas$pearson_r)
duct_text <- if (pc2_pass) {
  sprintf("PC-2 passed: Ductular injury (BICC1) adjusted fibrosis beta under the earlier set %.3f [%.3f, %.3f] (rebuilt %.3f [%.3f, %.3f]); it may be described as holding after adjustment in this arm, subject to the R12-01 qualification.",
          duct_e$composition_pc_beta_meta_fibrosis, duct_e$composition_pc_ci_lower_fibrosis, duct_e$composition_pc_ci_upper_fibrosis,
          duct_r$composition_pc_beta_meta_fibrosis, duct_r$composition_pc_ci_lower_fibrosis, duct_r$composition_pc_ci_upper_fibrosis)
} else {
  sprintf("PC-2 failed: Ductular injury (BICC1) adjusted fibrosis beta under the earlier set %.3f [%.3f, %.3f] (rebuilt %.3f [%.3f, %.3f]). The earlier-set adjustment removes a signal that the other analyses keep; Ductular injury is not described as holding after adjustment in this arm.",
          duct_e$composition_pc_beta_meta_fibrosis, duct_e$composition_pc_ci_lower_fibrosis, duct_e$composition_pc_ci_upper_fibrosis,
          duct_r$composition_pc_beta_meta_fibrosis, duct_r$composition_pc_ci_lower_fibrosis, duct_r$composition_pc_ci_upper_fibrosis)
}
ecm_text <- sprintf("PC-3 %s: Stromal ECM (IGFBP7) adjusted fibrosis beta earlier %.3f (q %.3g) vs rebuilt %.3f (q %.3g); described only with this range across sets.",
                    if (pc3_pass) "passed" else "failed",
                    ecm_e$composition_pc_beta_meta_fibrosis, ecm_e$composition_pc_q_value_fibrosis,
                    ecm_r$composition_pc_beta_meta_fibrosis, ecm_r$composition_pc_q_value_fibrosis)
fmt_calls <- function(x) if (length(x)) paste(sprintf("%s (%s)", names(x), ifelse(x > 0, "+", "-")), collapse = ", ") else "none"
md_table <- function(dt) {
  dt <- copy(dt)
  for (col in names(dt)) if (is.numeric(dt[[col]])) set(dt, j = col, value = signif(dt[[col]], 4))
  c(paste0("| ", paste(names(dt), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
    apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
d4_short <- d4[model_kind == "composition_pc" & axis == "fibrosis",
               .(program, set, beta_meta, ci_lower, ci_upper, q_value, tau2, q_heterogeneity, n_same_direction,
                 beta_GSE130970, beta_GSE135251, beta_GSE162694, beta_GSE174478)]
result_lines <- c(
  result_header(verdict),
  "## Primary output P1 (four-cohort arm, 457 records, composition-adjusted, BH over 224 tests, q < 0.05)", "",
  paste0("- Rebuilt fractions (same-run refit, equals stored): fibrosis ", fmt_calls(calls$rebuilt_fibrosis),
         "; NAS ", fmt_calls(calls$rebuilt_nas), ". Matches stored set: ", rebuilt_matches_stored, "."),
  paste0("- Earlier fractions: fibrosis ", n_e_fib, " calls: ", fmt_calls(calls$earlier_fibrosis),
         "; NAS ", n_e_nas, " calls: ", fmt_calls(calls$earlier_nas), "."),
  paste0("- Decision rule (identical three fibrosis programs with the same signs and an empty NAS set): ", decision, "."),
  paste0("- NC-1 precedence: NC-1 ", if (nc1_failed) "FAILED" else "PASSED", "; reading: ", verdict, "."), "",
  "## P2, co-reported (not a decision)", "", md_table(p2), "",
  "## Controls", "",
  paste0("- PC-1a (stored 457 arm reproduced, 12 columns, 112 programs): ", pc1a_ok,
         "; max abs difference ", signif(max(pc1a$max_abs_difference), 3), "; counts ", paste(pc1a_counts, collapse = "/"), "."),
  paste0("- PC-1b (stored 469-record map reproduced; adjusted 0/0): ", pc1b_ok,
         "; max abs difference ", signif(max(reproduction$max_abs_difference), 3), "."),
  paste0("- ", duct_text),
  paste0("- ", ecm_text),
  sprintf("- NC-1 (%d within-cohort shuffles of the earlier-set components, seed %d): median adjusted fibrosis calls %s (range %s to %s; predicted 12 to 18; fail if < 10); median Pearson r with unadjusted fibrosis betas %.4f (range %.4f to %.4f; fail if < 0.95); draws with 3 or fewer calls: %d. Outcome: %s.",
          nrow(nc1), h85_seed, nc1_median_count, min(nc1$n_q05_fibrosis), max(nc1$n_q05_fibrosis),
          nc1_median_r, min(nc1$pearson_r_fibrosis_vs_unadjusted), max(nc1$pearson_r_fibrosis_vs_unadjusted),
          sum(nc1$n_q05_fibrosis <= 3L), if (nc1_failed) "FAILED" else "PASSED"), "",
  "## D3, test variants (BH over 224 per set, model and test)", "", md_table(d3[, .(set, model, test, n_estimable, n_q05_fibrosis, n_q05_nas, min_q_fibrosis)]), "",
  paste0("Rebuilt knha counts in D3 (check against the producer's meta_analyze, expected 16/1/3/0): ",
         paste(d3_knha_check, collapse = "/"), "."), "",
  "## D1 and D2 (descriptive)", "",
  md_table(summary_rows[arm %in% c("D1_new774_plus_GSE174478_earlier_pcs_as_computed", "new774_three_cohorts",
                                   "new774_three_cohorts_earlierfrac") ,
                        .(arm, model, axis, n_units, n_eligible, bh_family_size, n_q05, n_q05_positive, n_q05_negative, min_q)]), "",
  "## D4, focal programs, composition-adjusted fibrosis", "", md_table(d4_short), "",
  "## Proposed text (returned as data; not installed)", "",
  "- File: docs/manuscript/05_figure_legends.md, lines 254-255", "",
  paste0("  old_text: ", legend_old), "",
  paste0("  new_text: ", legend_new), "",
  "- File: docs/RESULTS.md, new line after line 289 (fraction-set paragraph)", "",
  paste0("  new_text: ", results_line), "",
  "The R12-01 qualification (0/0 under truncated Knapp-Hartung) is unchanged by H85; see D3.", "",
  "## Files", "",
  paste0("- ", c(diagnostics_path, file.path(out_root, c("P1_adjusted_call_sets.tsv", "P2_between_set_agreement.tsv",
                                                       "nc1_draws.tsv", "nc1_permutation_indices.rds",
                                                       "D3_test_variant_counts.tsv", "D4_focal_programs.tsv",
                                                       "D5_per_program_both_sets.tsv", "diagnostic6_pc_comparison.tsv",
                                                       "diagnostic7_pc_r2_on_axes.tsv", "pc1a_stored_457_reproduction.tsv",
                                                       "input_hash_check.tsv", "producer_copy.diff")),
                 analysis_root))
)
writeLines(result_lines, result_path)
if (!h85_debug) {
  top_result <- file.path(h85_dir, "RESULT_H85.md")
  if (!file.exists(top_result)) file.copy(result_path, top_result)
}
cat(result_lines, sep = "\n")
cat("H85_COMPLETE\n")
