#!/usr/bin/env Rscript
# 07: liver proteome, and plasma coverage.
#
# This cohort carries the complete NASH-CRN grading that the bulk RNA substrate
# does not. NAS is the activity-family endpoint; Kleiner fibrosis is the
# fibrosis-family endpoint and is badly unbalanced at n = 58 (F0 13, F1 31,
# F2 11, F3 3, F4 0), so its cells will mostly print an upper bound rather than
# a decision. Both are fitted; step 15 picks the one matching each module's
# discovery family.
#
# NO IMPUTATION. A protein missing in a participant is missing. Scores average
# the measured members per participant, so the per-participant member coverage
# and its relation to the endpoint are recorded (a score that averages fewer
# members in sicker livers is a different score in sicker livers).
#
# PROTEIN GROUPS. A group listing several genes is used only if exactly one of
# them is in the universe, and is assigned to that gene; v1 took the first
# listed gene. 59 of 7,096 groups are multi-gene.
#
# THE COMPLEXITY SENSITIVITY. Ten per-run intensity-distribution descriptors
# are added as a refit; step 15 caps a proteome cell at `indeterminate` unless
# the association survives them. They are measurement descriptors but are not
# established as purely technical, and missing_fraction is a deterministic
# function of n_detected for a fixed feature roster; both facts are stated in
# the caption rather than hidden.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))

contract <- cam_contract()
out <- cam_dir("assays")
cam_assert(file.exists(file.path(cam_out_root(), "discovery", "READY")), "Run 05 first")

membership <- fread(file.path(cam_out_root(), "modules", "module_membership.tsv"))
modules <- split(membership$gene_symbol, membership$module_id)
module_ids <- names(modules)
universe <- fread(file.path(cam_out_root(), "universe", "universe.tsv"))$gene_symbol

prot <- fread(cam_input("proteome_liver", contract))
meta <- fread(cam_input("proteome_meta", contract))
value_cols <- names(prot)[4:ncol(prot)]
mat <- as.matrix(prot[, ..value_cols])
gene_lists <- lapply(strsplit(prot$Genes, ";", fixed = TRUE), function(g) trimws(g[nzchar(trimws(g))]))
gene <- vapply(gene_lists, function(g) {
  if (length(g) == 1L) return(g)
  inu <- g[g %in% universe]
  if (length(inu) == 1L) inu else NA_character_
}, character(1))
group_rule <- data.table(
  n_groups = nrow(prot),
  n_single_gene = sum(lengths(gene_lists) == 1L),
  n_multi_gene = sum(lengths(gene_lists) > 1L),
  n_multi_gene_assigned = sum(lengths(gene_lists) > 1L & !is.na(gene)),
  n_multi_gene_dropped = sum(lengths(gene_lists) > 1L & is.na(gene)))
cam_write_tsv(group_rule, file.path(out, "proteome_protein_group_rule.tsv"))

from_multi <- lengths(gene_lists) > 1L
detect <- rowMeans(!is.na(mat))
keep <- detect >= contract$universe$proteome_min_detection_fraction & !is.na(gene)
# Genes whose measurement is a multi-gene protein group: the intensity cannot
# be attributed to the named gene alone, and every module carrying one says so.
ambiguous_genes <- unique(gene[keep & from_multi])
mat <- log2(mat[keep, , drop = FALSE]); gene <- gene[keep]
present <- matrix(as.numeric(!is.na(mat)), nrow = nrow(mat), ncol = ncol(mat))
expr <- rowsum(replace(mat, is.na(mat), 0), group = gene) / rowsum(present, group = gene)
expr[!is.finite(expr)] <- NA_real_
cam_say("liver proteome: ", nrow(expr), " genes x ", ncol(expr), " participants")

pm <- meta[match(colnames(expr), liver_proteomics_filename)]
cam_assert(sum(is.na(pm$unique_identifier)) == 0L,
           "Some liver proteome columns do not join to the metadata")

complexity <- data.table(
  n_detected = colSums(!is.na(mat)),
  total_intensity = colSums(mat, na.rm = TRUE),
  median_intensity = apply(mat, 2L, stats::median, na.rm = TRUE),
  iqr_intensity = apply(mat, 2L, stats::IQR, na.rm = TRUE),
  sd_intensity = apply(mat, 2L, stats::sd, na.rm = TRUE),
  dynamic_range = apply(mat, 2L, function(v) diff(range(v, na.rm = TRUE))),
  top10_share = apply(mat, 2L, function(v) {
    v <- sort(2^v[!is.na(v)], decreasing = TRUE); sum(head(v, 10)) / sum(v) }),
  gini = apply(mat, 2L, function(v) {
    v <- sort(2^v[!is.na(v)]); n <- length(v)
    sum((2 * seq_len(n) - n - 1) * v) / (n * sum(v)) }),
  entropy = apply(mat, 2L, function(v) {
    p <- 2^v[!is.na(v)]; p <- p / sum(p); -sum(p * log(p)) }),
  missing_fraction = colMeans(is.na(mat))
)
complexity_mat <- scale(as.matrix(complexity))
complexity_mat[!is.finite(complexity_mat)] <- 0
cam_write_tsv(cbind(data.table(sample = colnames(expr)), complexity),
              file.path(out, "proteome_complexity_descriptors.tsv"))

z <- cam_zscore_rows(expr)
measured <- rownames(z)
pm[, kleiner_numeric := suppressWarnings(as.numeric(sub("^F", "", kleiner_fibrosis_grade)))]
pm[, nas_numeric := suppressWarnings(as.numeric(nafld_activity_score))]
pm[, age_numeric := suppressWarnings(as.numeric(alder))]
pm[, bmi_numeric := suppressWarnings(as.numeric(bmi))]
covar_df <- data.frame(age = pm$age_numeric, bmi = pm$bmi_numeric, sex = pm$gender)
endpoints <- list(
  nafld_activity_score = list(y = cam_endpoint_z(pm$nas_numeric, covar_df), family = "activity",
                              note = "NAS 0-8, z-scored"),
  kleiner_fibrosis_grade = list(y = cam_endpoint_z(pm$kleiner_numeric, covar_df), family = "fibrosis",
                                note = "Kleiner F0-F3 (no F4; F3 n=3), z-scored"))

rows <- list(); cov_rows <- list()
for (i in seq_along(module_ids)) {
  tb <- cam_testability(modules[[i]], measured, contract)
  idx <- intersect(modules[[i]], measured)
  score <- if (length(idx)) standardize_vector(colMeans(z[idx, , drop = FALSE], na.rm = TRUE))
           else rep(NA_real_, ncol(z))
  # Per-participant member coverage of the score actually averaged.
  member_frac <- if (length(idx)) colMeans(!is.na(z[idx, , drop = FALSE])) else rep(NA_real_, ncol(z))
  amb <- intersect(modules[[i]], ambiguous_genes)
  for (ep in names(endpoints)) {
    e <- endpoints[[ep]]
    d <- data.table(score = score, y = e$y, age = pm$age_numeric, bmi = pm$bmi_numeric, sex = pm$gender)
    fit <- if (tb$testable) cam_fit_score(d, "score", "y", c("age", "bmi", "sex")) else cam_fit_empty()
    dd <- cbind(d, as.data.table(complexity_mat))
    fit_cx <- if (tb$testable)
      cam_fit_score(dd, "score", "y", c("age", "bmi", "sex", colnames(complexity_mat))) else cam_fit_empty()
    miss_p <- if (tb$testable && stats::sd(member_frac, na.rm = TRUE) > 0)
      cam_fit_score(data.table(score = member_frac, y = e$y), "score", "y")$p_two_sided else NA_real_
    rows[[length(rows) + 1L]] <- data.table(
      assay = "proteome_liver_pxd051911", module_id = module_ids[i],
      endpoint = ep, endpoint_family = e$family, endpoint_note = e$note,
      unit = "participant", n_units = fit$n, df = fit$df,
      testable = tb$testable, n_members = tb$n_members,
      n_measured = tb$n_measured, fraction_measured = tb$fraction_measured,
      beta = fit$beta, se = fit$se, p_two_sided = fit$p_two_sided,
      complexity_adjusted_beta = fit_cx$beta, complexity_adjusted_p = fit_cx$p_two_sided,
      n_members_from_multigene_groups = length(amb),
      multigene_members = paste(amb, collapse = ";"),
      member_coverage_min = min(member_frac, na.rm = TRUE),
      member_coverage_median = stats::median(member_frac, na.rm = TRUE),
      member_missingness_vs_endpoint_p = miss_p)
  }
  cov_rows[[length(cov_rows) + 1L]] <- data.table(module_id = module_ids[i],
    sample = colnames(z), member_fraction_measured = member_frac)
}
res <- rbindlist(rows)
cam_assert_no_prohibited_columns(res, contract)
cam_write_tsv(res, file.path(out, "proteome_results.tsv"))
cam_write_tsv(rbindlist(cov_rows), file.path(out, "proteome_member_coverage_by_participant.tsv"))

# Plasma is coverage only: a different compartment and a different claim.
plasma <- fread(cam_input("proteome_plasma", contract), select = "ProteinAccessions")
uni <- fread(cam_input("uniprot_symbol_axis", contract),
             select = c("protein_accession", "gene_symbol"))
acc2sym <- setNames(uni$gene_symbol, uni$protein_accession)
plasma_sym <- unique(stats::na.omit(unname(acc2sym[
  vapply(strsplit(plasma$ProteinAccessions, ";", fixed = TRUE),
         function(a) { a <- trimws(a[nzchar(trimws(a))]); if (length(a)) a[1] else NA_character_ },
         character(1))])))
plasma_cov <- rbindlist(lapply(seq_along(module_ids), function(i) {
  tb <- cam_testability(modules[[i]], plasma_sym, contract)
  data.table(assay = "proteome_plasma_pxd051911", module_id = module_ids[i],
             role = "coverage_only", n_members = tb$n_members,
             n_measured = tb$n_measured, fraction_measured = tb$fraction_measured,
             testable = tb$testable)
}))
cam_write_tsv(plasma_cov, file.path(out, "plasma_coverage.tsv"))

cam_write_json(list(
  liver_participants = ncol(expr), liver_genes = nrow(expr),
  modules_testable = sum(res[endpoint == "nafld_activity_score", testable]),
  kleiner_distribution = as.list(table(pm$kleiner_fibrosis_grade)),
  nas_range = range(pm$nas_numeric, na.rm = TRUE),
  protein_group_rule = as.list(group_rule),
  ambiguous_genes_in_universe = length(ambiguous_genes),
  modules_with_ambiguous_member = length(unique(res[n_members_from_multigene_groups > 0, module_id])),
  plasma_symbols = length(plasma_sym),
  plasma_modules_testable = sum(plasma_cov$testable),
  n_complexity_descriptors = ncol(complexity_mat),
  direction_applied_at_scoring = FALSE
), file.path(out, "proteome_summary.json"))
writeLines("proteome scored", file.path(out, "READY_proteome"))
cam_say("07 complete")
