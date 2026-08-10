#!/usr/bin/env Rscript
# Export the selected SuSiE signal pair and variant-level shared posterior for
# one outcome-blind Plan 45 orientation-worklist row. This is a targeted
# reproduction of the corrected coloc model, not a new association scan.

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(susieR)
})

args <- commandArgs(trailingOnly = TRUE)
task_index <- as.integer(Sys.getenv(
  "PLAN45_ORIENTATION_INDEX",
  unset = if (length(args) >= 1L) args[[1L]] else "-1"
))
if (!is.finite(task_index) || task_index < 0L) {
  stop("PLAN45_ORIENTATION_INDEX must be a zero-based non-negative integer")
}

project_root <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
candidate_id <- Sys.getenv("PLAN45_CANDIDATE_ID", unset = "")
if (!grepl("^source-independent-risk-state-relay-[A-Za-z0-9._-]+$", candidate_id)) {
  stop("PLAN45_CANDIDATE_ID is missing or invalid")
}
candidate_root <- file.path(
  project_root,
  "Analysis/Multimodal_Program_Projection/candidates",
  candidate_id
)
worklist_path <- file.path(candidate_root, "orientation_worklist.tsv")
seal_path <- file.path(candidate_root, "ORIENTATION_WORKLIST_SEALED.json")
if (!file.exists(worklist_path) || !file.exists(seal_path)) {
  stop("Sealed orientation worklist is absent from candidate root")
}

worklist <- fread(worklist_path, sep = "\t", na.strings = character())
row_number <- task_index + 1L
if (row_number > nrow(worklist)) {
  stop(sprintf("Orientation index %d exceeds %d worklist rows", task_index, nrow(worklist)))
}
task <- worklist[row_number]
if (task$pair_export_status != "required_not_run" ||
    task$target_freeze_status != "prohibited_until_orientation_and_lineage_gates") {
  stop("Worklist row is not in the required fail-closed state")
}

relative_outputs <- c(
  task$required_pair_summary,
  task$required_variant_posterior,
  task$required_allele_audit
)
if (any(grepl("(^/|(^|/)\\.\\.(/|$))", relative_outputs))) {
  stop("Orientation output path is not a safe candidate-relative path")
}
output_paths <- file.path(candidate_root, relative_outputs)
if (any(file.exists(output_paths))) {
  stop("Refusing to overwrite an existing targeted pair export")
}
output_dir <- dirname(output_paths[[1L]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

atomic_fwrite <- function(x, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "")
  if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path)
}

fm_dir <- file.path(project_root, "GWAS/finemapping")
setwd(fm_dir)
source(file.path(fm_dir, "src/finemapping_functions.R"))

co_loc_p1 <- 1e-4
co_loc_p2 <- 1e-4
co_loc_p12 <- 5e-6
eqtl_n <- 1183
minimum_snps <- 100L
minimum_triple_snps <- 50L
ld_regularize <- 1e-3
susie_l <- 10L
eqtl_susie_dir <- file.path(fm_dir, "results/eqtl_susie_polyfun")
eqtl_dir <- file.path(project_root, "data/broadaway_eqtl")
panel_af_root <- file.path(fm_dir, "data/ld_ref/panel_af")

load_panel_af <- function(population, chromosome) {
  path <- file.path(
    panel_af_root,
    tolower(population),
    sprintf("chr%s.af.tsv.gz", chromosome)
  )
  if (!file.exists(path)) return(NULL)
  result <- fread(
    path,
    select = c("position", "bim_a1", "bim_a2", "af_a1"),
    showProgress = FALSE
  )
  result[, `:=`(bim_a1 = toupper(bim_a1), bim_a2 = toupper(bim_a2))]
  unique(result, by = "position")
}

panel_freq_of <- function(allele, other, bim_a1, bim_a2, af_a1) {
  pair_ok <- (allele == bim_a1 & other == bim_a2) |
    (allele == bim_a2 & other == bim_a1)
  fifelse(
    !pair_ok | is.na(pair_ok),
    NA_real_,
    fifelse(allele == bim_a1, af_a1, 1 - af_a1)
  )
}

strand_verdict_vs_panel <- function(observed_af, panel_af) {
  d_same <- abs(observed_af - panel_af)
  d_opposite <- abs(observed_af - (1 - panel_af))
  testable <- is.finite(observed_af) & is.finite(panel_af) &
    pmin(observed_af, 1 - observed_af) <= 0.40 &
    pmin(panel_af, 1 - panel_af) <= 0.40 &
    pmin(d_same, d_opposite) <= 0.15
  fifelse(
    !testable,
    NA_character_,
    fifelse(
      (d_opposite - d_same) >= 0.10,
      "same",
      fifelse((d_same - d_opposite) >= 0.10, "opposite", NA_character_)
    )
  )
}

registry <- fread(file.path(fm_dir, "config/gwas_registry.tsv"))
study <- registry[study_name == task$gwas_name]
if (nrow(study) != 1L) stop("Worklist GWAS does not map to exactly one registry row")
chromosome <- as.integer(task$chromosome)
chromosome_value <- chromosome
gene_id <- task$ensembl_id
if (!grepl("^ENSG[0-9]+$", gene_id)) stop("Worklist lacks one exact driving Ensembl ID")

gwas_path <- study$sumstats_path[[1L]]
af_sources_path <- file.path(fm_dir, "config/gwas_af_sources.tsv")
if (file.exists(af_sources_path)) {
  af_sources <- fread(af_sources_path)
  af_hit <- af_sources[study_name == task$gwas_name & file.exists(af_sumstats_path)]
  if (nrow(af_hit) == 1L) gwas_path <- af_hit$af_sumstats_path[[1L]]
}
gwas <- fread(gwas_path)[chromosome == chromosome_value]
if (nrow(gwas) == 0L) stop("No GWAS variants on target chromosome")
gwas[, merge_key := paste(chromosome, position, sep = ":")]

eqtl_path <- file.path(eqtl_dir, sprintf("chr%d_marginal_summary_results.tsv", chromosome))
eqtl_gene <- fread(eqtl_path)[ENSG == gene_id]
if (nrow(eqtl_gene) < minimum_snps) stop("Insufficient source eQTL variants")
if (uniqueN(eqtl_gene$GeneSymbol) != 1L || eqtl_gene$GeneSymbol[[1L]] != task$gene_symbol) {
  stop("Worklist symbol and source eQTL Ensembl mapping disagree")
}
eqtl_gene[, merge_key := paste(CHR, POS, sep = ":")]

merged <- merge(
  eqtl_gene[, .(
    merge_key,
    eqtl_pos = POS,
    eqtl_ea = toupper(EA),
    eqtl_nea = toupper(NEA),
    eqtl_beta_raw = Beta,
    eqtl_se = SE,
    eqtl_pval = PVAL,
    eqtl_eaf = EAF
  )],
  gwas[, .(
    merge_key,
    gwas_a1 = toupper(allele1),
    gwas_a2 = toupper(allele2),
    gwas_beta = beta,
    gwas_se = se,
    gwas_pval = pval,
    gwas_af = if ("af" %in% names(gwas)) af else NA_real_
  )],
  by = "merge_key"
)
merged <- merged[order(gwas_pval)][!duplicated(merge_key)]
if (nrow(merged) < minimum_snps) stop("Insufficient overlapping eQTL/GWAS variants")

n_merged_pre <- nrow(merged)
merged[, allele_match := eqtl_ea == gwas_a1 & eqtl_nea == gwas_a2]
merged[, allele_flip := eqtl_ea == gwas_a2 & eqtl_nea == gwas_a1]
n_match <- sum(merged$allele_match, na.rm = TRUE)
n_flip <- sum(merged$allele_flip, na.rm = TRUE)
n_unresolved <- n_merged_pre - n_match - n_flip
merged <- merged[allele_match | allele_flip]
merged[, eqtl_beta := fifelse(allele_flip, -eqtl_beta_raw, eqtl_beta_raw)]
merged[, is_palindromic :=
  (gwas_a1 %in% c("A", "T") & gwas_a2 %in% c("A", "T")) |
  (gwas_a1 %in% c("C", "G") & gwas_a2 %in% c("C", "G"))]
n_palindromic_present <- sum(merged$is_palindromic, na.rm = TRUE)

gwas_ancestry <- study$ancestry[[1L]]
pa_gwas <- load_panel_af(gwas_ancestry, chromosome)
pa_eqtl <- load_panel_af("EUR", chromosome)
merged[, `:=`(gwas_strand = NA_character_, eqtl_strand = NA_character_)]
if (!is.null(pa_gwas)) {
  merged[pa_gwas, on = .(eqtl_pos = position),
    `:=`(pg_a1 = i.bim_a1, pg_a2 = i.bim_a2, pg_af = i.af_a1)]
  merged[, gwas_strand := strand_verdict_vs_panel(
    gwas_af,
    panel_freq_of(gwas_a1, gwas_a2, pg_a1, pg_a2, pg_af)
  )]
}
if (!is.null(pa_eqtl)) {
  merged[pa_eqtl, on = .(eqtl_pos = position),
    `:=`(pe_a1 = i.bim_a1, pe_a2 = i.bim_a2, pe_af = i.af_a1)]
  merged[, eqtl_strand := strand_verdict_vs_panel(
    eqtl_eaf,
    panel_freq_of(eqtl_ea, eqtl_nea, pe_a1, pe_a2, pe_af)
  )]
}
merged[, af_testable := is_palindromic & !is.na(gwas_strand)]
merged[, af_inconsistent := af_testable & gwas_strand == "opposite"]
n_palindromic_af_tested <- sum(merged$af_testable, na.rm = TRUE)
n_palindromic_af_dropped <- sum(merged$af_inconsistent, na.rm = TRUE)
merged <- merged[!(af_inconsistent %in% TRUE)]
n_palindromic_kept <- sum(merged$is_palindromic, na.rm = TRUE)
if (nrow(merged) < minimum_snps) stop("Insufficient variants after harmonization")

eqtl_rds <- file.path(
  eqtl_susie_dir,
  sprintf("chr%d", chromosome),
  paste0(gene_id, "_susie.rds")
)
if (!file.exists(eqtl_rds)) stop("Pinned eQTL SuSiE fit is absent")
s_eqtl <- readRDS(eqtl_rds)
if (!isTRUE(s_eqtl$converged)) stop("Pinned eQTL SuSiE fit did not converge")
if (is.null(s_eqtl$ld_provenance) || s_eqtl$ld_provenance$ld_panel != "polyfun") {
  stop("Pinned eQTL SuSiE provenance is absent or not PolyFun EUR")
}

eqtl_snp_ids <- colnames(s_eqtl$lbf_variable)
eqtl_positions <- as.integer(sub("^[0-9]+:", "", eqtl_snp_ids))
gwas_sub <- merged[eqtl_pos %in% eqtl_positions]
if (nrow(gwas_sub) < minimum_triple_snps) stop("Insufficient GWAS/eQTL/LD overlap")
ss_for_ld <- data.frame(
  chromosome = chromosome,
  position = gwas_sub$eqtl_pos,
  allele1 = gwas_sub$gwas_a1,
  allele2 = gwas_sub$gwas_a2,
  beta = gwas_sub$gwas_beta,
  standard_error = gwas_sub$gwas_se
)
ld_result <- get_ld_per_locus(
  ss_per_locus = ss_for_ld,
  LOCUS = gene_id,
  CHR = chromosome,
  START = min(gwas_sub$eqtl_pos),
  END = max(gwas_sub$eqtl_pos),
  ancestry = gwas_ancestry
)
if (is.null(ld_result) || nrow(ld_result[[1L]]) < minimum_triple_snps) {
  stop("Targeted GWAS LD extraction failed")
}
ss_ld <- ld_result[[1L]]
r_matrix <- as.matrix(ld_result[[2L]])
r_regularized <- r_matrix + ld_regularize * diag(nrow(r_matrix))
gwas_z <- ss_ld$beta / ss_ld$standard_error
snp_ids <- paste0(chromosome, ":", ss_ld$position)
names(gwas_z) <- snp_ids
colnames(r_regularized) <- rownames(r_regularized) <- snp_ids

gwas_type <- study$trait_type[[1L]]
gwas_n <- study$N_tot[[1L]]
gwas_cases <- study$N_cases[[1L]]
gwas_controls <- if (!is.na(gwas_cases)) gwas_n - gwas_cases else NA_real_
gwas_n_eff <- if (
  gwas_type == "binary" && !is.na(gwas_cases) && gwas_cases > 0 &&
    !is.na(gwas_controls) && gwas_controls > 0
) {
  4 / (1 / gwas_cases + 1 / gwas_controls)
} else {
  gwas_n
}
set.seed(42)
s_gwas <- susie_rss(
  z = gwas_z,
  R = r_regularized,
  n = gwas_n_eff,
  L = susie_l,
  estimate_residual_variance = FALSE,
  residual_variance = 1,
  check_R = FALSE,
  min_abs_corr = 0.5,
  max_iter = 500
)
if (!isTRUE(s_gwas$converged)) stop("Targeted GWAS SuSiE fit did not converge")
common_snps <- intersect(snp_ids, colnames(s_eqtl$lbf_variable))
if (length(common_snps) < minimum_triple_snps) stop("Too few common SuSiE SNPs")
s_eqtl_sub <- s_eqtl
s_eqtl_sub$lbf_variable <- s_eqtl$lbf_variable[, common_snps, drop = FALSE]
s_gwas_sub <- s_gwas
s_gwas_sub$lbf_variable <- s_gwas$lbf_variable[, common_snps, drop = FALSE]
set.seed(42)
susie_result <- coloc.susie(
  s_gwas_sub,
  s_eqtl_sub,
  p1 = co_loc_p1,
  p2 = co_loc_p2,
  p12 = co_loc_p12
)
if (is.null(susie_result$summary) || nrow(susie_result$summary) == 0L ||
    is.null(susie_result$results)) {
  stop("Targeted coloc.susie returned no signal-pair result")
}
pair_index <- which.max(susie_result$summary$PP.H4.abf)
selected_pair <- as.data.table(susie_result$summary[pair_index, , drop = FALSE])
observed_pp4 <- selected_pair$PP.H4.abf[[1L]]
expected_pp4 <- as.numeric(task$susie_pp4)
pp4_delta <- observed_pp4 - expected_pp4
if (!is.finite(pp4_delta) || abs(pp4_delta) > 1e-8) {
  stop(sprintf(
    "Targeted pair PP.H4 failed portfolio reproduction: observed=%.12g expected=%.12g delta=%.3g",
    observed_pp4, expected_pp4, pp4_delta
  ))
}

posterior_column <- if (nrow(susie_result$summary) == 1L) {
  "SNP.PP.H4.abf"
} else {
  sprintf("SNP.PP.H4.row%d", pair_index)
}
if (!posterior_column %in% names(susie_result$results)) {
  stop("Selected pair lacks a variant-level SNP.PP.H4 column")
}
variant <- as.data.table(susie_result$results)[, .(
  snp,
  shared_posterior = get(posterior_column)
)]
variant <- merge(
  variant,
  merged[, .(
    snp = merge_key,
    position = eqtl_pos,
    gwas_effect_allele = gwas_a1,
    gwas_other_allele = gwas_a2,
    gwas_beta,
    gwas_se,
    eqtl_effect_allele_source = eqtl_ea,
    eqtl_other_allele_source = eqtl_nea,
    eqtl_beta_oriented_to_gwas_effect = eqtl_beta,
    eqtl_se,
    allele_match,
    allele_flip,
    is_palindromic,
    gwas_strand,
    eqtl_strand
  )],
  by = "snp",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(variant$gwas_beta) || anyNA(variant$eqtl_beta_oriented_to_gwas_effect)) {
  stop("Variant shared posterior could not be joined authoritatively to effects")
}
variant[, `:=`(
  risk_allele = fifelse(gwas_beta > 0, gwas_effect_allele, gwas_other_allele),
  risk_allele_effect_on_expression = fifelse(
    gwas_beta > 0,
    eqtl_beta_oriented_to_gwas_effect,
    -eqtl_beta_oriented_to_gwas_effect
  ),
  risk_to_expression_sign = sign(gwas_beta * eqtl_beta_oriented_to_gwas_effect)
)]
setorder(variant, -shared_posterior, snp)
variant[, posterior_rank := seq_len(.N)]
variant[, cumulative_shared_posterior := cumsum(shared_posterior)]
variant[, in_shared_95 := cumulative_shared_posterior - shared_posterior < 0.95]

positive_weight <- sum(variant$shared_posterior[variant$risk_to_expression_sign > 0])
negative_weight <- sum(variant$shared_posterior[variant$risk_to_expression_sign < 0])
dominant_direction <- if (positive_weight > negative_weight) "risk_increases_expression" else
  if (negative_weight > positive_weight) "risk_decreases_expression" else "unresolved"
consensus <- max(positive_weight, negative_weight)
leave_top <- variant[posterior_rank > 1L]
leave_top_total <- sum(leave_top$shared_posterior)
leave_top_positive <- if (leave_top_total > 0) {
  sum(leave_top$shared_posterior[leave_top$risk_to_expression_sign > 0]) / leave_top_total
} else {
  NA_real_
}
leave_top_negative <- if (leave_top_total > 0) {
  sum(leave_top$shared_posterior[leave_top$risk_to_expression_sign < 0]) / leave_top_total
} else {
  NA_real_
}
leave_top_direction <- if (!is.finite(leave_top_positive)) "unresolved" else
  if (leave_top_positive > leave_top_negative) "risk_increases_expression" else
    if (leave_top_negative > leave_top_positive) "risk_decreases_expression" else "unresolved"
leave_top_consensus <- suppressWarnings(max(leave_top_positive, leave_top_negative, na.rm = TRUE))
if (!is.finite(leave_top_consensus)) leave_top_consensus <- NA_real_

pair_summary <- data.table(
  orientation_uid = task$orientation_uid,
  gene_symbol = task$gene_symbol,
  ensembl_id = gene_id,
  gwas_name = task$gwas_name,
  chromosome = chromosome,
  selected_pair_row = pair_index,
  gwas_signal_index = selected_pair$idx1,
  eqtl_signal_index = selected_pair$idx2,
  gwas_signal_hit = selected_pair$hit1,
  eqtl_signal_hit = selected_pair$hit2,
  n_pair_snps = selected_pair$nsnps,
  portfolio_susie_pp4 = expected_pp4,
  reproduced_susie_pp4 = observed_pp4,
  pp4_delta = pp4_delta,
  positive_direction_posterior = positive_weight,
  negative_direction_posterior = negative_weight,
  orientation_consensus = consensus,
  oriented_risk_effect = dominant_direction,
  leave_top_orientation_consensus = leave_top_consensus,
  leave_top_oriented_risk_effect = leave_top_direction,
  consensus_gate_80 = consensus >= 0.80,
  leave_top_direction_agrees = leave_top_direction == dominant_direction,
  target_freeze_status = "prohibited_until_lineage_and_editability_gates"
)

allele_audit <- data.table(
  orientation_uid = task$orientation_uid,
  source_gwas_path = gwas_path,
  source_eqtl_path = eqtl_path,
  source_eqtl_susie_rds = eqtl_rds,
  eqtl_ld_panel = s_eqtl$ld_provenance$ld_panel,
  gwas_ld_panel = get_ld_base_dir(gwas_ancestry),
  n_merged_pre = n_merged_pre,
  n_allele_match = n_match,
  n_allele_flip = n_flip,
  n_allele_unresolved = n_unresolved,
  n_palindromic_present = n_palindromic_present,
  n_palindromic_af_tested = n_palindromic_af_tested,
  n_palindromic_af_dropped = n_palindromic_af_dropped,
  n_palindromic_kept = n_palindromic_kept,
  n_variant_posterior = nrow(variant),
  n_shared_95 = sum(variant$in_shared_95),
  coloc_version = as.character(packageVersion("coloc")),
  susieR_version = as.character(packageVersion("susieR")),
  data_table_version = as.character(packageVersion("data.table")),
  pair_pp4_reproduction_status = "pass",
  target_list_frozen = FALSE
)

atomic_fwrite(pair_summary, output_paths[[1L]])
atomic_fwrite(variant, output_paths[[2L]])
atomic_fwrite(allele_audit, output_paths[[3L]])
cat(sprintf(
  "Plan 45 pair export complete: %s; PP.H4=%.6f; orientation=%s; consensus=%.3f; target not frozen\n",
  task$orientation_uid,
  observed_pp4,
  dominant_direction,
  consensus
))
