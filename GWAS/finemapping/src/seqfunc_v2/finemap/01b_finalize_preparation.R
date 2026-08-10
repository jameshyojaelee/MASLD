#!/usr/bin/env Rscript

source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT", unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "GWAS", "finemapping", "src", "seqfunc_v2", "finemap", "common.R"
))

manifest <- fread(file.path(RUN_ROOT, "config", "study_manifest.tsv"))
fragment_dir <- file.path(RUN_ROOT, "work", "prep_fragments")
status_paths <- file.path(fragment_dir, sprintf("study_%02d_status.tsv", manifest$study_index))
locus_paths <- file.path(fragment_dir, sprintf("study_%02d_loci.tsv", manifest$study_index))
assert_true(all(file.exists(status_paths)), "Preparation status missing for study indices: %s",
            paste(manifest$study_index[!file.exists(status_paths)], collapse = ","))
assert_true(all(file.exists(locus_paths)), "Locus fragments missing for study indices: %s",
            paste(manifest$study_index[!file.exists(locus_paths)], collapse = ","))

study_status <- rbindlist(lapply(status_paths, fread), fill = TRUE)
assert_true(nrow(study_status) == 35L && uniqueN(study_status$study_name) == 35L,
            "Preparation did not produce one status row per frozen study")
immutable_fwrite(study_status, file.path(RUN_ROOT, "config", "preparation_study_status.tsv"))

# Strand certificates.  A study with no genome-wide-significant signal exits
# before this fragment is written and contributes no loci, so a missing fragment
# is expected and is simply absent from the table rather than an error.
strand_paths <- file.path(fragment_dir, sprintf("study_%02d_strand.tsv", manifest$study_index))
strand_paths <- strand_paths[file.exists(strand_paths)]
strand <- if (length(strand_paths)) rbindlist(lapply(strand_paths, fread), fill = TRUE) else data.table()
if (nrow(strand)) {
  # Recompute the label from the stored counts rather than trusting the one the
  # fragment was written with.  The counts are the measurement; the label is a
  # thresholded interpretation of it, so re-calibrating STRAND_CERT_MAX_OPP must
  # not require re-running the expensive per-study preparation.
  strand[, strand_certificate := mapply(strand_certificate_label, n_same, n_opp)]
  setorder(strand, study_index)
}
immutable_fwrite(strand, file.path(RUN_ROOT, "config", "strand_certificates.tsv"))
cat(sprintf("Strand certificates: %d studies (%s)\n", nrow(strand),
            if (nrow(strand)) paste(names(table(strand$strand_certificate)),
                                    table(strand$strand_certificate), sep = "=", collapse = " ") else "none"))

locus_parts <- lapply(locus_paths, function(p) {
  x <- fread(p)
  if (nrow(x) == 0L) return(NULL)
  x
})
loci <- rbindlist(locus_parts, fill = TRUE)
assert_true(nrow(loci) > 0L, "No GWS loci were prepared")
assert_true(!anyDuplicated(loci$locus_id), "Duplicate locus_id in prepared manifest")
assert_true(!anyDuplicated(loci[, .(study_name, chromosome, block_start, block_stop)]),
            "Duplicate study-block in prepared manifest")
setorder(loci, study_index, chromosome, block_start, block_stop)
loci[, locus_index := .I]
setcolorder(loci, c("locus_index", setdiff(names(loci), "locus_index")))
immutable_fwrite(loci, file.path(RUN_ROOT, "config", "locus_manifest.tsv"))

for (tier_name in c("small", "medium", "large", "xlarge", "unknown")) {
  tasks <- loci[memory_tier == tier_name, .(task_index = seq_len(.N), locus_index, locus_id)]
  immutable_fwrite(tasks, file.path(RUN_ROOT, "config", paste0("tasks_", tier_name, ".tsv")))
}

panel_provenance <- unique(loci[, .(
  ancestry, ld_panel_id, ld_panel_n, chromosome, block_start, block_stop,
  block_prefix, bim_path, ld_path, bim_exists, ld_exists, bim_bytes, ld_bytes,
  bim_sha256, ld_sha256, n_ld_variants, ld_min_position, ld_max_position,
  nominal_block_width, panel_covered_width, panel_coverage_fraction, panel_truncated
)])
setorder(panel_provenance, ancestry, chromosome, block_start)
immutable_fwrite(panel_provenance, file.path(RUN_ROOT, "config", "used_ld_blocks.tsv"))

summary <- list(
  run_id = RUN_ID, studies = nrow(study_status), loci = nrow(loci),
  gws_variants = sum(loci$n_gws_variants),
  loci_by_ancestry = as.list(table(loci$ancestry)),
  loci_by_memory_tier = as.list(table(loci$memory_tier)),
  missing_bim = sum(!loci$bim_exists), missing_ld = sum(!loci$ld_exists),
  truncated_panel_blocks = sum(loci$panel_truncated, na.rm = TRUE),
  gws_positionally_unrepresented = sum(loci$n_gws_positionally_unrepresented),
  strand_certificates = if (nrow(strand)) as.list(table(strand$strand_certificate)) else list(),
  studies_with_af = if (nrow(strand)) sum(strand$af_present) else 0L,
  studies_indeterminate_strand = if (nrow(strand)) sum(strand$strand_certificate == "indeterminate") else NA_integer_,
  canonical_outputs_mutated = FALSE
)
immutable_json(summary, file.path(RUN_ROOT, "config", "preparation_summary.json"))
cat(sprintf("Finalized preparation: %d studies, %d unique study-block loci\n", nrow(study_status), nrow(loci)))
