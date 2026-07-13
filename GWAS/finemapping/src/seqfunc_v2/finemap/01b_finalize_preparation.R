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
  canonical_outputs_mutated = FALSE
)
immutable_json(summary, file.path(RUN_ROOT, "config", "preparation_summary.json"))
cat(sprintf("Finalized preparation: %d studies, %d unique study-block loci\n", nrow(study_status), nrow(loci)))
