#!/usr/bin/env Rscript

source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT", unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "GWAS", "finemapping", "src", "seqfunc_v2", "finemap", "common.R"
))

tier_path <- file.path(FM_ROOT, "config", "gwas_trait_tier.tsv")
registry_path <- file.path(FM_ROOT, "config", "gwas_registry.tsv")
assert_true(file.exists(tier_path), "Missing tier table: %s", tier_path)
assert_true(file.exists(registry_path), "Missing registry: %s", registry_path)

tier <- fread(tier_path)
registry <- fread(registry_path)
required_tier <- c("study_name", "trait", "tier", "tier_label", "placement")
required_registry <- c(
  "study_name", "sumstats_path", "leadsnps_path", "ancestry", "trait_type",
  "N_tot", "N_cases", "ld_panel", "window_mb", "ld_panel_n"
)
assert_true(all(required_tier %in% names(tier)), "Tier table schema mismatch")
assert_true(all(required_registry %in% names(registry)), "Registry schema mismatch")

main <- tier[placement == "main"]
assert_true(nrow(main) == 35L, "Expected exactly 35 placement=main studies; found %d", nrow(main))
assert_true(all(main$tier %in% c(1L, 2L)), "Main manifest contains a non-Tier-1/2 study")
assert_true(uniqueN(main$study_name) == 35L, "Main tier table contains duplicate studies")
manifest <- merge(main, registry, by = "study_name", all.x = TRUE, sort = FALSE)
manifest <- manifest[match(main$study_name, study_name)]
assert_true(!anyNA(manifest$sumstats_path), "At least one main study is absent from the registry")

expected_ancestry <- c(EUR = 17L, AFR = 6L, EAS = 6L, AMR = 3L, SAS = 3L)
observed_ancestry <- table(manifest$ancestry)
assert_true(identical(as.integer(observed_ancestry[names(expected_ancestry)]), as.integer(expected_ancestry)),
            "Ancestry count mismatch; observed: %s",
            paste(names(observed_ancestry), observed_ancestry, collapse = ", "))
expected_tier <- c(`1` = 15L, `2` = 20L)
observed_tier <- table(manifest$tier)
assert_true(identical(as.integer(observed_tier[names(expected_tier)]), as.integer(expected_tier)),
            "Tier count mismatch; observed: %s", paste(names(observed_tier), observed_tier, collapse = ", "))

panel_dir <- c(
  EUR = "polyfun_eur", AFR = "1kg_afr", AMR = "1kg_amr",
  EAS = "1kg_eas", SAS = "1kg_sas"
)
panel_id <- c(
  EUR = "polyfun_ukbb_eur_337k", AFR = "1kg_phase3_afr",
  AMR = "1kg_phase3_amr", EAS = "1kg_phase3_eas", SAS = "1kg_phase3_sas"
)
panel_n <- c(EUR = 337000L, AFR = 661L, AMR = 347L, EAS = 504L, SAS = 489L)

manifest[, study_index := .I]
manifest[, sumstats_abs := normalizePath(file.path(FM_ROOT, sumstats_path), mustWork = TRUE)]
manifest[, source_leads_abs := normalizePath(file.path(FM_ROOT, leadsnps_path), mustWork = TRUE)]
manifest[, ld_panel_id_v2 := unname(panel_id[ancestry])]
manifest[, ld_panel_n_v2 := unname(panel_n[ancestry])]
manifest[, ld_panel_n_registry := as.integer(ld_panel_n)]
manifest[, ld_root := normalizePath(file.path(FM_ROOT, "data", "ld_ref", unname(panel_dir[ancestry])), mustWork = TRUE)]
manifest[, block_manifest := file.path(ld_root, "approx_LD_blocks.txt")]
assert_true(all(file.exists(manifest$block_manifest)), "At least one LD block manifest is missing")
manifest[, primary_locus_source := "sumstats_p_lt_5e-8"]
manifest[, genome_build := "GRCh37"]
manifest[, sumstats_sha256 := vapply(sumstats_abs, sha256_file, character(1))]
manifest[, source_leads_sha256 := vapply(source_leads_abs, sha256_file, character(1))]
manifest[, block_manifest_sha256 := vapply(block_manifest, sha256_file, character(1))]

manifest_cols <- c(
  "study_index", "study_name", "trait", "tier", "tier_label", "placement",
  "ancestry", "trait_type", "N_tot", "N_cases", "sumstats_abs", "sumstats_sha256",
  "source_leads_abs", "source_leads_sha256", "primary_locus_source", "genome_build",
  "ld_panel_id_v2", "ld_panel_n_v2", "ld_panel_n_registry", "ld_root", "block_manifest", "block_manifest_sha256"
)
manifest <- manifest[, ..manifest_cols]

ensure_dirs(
  file.path(RUN_ROOT, "config", "ld_blocks"), file.path(RUN_ROOT, "work", "prep_fragments"),
  file.path(RUN_ROOT, "work", "summary_stats"), file.path(RUN_ROOT, "results", "per_locus"),
  file.path(RUN_ROOT, "aggregate"), file.path(RUN_ROOT, "audit"), file.path(RUN_ROOT, "logs")
)
immutable_fwrite(manifest, file.path(RUN_ROOT, "config", "study_manifest.tsv"))
immutable_json(contract_list(), file.path(RUN_ROOT, "config", "run_contract.json"))

panel_rows <- unique(manifest[, .(ancestry, ld_panel_id_v2, ld_panel_n_v2, ld_root, block_manifest, block_manifest_sha256)])
setorder(panel_rows, ancestry)
for (i in seq_len(nrow(panel_rows))) {
  row <- panel_rows[i]
  blocks <- fread(row$block_manifest)
  assert_true(identical(names(blocks), c("chr", "start", "stop")),
              "Unexpected LD block schema: %s", row$block_manifest)
  blocks[, `:=`(
    ancestry = row$ancestry, ld_panel_id = row$ld_panel_id_v2,
    ld_panel_n = row$ld_panel_n_v2,
    block_prefix = file.path(row$ld_root, paste0("chr", chr), paste0(start, ".", stop), paste0(start, ".", stop))
  )]
  blocks[, `:=`(bim_exists = file.exists(paste0(block_prefix, ".bim")), ld_exists = file.exists(paste0(block_prefix, ".ld")))]
  immutable_fwrite(blocks, file.path(RUN_ROOT, "config", "ld_blocks", paste0(row$ancestry, ".tsv")))
}

provenance <- list(
  run_id = RUN_ID,
  frozen_release_date = "2026-07-13",
  project_root = normalizePath(PROJECT_ROOT),
  tier_table = list(path = normalizePath(tier_path), sha256 = sha256_file(tier_path)),
  registry = list(path = normalizePath(registry_path), sha256 = sha256_file(registry_path)),
  scripts = lapply(sort(list.files(SCRIPT_ROOT, full.names = TRUE)), function(p) {
    list(path = normalizePath(p), sha256 = sha256_file(p))
  }),
  counts = list(studies = nrow(manifest), tier1 = sum(manifest$tier == 1L), tier2 = sum(manifest$tier == 2L)),
  panel_size_note = "AFR registry metadata says 660; v2 freezes 661 observed .fam/sample rows. EUR 337000 is PolyFun/UKBB metadata.",
  note = "Standalone uniform-prior SuSiE sensitivity branch; no canonical promotion."
)
immutable_json(provenance, file.path(RUN_ROOT, "config", "provenance.json"))

cat(sprintf("Frozen %d studies at %s\n", nrow(manifest), file.path(RUN_ROOT, "config", "study_manifest.tsv")))
