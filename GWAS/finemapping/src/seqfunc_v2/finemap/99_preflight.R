#!/usr/bin/env Rscript

source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT", unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "GWAS", "finemapping", "src", "seqfunc_v2", "finemap", "common.R"
))

tier <- fread(file.path(FM_ROOT, "config", "gwas_trait_tier.tsv"))
registry <- fread(file.path(FM_ROOT, "config", "gwas_registry.tsv"))
main <- tier[placement == "main"]
x <- merge(main, registry, by = "study_name", all.x = TRUE)
assert_true(nrow(x) == 35L && uniqueN(x$study_name) == 35L, "Main study count is not exactly 35")
assert_true(sum(x$tier == 1L) == 15L && sum(x$tier == 2L) == 20L, "Tier counts changed")
expected_ancestry <- c(EUR = 17L, AFR = 6L, EAS = 6L, AMR = 3L, SAS = 3L)
observed <- table(x$ancestry)
assert_true(identical(as.integer(observed[names(expected_ancestry)]), as.integer(expected_ancestry)), "Ancestry counts changed")

required <- c("chromosome", "position", "allele1", "allele2", "beta", "se", "pval")
for (i in seq_len(nrow(x))) {
  p <- file.path(FM_ROOT, x$sumstats_path[i])
  assert_true(file.exists(p), "Missing sumstats: %s", p)
  header <- names(fread(p, nrows = 0L))
  assert_true(all(required %in% header), "Schema mismatch for %s", x$study_name[i])
}

panel_roots <- c(EUR = "polyfun_eur", AFR = "1kg_afr", AMR = "1kg_amr", EAS = "1kg_eas", SAS = "1kg_sas")
block_hashes <- character()
for (anc in names(panel_roots)) {
  p <- file.path(FM_ROOT, "data", "ld_ref", panel_roots[[anc]], "approx_LD_blocks.txt")
  blocks <- fread(p)
  assert_true(identical(names(blocks), c("chr", "start", "stop")), "Block schema mismatch for %s", anc)
  assert_true(nrow(blocks) == 1703L, "Expected 1,703 blocks for %s; found %d", anc, nrow(blocks))
  assert_true(all(blocks$stop > blocks$start), "Invalid block width for %s", anc)
  block_hashes[[anc]] <- sha256_file(p)
}
assert_true(uniqueN(block_hashes) == 1L, "Active panel block manifests are not identical")

stopifnot(
  identical(complement_allele(c("A", "C", "G", "T")), c("T", "G", "C", "A")),
  is_palindromic_pair("A", "T"), !is_palindromic_pair("A", "C"),
  memory_tier(5000L) == "small", memory_tier(5001L) == "medium",
  memory_tier(10001L) == "large", memory_tier(20001L) == "xlarge"
)

cat("Preflight PASS: exact 35-study contract, schemas, and five-panel block partition validated.\n")

