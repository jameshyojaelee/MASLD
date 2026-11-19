#!/usr/bin/env Rscript
# SuSiE-COLOC for the affected primary + 5 unaffected control genes on the PSD refit eQTL fits.
# Sources the hash-checked production 06_susie_coloc.R with patches:
#   - coloc_identity_guarded_runner (duplicate refusal + Cholesky PD check on PSD-projected R_reg)
#   - (a) compute allele_match/allele_flip and filter BEFORE deduplicating merge_key,
#         and drop any multiallelic duplicates on the same unordered alleles
#   - (b) GWAS susie_rss(..., max_iter = 5000)
#   - (c) get_ld_per_locus_unique() applies PSD correlation projection to GWAS LD (R_gwas)
# Env: GWAS_NAME, CHR_FILTER, EQTL_SUSIE_DIR, COLOC_OUT_ROOT, COLOC_GENE_FILTER,
#      *_LD_DIR, RUN_DIR, IDENTITY_LOSS_DIR, LD_DIAG_DIR.
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
run_dir <- Sys.getenv("RUN_DIR")
if (!nzchar(run_dir)) stop("RUN_DIR must be set")
original <- file.path(root, "GWAS/finemapping/src/06_susie_coloc.R")
hash <- strsplit(system2("sha256sum", shQuote(original), stdout = TRUE), "[[:space:]]+")[[1]][1]
stopifnot(hash == "e22b3ffbef48de9b96e220399275edc0113ec751fc765238ccefd076697114cc")
source(file.path(root, "GWAS/finemapping/src/coloc_variant_identity.R"))

replace_once <- function(lines, old, new) {
  hit <- which(lines == old)
  if (length(hit) != 1L) stop("patch anchor not found exactly once: ", old)
  c(lines[seq_len(hit - 1L)], new, lines[-seq_len(hit)])
}
lines <- readLines(original, warn = FALSE)
lines <- coloc_identity_guarded_runner(lines)            # duplicate guard + Cholesky PD check
# (a) Harmonize alleles BEFORE merge_key deduplication and drop multiallelic unordered-allele duplicates
lines <- replace_once(lines,
  '    by = "merge_key"',
  '    by = "merge_key", allow.cartesian = TRUE')
dedup_idx <- which(lines == '  merged <- merged[order(gwas_pval)][!duplicated(merge_key)]')
stopifnot(length(dedup_idx) == 1L,
          lines[dedup_idx + 1L] == '  if (nrow(merged) < MIN_SNPS) { n_skipped <- n_skipped + 1L; next }')
lines[dedup_idx] <- '  # Dedup moved after allele_match | allele_flip filter and unordered-allele duplicate removal'
lines[dedup_idx + 1L] <- '  if (nrow(merged) == 0L) { n_skipped <- n_skipped + 1L; next }'
lines <- replace_once(lines,
  '  merged[allele_flip == TRUE, eqtl_beta := -eqtl_beta]',
  c('  merged[allele_flip == TRUE, eqtl_beta := -eqtl_beta]',
    '  if (nrow(merged) > 0L) {',
    '    merged[, .pair_key := paste(merge_key, pmin(gwas_a1, gwas_a2), pmax(gwas_a1, gwas_a2), sep = ":")]',
    '    merged <- merged[!(duplicated(.pair_key) | duplicated(.pair_key, fromLast = TRUE))]',
    '    merged[, .pair_key := NULL]',
    '    merged <- merged[order(gwas_pval)][!duplicated(merge_key)]',
    '  }',
    '  if (nrow(merged) < MIN_SNPS) { n_skipped <- n_skipped + 1L; next }'))
# (b) max_iter = 5000 for GWAS susie_rss
lines <- replace_once(lines,
  '              max_iter                   = 500',
  '              max_iter                   = 5000')
# (c) Use get_ld_per_locus_unique() (includes PSD correlation projection for R_gwas)
lines <- replace_once(lines, '          snp_ids <- paste0(chr_num, ":", ss_ld$position)',
                      '          snp_ids <- ss_ld$event_id')
lines <- replace_once(lines, '        ld_result <- get_ld_per_locus(',
                      '        ld_result <- get_ld_per_locus_unique(')
lines <- replace_once(lines, '      eqtl_positions <- as.integer(sub("^[0-9]+:", "", eqtl_snp_ids))',
                      '      eqtl_positions <- as.integer(sub("^[0-9]+:([0-9]+):.*$", "\\\\1", eqtl_snp_ids))')
lines <- replace_once(lines, '            if (!all(grepl("^[0-9]+:[0-9]+$", eqtl_snp_ids))) {',
                      '            if (!all(grepl("^[0-9]+:[0-9]+:[A-Za-z]+:[A-Za-z]+$", eqtl_snp_ids))) {')
lines <- replace_once(lines, 'source(file.path(FM_DIR, "src/coloc_signal_eligibility.R"))',
  c('source(file.path(FM_DIR, "src/coloc_signal_eligibility.R"))',
    sprintf('source("%s")', file.path(run_dir, "code/unique_ld_match.R"))))
patched_path <- file.path(run_dir, "code/patched_06_susie_coloc.R")
if (!file.exists(patched_path)) writeLines(lines, patched_path)
invisible(parse(text = lines))
if (nzchar(Sys.getenv("PATCH_ONLY"))) quit(save = "no", status = 0)
connection <- textConnection(lines)
tryCatch(source(connection, local = globalenv()), finally = close(connection))
