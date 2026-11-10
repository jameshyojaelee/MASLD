#!/usr/bin/env Rscript
# SuSiE-COLOC for the affected primary genes on the refit eQTL fits.
# Sources the hash-checked production 06_susie_coloc.R with five text patches.
# Env as in r2 (GWAS_NAME, CHR_FILTER, EQTL_SUSIE_DIR, COLOC_OUT_ROOT,
# COLOC_GENE_FILTER, *_LD_DIR) plus RUN_DIR and IDENTITY_LOSS_DIR.
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
run_dir <- Sys.getenv("RUN_DIR")
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
lines <- coloc_identity_guarded_runner(lines)            # existing guard
# Parity: drop the Nov 7 Cholesky LD gate; r2 production fits never applied it (check_R=FALSE).
gate <- grep("if (!coloc_reference_ld_valid(R_reg))", lines, fixed = TRUE)
stopifnot(length(gate) == 1L, trimws(lines[gate + 3L]) == "}")
lines <- lines[-(gate:(gate + 3L))]
stopifnot(!any(grepl("coloc_reference_ld_valid", lines, fixed = TRUE)))
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
writeLines(lines, file.path(run_dir, "code/patched_06_susie_coloc.R"))
invisible(parse(text = lines))
if (nzchar(Sys.getenv("PATCH_ONLY"))) quit(save = "no", status = 0)
connection <- textConnection(lines)
tryCatch(source(connection, local = globalenv()), finally = close(connection))
