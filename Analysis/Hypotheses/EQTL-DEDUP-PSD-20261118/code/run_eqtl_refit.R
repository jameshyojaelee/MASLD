#!/usr/bin/env Rscript
# eQTL SuSiE refits (affected primary + 5 unaffected control genes), one-to-one allelic matching + PSD LD.
# Sources the hash-checked production 08_run_eqtl_susie.R with patches:
#   - coloc_identity_guarded_eqtl_runner (duplicate refusal, Cholesky LD check, set.seed(42))
#   - max_iter = 5000 (matching src/perf/08_refit_nonconverged.R)
#   - set.seed(42) preserved immediately before susie_rss()
#   - get_ld_per_locus_unique() with PSD correlation projection + pre_min_eig logging
#   - OUTPUT_DIR defaults to EQTL-DEDUP-PSD-20261118/eqtl_susie/
# Env: CHR, EQTL_SUSIE_DIR, TARGET_GENES, IDENTITY_LOSS_DIR, LD_DIAG_DIR, RUN_DIR.
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
run_dir <- Sys.getenv("RUN_DIR")
if (!nzchar(run_dir)) stop("RUN_DIR must be set")
original <- file.path(root, "GWAS/finemapping/src/08_run_eqtl_susie.R")
hash <- strsplit(system2("sha256sum", shQuote(original), stdout = TRUE), "[[:space:]]+")[[1]][1]
stopifnot(hash == "c6f24b8e7716f54e4bbf11d87bb0eb0ba1174d7f2e4fb78742ef4d10df0bf1cd")
source(file.path(root, "GWAS/finemapping/src/coloc_variant_identity.R"))

replace_once <- function(lines, old, new) {
  hit <- which(lines == old)
  if (length(hit) != 1L) stop("patch anchor not found exactly once: ", old)
  c(lines[seq_len(hit - 1L)], new, lines[-seq_len(hit)])
}
lines <- readLines(original, warn = FALSE)
lines <- coloc_identity_guarded_eqtl_runner(lines)      # guard + Cholesky check + set.seed(42)
# Default EQTL_SUSIE_DIR to EQTL-DEDUP-PSD-20261118/eqtl_susie/
lines <- replace_once(lines,
  'OUTPUT_DIR <- Sys.getenv("EQTL_SUSIE_DIR", unset = file.path(FM_DIR, "results/eqtl_susie"))',
  sprintf('OUTPUT_DIR <- Sys.getenv("EQTL_SUSIE_DIR", unset = "%s")',
          file.path(root, "Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/eqtl_susie")))
# Place set.seed(42) immediately before susie_rss() as well (no RNG calls intervene)
lines <- replace_once(lines,
  '  s <- tryCatch(',
  c('  set.seed(42)',
    '  s <- tryCatch('))
lines <- replace_once(lines,
  '      max_iter                   = 500',
  '      max_iter                   = 5000')
lines <- replace_once(lines, '  snp_ids <- paste0(CHR, ":", ss_filtered$position)',
                      '  snp_ids <- ss_filtered$event_id')
lines <- replace_once(lines, '    get_ld_per_locus(ss, LOCUS = ensg_id, CHR = CHR, START = START, END = END, ancestry = "EUR"),',
                      '    get_ld_per_locus_unique(ss, LOCUS = ensg_id, CHR = CHR, START = START, END = END, ancestry = "EUR"),')
lines <- replace_once(lines, 'source(file.path(FM_DIR, "src/finemapping_functions.R"))',
  c('source(file.path(FM_DIR, "src/finemapping_functions.R"))',
    sprintf('source("%s")', file.path(run_dir, "code/unique_ld_match.R"))))
lines <- replace_once(lines, 'ensg_list  <- dt[, unique(ENSG)]',
  c('ensg_list  <- dt[, unique(ENSG)]',
    'target_tab <- fread(Sys.getenv("TARGET_GENES"))',
    'ensg_list  <- intersect(ensg_list, target_tab[in_refit == TRUE & chr == CHR, ensembl])',
    'stopifnot(length(ensg_list) == target_tab[in_refit == TRUE & chr == CHR, .N])'))
writeLines(lines, file.path(run_dir, "code/patched_08_run_eqtl_susie.R"))
invisible(parse(text = lines))
if (nzchar(Sys.getenv("PATCH_ONLY"))) quit(save = "no", status = 0)
connection <- textConnection(lines)
tryCatch(source(connection, local = globalenv()), finally = close(connection))
