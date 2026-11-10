#!/usr/bin/env Rscript
# eQTL SuSiE refits for the affected primary genes, one-to-one allelic matching.
# Sources the hash-checked production 08_run_eqtl_susie.R with four text patches.
# Env: CHR, EQTL_SUSIE_DIR (run-owned), TARGET_GENES (TSV with ensembl, in_refit),
#      IDENTITY_LOSS_DIR, RUN_DIR.
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
run_dir <- Sys.getenv("RUN_DIR")
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
lines <- coloc_identity_guarded_eqtl_runner(lines)      # existing guard + seed 42
# Parity: drop the Nov 7 Cholesky LD gate; r2 production fits never applied it (check_R=FALSE).
gate <- grep("coloc_reference_ld_valid(R_mat", lines, fixed = TRUE)
stopifnot(length(gate) == 1L, trimws(lines[gate + 6L]) == "}")
lines <- lines[-(gate:(gate + 6L))]
stopifnot(!any(grepl("coloc_reference_ld_valid", lines, fixed = TRUE)))
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
