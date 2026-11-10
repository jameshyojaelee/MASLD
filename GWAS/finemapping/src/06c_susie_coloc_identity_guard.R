#!/usr/bin/env Rscript
# Successor entrypoint: preserve the read-only r2 producer while preventing
# duplicate position IDs from reaching SuSiE/coloc in new isolated runs.
# Usage is unchanged: GWAS_NAME / CHR_FILTER or the two original arguments.
root <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
original <- file.path(root, "GWAS/finemapping/src/06_susie_coloc.R")
hash <- strsplit(system2("sha256sum", shQuote(original), stdout = TRUE), "[[:space:]]+")[[1]][1]
stopifnot(hash == "e22b3ffbef48de9b96e220399275edc0113ec751fc765238ccefd076697114cc")
source(file.path(root, "GWAS/finemapping/src/coloc_variant_identity.R"))
lines <- coloc_identity_guarded_runner(readLines(original, warn = FALSE))
connection <- textConnection(lines)
tryCatch(source(connection, local = globalenv()), finally = close(connection))
