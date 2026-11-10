#!/usr/bin/env Rscript
# Use only a new explicit output directory; historical eQTL fits are preserved.
root <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
destination <- Sys.getenv("EQTL_SUSIE_DIR")
if (!nzchar(destination) || dir.exists(destination)) {
  stop("EQTL_SUSIE_DIR must name a new output directory; do not reuse historical fits")
}
original <- file.path(root,"GWAS/finemapping/src/08_run_eqtl_susie.R")
hash <- strsplit(system2("sha256sum",shQuote(original),stdout=TRUE),"[[:space:]]+")[[1]][1]
stopifnot(hash=="c6f24b8e7716f54e4bbf11d87bb0eb0ba1174d7f2e4fb78742ef4d10df0bf1cd")
source(file.path(root,"GWAS/finemapping/src/coloc_variant_identity.R"))
lines <- coloc_identity_guarded_eqtl_runner(readLines(original,warn=FALSE))
connection <- textConnection(lines)
tryCatch(source(connection,local=globalenv()),finally=close(connection))
