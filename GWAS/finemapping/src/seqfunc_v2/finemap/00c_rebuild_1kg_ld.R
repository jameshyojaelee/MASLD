#!/usr/bin/env Rscript

# Rebuild the 1000 Genomes per-block LD matrices as mean-imputed Gram matrices.
#
# WHY. The shipped `.ld` files were produced with `plink --r square`, whose
# estimator is not the genotype dosage correlation and is not the Gram matrix of
# any single dataset -- so it can be, and is, indefinite. Measured on
# 1kg_afr/chr8/9154694.9640787 (661 samples, 3,361 variants, 0.0043% missing
# calls, 6-significant-figure storage, so neither missingness nor rounding
# explains it):
#
#     shipped  min eigenvalue = -0.189664   -> needs ridge 0.300
#     rebuilt  min eigenvalue = -0.000000   -> needs ridge 0.001
#
# -0.000000 is exactly what rank deficiency should give: with p >> n the true
# correlation matrix is singular but positive SEMI-definite. Small sample size
# cannot produce negative eigenvalues; only a non-Gram estimator can.
#
# This matters because the ridge required to factorise the shipped matrices
# (median 0.300 for AFR/AMR/SAS) exceeds RIDGE_LAMBDA_PRIMARY_MAX, so every
# African, Admixed-American and South-Asian locus was excluded from the anchor
# pool. That exclusion was an artifact of panel construction, not of ancestry or
# of GWAS power.
#
# SAFETY. The originals are SHARED with the canonical COLOC pipeline
# (finemapping_functions.R get_ld_base_dir()). This writes to a parallel root
# `1kg_<pop>_gram/` and never touches `1kg_<pop>/`.
#
# Usage: Rscript 00c_rebuild_1kg_ld.R <block_row_index>
#        (row of config/used_ld_blocks.tsv restricted to 1kg panels)

source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT", unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "GWAS", "finemapping", "src", "seqfunc_v2", "finemap", "common.R"
))

PLINK2 <- Sys.getenv("PLINK2_BIN", unset = "/nfs/sw/plink/plink-2.0a5.13/plink")
LD_REF <- file.path(FM_ROOT, "data", "ld_ref")
# Six decimals matches the storage precision of the shipped files, so the
# rebuild is not penalised or flattered by a format change.
LD_DIGITS <- as.integer(Sys.getenv("LD_REBUILD_DIGITS", unset = "6"))

args <- commandArgs(trailingOnly = TRUE)
row_index <- suppressWarnings(as.integer(args[[1L]]))
assert_true(!is.na(row_index), "Usage: 00c_rebuild_1kg_ld.R <block_row_index>")

# Two modes.  LD_REBUILD_POP=<pop> rebuilds that panel's FULL block manifest --
# needed by canonical COLOC, which attempts SuSiE for every gene with a converged
# eQTL fit and therefore touches essentially every block genome-wide.  Without it
# the index addresses only the blocks a given seqfunc run consumed.
rebuild_pop <- Sys.getenv("LD_REBUILD_POP", unset = "")
if (nzchar(rebuild_pop)) {
  assert_true(rebuild_pop %in% c("afr", "amr", "eas", "sas"),
              "LD_REBUILD_POP must be one of afr/amr/eas/sas (got '%s')", rebuild_pop)
  man <- file.path(LD_REF, paste0("1kg_", rebuild_pop), "approx_LD_blocks.txt")
  assert_true(file.exists(man), "Missing block manifest: %s", man)
  bl <- fread(man)
  assert_true(identical(names(bl), c("chr", "start", "stop")), "Unexpected block schema: %s", man)
  setorder(bl, chr, start)
  assert_true(row_index >= 1L && row_index <= nrow(bl),
              "block_row_index out of range: %d (have %d)", row_index, nrow(bl))
  r <- bl[row_index]
  blk_prefix <- file.path(LD_REF, paste0("1kg_", rebuild_pop), paste0("chr", r$chr),
                          paste0(r$start, ".", r$stop), paste0(r$start, ".", r$stop))
  b <- data.table(ld_panel_id = paste0("1kg_phase3_", rebuild_pop), chromosome = as.integer(r$chr),
                  block_start = as.integer(r$start), block_stop = as.integer(r$stop),
                  block_prefix = blk_prefix)
} else {
  blocks <- fread(file.path(RUN_ROOT, "config", "used_ld_blocks.tsv"))
  blocks <- blocks[grepl("^1kg_", ld_panel_id)]
  setorder(blocks, ld_panel_id, chromosome, block_start)
  assert_true(row_index >= 1L && row_index <= nrow(blocks),
              "block_row_index out of range: %d (have %d)", row_index, nrow(blocks))
  b <- blocks[row_index]
}

pop <- sub("^1kg_phase3_", "", b$ld_panel_id)
src_bim <- paste0(b$block_prefix, ".bim")
if (!file.exists(src_bim) || !file.exists(paste0(b$block_prefix, ".ld"))) {
  cat(sprintf("Skip: no source block at %s\n", b$block_prefix)); quit(save = "no", status = 0L)
}
chr_prefix <- file.path(LD_REF, paste0("1kg_", pop), sprintf("chr%d_%s", b$chromosome, pop))
for (ext in c(".bed", ".bim", ".fam")) {
  assert_true(file.exists(paste0(chr_prefix, ext)), "Missing panel file: %s", paste0(chr_prefix, ext))
}

out_dir <- file.path(LD_REF, paste0("1kg_", pop, "_gram"),
                     sprintf("chr%d", b$chromosome),
                     sprintf("%d.%d", b$block_start, b$block_stop))
ensure_dirs(out_dir)
out_stem <- file.path(out_dir, sprintf("%d.%d", b$block_start, b$block_stop))
out_ld <- paste0(out_stem, ".ld")
out_bim <- paste0(out_stem, ".bim")
qc_path <- paste0(out_stem, ".rebuild_qc.json")
if (file.exists(qc_path) && file.exists(out_ld)) {
  cat(sprintf("Idempotent skip: %s\n", out_stem)); quit(save = "no", status = 0L)
}

bim <- fread(src_bim, header = FALSE)
assert_true(ncol(bim) >= 6L, "BIM has fewer than six columns: %s", src_bim)
ids <- as.character(bim$V2)
cat(sprintf("[%s chr%d %d-%d] %d variants\n", pop, b$chromosome, b$block_start, b$block_stop, length(ids)))

tmp <- file.path(tempdir(), sprintf("ldrb_%s_%d_%d", pop, b$chromosome, Sys.getpid()))
on.exit(unlink(paste0(tmp, c(".snps", ".raw", ".log"))), add = TRUE)
writeLines(ids, paste0(tmp, ".snps"))
st <- system2(PLINK2, c("--bfile", chr_prefix, "--extract", paste0(tmp, ".snps"),
                        "--export", "A", "--out", tmp), stdout = FALSE, stderr = FALSE)
assert_true(identical(as.integer(st), 0L), "plink2 --export A failed for %s", out_stem)

raw <- fread(paste0(tmp, ".raw"))
G <- as.matrix(raw[, -(1:6)])
# plink2 appends _<countedallele> to each exported column name.
colnames(G) <- sub("_[^_]*$", "", colnames(G))
# Reindex to the BLOCK BIM order: 02_run_susie_locus.R addresses the matrix by
# BIM row (ld_index), so any reordering here would silently mis-map every variant.
idx <- match(ids, colnames(G))
assert_true(!anyNA(idx), "%d block variants absent from the exported genotypes for %s",
            sum(is.na(idx)), out_stem)
G <- G[, idx, drop = FALSE]
n_samples <- nrow(G)
frac_missing <- mean(is.na(G))

# Mean-impute so the covariance is the Gram matrix of a complete dataset, which
# is positive semi-definite by construction.
for (j in seq_len(ncol(G))) {
  m <- is.na(G[, j])
  if (any(m)) G[m, j] <- mean(G[, j], na.rm = TRUE)
}
sdv <- apply(G, 2, sd)
poly <- sdv > 0
R <- matrix(NA_real_, ncol(G), ncol(G))
R[poly, poly] <- cor(G[, poly, drop = FALSE])
R <- round(R, LD_DIGITS)
diag(R)[poly] <- 1

# QC: how much ridge does the rebuilt matrix actually need?
Rf <- R[poly, poly, drop = FALSE]
needed <- NA_real_
for (l in SUSIE_RIDGE_LADDER) {
  if (tryCatch({ invisible(chol(Rf + l * diag(nrow(Rf)))); TRUE }, error = function(e) FALSE)) {
    needed <- l; break
  }
}

tmp_ld <- paste0(out_ld, ".tmp.", Sys.getpid())
fwrite(as.data.table(R), tmp_ld, sep = " ", col.names = FALSE, na = "nan", quote = FALSE)
assert_true(file.rename(tmp_ld, out_ld), "Atomic rename failed: %s", out_ld)
invisible(file.copy(src_bim, out_bim, overwrite = TRUE))

atomic_json(list(
  panel_pop = pop, chromosome = b$chromosome,
  block_start = b$block_start, block_stop = b$block_stop,
  n_variants = ncol(G), n_polymorphic = sum(poly), n_samples = n_samples,
  fraction_missing_calls = frac_missing,
  ridge_needed_rebuilt = needed,
  estimator = "mean_imputed_genotype_dosage_correlation",
  rationale = "plink --r square is not a Gram matrix and is indefinite; this is PSD by construction",
  digits = LD_DIGITS,
  source_bim = src_bim, source_bim_sha256 = sha256_file(src_bim),
  out_ld_sha256 = sha256_file(out_ld),
  canonical_outputs_mutated = FALSE
), qc_path)

cat(sprintf("  n=%d  polymorphic=%d/%d  missing=%.4f%%  ridge_needed=%s  -> %s\n",
            n_samples, sum(poly), ncol(G), 100 * frac_missing, format(needed), out_ld))
