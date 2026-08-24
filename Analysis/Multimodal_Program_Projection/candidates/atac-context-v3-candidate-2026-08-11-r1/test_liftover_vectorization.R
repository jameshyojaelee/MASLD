#!/usr/bin/env Rscript
# Equivalence fixture: sealed per-element loop vs vectorized replacement.
# Both implementations are reproduced verbatim from their source files.
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer)
})
OUTDIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/atac-v3-genetics-upstream-manifest-20260819T144901Z/liftover_vectorization_draft"
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CAND <- file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1")

sealed_impl <- function(lifted, variants) {
  rbindlist(lapply(seq_along(lifted), function(index) {
    value <- lifted[[index]]
    if (length(value) == 0L) {
      return(data.table(
        variant_key = variants$variant_key[index], n_liftover_mappings = 0L,
        hg38_chrom = NA_character_, hg38_position_1based = NA_integer_
      ))
    }
    data.table(
      variant_key = variants$variant_key[index],
      n_liftover_mappings = length(value),
      hg38_chrom = as.character(seqnames(value))[1L],
      hg38_position_1based = start(value)[1L]
    )
  }))
}

vectorized_impl <- function(lifted, variants) {
  n_liftover_mappings <- as.integer(elementNROWS(lifted))
  flat <- unlist(lifted, use.names = FALSE)
  first_index <- head(cumsum(c(1L, n_liftover_mappings)), -1L)
  mapped <- n_liftover_mappings > 0L
  hg38_chrom <- rep(NA_character_, length(lifted))
  hg38_position_1based <- rep(NA_integer_, length(lifted))
  if (any(mapped)) {
    flat_chrom <- as.character(seqnames(flat))
    flat_start <- start(flat)
    hg38_chrom[mapped] <- flat_chrom[first_index[mapped]]
    hg38_position_1based[mapped] <- flat_start[first_index[mapped]]
  }
  data.table(
    variant_key = variants$variant_key,
    n_liftover_mappings = n_liftover_mappings,
    hg38_chrom = hg38_chrom,
    hg38_position_1based = hg38_position_1based
  )
}

fail <- 0L
check <- function(label, ok) {
  cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) fail <<- fail + 1L
}

# ---- Test A: synthetic, forced edge cases -------------------------------
cat("Test A: synthetic GRangesList with 0 / 1 / many mappings\n")
set.seed(17)
N <- 4000L
counts <- sample(c(0L, 1L, 1L, 1L, 2L, 3L), N, replace = TRUE)
total <- sum(counts)
flat_gr <- GRanges(
  seqnames = sample(paste0("chr", 1:22), total, TRUE),
  ranges = IRanges(start = sample(1e6:9e7, total, TRUE), width = 1L)
)
grl <- relist(flat_gr, PartitioningByEnd(cumsum(counts)))
variants <- data.table(variant_key = sprintf("v%06d", seq_len(N)))
a_old <- sealed_impl(grl, variants); a_new <- vectorized_impl(grl, variants)
check("column names identical", identical(names(a_old), names(a_new)))
check("row count identical", identical(nrow(a_old), nrow(a_new)))
check("all.equal exact", isTRUE(all.equal(a_old, a_new, check.attributes = FALSE)))
check("identical()", identical(a_old, a_new))
cat(sprintf("    (%d unmapped, %d single, %d multi)\n",
            sum(counts == 0L), sum(counts == 1L), sum(counts > 1L)))

# ---- Test B: real chain, real replayed variants -------------------------
cat("Test B: frozen hg19->hg38 chain on real replayed variants\n")
chain_path <- file.path(BASE, "data/broadaway_eqtl/hg19ToHg38.over.chain")
files <- Sys.glob(file.path(CAND, "genetics/replay_execution/batch_*/exports/*/chr*/*/variant_posteriors.tsv.gz"))
stopifnot(length(files) > 0L)
set.seed(3)
sample_files <- sample(files, min(40L, length(files)))
post <- rbindlist(lapply(sample_files, fread), fill = TRUE)
post[, hg19_chrom := paste0("chr", chr)]
post[, variant_key := paste(hg19_chrom, hg19_position, allele1, allele2, sep = ":")]
variants <- unique(post[, .(variant_key, hg19_chrom, hg19_position)], by = "variant_key")
ranges <- GRanges(seqnames = variants$hg19_chrom,
                  ranges = IRanges(start = variants$hg19_position, width = 1L))
names(ranges) <- variants$variant_key
lifted <- liftOver(ranges, import.chain(chain_path))
t_old <- system.time(b_old <- sealed_impl(lifted, variants))["elapsed"]
t_new <- system.time(b_new <- vectorized_impl(lifted, variants))["elapsed"]
check("identical() on real data", identical(b_old, b_new))
cat(sprintf("    variants=%d  unmapped=%d  multi=%d\n", nrow(variants),
            sum(b_new$n_liftover_mappings == 0L), sum(b_new$n_liftover_mappings > 1L)))
cat(sprintf("    sealed %.2fs vs vectorized %.3fs  =>  %.0fx faster\n",
            t_old, t_new, t_old / max(t_new, 1e-6)))

cat(if (fail == 0L) "\nLIFTOVER_VECTORIZATION_FIXTURE\tPASS\n" else
    sprintf("\nLIFTOVER_VECTORIZATION_FIXTURE\tFAIL\t%d\n", fail))
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
quit(status = if (fail == 0L) 0L else 1L)
