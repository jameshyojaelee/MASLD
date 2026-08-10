# get_ld_per_locus_cached.R -- performance override for get_ld_per_locus().
#
# Sourced AFTER src/finemapping_functions.R to shadow that definition. The
# canonical file is not modified, so nothing outside a run that explicitly
# sources this file can be affected.
#
# WHY (all measured on logs/coloc_rerun/coloc_19648525_66.out, chr1, 1,671 genes):
#
#   * Each gene re-reads its LD block(s) from GPFS. 171 distinct block-sets serve
#     1,671 gene-loads -- consecutive genes share a block-set (run lengths:
#     median 6, mean 9.8, max 62) because eGenes arrive in genomic order. A
#     one-slot cache therefore captures essentially all reuse: 1,671 loads -> 171.
#
#   * bdiag() stores the result sparse. The matrices are 0.0% sparse (measured:
#     66/66 fields nonzero in a sampled .ld row), so sparse costs 12 bytes/element
#     against 8 dense, and makes the two per-gene subsets far slower. Dense
#     block-diagonal assembly produces identical VALUES.
#
#   * Scale of the waste: median LD matrix loaded is 23,644 x 23,644 (4.5 GB
#     dense) while SuSiE fits on a median of 3,567 matched variants -- 49x more
#     matrix than the model uses.
#
# IDENTITY: this changes only (a) how many times the same bytes are read and
# (b) dense vs sparse storage of the same numbers. Everything from the
# summary-stat matching onward is byte-for-byte the original code. That claim is
# not asserted, it is tested -- see src/perf/verify_cache_identity.sbatch, which
# reruns a COMPLETED task under both the original and this override and diffs
# the outputs, with an unpatched arm as the determinism control.
#
# MEMORY: the cache holds ONE block-set. Median 4.5 GB, worst observed M=47,064
# -> 17.7 GB. Peak during assembly is ~2-3x that for the largest blocks, which
# is why the verification runs at the same 96G as production.

.LD_CACHE <- new.env(parent = emptyenv())

get_ld_per_locus <- function(ss_per_locus, LOCUS, CHR, START, END, ancestry = "EUR") {
  path_to_LD <- find_LD_block(CHR, START, END, ancestry)
  if (is.null(path_to_LD)) {
    return(NULL)
  }

  # ---- cached region: everything here depends only on the block paths --------
  key <- paste(c(ancestry, path_to_LD), collapse = "|")
  ent <- .LD_CACHE$entry
  if (is.null(ent) || !identical(ent$key, key)) {
    .LD_CACHE$entry <- NULL   # drop the old block BEFORE allocating the new one
    gc(FALSE)

    bim <- data.frame()
    mats <- list()
    for (file in path_to_LD) {
      bim_file <- paste0(file, ".bim")
      ld_file <- paste0(file, ".ld")
      if (!file.exists(bim_file) || !file.exists(ld_file)) {
        warning(paste("LD files not found:", file))
        next
      }
      BIM <- fread(bim_file)
      LD <- fread(ld_file)
      bim <- rbind(bim, BIM)
      mats[[length(mats) + 1L]] <- as.matrix(LD)
    }

    if (nrow(bim) == 0) {
      warning(paste("No LD data loaded for locus", LOCUS))
      return(NULL)
    }

    # Dense block-diagonal. Identical values to bdiag(); see header for why the
    # sparse representation is the wrong one for a 0%-sparse matrix.
    if (length(mats) == 1L) {
      ld <- mats[[1L]]
    } else {
      tot <- sum(vapply(mats, nrow, integer(1)))
      ld <- matrix(0, nrow = tot, ncol = tot)
      off <- 0L
      for (m in mats) {
        ix <- off + seq_len(nrow(m))
        ld[ix, ix] <- m
        off <- off + nrow(m)
      }
    }
    rm(mats)

    colnames(bim) <- c("chr", "rsid", "dk", "pos", "alt", "ref")
    bim$SNP <- paste(bim$chr, bim$pos, bim$alt, bim$ref, sep = ":")
    colnames(ld) <- rownames(ld) <- bim$SNP

    # Original always paid this subset. Skip it when it is a no-op -- that is a
    # full M^2 copy avoided, and when it is NOT a no-op the behaviour is identical.
    keep_uniq <- !duplicated(bim$SNP)
    if (!all(keep_uniq)) {
      ld <- ld[keep_uniq, keep_uniq]
      bim <- bim[keep_uniq]
    }

    .LD_CACHE$entry <- list(key = key, bim = bim, ld = ld)
    ent <- .LD_CACHE$entry
  }

  # bim is MUTATED below (ss_matched column), so hand out a copy or the cache is
  # poisoned for the next gene. ld is only ever read via ld[i, j], which
  # allocates a new matrix, so it is safe to share.
  bim <- data.table::copy(ent$bim)
  ld <- ent$ld
  # ---- end cached region; everything below is the original code verbatim -----

  # Match summary stats to LD variants (handle allele flips)
  ss_per_locus$SNP1 <- paste(ss_per_locus$chromosome,
                             ss_per_locus$position,
                             ss_per_locus$allele1,
                             ss_per_locus$allele2, sep = ":")
  ss_per_locus$SNP2 <- paste(ss_per_locus$chromosome,
                             ss_per_locus$position,
                             ss_per_locus$allele2,
                             ss_per_locus$allele1, sep = ":")

  bim$ss_matched <- bim$SNP %in% c(ss_per_locus$SNP2, ss_per_locus$SNP1)
  n_matched <- sum(bim$ss_matched)
  cat(paste0("  Matched ", n_matched, " / ", nrow(bim), " LD variants to summary stats\n"))

  if (n_matched < 10) {
    warning(paste("Too few matched variants for locus", LOCUS, ":", n_matched))
    return(NULL)
  }

  ld_filtered <- ld[bim$ss_matched, bim$ss_matched]
  bim_filtered <- bim[bim$ss_matched, ]

  # NOTE: rbind puts direct matches (SNP1) first, so !duplicated keeps direct match
  # over flipped match. This is intentional — direct allele orientation is preferred.
  ss_filtered <- rbind(
    ss_per_locus %>%
      filter(SNP1 %in% bim_filtered$SNP) %>%
      rename(SNP = SNP1) %>%
      select(!c(SNP2)),
    ss_per_locus %>%
      filter(SNP2 %in% bim_filtered$SNP) %>%
      mutate(beta = -1 * beta) %>%
      rename(SNP = SNP2, allele1 = allele2, allele2 = allele1) %>%
      select(!c(SNP1))
  )

  # Deduplicate (keep first match per SNP) and align to bim_filtered order
  ss_filtered <- ss_filtered[!duplicated(ss_filtered$SNP), ]
  ss_filtered <- ss_filtered[match(bim_filtered$SNP, ss_filtered$SNP), ]
  ss_filtered <- ss_filtered[!is.na(ss_filtered$SNP), ]

  # Subset LD to only variants present in both SS and LD
  keep <- bim_filtered$SNP %in% ss_filtered$SNP
  if (sum(keep) < 10) {
    warning(paste("Too few variants after alignment for locus", LOCUS))
    return(NULL)
  }
  bim_filtered <- bim_filtered[keep, ]
  ld_filtered <- as.matrix(ld_filtered)[keep, keep]
  ss_filtered <- ss_filtered[match(bim_filtered$SNP, ss_filtered$SNP), ]

  # Handle NAs in LD matrix — iteratively remove high-NA variants
  if (sum(is.na(ld_filtered)) != 0) {
    cat("  LD matrix contains NAs —")
    ld_mat <- as.matrix(ld_filtered)
    n_before <- nrow(ld_mat)
    for (iter in seq_len(100)) {
      na_per_var <- colSums(is.na(ld_mat))
      if (max(na_per_var) == 0) break
      # Remove variants with >50% NAs first (monomorphic/low-quality)
      bad <- na_per_var > (nrow(ld_mat) * 0.5)
      if (sum(bad) == 0) bad <- na_per_var > 0  # then any remaining NA variants
      ld_mat <- ld_mat[!bad, !bad, drop = FALSE]
      ss_filtered <- ss_filtered[!bad, ]
      if (nrow(ld_mat) < 10) break
    }
    cat(" removed", n_before - nrow(ld_mat), "bad variants,", nrow(ld_mat), "remaining\n")
    ld_filtered <- ld_mat
    if (nrow(ld_filtered) < 10) {
      warning(paste("Too few variants after NA removal for locus", LOCUS))
      return(NULL)
    }
  }

  # Return LD as dense matrix (avoid extra data.frame copy — saves ~800MB per large locus)
  if (!is.matrix(ld_filtered)) ld_filtered <- as.matrix(ld_filtered)
  return(list(ss_filtered, ld_filtered))
}
