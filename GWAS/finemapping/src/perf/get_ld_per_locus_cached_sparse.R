# get_ld_per_locus_cached_sparse.R -- CACHE ONLY. Keeps bdiag()/sparse exactly as
# the canonical function has it. The single change is memoization.
#
# WHY THIS EXISTS: the first attempt (get_ld_per_locus_cached.R) changed TWO
# things at once -- memoization AND sparse->dense assembly -- and then measured
# one number, 3.04x. That is a confounded experiment: it cannot say which change
# bought the speed, and it came back NOT bit-identical (lambda_s_locus moved by
# up to 4.1e-14, because estimate_s_rss runs a Brent optimiser whose finite
# tolerance amplifies last-bit differences from a different linear-algebra path).
#
# This arm isolates the memoization. Prediction on record, before it runs:
#   * bit-identical to the reference, because nothing about the numbers or their
#     storage changes -- only how many times the same bytes are read;
#   * captures most of the 3.04x, because the cached region contains BOTH the
#     fread and the bdiag, which are the two expensive steps. What it will NOT
#     recover is faster per-gene subsetting of a dense vs sparse matrix.
# If it is bit-identical and close to 3x, it is strictly the better option and
# the dense change is not worth its provenance cost.

.LD_CACHE_SP <- new.env(parent = emptyenv())

get_ld_per_locus <- function(ss_per_locus, LOCUS, CHR, START, END, ancestry = "EUR") {
  path_to_LD <- find_LD_block(CHR, START, END, ancestry)
  if (is.null(path_to_LD)) {
    return(NULL)
  }

  key <- paste(c(ancestry, path_to_LD), collapse = "|")
  ent <- .LD_CACHE_SP$entry
  if (is.null(ent) || !identical(ent$key, key)) {
    .LD_CACHE_SP$entry <- NULL
    gc(FALSE)

    # ---- verbatim from finemapping_functions.R:159-184 ----------------------
    bim <- data.frame()
    ld <- matrix(ncol = 0, nrow = 0)
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
      ld <- bdiag(ld, as.matrix(LD))
    }

    if (nrow(bim) == 0) {
      warning(paste("No LD data loaded for locus", LOCUS))
      return(NULL)
    }

    colnames(bim) <- c("chr", "rsid", "dk", "pos", "alt", "ref")
    bim$SNP <- paste(bim$chr, bim$pos, bim$alt, bim$ref, sep = ":")
    colnames(ld) <- rownames(ld) <- bim$SNP

    ld <- ld[!duplicated(bim$SNP), !duplicated(bim$SNP)]
    bim <- bim[!duplicated(bim$SNP)]
    # ---- end verbatim -------------------------------------------------------

    .LD_CACHE_SP$entry <- list(key = key, bim = bim, ld = ld)
    ent <- .LD_CACHE_SP$entry
  }

  # bim gains an ss_matched column below, so hand out a copy. ld is only read via
  # ld[i, j], which allocates, so sharing it is safe.
  bim <- data.table::copy(ent$bim)
  ld <- ent$ld

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
