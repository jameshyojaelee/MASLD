#!/usr/bin/env Rscript
# 04d_reexport_carma_txt.R
# Re-exports CARMA .txt.gz files from .rds backups
#
# PROBLEM: data.table::fwrite(compress="gzip") produced 0-byte files on this
# system because zlib headers were not found at compile time. The .rds files
# contain the actual CARMA results. This script reconstructs the tabular output
# from .rds + SuSiE tsv (for LD-matched variant metadata).
#
# STRATEGY:
#   1. Scan all CARMA directories for .txt.gz files that are 0 bytes or unreadable
#   2. For each, load the .rds (CARMA result object) and the SuSiE .tsv
#      (which has the same LD-matched variant set in the same order)
#   3. Reconstruct the CARMA output table (SNP, chr, pos, a1, a2, beta, se, pval,
#      locus, Z, PIP, CS, outlier) and write via gzfile() + write.table()
#
# Usage: Rscript src/04d_reexport_carma_txt.R

suppressPackageStartupMessages({
  library(data.table)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
setwd(FM_DIR)

cat("============================================================\n")
cat("Re-exporting CARMA .txt.gz from .rds backups\n")
cat("============================================================\n")

# ---- Step 1: Find all CARMA .txt.gz files ----
all_txt_gz <- list.files(file.path(FM_DIR, "output"),
                         pattern = "\\.txt\\.gz$",
                         recursive = TRUE, full.names = TRUE)
# Filter to CARMA directories only
all_txt_gz <- all_txt_gz[grepl("/CARMA/", all_txt_gz)]
cat("Total CARMA .txt.gz files found:", length(all_txt_gz), "\n")

# ---- Step 2: Identify broken files (0 bytes or fread fails) ----
is_broken <- logical(length(all_txt_gz))

for (i in seq_along(all_txt_gz)) {
  f <- all_txt_gz[i]
  if (file.size(f) == 0) {
    is_broken[i] <- TRUE
    next
  }
  # Try reading — catch corrupt gzip
  tryCatch({
    dt <- fread(f, nrows = 1)
    if (ncol(dt) == 0 || nrow(dt) == 0) is_broken[i] <- TRUE
  }, error = function(e) {
    is_broken[i] <<- TRUE
  })
}

broken_files <- all_txt_gz[is_broken]
cat("Broken/empty CARMA .txt.gz files:", length(broken_files), "\n")

if (length(broken_files) == 0) {
  cat("Nothing to re-export. Exiting.\n")
  quit(save = "no", status = 0)
}

# ---- Step 3: Re-export each broken file ----
n_success <- 0L
n_fail <- 0L
failures <- character(0)

for (f in broken_files) {
  tryCatch({
    carma_dir <- dirname(f)
    base <- sub("\\.txt\\.gz$", "", basename(f))

    # Parse the filename: {study}_{LDpanel}_{window}Mb_{locus}
    # LDpanel is UKBB (EUR) or 1KG_EAS (EAS)
    # Anchor to the trailing pattern (regmatches) so studies whose names
    # already contain "_UKBB_" (e.g., 2023_36280732_NAFLD_UKBB_EUR) split on
    # the *last* delimiter, not the first.
    m <- regmatches(base, regexec(
      "^(.+)_(UKBB|1KG_EAS)_([0-9.]+)Mb_(.+)$", base))[[1]]
    if (length(m) == 5) {
      study   <- m[2]
      ldpanel <- m[3]
      rest    <- paste0(m[4], "Mb_", m[5])
    } else if (grepl("_UKBB_", base)) {
      # Fallback for legacy parser path
      ldpanel <- "UKBB"
      parts <- strsplit(base, "_UKBB_")[[1]]
      study <- parts[1]
      rest <- parts[2]  # e.g., "0.5Mb_10.113949664"
    } else if (grepl("_1KG_EAS_", base)) {
      ldpanel <- "1KG_EAS"
      parts <- strsplit(base, "_1KG_EAS_")[[1]]
      study <- parts[1]
      rest <- parts[2]
    } else {
      stop("Cannot parse LDpanel from filename: ", base)
    }

    # rest = "{window}Mb_{locus}" -> split on "Mb_"
    rest_parts <- strsplit(rest, "Mb_")[[1]]
    window <- rest_parts[1]  # e.g., "0.5"
    locus <- rest_parts[2]   # e.g., "10.113949664"

    cat("\n--- Re-exporting:", base, "---\n")
    cat("  Study:", study, " LDpanel:", ldpanel, " Window:", window,
        "Mb  Locus:", locus, "\n")

    # Locate the .rds file (same directory, same name)
    rds_path <- file.path(carma_dir, paste0(base, ".rds"))
    if (!file.exists(rds_path)) {
      stop("RDS not found: ", rds_path)
    }

    # Locate the SuSiE .tsv (variant metadata in LD-matched order)
    susie_dir <- sub("/CARMA$", "/susie", carma_dir)
    susie_pattern <- paste0("^", ldpanel, "_", locus, "_cov0\\.95_")
    susie_candidates <- list.files(susie_dir, pattern = susie_pattern,
                                   full.names = TRUE)
    susie_candidates <- susie_candidates[grepl("\\.tsv$", susie_candidates)]

    if (length(susie_candidates) == 0) {
      stop("No SuSiE .tsv found in ", susie_dir, " matching ", susie_pattern)
    }
    susie_path <- susie_candidates[1]

    # Load data
    rds <- readRDS(rds_path)
    susie_dt <- fread(susie_path)

    # Validate dimensions
    n_pips <- length(rds[[1]]$PIPs)
    n_susie <- nrow(susie_dt)

    if (n_pips == n_susie) {
      # Primary path: SuSiE tsv has matching variant count
      result <- data.table(
        SNP        = susie_dt$SNP,
        chromosome = susie_dt$chromosome,
        position   = susie_dt$position,
        allele1    = susie_dt$allele1,
        allele2    = susie_dt$allele2,
        beta       = susie_dt$beta,
        se         = susie_dt$se,
        pval       = susie_dt$pval,
        locus      = susie_dt$locus,
        Z          = susie_dt$beta / susie_dt$se,
        PIP        = rds[[1]]$PIPs,
        CS         = 0L,
        outlier    = 0L
      )
    } else {
      # Fallback: SuSiE tsv has different variant set (e.g., different LD filter).
      # Use SuSiE RDS pip names (chr:pos:a1:a2) to reconstruct variant metadata,
      # then join with SS file for beta/se/pval.
      cat("  Dimension mismatch (PIPs=", n_pips, ", SuSiE tsv=", n_susie,
          ") — using SuSiE RDS fallback\n")

      susie_rds_path <- file.path(susie_dir, paste0(base, ".rds"))
      if (!file.exists(susie_rds_path)) {
        stop(sprintf("SuSiE RDS not found for fallback: %s", susie_rds_path))
      }
      susie_rds <- readRDS(susie_rds_path)

      # Variant names from SuSiE RDS pip vector
      var_names <- names(susie_rds$pip)
      if (is.null(var_names) || length(var_names) != n_pips) {
        stop(sprintf("SuSiE RDS pip names unusable: %d names for %d PIPs",
                      length(var_names), n_pips))
      }

      # Parse chr:pos:a1:a2 into columns
      parts_mat <- do.call(rbind, strsplit(var_names, ":"))
      snp_chr <- as.integer(parts_mat[, 1])
      snp_pos <- as.integer(parts_mat[, 2])
      snp_a1  <- parts_mat[, 3]
      snp_a2  <- parts_mat[, 4]

      # Read the original SS file to get beta/se/pval
      ss_dir <- sub("/CARMA$", "/ss", carma_dir)
      ss_path <- file.path(ss_dir, paste0(study, "_", window, "Mb_", locus, ".txt"))
      if (!file.exists(ss_path)) {
        stop("SS file not found: ", ss_path)
      }
      ss <- fread(ss_path)

      # Build a lookup key: chr:pos:a1:a2 for both orientations
      ss$key_fwd <- paste(ss$chromosome, ss$position, ss$allele1, ss$allele2, sep = ":")
      ss$key_rev <- paste(ss$chromosome, ss$position, ss$allele2, ss$allele1, sep = ":")

      result <- data.table(
        SNP        = var_names,
        chromosome = snp_chr,
        position   = snp_pos,
        allele1    = snp_a1,
        allele2    = snp_a2,
        beta       = NA_real_,
        se         = NA_real_,
        pval       = NA_real_,
        locus      = NA_real_,
        Z          = NA_real_,
        PIP        = rds[[1]]$PIPs,
        CS         = 0L,
        outlier    = 0L
      )

      # Join forward match
      fwd_idx <- match(var_names, ss$key_fwd)
      matched_fwd <- !is.na(fwd_idx)
      result$beta[matched_fwd]  <- ss$beta[fwd_idx[matched_fwd]]
      result$se[matched_fwd]    <- ss$se[fwd_idx[matched_fwd]]
      result$pval[matched_fwd]  <- ss$pval[fwd_idx[matched_fwd]]
      result$locus[matched_fwd] <- ss$locus[fwd_idx[matched_fwd]]

      # Join reverse match (allele flip) for remaining
      remaining <- !matched_fwd
      rev_idx <- match(var_names[remaining], ss$key_rev)
      matched_rev <- !is.na(rev_idx)
      result$beta[which(remaining)[matched_rev]]  <- -ss$beta[rev_idx[matched_rev]]
      result$se[which(remaining)[matched_rev]]    <- ss$se[rev_idx[matched_rev]]
      result$pval[which(remaining)[matched_rev]]  <- ss$pval[rev_idx[matched_rev]]
      result$locus[which(remaining)[matched_rev]] <- ss$locus[rev_idx[matched_rev]]

      # Compute Z from matched beta/se
      result$Z <- result$beta / result$se

      n_matched <- sum(!is.na(result$beta))
      cat("  Matched", n_matched, "of", n_pips, "variants to SS file\n")
      if (n_matched < n_pips * 0.5) {
        warning(sprintf("Low match rate: %d/%d (%.1f%%)", n_matched, n_pips,
                        100 * n_matched / n_pips))
      }
    }

    # Assign credible sets
    cs_obj <- rds[[1]]$`Credible set`
    if (length(cs_obj) >= 2 && length(cs_obj[[2]]) > 0) {
      for (l in seq_along(cs_obj[[2]])) {
        idx <- cs_obj[[2]][[l]]
        # Guard against out-of-range indices
        idx <- idx[idx >= 1 & idx <= nrow(result)]
        result$CS[idx] <- l
      }
    }

    # Assign outliers
    outlier_idx <- rds[[1]]$Outliers$Index
    if (length(outlier_idx) > 0) {
      outlier_idx <- outlier_idx[outlier_idx >= 1 & outlier_idx <= nrow(result)]
      result$outlier[outlier_idx] <- 1L
    }

    # Write via gzfile() to bypass the fwrite compression bug
    con <- gzfile(f, "wt")
    write.table(result, con, sep = "\t", row.names = FALSE, quote = FALSE,
                na = "NA")
    close(con)

    new_size <- file.size(f)
    cat("  Written:", nrow(result), "variants,", new_size, "bytes\n")

    # Verify the file is readable
    check <- fread(f, nrows = 2)
    if (ncol(check) < 10) {
      stop("Verification failed: re-exported file has too few columns")
    }

    n_success <- n_success + 1L

  }, error = function(e) {
    cat("  FAILED:", conditionMessage(e), "\n")
    n_fail <<- n_fail + 1L
    failures <<- c(failures, paste0(basename(f), ": ", conditionMessage(e)))
  })
}

# ---- Step 4: Report ----
cat("\n============================================================\n")
cat("Re-export complete\n")
cat("  Successfully re-exported:", n_success, "files\n")
cat("  Failed:", n_fail, "files\n")
if (length(failures) > 0) {
  cat("\nFailure details:\n")
  for (msg in failures) {
    cat("  -", msg, "\n")
  }
}
cat("============================================================\n")
