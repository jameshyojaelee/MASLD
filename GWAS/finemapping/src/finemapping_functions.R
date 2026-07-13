#!/usr/bin/env Rscript
# finemapping_functions.R
# Adapted from Lappalainen-Singh lab pipeline (finemapping_autoimmune)
# Changes: configurable paths (no hardcoded lab paths), removed GCS functions

library(data.table)
library(dplyr)
library(Matrix)
library(readr)
library(tidyr)
library(stringr)

# ---------------------------------------------------------------------------
# Configuration — set via environment variable or default
# ---------------------------------------------------------------------------
# LD panel dispatch:
#   LD_PANEL env var ∈ {"1kg", "topld", "ukbb", "polyfun" (default)}
#     1kg     → 1000 Genomes Phase 3 per ancestry (1kg_eur, 1kg_eas, 1kg_afr, 1kg_sas)
#     topld   → TOP-LD (Huang 2022 AJHG) per ancestry (topld_eur, topld_eas, topld_afr, topld_sas)
#     ukbb    → sghatan UKBB EUR (EUR only; AFR/EAS/SAS fall back to 1kg). DO NOT USE
#               for new computation — symlink is preserved only for v1 snapshot reads.
#     polyfun → PolyFun precomputed UKBB EUR LD (Weissbrod 2020; ~337K UKBB EUR;
#               polyfun_eur). EUR only; AFR/EAS/SAS fall back to 1kg with a
#               one-time warning.
# Per-ancestry env vars (UKBB_LD_DIR, EAS_LD_DIR, AFR_LD_DIR, SAS_LD_DIR) override
# the dispatch entirely if set — use only for one-off custom panel paths.
get_ld_base_dir <- function(ancestry = "EUR") {
  # Default = polyfun (review/IND-2): matches the documented production EUR LD reference
  # (Phase 10, 2026-05-06). Non-EUR ancestries fall back to 1kg below, so this is safe for
  # all callers and prevents an un-set run from silently using the 1kg sensitivity panel.
  panel <- tolower(Sys.getenv("LD_PANEL", unset = "polyfun"))
  if (!panel %in% c("1kg", "topld", "ukbb", "polyfun")) {
    stop(paste0("LD_PANEL='", panel, "' not recognized. Use one of: 1kg, topld, ukbb, polyfun."))
  }
  proj_root <- Sys.getenv("MASLD_PROJECT_ROOT",
    unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
  ld_ref <- file.path(proj_root, "GWAS/finemapping/data/ld_ref")

  # Per-ancestry env override beats panel dispatch
  override_env <- switch(ancestry,
    "EUR" = "UKBB_LD_DIR",
    "EAS" = "EAS_LD_DIR",
    "AFR" = "AFR_LD_DIR",
    "AMR" = "AMR_LD_DIR",
    "SAS" = "SAS_LD_DIR",
    stop(paste("Unsupported ancestry:", ancestry)))
  override_path <- Sys.getenv(override_env, unset = "")
  if (nzchar(override_path)) return(override_path)

  # LD_PANEL=ukbb deprecation (Phase 10, 2026-05-06): PolyFun replaced sghatan
  # UKBB v1 as the production EUR LD reference (concordance r=0.986). The v1
  # snapshot has been moved to archive/ld_panel_v1_sghatan_2026-05-06/ and is
  # read-only / sensitivity-comparison only. Any new computation that requests
  # LD_PANEL=ukbb is almost certainly a mistake — use polyfun (default-equivalent
  # for EUR) or 1kg.
  if (panel == "ukbb" && ancestry == "EUR") {
    if (!isTRUE(getOption(".ukbb_eur_deprecated_warned", FALSE))) {
      warning(paste0(
        "LD_PANEL='ukbb' is DEPRECATED for new computation. The sghatan UKBB v1 ",
        "EUR LD panel was retired 2026-05-06; PolyFun replaced it (concordance ",
        "r=0.986). The v1 snapshot is preserved at ",
        "archive/ld_panel_v1_sghatan_2026-05-06/ for sensitivity comparisons ",
        "only. Use LD_PANEL=polyfun (production) or LD_PANEL=1kg (fallback)."),
        call. = FALSE)
      options(.ukbb_eur_deprecated_warned = TRUE)
    }
  }

  # Panel × ancestry dispatch
  if (panel == "polyfun" && ancestry != "EUR") {
    if (!isTRUE(getOption(".polyfun_nonEUR_warned", FALSE))) {
      message(sprintf(
        "[get_ld_base_dir] LD_PANEL=polyfun has no %s LD; falling back to 1kg_%s.",
        ancestry, tolower(ancestry)))
      options(.polyfun_nonEUR_warned = TRUE)
    }
    return(file.path(ld_ref, paste0("1kg_", tolower(ancestry))))
  }
  subdir <- switch(panel,
    "1kg"     = paste0("1kg_",   tolower(ancestry)),
    "topld"   = paste0("topld_", tolower(ancestry)),
    "ukbb"    = if (ancestry == "EUR") "ukbb_eur" else paste0("1kg_", tolower(ancestry)),
    "polyfun" = "polyfun_eur")  # only reaches here when ancestry == "EUR"
  return(file.path(ld_ref, subdir))
}

get_fm_base_dir <- function() {
  file.path(Sys.getenv("MASLD_PROJECT_ROOT",
    unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
    "GWAS/finemapping")
}

# ---------------------------------------------------------------------------
# Stage 1: Extract per-locus summary stats
# ---------------------------------------------------------------------------

get_ss_per_locus <- function(df.sumstats, LOCUS, CHR, START, END) {
  df_chrom <- df.sumstats %>%
    filter(chromosome == CHR) %>%
    arrange(position)
  ROW_START <- which(df_chrom[, "position"] > START)[1]
  ROW_END <- which(df_chrom[, "position"] < END)[length(which(df_chrom[, "position"] < END))]
  if (is.na(ROW_START) || is.na(ROW_END) || ROW_START > ROW_END) {
    warning(paste("No variants found for locus", LOCUS, "chr", CHR, START, "-", END))
    return(NULL)
  }
  df_chrom[ROW_START:ROW_END, "locus"] <- as.character(LOCUS)
  return(df_chrom[ROW_START:ROW_END, ])
}

save_ss_per_locus <- function(sumstats_name, df.sumstats, LOCUS, CHR, START, END,
                              window_mb, ld_pop, base_dir = get_fm_base_dir()) {
  ss <- get_ss_per_locus(df.sumstats, as.character(LOCUS), CHR, START, END)
  if (is.null(ss) || nrow(ss) == 0) {
    warning(paste("Skipping locus", LOCUS, "- no variants in region"))
    return(invisible(NULL))
  }
  ss_dir <- file.path(base_dir, "output", sumstats_name,
                      paste0(ld_pop, "_", window_mb, "Mb"), "ss")
  dir.create(ss_dir, recursive = TRUE, showWarnings = FALSE)
  ss_path <- file.path(ss_dir, paste0(sumstats_name, "_", window_mb, "Mb_",
                                      as.character(LOCUS), ".txt"))
  write_delim(ss, ss_path, delim = "\t")
  cat("  Saved", nrow(ss), "variants for locus", LOCUS, "\n")
}

# ---------------------------------------------------------------------------
# Stage 2: LD matrix extraction (local UKBB or 1000G)
# ---------------------------------------------------------------------------

find_LD_block <- function(CHR, START, END, ancestry = "EUR") {
  ld_base <- get_ld_base_dir(ancestry)
  blocks_file <- file.path(ld_base, "approx_LD_blocks.txt")
  if (!file.exists(blocks_file)) {
    stop(paste("LD blocks file not found:", blocks_file))
  }
  approx_blocks <- fread(blocks_file)
  matching <- approx_blocks %>% filter(chr == CHR, start < END, stop > START)
  num_blocks <- nrow(matching)
  cat(paste0("  Number of LD blocks: ", num_blocks, "\n"))
  if (num_blocks == 0) {
    warning(paste("No LD blocks found for chr", CHR, START, "-", END))
    return(NULL)
  }
  LD_start <- matching %>% pull(start)
  LD_stop <- matching %>% pull(stop)
  filepath <- file.path(ld_base, paste0("chr", CHR),
                        paste0(LD_start, ".", LD_stop),
                        paste0(LD_start, ".", LD_stop))
  return(filepath)
}

get_ld_per_locus <- function(ss_per_locus, LOCUS, CHR, START, END, ancestry = "EUR") {
  path_to_LD <- find_LD_block(CHR, START, END, ancestry)
  if (is.null(path_to_LD)) {
    return(NULL)
  }

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

# ---------------------------------------------------------------------------
# Stage 3a: Run SuSiE
# ---------------------------------------------------------------------------

run_susie <- function(sumstats, LDmat, N_tot, N_cases, sumstats_name, ld_pop,
                      window_mb, locus, LDpanel, coverage = 0.95,
                      base_dir = get_fm_base_dir()) {
  if (nrow(sumstats) != nrow(LDmat)) {
    warning(paste("SS dimensions", nrow(sumstats), "do not match LD", nrow(LDmat)))
  }
  colnames(LDmat) <- rownames(LDmat) <- paste(sumstats$chromosome, sumstats$position,
                                               sumstats$allele1, sumstats$allele2, sep = ":")

  out_dir <- file.path(base_dir, "output", sumstats_name,
                       paste0(ld_pop, "_", window_mb, "Mb"), "susie")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    # LD ridge regularization to ensure positive-definiteness
    R_mat <- as.matrix(LDmat) + 1e-6 * diag(nrow(LDmat))

    # C6 (2026-07-05): effective sample size for case-control (log-OR) strata.
    #   N_eff = 4 / (1/N_cases + 1/N_ctrl)  (Willer 2010; Kanai 2022 SuSiE-RSS).
    # Quantitative strata (N_cases is NA) -> N_eff == N_tot. Matches the GWAS-side
    # SuSiE-COLOC arm in 06_susie_coloc.R (gwas_n_eff). Keeping N_tot over-states the
    # information under case:control imbalance (e.g. deCODE NAFLD N_eff/N_tot=0.009).
    N_eff <- if (!is.na(N_cases) && N_cases > 0 && (N_tot - N_cases) > 0) {
      4 / (1/N_cases + 1/(N_tot - N_cases))
    } else {
      N_tot
    }

    # LD-consistency diagnostic (Zou 2022 SuSiE-RSS; Kanai 2022): estimate the
    # inconsistency parameter s between the z-scores and the (out-of-sample /
    # small-panel) LD reference. High s flags loci where the reference LD does not
    # match the GWAS -> credible sets are unreliable (the false-credible-set risk
    # for small non-EUR 1000G panels: AMR n~347, EAS 504, AFR ~660). Recorded on
    # every output row so aggregation (04) can gate/label exploratory loci; it
    # never blocks the run.
    lambda_s <- tryCatch(
      as.numeric(susieR::estimate_s_rss(z = sumstats$beta / sumstats$se,
                                        R = R_mat, n = N_eff)),
      error = function(e) NA_real_
    )
    cat("  LD-consistency lambda_s:",
        ifelse(is.na(lambda_s), "NA", round(lambda_s, 4)), "\n")
    sumstats$lambda_s <- lambda_s

    if (is.na(N_cases)) {
      # quantitative: N_eff == N_tot; bhat/shat mode unchanged.
      set.seed(42)  # C7 reproducibility: deterministic FM-master susie_rss (quant branch)
      susie_output <- susieR::susie_rss(R = R_mat, n = N_eff,
                                        bhat = sumstats$beta, shat = sumstats$se,
                                        coverage = coverage)
    } else {
      # binary (log-OR): standardized z-score representation scaled by N_eff (C6,
      # 2026-07-05). Supersedes the linear-on-0/1 var_y = phi*(1-phi) approximation,
      # which is mis-specified for logistic (SAIGE/REGENIE log-OR) effect sizes and
      # now matches the GWAS-side SuSiE-COLOC arm.
      set.seed(42)  # C7 reproducibility: deterministic FM-master susie_rss (case-control branch)
      susie_output <- susieR::susie_rss(R = R_mat, n = N_eff,
                                        z = sumstats$beta / sumstats$se,
                                        coverage = coverage)
    }

    saveRDS(susie_output, file = file.path(out_dir, paste0(sumstats_name, "_",
            LDpanel, "_", window_mb, "Mb_", locus, ".rds")))

    if (susie_output$converged) {
      cat("  SuSiE converged for locus", locus, "\n")
      sumstats$converged <- TRUE
      sumstats$PIP <- susie_output$pip
      cs_out <- susie_output$sets
      sumstats$CS <- 0L
      if (!is.null(cs_out$cs)) {
        for (j in seq_along(cs_out$cs)) {
          sumstats$CS[cs_out$cs[[j]]] <- j
        }
        cat("    Found", length(cs_out$cs), "credible sets\n")
      }
      write_tsv(sumstats, file.path(out_dir, paste0(LDpanel, "_", locus,
                "_cov", coverage, "_adjustedvar.tsv")))
    } else {
      cat("  SuSiE did NOT converge for locus", locus, "\n")
      sumstats$converged <- FALSE
      sumstats$PIP <- susie_output$pip
      sumstats$CS <- 0L
      write_tsv(sumstats, file.path(out_dir, paste0(LDpanel, "_", locus,
                "_cov", coverage, "_notconverged.tsv")))
    }
    return(invisible(susie_output))
  }, error = function(e) {
    cat("  SuSiE ERROR for locus", locus, ":", conditionMessage(e), "\n")
    return(invisible(NULL))
  })
}

# ---------------------------------------------------------------------------
# Stage 3b: Run CARMA
# ---------------------------------------------------------------------------

run_CARMA <- function(sumstat, ld, sumstats_name, ld_pop, window_mb, locus,
                      LDpanel, base_dir = get_fm_base_dir()) {
  out_dir <- file.path(base_dir, "output", sumstats_name,
                       paste0(ld_pop, "_", window_mb, "Mb"), "CARMA")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    z.list <- list()
    ld.list <- list()
    lambda.list <- list()
    sumstat$Z <- sumstat$beta / sumstat$se
    z.list[[1]] <- sumstat$Z
    ld.list[[1]] <- as.matrix(ld)
    lambda.list[[1]] <- 1

    set.seed(42)  # C7 reproducibility: deterministic FM-master CARMA fine-mapping
    CARMA.results <- CARMA::CARMA(z.list, ld.list, lambda.list = lambda.list,
                                   outlier.switch = TRUE, rho.index = 0.95)

    saveRDS(CARMA.results, file = file.path(out_dir, paste0(sumstats_name, "_",
            LDpanel, "_", window_mb, "Mb_", locus, ".rds")))

    cat("  CARMA completed for locus", locus, "\n")
    if (length(CARMA.results[[1]]$Outliers$Index) > 0) {
      cat("    Outliers:", length(CARMA.results[[1]]$Outliers$Index), "variants\n")
    }

    sumstat.result <- sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0L)

    if (length(CARMA.results[[1]]$`Credible set`[[2]]) != 0) {
      for (l in seq_along(CARMA.results[[1]]$`Credible set`[[2]])) {
        sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]] <- l
      }
      cat("    Found", length(CARMA.results[[1]]$`Credible set`[[2]]), "credible sets\n")
    }

    sumstat.result$outlier <- 0L
    sumstat.result$outlier[CARMA.results[[1]]$Outliers$Index] <- 1L

    fwrite(x = sumstat.result,
           file = file.path(out_dir, paste0(sumstats_name, "_", LDpanel, "_",
                  window_mb, "Mb_", locus, ".txt.gz")),
           sep = "\t", quote = FALSE, na = "NA", row.names = FALSE,
           col.names = TRUE, compress = "gzip")

    return(invisible(CARMA.results))
  }, error = function(e) {
    cat("  CARMA ERROR for locus", locus, ":", conditionMessage(e), "\n")
    return(invisible(NULL))
  })
}
