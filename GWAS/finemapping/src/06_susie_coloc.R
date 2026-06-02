#!/usr/bin/env Rscript
# 06_susie_coloc.R — SuSiE-COLOC (+ ABF) : GWAS × Broadaway liver eQTLs
# Runs coloc.abf() for all eGenes on one chromosome, AND attempts SuSiE-COLOC using
# pre-computed eQTL SuSiE fits + ancestry-matched GWAS LD (falls back to ABF when LD/SuSiE
# is unavailable). (Header corrected — this script DOES use an LD matrix for the SuSiE arm.)
#
# Usage: Rscript 06_susie_coloc.R <gwas_name> <chr_num>
# Or:    GWAS_NAME=xxx CHR_FILTER=N Rscript 06_susie_coloc.R

args <- commandArgs(trailingOnly = TRUE)
gwas_name <- Sys.getenv("GWAS_NAME", unset = if (length(args) >= 1) args[1] else "")
chr_num <- as.integer(Sys.getenv("CHR_FILTER", unset = if (length(args) >= 2) args[2] else "0"))

if (gwas_name == "" || chr_num == 0) {
  stop("Usage: GWAS_NAME=xxx CHR_FILTER=N Rscript 06_susie_coloc.R\n  Or: Rscript 06_susie_coloc.R <gwas_name> <chr_num>")
}

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
setwd(FM_DIR)

library(data.table)
library(coloc)
library(susieR)
source(file.path(FM_DIR, "src/finemapping_functions.R"))

cat("Versions: coloc", as.character(packageVersion("coloc")),
    "| susieR", as.character(packageVersion("susieR")),
    "| data.table", as.character(packageVersion("data.table")), "\n")

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
COLOC_P1 <- 1e-4
COLOC_P2 <- 1e-4
COLOC_P12 <- 5e-6
EQTL_N <- 1183
MIN_SNPS <- 100L  # Minimum overlapping SNPs for COLOC (raised from 10; see threshold audit)

EQTL_SUSIE_DIR  <- file.path(FM_DIR, "results/eqtl_susie")
MIN_TRIPLE_SNPS <- 50L
LD_REGULARIZE   <- 1e-3
SUSIE_L         <- 10L
SKIP_SUSIE      <- Sys.getenv("SKIP_SUSIE", "0")  # "1" = fast ABF-only pass (skip per-locus LD + SuSiE)

EQTL_DIR <- file.path(BASE_DIR, "data/broadaway_eqtl")
# Output dir suffix lets parallel LD panels write side-by-side without overwriting:
#   COLOC_OUT_SUFFIX="_1kg" → results/susie_coloc_1kg/<gwas>/
#   default ""              → results/susie_coloc/<gwas>/  (back-compat)
OUT_SUFFIX <- Sys.getenv("COLOC_OUT_SUFFIX", unset = "")
OUT_DIR <- file.path(FM_DIR, paste0("results/susie_coloc", OUT_SUFFIX), gwas_name)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(OUT_DIR, paste0("susie_coloc_chr", chr_num, ".csv"))

cat("============================================================\n")
cat("ABF COLOC: GWAS =", gwas_name, " | chr =", chr_num, "\n")
cat("============================================================\n")

# ---------------------------------------------------------------------------
# Load GWAS registry
# ---------------------------------------------------------------------------
registry <- read.delim("config/gwas_registry.tsv", stringsAsFactors = FALSE)
study_row <- registry[registry$study_name == gwas_name, ]
if (nrow(study_row) == 0) stop(paste("GWAS", gwas_name, "not found in registry"))

gwas_path <- study_row$sumstats_path
gwas_n <- study_row$N_tot
gwas_n_cases <- study_row$N_cases
gwas_type <- study_row$trait_type

gwas_ancestry <- study_row$ancestry
if (length(gwas_ancestry) == 0 || is.na(gwas_ancestry) || gwas_ancestry == "") {
  gwas_ancestry <- "EUR"
}
cat("GWAS ancestry:", gwas_ancestry, "\n")
# NOTE: Broadaway eQTLs are European; ABF COLOC is valid across ancestries.
# SuSiE-COLOC uses ancestry-matched LD for GWAS fine-mapping; eQTL SuSiE fits
# always use EUR LD (cross-ancestry design — eQTL cohort is European).
if (gwas_ancestry != "EUR") {
  cat("Non-EUR GWAS detected: SuSiE-COLOC will use", gwas_ancestry,
      "LD for GWAS; ABF COLOC proceeds normally for all genes.\n")
}

gwas_coloc_type <- ifelse(gwas_type == "binary", "cc", "quant")
cat("GWAS:", gwas_path, " N:", gwas_n, " Type:", gwas_coloc_type, "\n")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cat("\n--- Loading GWAS ---\n")
gwas <- fread(gwas_path)
gwas <- gwas[chromosome == chr_num]
cat("  GWAS variants on chr", chr_num, ":", nrow(gwas), "\n")

cat("--- Loading Broadaway eQTLs ---\n")
eqtl_file <- file.path(EQTL_DIR, paste0("chr", chr_num, "_marginal_summary_results.tsv"))
if (!file.exists(eqtl_file)) stop(paste("eQTL file not found:", eqtl_file))
eqtl_all <- fread(eqtl_file)
egenes <- unique(eqtl_all$ENSG)
cat("  eQTL variants:", nrow(eqtl_all), " | eGenes:", length(egenes), "\n")

# Pre-compute GWAS merge keys once
gwas[, merge_key := paste(chromosome, position, sep = ":")]

# ---------------------------------------------------------------------------
# Per-eGene ABF COLOC
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Resume from checkpoint (skip already-processed genes)
# Only resume from CHECKPOINT files (not final output) to avoid treating stale
# ABF-only results as "done". Checkpoints are created by this script during the
# current SuSiE-COLOC run and always have the PP.H4.susie column.
# ---------------------------------------------------------------------------
chk_file <- file.path(OUT_DIR, paste0("checkpoint_chr", chr_num, ".csv"))
done_genes <- character(0)
if (file.exists(chk_file)) {
  chk <- tryCatch(fread(chk_file), error = function(e) NULL)
  if (!is.null(chk) && nrow(chk) > 0 && "ensembl" %in% names(chk) &&
      "PP.H4.susie" %in% names(chk)) {
    done_genes <- unique(chk$ensembl)
    cat("  Resuming from checkpoint:", length(done_genes), "genes already done\n")
  } else if (!is.null(chk)) {
    cat("  Checkpoint exists but lacks SuSiE columns — ignoring (stale ABF-only)\n")
    file.remove(chk_file)
  }
}
# Also check per-worker checkpoint files
worker_chk <- file.path(OUT_DIR, paste0("checkpoint_chr", chr_num, "_w",
                         Sys.getenv("WORKER_ID", unset = "0"), ".csv"))
if (file.exists(worker_chk) && length(done_genes) == 0) {
  wchk <- tryCatch(fread(worker_chk), error = function(e) NULL)
  if (!is.null(wchk) && nrow(wchk) > 0 && "ensembl" %in% names(wchk) &&
      "PP.H4.susie" %in% names(wchk)) {
    done_genes <- unique(wchk$ensembl)
    cat("  Resuming from worker checkpoint:", length(done_genes), "genes already done\n")
  }
}

# ---------------------------------------------------------------------------
# Worker parallelization: WORKER_ID (0-indexed) and N_WORKERS for stride
# Each worker processes genes where (i %% N_WORKERS == WORKER_ID)
# ---------------------------------------------------------------------------
WORKER_ID  <- as.integer(Sys.getenv("WORKER_ID", unset = "0"))
N_WORKERS  <- as.integer(Sys.getenv("N_WORKERS", unset = "1"))
if (N_WORKERS > 1) {
  cat("  Parallel mode: worker", WORKER_ID, "of", N_WORKERS, "(stride pattern)\n")
}

cat("Processing", length(egenes), "eGenes (SuSiE-COLOC with ABF fallback)\n")
n_susie_ok      <- 0L
n_abf_fallback  <- 0L
n_ld_fail       <- 0L

results <- vector("list", length(egenes))
n_tested <- 0L
n_skipped <- 0L
n_resumed <- 0L

for (i in seq_along(egenes)) {
  # Stride: skip genes not assigned to this worker
  if (N_WORKERS > 1 && ((i - 1L) %% N_WORKERS != WORKER_ID)) next

  gene_id <- egenes[i]

  # Resume: skip already-processed genes
  if (gene_id %in% done_genes) { n_resumed <- n_resumed + 1L; next }

  eqtl_gene <- eqtl_all[ENSG == gene_id]
  gene_symbol <- eqtl_gene$GeneSymbol[1]

  if (nrow(eqtl_gene) < MIN_SNPS) { n_skipped <- n_skipped + 1L; next }

  # Merge eQTL + GWAS on position (both hg19)
  eqtl_gene[, merge_key := paste(CHR, POS, sep = ":")]
  merged <- merge(
    eqtl_gene[, .(merge_key, eqtl_pos = POS, eqtl_ea = EA, eqtl_nea = NEA,
                   eqtl_beta = Beta, eqtl_se = SE, eqtl_pval = PVAL)],
    gwas[, .(merge_key, gwas_a1 = allele1, gwas_a2 = allele2,
             gwas_beta = beta, gwas_se = se, gwas_pval = pval)],
    by = "merge_key"
  )
  merged <- merged[order(gwas_pval)][!duplicated(merge_key)]
  if (nrow(merged) < MIN_SNPS) { n_skipped <- n_skipped + 1L; next }

  # Allele harmonization
  merged[, allele_match := (eqtl_ea == gwas_a1 & eqtl_nea == gwas_a2)]
  merged[, allele_flip := (eqtl_ea == gwas_a2 & eqtl_nea == gwas_a1)]
  merged <- merged[allele_match | allele_flip]
  merged[allele_flip == TRUE, eqtl_beta := -eqtl_beta]

  # Remove strand-ambiguous SNPs (A/T and C/G pairs)
  merged <- merged[!((gwas_a1 %in% c("A","T") & gwas_a2 %in% c("A","T")) |
                     (gwas_a1 %in% c("C","G") & gwas_a2 %in% c("C","G")))]

  if (nrow(merged) < MIN_SNPS) { n_skipped <- n_skipped + 1L; next }

  # --- Run ABF COLOC ---
  pp_abf <- rep(NA_real_, 5); names(pp_abf) <- paste0("PP.H", 0:4)
  top_snp <- NA_character_
  top_snp_pp <- NA_real_

  tryCatch({
    d1 <- list(beta = merged$gwas_beta, varbeta = merged$gwas_se^2,
               N = gwas_n, type = gwas_coloc_type, snp = merged$merge_key)
    if (gwas_coloc_type == "cc" && gwas_n_cases > 0) d1$s <- gwas_n_cases / gwas_n
    if (gwas_coloc_type == "quant") {
      d1$sdY <- 1
      # NOTE: sdY=1 assumes betas are on a standardized scale. If betas are on
      # the raw phenotype scale, ABF calibration may be slightly off. Verify
      # that GWAS summary stats were standardized before harmonization.
    }

    d2 <- list(beta = merged$eqtl_beta, varbeta = merged$eqtl_se^2,
               N = EQTL_N, type = "quant", sdY = 1, snp = merged$merge_key)

    abf_res <- suppressMessages(suppressWarnings(
      coloc.abf(d1, d2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12)
    ))
    pp_abf <- abf_res$summary[paste0("PP.H", 0:4, ".abf")]
    names(pp_abf) <- paste0("PP.H", 0:4)

    if (!is.null(abf_res$results)) {
      top_idx <- which.max(abf_res$results$SNP.PP.H4)
      top_snp <- abf_res$results$snp[top_idx]
      top_snp_pp <- abf_res$results$SNP.PP.H4[top_idx]
    }
  }, error = function(e) {})

  # ---- SuSiE-COLOC attempt ----
  # Uses pre-computed eQTL SuSiE fits if available; falls back gracefully.
  # GWAS SuSiE uses ancestry-matched LD (EUR or EAS); eQTL fits always EUR.
  susie_pp4    <- NA_real_
  susie_pp3    <- NA_real_
  susie_method <- "abf_only"
  n_cs_pairs   <- NA_integer_

  eqtl_rds <- file.path(EQTL_SUSIE_DIR, paste0("chr", chr_num),
                         paste0(gene_id, "_susie.rds"))
  if (SKIP_SUSIE != "1" && file.exists(eqtl_rds)) {
    tryCatch({
      s_eqtl <- readRDS(eqtl_rds)

      # SNP IDs in the pre-computed eQTL fit are "CHR:POS"
      eqtl_snp_ids  <- colnames(s_eqtl$lbf_variable)
      eqtl_positions <- as.integer(sub("^[0-9]+:", "", eqtl_snp_ids))

      # Restrict to GWAS variants at those positions (already in merged table)
      gwas_sub <- merged[eqtl_pos %in% eqtl_positions]

      if (nrow(gwas_sub) >= MIN_TRIPLE_SNPS) {
        # Build summary-stats data.frame for get_ld_per_locus
        # The function expects: chromosome, position, allele1, allele2, beta, standard_error
        ss_for_ld <- data.frame(
          chromosome     = chr_num,
          position       = gwas_sub$eqtl_pos,
          allele1        = gwas_sub$gwas_a1,
          allele2        = gwas_sub$gwas_a2,
          beta           = gwas_sub$gwas_beta,
          standard_error = gwas_sub$gwas_se
        )

        ld_result <- get_ld_per_locus(
          ss_per_locus = ss_for_ld,
          LOCUS        = gene_id,
          CHR          = chr_num,
          START        = min(gwas_sub$eqtl_pos),
          END          = max(gwas_sub$eqtl_pos),
          ancestry     = gwas_ancestry
        )

        if (!is.null(ld_result) && nrow(ld_result[[1]]) >= MIN_TRIPLE_SNPS) {
          ss_ld <- ld_result[[1]]   # summary stats aligned to LD variants
          R_mat <- as.matrix(ld_result[[2]])

          R_reg <- R_mat + LD_REGULARIZE * diag(nrow(R_mat))

          gwas_z  <- ss_ld$beta / ss_ld$standard_error
          snp_ids <- paste0(chr_num, ":", ss_ld$position)
          names(gwas_z) <- snp_ids
          colnames(R_reg) <- rownames(R_reg) <- snp_ids

          # Fine-map GWAS with SuSiE (out-of-sample LD: fix residual_variance=1)
          s_gwas <- tryCatch(
            susie_rss(
              z                          = gwas_z,
              R                          = R_reg,
              n                          = gwas_n,
              L                          = SUSIE_L,
              estimate_residual_variance = FALSE,
              residual_variance          = 1,
              check_R                    = FALSE,
              max_iter                   = 500
            ),
            error = function(e) NULL
          )

          if (!is.null(s_gwas)) {
            # Subset eQTL SuSiE to SNPs shared with GWAS LD set
            common_snps <- intersect(snp_ids, colnames(s_eqtl$lbf_variable))

            if (length(common_snps) >= MIN_TRIPLE_SNPS) {
              s_eqtl_sub <- s_eqtl
              s_eqtl_sub$lbf_variable <- s_eqtl$lbf_variable[
                , common_snps, drop = FALSE]

              s_gwas_sub <- s_gwas
              s_gwas_sub$lbf_variable <- s_gwas$lbf_variable[
                , common_snps, drop = FALSE]

              susie_res <- tryCatch(
                coloc.susie(s_gwas_sub, s_eqtl_sub),
                error = function(e) NULL
              )

              if (!is.null(susie_res) && !is.null(susie_res$summary) &&
                  nrow(susie_res$summary) > 0 &&
                  "PP.H4.abf" %in% names(susie_res$summary)) {
                best_row   <- susie_res$summary[which.max(susie_res$summary$PP.H4.abf), ]
                susie_pp4  <- best_row$PP.H4.abf
                susie_pp3  <- best_row$PP.H3.abf
                n_cs_pairs <- nrow(susie_res$summary)
                susie_method <- "susie"
              }
            }
          }
        } else {
          n_ld_fail <- n_ld_fail + 1L
        }
      }
    }, error = function(e) {
      cat("  SuSiE-COLOC failed for", gene_id, ":", conditionMessage(e), "\n")
    })

    # Free LD/SuSiE objects from this iteration to keep memory bounded
    rm(list = intersect(ls(), c("ld_result", "ss_ld", "R_mat", "R_reg",
       "s_gwas", "s_eqtl", "s_eqtl_sub", "s_gwas_sub", "susie_res",
       "ss_for_ld", "gwas_sub")))
    if (i %% 100 == 0) gc(verbose = FALSE)

    if (susie_method == "abf_only") {
      susie_method <- "abf_fallback"
      n_abf_fallback <- n_abf_fallback + 1L
    } else {
      n_susie_ok <- n_susie_ok + 1L
    }
  }
  # susie_method stays "abf_only" when no eQTL .rds exists (no attempt made)

  n_tested <- n_tested + 1L
  results[[i]] <- data.table(
    gene        = gene_symbol,
    ensembl     = gene_id,
    chr         = chr_num,
    gwas_name   = gwas_name,
    PP.H0.abf   = pp_abf[1],
    PP.H1.abf   = pp_abf[2],
    PP.H2.abf   = pp_abf[3],
    PP.H3.abf   = pp_abf[4],
    PP.H4.abf   = pp_abf[5],
    PP.H3.susie = susie_pp3,
    PP.H4.susie = susie_pp4,
    n_cs_pairs  = n_cs_pairs,
    n_snps      = nrow(merged),
    method      = susie_method,
    top_snp     = top_snp,
    top_snp_PP  = top_snp_pp
  )

  # Progress + checkpoint every 200 genes
  if (n_tested %% 200 == 0 && n_tested > 0) {
    new_results <- results[!sapply(results, is.null)]
    if (length(new_results) > 0) {
      # Merge new results with any resumed checkpoint
      chk_out <- if (N_WORKERS > 1) {
        file.path(OUT_DIR, paste0("checkpoint_chr", chr_num, "_w", WORKER_ID, ".csv"))
      } else {
        file.path(OUT_DIR, paste0("checkpoint_chr", chr_num, ".csv"))
      }
      new_dt <- rbindlist(new_results, fill = TRUE)
      # Merge with resumed checkpoint (shared or per-worker)
      chk_resume <- if (length(done_genes) > 0 && file.exists(chk_file)) {
        chk_file
      } else if (length(done_genes) > 0 && exists("worker_chk") && file.exists(worker_chk)) {
        worker_chk
      } else {
        NULL
      }
      if (!is.null(chk_resume)) {
        old_dt <- tryCatch(fread(chk_resume), error = function(e) NULL)
        if (!is.null(old_dt)) new_dt <- rbindlist(list(old_dt, new_dt), fill = TRUE)
        new_dt <- new_dt[!duplicated(ensembl)]
      }
      fwrite(new_dt, chk_out)
    }
    cat(sprintf("  Progress: %d/%d (tested: %d, resumed: %d, skipped: %d | susie_ok: %d, abf_fallback: %d, ld_fail: %d)\n",
                i, length(egenes), n_tested, n_resumed, n_skipped,
                n_susie_ok, n_abf_fallback, n_ld_fail))
  }
}

# ---------------------------------------------------------------------------
# Save final results
# ---------------------------------------------------------------------------
results <- results[!sapply(results, is.null)]

if (length(results) > 0) {
  final <- rbindlist(results, fill = TRUE)

  # Merge with resumed checkpoint results (shared OR per-worker checkpoint)
  resume_file <- if (length(done_genes) > 0 && file.exists(chk_file)) {
    chk_file
  } else if (length(done_genes) > 0 && exists("worker_chk") && file.exists(worker_chk)) {
    worker_chk
  } else {
    NULL
  }
  if (!is.null(resume_file)) {
    old_dt <- tryCatch(fread(resume_file), error = function(e) NULL)
    if (!is.null(old_dt)) final <- rbindlist(list(old_dt, final), fill = TRUE)
    final <- final[!duplicated(ensembl)]
  }

  # Workers write to separate files; final merge done by worker 0 or wrapper
  if (N_WORKERS > 1) {
    worker_file <- file.path(OUT_DIR, paste0("worker_chr", chr_num, "_w", WORKER_ID, ".csv"))
    fwrite(final, worker_file)
    cat("Worker", WORKER_ID, "wrote", nrow(final), "genes to", worker_file, "\n")
  } else {
    fwrite(final, out_file)
  }

  cat("\n============================================================\n")
  cat("Results:", out_file, "\n")
  cat("Genes tested:", nrow(final), " | Skipped:", n_skipped, "\n")
  cat("Method breakdown:",
      "susie =", n_susie_ok, "|",
      "abf_fallback =", n_abf_fallback, "|",
      "ld_fail (within susie attempt) =", n_ld_fail, "|",
      "abf_only (no eQTL .rds) =",
      sum(final$method == "abf_only", na.rm = TRUE), "\n")
  cat("--- ABF COLOC ---\n")
  cat("PP.H4.abf > 0.8:", sum(final$PP.H4.abf > 0.8, na.rm = TRUE), "\n")
  cat("PP.H4.abf > 0.5:", sum(final$PP.H4.abf > 0.5, na.rm = TRUE), "\n")
  cat("PP.H4.abf > 0.3:", sum(final$PP.H4.abf > 0.3, na.rm = TRUE), "\n")
  cat("--- SuSiE COLOC (method == 'susie' only) ---\n")
  cat("PP.H4.susie > 0.8:", sum(final$PP.H4.susie > 0.8, na.rm = TRUE), "\n")
  cat("PP.H4.susie > 0.5:", sum(final$PP.H4.susie > 0.5, na.rm = TRUE), "\n")
  cat("PP.H4.susie > 0.3:", sum(final$PP.H4.susie > 0.3, na.rm = TRUE), "\n")
  cat("============================================================\n")
} else {
  cat("WARNING: No results produced for chr", chr_num, "\n")
}

# Clean up checkpoint
chk <- file.path(OUT_DIR, paste0("checkpoint_chr", chr_num, ".csv"))
if (file.exists(chk)) file.remove(chk)
