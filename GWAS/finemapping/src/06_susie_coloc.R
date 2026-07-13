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
# NOTE: OUT_DIR is created AFTER the registry validation (below), so a stratum not
# in the registry never leaves an empty results/susie_coloc/<gwas>/ dir that 07's
# list.dirs() glob would later ingest into gene_level_coloc.csv (hardening 2026-07-04).

cat("============================================================\n")
cat("ABF COLOC: GWAS =", gwas_name, " | chr =", chr_num, "\n")
cat("============================================================\n")

# ---------------------------------------------------------------------------
# Load GWAS registry
# ---------------------------------------------------------------------------
registry <- read.delim("config/gwas_registry.tsv", stringsAsFactors = FALSE)
study_row <- registry[registry$study_name == gwas_name, ]
if (nrow(study_row) == 0) stop(paste("GWAS", gwas_name, "not found in registry"))

# registry validated -> now safe to create the output dir + resolve the output path
OUT_DIR <- file.path(FM_DIR, paste0("results/susie_coloc", OUT_SUFFIX), gwas_name)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(OUT_DIR, paste0("susie_coloc_chr", chr_num, ".csv"))

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
# C3 (2026-07-05): registry -> coloc statistical-spec assertion guard.
# Validate that the coloc spec derived from this registry row is internally
# consistent BEFORE any coloc.abf()/coloc.susie() call, so a future registry
# edit can never silently (a) run a binary trait as quant / a quant as cc via an
# unmapped trait_type, (b) set an out-of-range case fraction s = N_cases/N_tot
# (needs 0 < s < 1), or (c) attach a case fraction to a quantitative trait.
# This asserts EXACTLY the type/s/sdY mapping the coloc datasets use below
# (d1$type = gwas_coloc_type; d1$s = N_cases/N_tot only when cc & N_cases>0,
#  line ~262; d1$sdY = 1 for quant; d2 eQTL always type="quant", sdY=1) and the
# per-row spec documented in results/qa_campaign/coloc_spec_by_registry_row.tsv.
# Inert on the current spec-clean registry (all 50 rows pass); a STANDING guard
# for future runs. Does NOT touch the coloc engine, the ABF PP4, or on-disk
# results (mirrors the C6 N_eff / C8 anchor standing-guard pattern).
if (!gwas_type %in% c("binary", "quantitative")) {
  stop(sprintf(paste0("C3 spec guard: unmapped trait_type '%s' for %s ",
                      "(must be binary|quantitative -> cc|quant)"),
               gwas_type, gwas_name))
}
if (gwas_coloc_type == "cc") {
  if (is.na(gwas_n_cases) || gwas_n_cases <= 0) {
    stop(sprintf(paste0("C3 spec guard: binary trait %s has N_cases=%s (<=0) ",
                        "-> coloc case fraction 's' would be unset/invalid"),
                 gwas_name, ifelse(is.na(gwas_n_cases), "NA", gwas_n_cases)))
  }
  if (gwas_n_cases >= gwas_n) {
    stop(sprintf(paste0("C3 spec guard: binary trait %s has N_cases(%s) >= ",
                        "N_tot(%s) -> s = N_cases/N_tot >= 1 (invalid)"),
                 gwas_name, gwas_n_cases, gwas_n))
  }
  cat(sprintf("  C3 spec guard OK: cc, s = N_cases/N_tot = %d/%d = %.4f (0<s<1)\n",
              gwas_n_cases, gwas_n, gwas_n_cases / gwas_n))
} else {
  # quantitative -> quant, sdY = 1, no case fraction (d1$s must NOT be set)
  if (!is.na(gwas_n_cases) && gwas_n_cases > 0) {
    stop(sprintf(paste0("C3 spec guard: quantitative trait %s has N_cases=%s ",
                        "(>0) -> must be 0/NA for a quant (sdY=1) spec"),
                 gwas_name, gwas_n_cases))
  }
  cat(sprintf("  C3 spec guard OK: quant, sdY=1, N_cases=%s (no case fraction)\n",
              ifelse(is.na(gwas_n_cases), "NA", gwas_n_cases)))
}

# ---------------------------------------------------------------------------
# C6 (2026-07-05): effective sample size for the GWAS-side SuSiE fine-map arm.
# For a binary (case-control) trait whose summary stats are log-odds (logistic
# SAIGE/REGENIE -- every NAFLD/NASH/cirrhosis/chronic-liver stratum here), the
# z-score-mode susie_rss()/estimate_s_rss() "n" must be the EFFECTIVE sample size
#   N_eff = 4 / (1/N_cases + 1/N_ctrl)     (Willer 2010; Yang 2012; Kanai 2022
# SuSiE-RSS for case-control), NOT N_tot. Passing N_tot over-states the
# information under case:control imbalance (e.g. NAFLD_2021 N_tot=778,614 but
# N_eff=33,371) and can inflate GWAS-side credible-set precision. Quantitative
# strata carry N_cases=0 -> N_eff falls through to N_tot (N_eff==N by definition).
#
# SCOPE / BLAST RADIUS (see results/qa_campaign/audit/C6_neff_audit.R):
#   * affects ONLY the SuSiE-COLOC arm  = 22,072 of 927,225 coloc rows (~2.4%);
#     the ABF arm = 905,153 rows (835,143 abf_fallback + 70,010 abf_only).
#   * and only the 19 binary strata (31 quantitative strata are N-invariant).
#   * the PRIMARY/headline arm is coloc.abf() below, which is a Wakefield ABF on
#     beta+varbeta -- N enters coloc.abf ONLY when varbeta is absent, and varbeta
#     is ALWAYS supplied here (merged$gwas_se^2), so the ABF arm is N-INVARIANT
#     and is deliberately left unchanged. Current on-disk results are unaffected
#     until 06 is re-run by the FM pipeline; this is a standing code guard.
gwas_n_ctrl <- if (!is.na(gwas_n_cases)) gwas_n - gwas_n_cases else NA_real_
gwas_n_eff <- if (gwas_coloc_type == "cc" && !is.na(gwas_n_cases) &&
                  gwas_n_cases > 0 && !is.na(gwas_n_ctrl) && gwas_n_ctrl > 0) {
  4 / (1/gwas_n_cases + 1/gwas_n_ctrl)   # N_eff (Willer 2010) for case-control
} else {
  gwas_n                                 # quantitative: N_eff == N_tot
}
cat(sprintf("  GWAS-side SuSiE effective N: N_tot=%s N_cases=%s -> N_eff=%s [%s]\n",
            gwas_n, ifelse(is.na(gwas_n_cases), "NA", gwas_n_cases),
            round(gwas_n_eff),
            ifelse(gwas_coloc_type == "cc", "binary:N_eff", "quant:N_tot")))

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
  # C1 harmonization tally (per-locus): tally match / flip / unresolved BEFORE the
  # match|flip filter, and strand-ambiguous palindromes dropped AFTER the flip. These
  # per-locus counts are emitted in the results table below so every future coloc run
  # carries an auditable allele-harmonization trail (n_match/n_flip/n_unresolved/n_ambiguous).
  # The flip/palindrome LOGIC below is unchanged — only the counters are new.
  n_merged_pre <- nrow(merged)
  merged[, allele_match := (eqtl_ea == gwas_a1 & eqtl_nea == gwas_a2)]
  merged[, allele_flip := (eqtl_ea == gwas_a2 & eqtl_nea == gwas_a1)]
  n_match_locus      <- sum(merged$allele_match, na.rm = TRUE)
  n_flip_locus       <- sum(merged$allele_flip,  na.rm = TRUE)
  n_unresolved_locus <- n_merged_pre - n_match_locus - n_flip_locus  # alleles reconcile in neither orientation -> dropped
  merged <- merged[allele_match | allele_flip]
  merged[allele_flip == TRUE, eqtl_beta := -eqtl_beta]

  # Remove strand-ambiguous SNPs (A/T and C/G pairs)
  n_pre_ambiguous <- nrow(merged)
  merged <- merged[!((gwas_a1 %in% c("A","T") & gwas_a2 %in% c("A","T")) |
                     (gwas_a1 %in% c("C","G") & gwas_a2 %in% c("C","G")))]
  n_ambiguous_locus <- n_pre_ambiguous - nrow(merged)  # strand-ambiguous palindromes dropped after the flip

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
  susie_lambda_s <- NA_real_   # C4: per-gene GWAS-side LD-consistency (estimate_s_rss)

  eqtl_rds <- file.path(EQTL_SUSIE_DIR, paste0("chr", chr_num),
                         paste0(gene_id, "_susie.rds"))
  if (SKIP_SUSIE != "1" && file.exists(eqtl_rds)) {
    tryCatch({
      s_eqtl <- readRDS(eqtl_rds)

      # Guard: a non-converged eQTL SuSiE fit can report PIP=1.0 for multiple
      # variants at once (physically impossible; PIPs within a credible set
      # should sum to ~1 for a single causal signal). Feeding such a fit into
      # coloc.susie() emits a spurious high PP.H4 (the PNPLA3/TM6SF2 PIP-myth
      # producer). Skip the SuSiE arm and fall through to ABF when the eQTL fit
      # did not converge (or lacks the converged flag).
      eqtl_converged <- isTRUE(s_eqtl$converged)
      if (!eqtl_converged) {
        stop("eQTL SuSiE fit not converged — skipping SuSiE-COLOC (ABF fallback)")
      }

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

          # C4: LD-consistency diagnostic for THIS gene's cis-region GWAS fine-map
          # (susieR::estimate_s_rss; Zou 2022 SuSiE-RSS). High lambda_s => the
          # (small non-EUR 1000G) LD reference does not match the GWAS sample, so
          # the coloc.susie() PP.H4 for this gene is LD-unreliable. Emitted per gene
          # (column lambda_s_locus) so 07 carries a genuine per-locus lambda_s
          # without joining the FM master. Cheap relative to the SuSiE fit below.
          susie_lambda_s <- tryCatch(
            # C6: N_eff (binary) / N_tot (quant) -- see the gwas_n_eff block above.
            as.numeric(susieR::estimate_s_rss(z = gwas_z, R = R_reg, n = gwas_n_eff)),
            error = function(e) NA_real_
          )

          # Fine-map GWAS with SuSiE (out-of-sample LD: fix residual_variance=1)
          # C6: n = gwas_n_eff -> N_eff for binary (case-control) strata, N_tot for
          # quantitative (see the gwas_n_eff block above). SuSiE-arm-only; the ABF
          # primary arm is N-invariant and unchanged.
          # PURITY: susie_rss() uses the susieR default min_abs_corr = 0.5, so any
          # credible set reported into coloc.susie() already meets the standard 0.5
          # within-CS min-|r| purity floor (Wang 2020) -- low-purity (LD-diffuse)
          # sets are pruned by susieR before they can drive a spurious PP.H4.
          set.seed(42)  # C7 reproducibility (2026-07-05): deterministic GWAS-side susie_rss
          s_gwas <- tryCatch(
            susie_rss(
              z                          = gwas_z,
              R                          = R_reg,
              n                          = gwas_n_eff,
              L                          = SUSIE_L,
              estimate_residual_variance = FALSE,
              residual_variance          = 1,
              check_R                    = FALSE,
              min_abs_corr               = 0.5,   # C6: explicit CS purity floor (susieR default)
              max_iter                   = 500
            ),
            error = function(e) NULL
          )

          # Same convergence guard on the GWAS-side SuSiE fit: a non-converged
          # susie_rss() solution yields unreliable PIPs that propagate into a
          # spurious coloc.susie() PP.H4. Require convergence before colocalizing.
          if (!is.null(s_gwas) && isTRUE(s_gwas$converged)) {
            # Subset eQTL SuSiE to SNPs shared with GWAS LD set
            common_snps <- intersect(snp_ids, colnames(s_eqtl$lbf_variable))

            if (length(common_snps) >= MIN_TRIPLE_SNPS) {
              s_eqtl_sub <- s_eqtl
              s_eqtl_sub$lbf_variable <- s_eqtl$lbf_variable[
                , common_snps, drop = FALSE]

              s_gwas_sub <- s_gwas
              s_gwas_sub$lbf_variable <- s_gwas$lbf_variable[
                , common_snps, drop = FALSE]

              set.seed(42)  # C7 reproducibility (2026-07-05): deterministic coloc.susie
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
    n_match     = n_match_locus,       # C1 harmonization tally: eqtl EA/NEA already aligned to GWAS a1/a2
    n_flip      = n_flip_locus,        # C1 harmonization tally: eqtl EA/NEA reversed -> eqtl_beta negated
    n_unresolved = n_unresolved_locus, # C1 harmonization tally: alleles matched neither orientation (dropped)
    n_ambiguous = n_ambiguous_locus,   # C1 harmonization tally: strand-ambiguous A/T or C/G palindromes (dropped)
    method      = susie_method,
    lambda_s_locus = susie_lambda_s,  # C4: GWAS-side LD-consistency for this gene's cis-region
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
  # C1 harmonization tally (this run's newly-tested loci; resumed rows lack these cols -> na.rm)
  if (all(c("n_match","n_flip","n_unresolved","n_ambiguous") %in% names(final))) {
    cat("--- Allele harmonization tally (SNP-instances across tested loci) ---\n")
    cat("  matched (a1/a2 aligned) :", sum(final$n_match, na.rm = TRUE), "\n")
    cat("  flipped (eqtl_beta neg.):", sum(final$n_flip, na.rm = TRUE), "\n")
    cat("  unresolved (dropped)    :", sum(final$n_unresolved, na.rm = TRUE), "\n")
    cat("  strand-ambiguous dropped:", sum(final$n_ambiguous, na.rm = TRUE), "\n")
  }
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
