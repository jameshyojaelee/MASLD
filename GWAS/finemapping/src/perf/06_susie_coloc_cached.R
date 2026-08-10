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
source(file.path(FM_DIR, "src/perf/get_ld_per_locus_cached.R"))  # PERF OVERRIDE (verification)

# --- 1000G panel allele frequencies, for the palindrome strand check ---------
# Built by src/seqfunc_v2/finemap/00a_build_panel_af.R (110 sidecars,
# 49,978,712 variants).  Schema: chromosome position bim_a1 bim_a2 af_a1
# n_obs_alleles panel_pop, where af_a1 is the frequency of bim_a1.
PANEL_AF_ROOT <- Sys.getenv("PANEL_AF_ROOT",
                            unset = file.path(FM_DIR, "data/ld_ref/panel_af"))
load_panel_af <- function(pop, chr) {
  p <- file.path(PANEL_AF_ROOT, tolower(pop), sprintf("chr%s.af.tsv.gz", chr))
  if (!file.exists(p)) return(NULL)
  a <- fread(p, select = c("position", "bim_a1", "bim_a2", "af_a1"),
             showProgress = FALSE)
  a[, `:=`(bim_a1 = toupper(bim_a1), bim_a2 = toupper(bim_a2))]
  unique(a, by = "position")
}
# Frequency of `allele` in the panel, read under the LETTERS assumption (i.e.
# assuming the dataset and the panel are on the same strand).  For a palindrome
# both panel alleles are the same two letters, so this is exactly the reading
# whose correctness we are testing -- if the dataset is actually reverse-strand
# the returned value will be ~1-x rather than ~x, which is what the calibrated
# rule below detects.  NA where the variant is absent or the alleles disagree.
#
# Requires the full allele PAIR to match the panel record, not just the one
# allele being looked up: a one-allele match would accept a record whose other
# allele differs (a different variant at the same position).  Measured at
# 878 of 14,213,794 GWAS-side verdicts (0.0062%) before this was tightened --
# immaterial, but there is no reason to carry it.
panel_freq_of <- function(allele, other, bim_a1, bim_a2, af_a1) {
  pair_ok <- (allele == bim_a1 & other == bim_a2) |
             (allele == bim_a2 & other == bim_a1)
  fifelse(!pair_ok | is.na(pair_ok), NA_real_,
          fifelse(allele == bim_a1, af_a1, 1 - af_a1))
}
# Calibration of record: resolve_palindrome_af() in
# src/seqfunc_v2/finemap/common.R.  Keep the three constants in sync.
AF_PAL_MAF_MAX  <- 0.40   # near 0.5 the two hypotheses are indistinguishable
AF_PAL_MAX_DIST <- 0.15   # neither hypothesis fits -> abstain, not evidence
AF_PAL_MARGIN   <- 0.10   # required separation before calling a strand verdict
# Returns "same" / "opposite" / NA (abstain) per variant.
strand_verdict_vs_panel <- function(obs_af, panel_af) {
  d_same <- abs(obs_af - panel_af)
  d_opp  <- abs(obs_af - (1 - panel_af))
  testable <- is.finite(obs_af) & is.finite(panel_af) &
    pmin(obs_af, 1 - obs_af)     <= AF_PAL_MAF_MAX &
    pmin(panel_af, 1 - panel_af) <= AF_PAL_MAF_MAX &
    pmin(d_same, d_opp)          <= AF_PAL_MAX_DIST
  fifelse(!testable, NA_character_,
    fifelse((d_opp - d_same) >= AF_PAL_MARGIN, "same",
      fifelse((d_same - d_opp) >= AF_PAL_MARGIN, "opposite", NA_character_)))
}

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

# Env-overridable so COLOC can be pointed at a regenerated, LD-pinned eQTL set.
EQTL_SUSIE_DIR  <- Sys.getenv("EQTL_SUSIE_DIR", unset = file.path(FM_DIR, "results/eqtl_susie"))
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
# Prefer the allele-frequency-augmented copy when available. It is a LEFT-JOIN of
# `af` onto this same file, so the seven core columns are byte-identical; the
# extra column enables the palindrome frequency guard below. `data/sumstats/` is
# never modified.
.af_src <- file.path(FM_DIR, "config", "gwas_af_sources.tsv")
if (file.exists(.af_src) && !identical(Sys.getenv("COLOC_USE_AF_SUMSTATS", unset = "1"), "0")) {
  .af_tab <- tryCatch(data.table::fread(.af_src), error = function(e) NULL)
  if (!is.null(.af_tab) && "study_name" %in% names(.af_tab)) {
    .hit <- .af_tab[study_name == gwas_name & file.exists(af_sumstats_path)]
    if (nrow(.hit) == 1L) {
      gwas_path <- .hit$af_sumstats_path[[1L]]
      cat("Using AF-augmented sumstats:", gwas_path, "\n")
    }
  }
}
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

# 1000G panel AF for the palindrome strand check.  Loaded ONCE per run, not per
# gene -- chr_num is a single chromosome (line 14 rejects chr_num == 0), so
# these two tables serve every gene in the loop below.
pa_gwas <- load_panel_af(gwas_ancestry, chr_num)
pa_eqtl <- load_panel_af("eur",         chr_num)
cat("Panel AF for palindrome check: ", gwas_ancestry, " ",
    if (is.null(pa_gwas)) "MISSING" else paste0(nrow(pa_gwas), " variants"),
    " | eQTL-side eur ",
    if (is.null(pa_eqtl)) "MISSING" else paste0(nrow(pa_eqtl), " variants"), "\n", sep = "")

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
                   eqtl_beta = Beta, eqtl_se = SE, eqtl_pval = PVAL,
                   eqtl_eaf = EAF)],
    gwas[, .(merge_key, gwas_a1 = allele1, gwas_a2 = allele2,
             gwas_beta = beta, gwas_se = se, gwas_pval = pval,
             gwas_af = if ("af" %in% names(gwas)) af else NA_real_)],
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

  # --- Strand-ambiguous (palindromic A/T, C/G) SNPs -------------------------
  # These were previously dropped UNCONDITIONALLY (~12.6% of variants per locus).
  # The cost was concrete: PNPLA3 rs738409 is itself a C/G palindrome, so it
  # appeared 0 times in 29,262 chr22 top_snp values and neighbours 3 bp away
  # stood in as the reported top SNP.
  #
  # The drop was justified because strand could not be inferred. That constraint
  # no longer holds. BOTH datasets are certified forward-strand against the same
  # 1000G reference, measured on NON-palindromic SNVs where orientation is fixed
  # by letters alone:
  #   Broadaway eQTL : 6,022 same / 0 opposite   (opp_fraction 0.00e+00)
  #   all 32 GWAS    : opp_fraction <= 1.3e-03   (config/strand_certificates.tsv)
  # With both sides on the same strand, `eqtl_ea == gwas_a1` is a genuine match
  # rather than an ambiguity, so palindromes are alignable by letters.
  #
  # NOTE the certificate is an OFFLINE argument, not a runtime check: nothing in
  # this file reads config/strand_certificates.tsv, and 18 of the 50 registry
  # studies have no certificate at all (the table covers the 35-study finemap
  # portfolio).  The per-variant guard below is therefore the only strand
  # evidence this script actually evaluates.
  #
  # It compares each dataset against 1000G, ANCESTRY-MATCHED ON BOTH SIDES --
  # never GWAS against the eQTL EAF, which is what an earlier version did and
  # which is confounded by construction (Broadaway is European, so against an
  # AFR/EAS GWAS "strand-flipped" and "population-drifted" have an identical
  # frequency signature).  Measured control false-positive rates, using
  # NON-palindromic variants whose orientation is fixed by letters so every flag
  # there is by definition wrong:
  #
  #   comparison                      AFR      AMR      EAS      SAS      EUR
  #   vs eQTL EAF (cross-ancestry)  3.329%   0.020%   3.771%      --    0.005%
  #   vs own-ancestry 1000G         0.004%   0.005%   0.001%   0.001%   0.034%
  #
  # a 300-3,000x reduction, with cor(gwas_af, panel) rising from 0.76-0.79 to
  # 0.986-0.998.  Low FPR alone would be satisfied by a rule that never fires,
  # so recall was measured separately by artificially strand-flipping half of a
  # real palindrome set: recall 0.84-0.92 at precision ~1.0000, with the misses
  # being deliberate abstentions (MAF>0.40 or neither hypothesis fits) rather
  # than failures.  Independent letters-only corroboration: across ~7.6M GWAS
  # and ~424K eQTL variants, complement-only matches = 0.
  #
  # Set COLOC_DROP_PALINDROMIC=1 to restore the old unconditional behaviour.
  n_pre_ambiguous <- nrow(merged)
  merged[, is_palindromic := (gwas_a1 %in% c("A","T") & gwas_a2 %in% c("A","T")) |
                             (gwas_a1 %in% c("C","G") & gwas_a2 %in% c("C","G"))]
  n_palindromic_present <- sum(merged$is_palindromic, na.rm = TRUE)
  if (identical(Sys.getenv("COLOC_DROP_PALINDROMIC", unset = "0"), "1")) {
    # `merged[!is_palindromic]` would be a data.table NOT-JOIN (DT[!i]), not a
    # logical negation, and errors with "not found in calling scope". The
    # pre-refactor code passed a full expression, which evaluates normally; a
    # bare `!symbol` does not.
    merged <- merged[is_palindromic == FALSE]
    n_palindromic_kept <- 0L; n_palindromic_af_tested <- 0L; n_palindromic_af_dropped <- 0L
    pairmiss <- NA_real_
  } else {
    # `eqtl_eaf` is the frequency of `eqtl_ea`, the ORIGINAL eQTL effect allele.
    # The flip above negates `eqtl_beta` only -- it does not touch `eqtl_eaf` --
    # so on flipped rows eqtl_ea == gwas_a2 and freq(gwas_a1) = 1 - eqtl_eaf.
    # After this transform both quantities refer to gwas_a1.  Measured on
    # MVP_NAFLD_EUR chr22 (18,561 flipped rows): cor(gwas_af, eqtl_eaf) = -0.982
    # raw and +0.982 transformed, so the flipped subset is what identifies the
    # direction and only this orientation is consistent with it.
    merged[, eqtl_eaf_oriented := fifelse(allele_flip, 1 - eqtl_eaf, eqtl_eaf)]

    # Two ancestry-concordant checks, never one cross-ancestry check.  If each
    # dataset is independently forward-strand against 1000G, then the two are on
    # the same strand as each other and the letters match above is genuine --
    # and neither comparison crosses an ancestry boundary, so drift cannot enter.
    # The eQTL side uses the EUR panel because Broadaway is European.
    merged[, `:=`(gwas_strand = NA_character_, eqtl_strand = NA_character_)]
    if (!is.null(pa_gwas)) {
      merged[pa_gwas, on = .(eqtl_pos = position),
             `:=`(pg_a1 = i.bim_a1, pg_a2 = i.bim_a2, pg_af = i.af_a1)]
      merged[, gwas_strand := strand_verdict_vs_panel(
        gwas_af, panel_freq_of(gwas_a1, gwas_a2, pg_a1, pg_a2, pg_af))]
    }
    if (!is.null(pa_eqtl)) {
      merged[pa_eqtl, on = .(eqtl_pos = position),
             `:=`(pe_a1 = i.bim_a1, pe_a2 = i.bim_a2, pe_af = i.af_a1)]
      merged[, eqtl_strand := strand_verdict_vs_panel(
        eqtl_eaf, panel_freq_of(eqtl_ea, eqtl_nea, pe_a1, pe_a2, pe_af))]
    }
    # ONLY the GWAS side vetoes.  `eqtl_strand` is computed and recorded because
    # it is a genuine QC statistic that would catch a future eQTL-panel
    # regression, but it does not drop anything.  The two sides are not
    # comparable in quality, measured against the same control (non-palindromic
    # variants, where letters fix orientation so every "opposite" is wrong):
    #
    #   GWAS side : palindromic/control ratio 5.7-10.5 at p 2e-42 to 9e-66 in
    #               five studies -- real, per-study, overwhelming.
    #   eQTL side : 30/79,507 palindromic vs 164/508,313 control = ratio 1.17,
    #               Fisher p = 0.40, 95% CI [0.765, 1.734] SPANS 1.
    #
    # The eQTL side fires at the same rate where it must be wrong as where it
    # would be acted on, i.e. it is at its own noise floor: expected 25.7,
    # observed 30, excess ~4.3 genuine calls.  Letting it veto cost 985 of 2,118
    # total drops (46.5%) to catch those ~4 -- and because the eQTL comparison is
    # identical for all 50 studies, one verdict removed a variant from every
    # study at once, including four AF-less studies with no way to corroborate
    # it.  The calls are also a reference artifact rather than per-variant error:
    # 62.4% cluster in 100kb windows holding >1 (chr11 alone has 137 of 194),
    # they concentrate near MAF 0.5 (median |eaf-0.5| = 0.200 vs 0.319 baseline),
    # and none has an allele-pair mismatch or a multi-allelic position.  At
    # n=379 the panel SE is 0.018, so the 0.10 margin is 5.5 SE and sampling
    # noise cannot cross it -- the cause is regional divergence between the
    # Broadaway cohort and 1kg EUR, not small-n.
    # POSITION-INTEGRITY GUARD. The panel join is by POSITION, so if a study's
    # coordinates are off (liftover or a mis-documented source build) variants
    # join to the WRONG panel record. That error is invisible here for
    # palindromes -- any palindrome pair-matches any palindrome, {C,G} or {A,T}
    # on both sides -- while NON-palindromic variants fail the pair check and
    # drop out silently. So a position error manufactures spurious "opposite"
    # calls in the palindromic stratum against a control that has been quietly
    # cleaned. Measured across 50 studies x 22 chromosomes: palindromic
    # "opposite" calls sit in windows with 79x the non-palindromic pair-mismatch
    # rate (pooled Wilcoxon p = 0), and in hg19-native studies with no
    # coordinate transformation there is NO elevation at all (0.00038 vs
    # 0.00049) -- the fingerprint is present exactly where a position transform
    # happened and absent exactly where none did.
    #
    # Non-palindromic pair-mismatch is therefore a direct, run-time detector of
    # the failure, and it is free: panel_freq_of() already returns NA on a pair
    # mismatch. When it is elevated the AF verdicts cannot be trusted, so the
    # guard abstains wholesale and palindromes fall through to the letters match.
    # The signal is LOCAL, so this must be measured per WINDOW, not per run.
    # Chromosome-wide the rate is 0.0003-0.004 for every study regardless of
    # mode -- position errors are concentrated, and averaging over a chromosome
    # dilutes a 0.09 window into a 0.004 chromosome and detects nothing.
    # Measured per 100kb window: palindromic "opposite" calls sit at 0.0876
    # (join_lift) and 0.4607 (existing) against same-call baselines of 0.00079
    # and 0.03174, while hg19-native studies show 0.00038 vs 0.00049 -- no
    # elevation at all.
    # MEASURED AND RECORDED, NOT ACTED ON. `nonpal_pairmiss` is emitted per gene
    # so the rerun's own output carries the diagnostic and the threshold can be
    # calibrated at portfolio scale afterwards. It deliberately does NOT gate
    # anything yet: the effect is local (a 0.09 window averages to 0.004 across a
    # chromosome, so run-level gating detects nothing), and a windowed threshold
    # could not be validated here -- chr22 carries only 2-8 "opposite" calls per
    # study against the ~11,000 the portfolio-scale measurement pooled, and a
    # trial 0.02 cut fired on 3 of 8 calls in an hg19-native study that has no
    # position problem. Gating canonical code on an unvalidated cut is what
    # produced the af_over_certificate defect; measure first.
    PAIRMISS_MIN_N <- 20L
    WINDOW_BP      <- 1e5
    pairmiss <- NA_real_
    # `pg_*` exist only when a panel sidecar was found for this ancestry.
    if (all(c("pg_a1","pg_a2","pg_af") %in% names(merged))) {
      merged[, np_miss := ifelse(is_palindromic == FALSE & is.finite(pg_af),
                                 as.numeric(is.na(panel_freq_of(gwas_a1, gwas_a2,
                                                                pg_a1, pg_a2, pg_af))),
                                 NA_real_)]
      merged[, af_window := floor(eqtl_pos / WINDOW_BP)]
      merged[, win_pairmiss := {
        n <- sum(!is.na(np_miss))
        if (n >= PAIRMISS_MIN_N) rep(mean(np_miss, na.rm = TRUE), .N) else rep(NA_real_, .N)
      }, by = af_window]
      # Reported as the rate in the windows where the guard actually fired, which
      # is the quantity the portfolio measurement found elevated (79x).
      opp <- merged[is_palindromic == TRUE & gwas_strand %in% "opposite"]
      pairmiss <- if (nrow(opp)) suppressWarnings(mean(opp$win_pairmiss, na.rm = TRUE)) else NA_real_
    }
    merged[, af_testable := is_palindromic & !is.na(gwas_strand)]
    merged[, af_inconsistent := af_testable %in% TRUE & gwas_strand %in% "opposite"]
    n_palindromic_af_tested  <- sum(merged$af_testable, na.rm = TRUE)
    n_palindromic_af_dropped <- sum(merged$af_inconsistent, na.rm = TRUE)
    merged <- merged[!(af_inconsistent %in% TRUE)]
    n_palindromic_kept <- sum(merged$is_palindromic, na.rm = TRUE)
  }
  # Kept as "palindromes PRESENT", which is what it meant when every palindrome
  # was dropped, and what results/qa_campaign/audit/C1_harmonization_audit.R
  # independently recomputes from raw data.  Defining it as rows-lost instead
  # would silently break comparability with both.  The retained/tested/dropped
  # breakdown lives in the three companion columns.
  n_ambiguous_locus <- n_palindromic_present

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
    n_ambiguous = n_ambiguous_locus,   # C1 harmonization tally: strand-ambiguous A/T or C/G palindromes PRESENT
    n_palindromic_kept = n_palindromic_kept,             # retained via the dual forward-strand certificate
    n_palindromic_af_tested = n_palindromic_af_tested,   # palindromes the AF guard could actually adjudicate
    nonpal_pairmiss = pairmiss,                          # position-integrity diagnostic, RECORDED not acted on
    n_palindromic_af_dropped = n_palindromic_af_dropped, # of those, dropped as genuinely opposite-strand
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
    cat("  strand-ambiguous present:", sum(final$n_ambiguous, na.rm = TRUE), "\n")
    if ("n_palindromic_kept" %in% names(final)) {
      cat("    ... retained          :", sum(final$n_palindromic_kept, na.rm = TRUE), "\n")
      cat("    ... AF-adjudicable    :", sum(final$n_palindromic_af_tested, na.rm = TRUE), "\n")
      cat("    ... dropped (AF-flip) :", sum(final$n_palindromic_af_dropped, na.rm = TRUE), "\n")
    }
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
