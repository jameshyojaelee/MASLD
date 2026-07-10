#!/usr/bin/env Rscript
# qc_positive_controls.R
# ----------------------------------------------------------------------------
# Positive-control QC harness for the harmonized MVP (Million Veteran Program)
# GWAS sumstats produced by format_mvp_for_coloc.R.
#
# Its job is to CATCH harmonization bugs, not to bless the output. It is
# deliberately adversarial. For every existing
#   GWAS/finemapping/data/sumstats/MVP_<name>_reformatted_hg19.tsv
# it checks:
#   1. Schema / integrity  : exact header, chr 1..22, finite beta, se>0,
#                            pval in [0,1], exact-variant duplicates, row count.
#   2. liftOver / build     : the harmonizer lifts hg38 -> hg19. We confirm the
#                            file is hg19 by testing whether the genome-wide
#                            signal at known liver loci sits at the *hg19*
#                            coordinate and NOT at the *hg38* coordinate
#                            (a silent liftOver failure leaves hg38 coords ->
#                            the peak appears in the hg38 window = BUILD_MISMATCH).
#   3. Sign concordance     : (a) internal cross-ancestry consistency -- for a
#                            given locus+trait the aligned-allele beta must have
#                            the SAME sign across ancestries (a flip = allele
#                            orientation / OR-inversion bug); (b) biological
#                            anchor -- PNPLA3-G / TM6SF2-T raise disease & enzymes
#                            (expected beta > 0).
#   4. Power-aware presence : per-stratum count of genome-wide-sig (p<5e-8)
#                            variants; expect many (EUR) -> few (AFR/AMR) -> ~0
#                            (EAS NAFLD, 505 cases).
#   5. FIB-4 components      : confirm ALT / AST / Platelet strata exist.
#
# Handles missing files gracefully (the harmonization array may still be running)
# and reports which of the 32 strata currently exist. Re-run once all 32 are done.
#
# Usage:  micromamba activate rnaseq && Rscript qc_positive_controls.R
# ----------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SUMSTATS_DIR <- file.path(ROOT, "GWAS/finemapping/data/sumstats")
MANIFEST     <- file.path(ROOT, "GWAS/MR_Data/MVP/mvp_manifest.tsv")
OUT_CSV      <- file.path(ROOT, "GWAS/MR_Data/MVP/qc_positive_controls.csv")

EXPECTED_HEADER <- c("chromosome","position","allele1","allele2","beta","se","pval")
WINDOW   <- 5000L      # +/- bp around a positive-control locus
GWSIG    <- 5e-8
SIGNAL_P <- 1e-6       # a window with min-p below this is "has signal" (for build test)
NOMINAL  <- 0.05       # only strata with locus p<NOMINAL enter the sign-concordance vote

setDTthreads(min(4L, parallel::detectCores()))

# --- positive-control liver loci (HG19 + HG38 positions) --------------------
# a_eff = disease / enzyme-INCREASING allele (the allele we orient beta to).
# exp_sign_disease = expected sign of the a_eff beta for disease/enzyme traits.
#   PNPLA3-G, TM6SF2-T : rock-solid risk-raising at NAFLD/ALT/AST.
#   GCKR-T            : raises ALT / triglycerides / NAFLD (solid).
#   MARC1, HSD17B13, TRIB1: sign NOT asserted (indel / region / direction debated)
#     -> presence + internal cross-ancestry consistency only.
# sign_checkable = SNP with unambiguous alleles -> participates in sign vote.
loci <- data.table(
  locus            = c("PNPLA3",   "TM6SF2",     "GCKR",      "HSD17B13",   "MARC1",     "TRIB1"),
  rsid             = c("rs738409", "rs58542926", "rs1260326", "rs72613567", "rs2642438", "TRIB1_region"),
  chr              = c(22L,        19L,          2L,          4L,           1L,          8L),
  hg19             = c(44324727L,  19379549L,    27730940L,   88231392L,    220970028L,  126490972L),
  hg38             = c(43928847L,  19268740L,    27508073L,   87310241L,    220796686L,  125478726L),
  a_eff            = c("G",        "T",          "T",         NA,           "A",         NA),
  a_oth            = c("C",        "C",          "C",         NA,           "G",         NA),
  exp_sign_disease = c( 1L,         1L,           1L,         NA_integer_,  NA_integer_, NA_integer_),
  sign_checkable   = c(TRUE,       TRUE,         TRUE,        FALSE,        TRUE,        FALSE)
)

# ---------------------------------------------------------------------------
manifest <- fread(MANIFEST)
setnames(manifest, 1, "gwas_name")
# trait label = gwas_name minus "MVP_" prefix and "_<ancestry>" suffix
manifest[, trait := sub("^MVP_", "", gwas_name)]
manifest[, trait := mapply(function(t, a) sub(paste0("_", a, "$"), "", t), trait, ancestry)]
all_names <- manifest$gwas_name
fpath <- function(nm) file.path(SUMSTATS_DIR, paste0(nm, "_reformatted_hg19.tsv"))
manifest[, file := fpath(gwas_name)]
manifest[, exists := file.exists(file)]

cat(sprintf("\n=== MVP positive-control QC harness ===\n"))
cat(sprintf("strata in manifest : %d\n", nrow(manifest)))
cat(sprintf("files present now  : %d\n", sum(manifest$exists)))
cat(sprintf("files missing      : %d  (harmonization array may still be running)\n\n",
            sum(!manifest$exists)))

# ---------------------------------------------------------------------------
# per-file processor: returns list(stratum = 1-row dt, loci = per-locus dt)
process_one <- function(nm, ancestry, trait, trait_type, file) {
  res <- list(stratum = NULL, loci = NULL, error = NA_character_)
  hdr <- tryCatch(names(fread(file, nrows = 0)), error = function(e) NULL)
  schema_ok <- !is.null(hdr) && identical(hdr, EXPECTED_HEADER)
  if (!schema_ok) {
    res$stratum <- data.table(gwas_name = nm, ancestry = ancestry, trait = trait,
      trait_type = trait_type, n_var = NA_integer_, n_dup_exact = NA_integer_,
      schema_ok = FALSE, header = if (is.null(hdr)) "UNREADABLE" else paste(hdr, collapse=","),
      chr_ok = NA, beta_finite_ok = NA, se_pos_ok = NA, pval_range_ok = NA,
      min_pval = NA_real_, n_gwsig = NA_integer_, integrity_ok = FALSE)
    res$error <- "schema_fail"
    return(res)
  }
  d <- tryCatch(
    fread(file, colClasses = list(integer = "chromosome", character = c("allele1","allele2"))),
    error = function(e) NULL)
  if (is.null(d) || !nrow(d)) {
    res$stratum <- data.table(gwas_name = nm, ancestry = ancestry, trait = trait,
      trait_type = trait_type, n_var = 0L, n_dup_exact = NA_integer_, schema_ok = TRUE,
      header = paste(hdr, collapse=","), chr_ok = NA, beta_finite_ok = NA, se_pos_ok = NA,
      pval_range_ok = NA, min_pval = NA_real_, n_gwsig = 0L, integrity_ok = FALSE)
    res$error <- "empty"
    return(res)
  }

  # --- integrity ---
  chr_ok        <- all(d$chromosome >= 1L & d$chromosome <= 22L)
  beta_finite   <- all(is.finite(d$beta))
  se_pos        <- all(is.finite(d$se) & d$se > 0)
  pval_range    <- all(is.finite(d$pval) & d$pval >= 0 & d$pval <= 1)
  key <- paste(d$chromosome, d$position, d$allele1, d$allele2, sep=":")
  n_dup <- sum(duplicated(key))
  n_gwsig <- sum(d$pval < GWSIG, na.rm = TRUE)
  integrity_ok <- chr_ok && beta_finite && se_pos && pval_range && n_dup == 0L

  stratum <- data.table(gwas_name = nm, ancestry = ancestry, trait = trait,
    trait_type = trait_type, n_var = nrow(d), n_dup_exact = n_dup, schema_ok = TRUE,
    header = paste(hdr, collapse=","), chr_ok = chr_ok, beta_finite_ok = beta_finite,
    se_pos_ok = se_pos, pval_range_ok = pval_range, min_pval = min(d$pval, na.rm = TRUE),
    n_gwsig = n_gwsig, integrity_ok = integrity_ok)

  # --- positive-control loci ---
  setkey(d, chromosome, position)
  lrows <- vector("list", nrow(loci))
  for (i in seq_len(nrow(loci))) {
    L <- loci[i]
    w19 <- d[chromosome == L$chr & position >= L$hg19 - WINDOW & position <= L$hg19 + WINDOW]
    w38 <- d[chromosome == L$chr & position >= L$hg38 - WINDOW & position <= L$hg38 + WINDOW]
    hg19_n   <- nrow(w19);  hg38_n <- nrow(w38)
    hg19_minp <- if (hg19_n) min(w19$pval) else NA_real_
    hg38_minp <- if (hg38_n) min(w38$pval) else NA_real_

    lead_pos <- lead_a1 <- lead_a2 <- NA; lead_beta <- aligned_beta <- NA_real_
    lead_dist <- NA_integer_; allele_match <- NA
    if (hg19_n) {
      lead <- w19[which.min(pval)]
      lead_pos <- lead$position; lead_a1 <- lead$allele1; lead_a2 <- lead$allele2
      lead_beta <- lead$beta; lead_dist <- as.integer(lead$position - L$hg19)
      if (L$sign_checkable) {
        # prefer the exact bi-allelic SNP matching the locus alleles for the sign call
        m <- w19[(allele1 == L$a_eff & allele2 == L$a_oth) |
                 (allele1 == L$a_oth & allele2 == L$a_eff)]
        if (nrow(m)) {
          mm <- m[which.min(pval)]
          allele_match <- TRUE
          lead_pos <- mm$position; lead_a1 <- mm$allele1; lead_a2 <- mm$allele2
          lead_beta <- mm$beta; lead_dist <- as.integer(mm$position - L$hg19)
          aligned_beta <- if (mm$allele1 == L$a_eff) mm$beta else -mm$beta
        } else {
          allele_match <- FALSE
        }
      }
    }
    # build verdict (only informative where a real signal exists)
    build <- "uninformative_no_signal"
    if (!is.na(hg19_minp) && hg19_minp < SIGNAL_P &&
        (is.na(hg38_minp) || hg19_minp <= hg38_minp)) {
      build <- "hg19_confirmed"
    } else if (!is.na(hg38_minp) && hg38_minp < SIGNAL_P &&
               (is.na(hg19_minp) || hg38_minp < hg19_minp)) {
      build <- "BUILD_MISMATCH_signal_at_hg38"
    }

    lrows[[i]] <- data.table(
      gwas_name = nm, ancestry = ancestry, trait = trait, trait_type = trait_type,
      locus = L$locus, rsid = L$rsid, chr = L$chr,
      hg19_pos = L$hg19, hg38_pos = L$hg38,
      n_in_hg19_window = hg19_n, hg19_min_p = hg19_minp,
      n_in_hg38_window = hg38_n, hg38_min_p = hg38_minp,
      present_5kb = hg19_n > 0L, gwsig_at_locus = !is.na(hg19_minp) & hg19_minp < GWSIG,
      lead_pos = lead_pos, lead_dist_bp = lead_dist,
      lead_allele1 = lead_a1, lead_allele2 = lead_a2, lead_beta = lead_beta,
      a_eff = L$a_eff, aligned_beta = aligned_beta, allele_match = allele_match,
      exp_sign_disease = L$exp_sign_disease, sign_checkable = L$sign_checkable,
      build_verdict = build)
  }
  res$stratum <- stratum
  res$loci <- rbindlist(lrows)
  res
}

# ---------------------------------------------------------------------------
present <- manifest[exists == TRUE]
strat_list <- list(); loci_list <- list()
for (i in seq_len(nrow(present))) {
  r <- present[i]
  cat(sprintf("[%2d/%2d] %-24s ", i, nrow(present), r$gwas_name)); flush.console()
  out <- tryCatch(process_one(r$gwas_name, r$ancestry, r$trait, r$trait_type, r$file),
                  error = function(e) list(stratum = NULL, loci = NULL, error = conditionMessage(e)))
  if (!is.null(out$stratum)) {
    strat_list[[length(strat_list)+1]] <- out$stratum
    if (!is.null(out$loci)) loci_list[[length(loci_list)+1]] <- out$loci
    s <- out$stratum
    cat(sprintf("n=%s integ=%s gwsig=%s\n",
        format(s$n_var, big.mark=","), s$integrity_ok, s$n_gwsig))
  } else {
    cat(sprintf("ERROR: %s\n", out$error))
  }
  gc(verbose = FALSE)
}

if (!length(strat_list)) { cat("\nNo readable MVP strata yet. Re-run when files exist.\n"); quit(status = 0) }

strata <- rbindlist(strat_list, fill = TRUE)
loci_dt <- if (length(loci_list)) rbindlist(loci_list, fill = TRUE) else data.table()

# --- sign-concordance vote: per (locus, trait) across ancestries -----------
flips <- data.table()
if (nrow(loci_dt)) {
  votable <- loci_dt[sign_checkable == TRUE & allele_match == TRUE &
                     !is.na(aligned_beta) & !is.na(hg19_min_p) & hg19_min_p < NOMINAL]
  if (nrow(votable)) {
    flips <- votable[, .(
        n_strata = .N,
        ancestries = paste(ancestry, collapse=","),
        signs = paste(sprintf("%s%+0.3f", ancestry, aligned_beta), collapse="; "),
        n_pos = sum(aligned_beta > 0), n_neg = sum(aligned_beta < 0),
        flip = (sum(aligned_beta > 0) > 0 & sum(aligned_beta < 0) > 0)
      ), by = .(locus, trait, exp_sign_disease)]
  }
}

# --- write report ----------------------------------------------------------
# Long CSV at (stratum x locus) grain with stratum-level integrity denormalized.
report <- if (nrow(loci_dt)) {
  merge(loci_dt,
        strata[, .(gwas_name, n_var, n_dup_exact, schema_ok, integrity_ok,
                   stratum_min_pval = min_pval, stratum_n_gwsig = n_gwsig)],
        by = "gwas_name", all.x = TRUE)
} else strata
fwrite(report, OUT_CSV)

# ===========================================================================
# STDOUT SUMMARY
# ===========================================================================
cat("\n\n================ SUMMARY ================\n")

cat(sprintf("\n[1] FILE PRESENCE: %d/%d strata harmonized\n", sum(manifest$exists), nrow(manifest)))
miss <- manifest[exists == FALSE, gwas_name]
if (length(miss)) cat("    MISSING:", paste(miss, collapse=", "), "\n")

cat("\n[2] SCHEMA / INTEGRITY (per stratum):\n")
print(strata[, .(gwas_name, ancestry, trait, n_var, n_dup_exact,
                 schema_ok, integrity_ok, n_gwsig)], nrow = 200)
bad <- strata[integrity_ok == FALSE]
if (nrow(bad)) {
  cat("\n    *** INTEGRITY FAILURES ***\n")
  print(bad[, .(gwas_name, schema_ok, chr_ok, beta_finite_ok, se_pos_ok,
                pval_range_ok, n_dup_exact)])
} else cat("    all readable strata PASS schema + integrity.\n")

cat("\n[3] POWER-AWARE PRESENCE (genome-wide-sig p<5e-8 per stratum):\n")
pw <- strata[order(trait, ancestry), .(gwas_name, ancestry, trait, trait_type, n_gwsig)]
print(pw, nrow = 200)

if (nrow(loci_dt)) {
  cat("\n[4] POSITIVE-CONTROL LOCI -- min p-value at hg19 position (+/-5kb):\n")
  pc <- dcast(loci_dt, gwas_name + ancestry + trait ~ locus, value.var = "hg19_min_p")
  print(pc, nrow = 200)

  cat("\n[5] BUILD CHECK (liftOver hg38->hg19): per-locus verdicts where signal exists\n")
  cat("    'hg19_confirmed' = peak at hg19 coord (GOOD). 'BUILD_MISMATCH_signal_at_hg38' = RED FLAG.\n")
  bv <- loci_dt[build_verdict != "uninformative_no_signal",
                .N, by = .(locus, build_verdict)][order(locus)]
  print(bv)
  mismatch <- loci_dt[build_verdict == "BUILD_MISMATCH_signal_at_hg38"]
  if (nrow(mismatch)) {
    cat("\n    *** BUILD MISMATCH DETECTED ***\n")
    print(mismatch[, .(gwas_name, locus, hg19_min_p, hg38_min_p)])
  } else cat("    no build mismatches among strata with signal.\n")

  cat("\n[6] SIGN CONCORDANCE (aligned-allele beta, per locus x trait across ancestry):\n")
  cat("    a_eff = disease/enzyme-raising allele; aligned_beta>0 means that allele raises the trait.\n")
  if (nrow(flips)) {
    print(flips[order(locus, trait)], nrow = 200)
    fl <- flips[flip == TRUE]
    if (nrow(fl)) {
      cat("\n    *** SIGN FLIPS (allele-orientation bug suspected) ***\n")
      print(fl[, .(locus, trait, ancestries, signs)])
    } else cat("    no cross-ancestry sign flips among checkable loci.\n")
    # biological anchor: PNPLA3/TM6SF2/GCKR a_eff should be POSITIVE on disease/enzyme traits
    anchor <- loci_dt[sign_checkable == TRUE & allele_match == TRUE & !is.na(aligned_beta) &
                      !is.na(exp_sign_disease) & hg19_min_p < NOMINAL &
                      trait %in% c("NAFLD","Cirrhosis","ChronLiver","Phe571_8","Phe571_81",
                                   "ALT","AST")]
    if (nrow(anchor)) {
      anchor[, dir_ok := sign(aligned_beta) == exp_sign_disease]
      cat(sprintf("\n    Biological anchor (PNPLA3-G/TM6SF2-T/GCKR-T raise disease+enzymes):\n"))
      cat(sprintf("      %d/%d anchor calls match expected risk direction.\n",
                  sum(anchor$dir_ok), nrow(anchor)))
      wrong <- anchor[dir_ok == FALSE]
      if (nrow(wrong)) {
        cat("      *** WRONG-DIRECTION anchors (possible OR-inversion / beta-sign bug) ***\n")
        print(wrong[, .(gwas_name, locus, a_eff, lead_allele1, lead_allele2,
                        aligned_beta, hg19_min_p)])
      }
    }
  } else cat("    not enough significant locus hits yet for a sign vote.\n")
}

cat("\n[7] FIB-4 COMPONENTS (ALT, AST, Platelet must be present):\n")
for (comp in c("ALT","AST","Platelet")) {
  have <- manifest[trait == comp & exists == TRUE, ancestry]
  want <- manifest[trait == comp, ancestry]
  cat(sprintf("    %-9s : %d/%d ancestries present (%s)\n",
      comp, length(have), length(want),
      if (length(have)) paste(have, collapse=",") else "NONE YET"))
}

cat(sprintf("\nReport written: %s\n", OUT_CSV))
cat(sprintf("(%d/%d strata present; RE-RUN after all 32 are harmonized.)\n\n",
            sum(manifest$exists), nrow(manifest)))
