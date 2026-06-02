#!/usr/bin/env Rscript
# 09_identify_shared_loci.R
# Identify loci that are genome-wide significant (p < 5e-8) in BOTH a EUR
# and a matched EAS GWAS, then define merged windows for SuSiEX joint fine-mapping.
#
# Input:
#   gwas_registry.tsv            — GWAS registry with ancestry, sumstats paths, N
#   data/sumstats/*.tsv          — Preprocessed GWAS summary statistics (hg19)
# Output:
#   results/susiex/shared_loci.csv
#
# Usage: Rscript 09_identify_shared_loci.R
#   Override root: MASLD_PROJECT_ROOT=/path Rscript 09_identify_shared_loci.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
FM_DIR <- Sys.getenv("FM_DIR",
  unset = file.path(Sys.getenv("MASLD_PROJECT_ROOT",
    unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
    "GWAS/finemapping"))

REGISTRY_FILE  <- file.path(FM_DIR, "config/gwas_registry.tsv")
OUTPUT_DIR     <- file.path(FM_DIR, "results/susiex")
GW_PVAL        <- 5e-8       # genome-wide significance threshold
CLUMP_KB       <- 500        # clumping window (each side, kb)
LOCUS_MATCH_KB <- 500        # max distance to call EUR/EAS loci "shared"
FLANK_KB       <- 100        # extra flanking added to merged window

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Define EUR <-> EAS trait pairs
# Registry uses these study_name values for EUR liver enzymes and EAS BBJ.
# Each pair is (eur_study_name, eas_study_name, trait_label).
# ---------------------------------------------------------------------------
TRAIT_PAIRS <- data.frame(
  eur_study  = c("UKBB_ALT",  "UKBB_AST",  "UKBB_GGT"),
  eas_study  = c("BBJ_ALT",   "BBJ_AST",   "BBJ_GGT"),
  trait_pair = c("ALT",       "AST",       "GGT"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Load and parse GWAS registry (skip comment lines starting with #)
# ---------------------------------------------------------------------------
load_registry <- function(path) {
  lines <- readLines(path)
  active <- lines[!grepl("^\\s*#", lines) & nchar(trimws(lines)) > 0]
  fread(text = paste(active, collapse = "\n"), sep = "\t")
}

registry <- load_registry(REGISTRY_FILE)
cat("Registry loaded:", nrow(registry), "active studies\n")

# ---------------------------------------------------------------------------
# Helper: load sumstats and standardise column names
# Preprocessed files have columns: chromosome, position, allele1, allele2,
# beta, se, pval  (EUR UKBB format)
# BBJ files: chromosome, base_pair_location, pos_hg19, effect_allele,
#             other_allele, beta, standard_error, p_value, EAF
# ---------------------------------------------------------------------------
load_sumstats <- function(study_name_arg, registry, fm_dir) {
  row <- registry[registry$study_name == study_name_arg, ]
  if (nrow(row) == 0) stop(paste("Study not found in registry:", study_name_arg))
  path <- file.path(fm_dir, row$sumstats_path)
  if (!file.exists(path)) stop(paste("Sumstats file not found:", path))
  cat("  Loading", study_name_arg, "from", basename(path), "...\n")
  ss <- fread(path)

  # Detect column layout and rename to canonical names:
  #   CHR, POS, A1 (effect), A2 (other), BETA, SE, PVAL
  nms <- tolower(names(ss))
  names(ss) <- nms

  # BBJ format: chromosome, base_pair_location, pos_hg19, effect_allele, other_allele
  if ("pos_hg19" %in% nms) {
    # BBJ files: use pos_hg19 as the hg19 coordinate
    setnames(ss, "pos_hg19",           "POS",  skip_absent = TRUE)
    setnames(ss, "chromosome",         "CHR",  skip_absent = TRUE)
    setnames(ss, "effect_allele",      "A1",   skip_absent = TRUE)
    setnames(ss, "other_allele",       "A2",   skip_absent = TRUE)
    setnames(ss, "beta",               "BETA", skip_absent = TRUE)
    setnames(ss, "standard_error",     "SE",   skip_absent = TRUE)
    setnames(ss, "p_value",            "PVAL", skip_absent = TRUE)
  } else if ("position" %in% nms) {
    # EUR (UKBB / preprocessed) format
    setnames(ss, "chromosome", "CHR",  skip_absent = TRUE)
    setnames(ss, "position",   "POS",  skip_absent = TRUE)
    setnames(ss, "allele1",    "A1",   skip_absent = TRUE)
    setnames(ss, "allele2",    "A2",   skip_absent = TRUE)
    setnames(ss, "beta",       "BETA", skip_absent = TRUE)
    setnames(ss, "se",         "SE",   skip_absent = TRUE)
    setnames(ss, "pval",       "PVAL", skip_absent = TRUE)
  } else {
    stop(paste("Unrecognised sumstats format for study:", study_name_arg,
               "| columns:", paste(names(ss), collapse = ", ")))
  }

  # Keep only the canonical columns; drop rows with NA POS or PVAL
  ss <- ss[, c("CHR", "POS", "A1", "A2", "BETA", "SE", "PVAL"), with = FALSE]
  # Force numeric: fread may read pval as character if file has mixed formatting
  ss$PVAL <- as.numeric(ss$PVAL)
  ss$BETA <- as.numeric(ss$BETA)
  ss$SE   <- as.numeric(ss$SE)
  ss <- ss[!is.na(POS) & !is.na(PVAL)]
  ss$CHR <- as.integer(ss$CHR)
  ss$POS <- as.integer(ss$POS)
  cat("    Loaded", nrow(ss), "variants\n")
  ss
}

# ---------------------------------------------------------------------------
# Helper: clump GW-sig variants into loci by sliding 500kb windows
# Sorts by p-value; iteratively picks top hit, masks ±CLUMP_KB variants.
# Returns data.frame of lead SNPs with columns: chr, pos, pval.
# ---------------------------------------------------------------------------
clump_loci <- function(ss, pval_thresh = GW_PVAL, clump_kb = CLUMP_KB) {
  sig <- ss[PVAL < pval_thresh][order(PVAL)]
  if (nrow(sig) == 0) return(data.frame(chr = integer(), pos = integer(),
                                         pval = numeric()))
  leads <- vector("list", 0)
  used  <- rep(FALSE, nrow(sig))
  for (i in seq_len(nrow(sig))) {
    if (used[i]) next
    lead <- sig[i]
    leads <- c(leads, list(data.frame(chr = lead$CHR, pos = lead$POS,
                                      pval = lead$PVAL)))
    # Mask all variants within ±clump_kb on the same chromosome
    within_window <- sig$CHR == lead$CHR &
                     abs(sig$POS - lead$POS) <= clump_kb * 1000
    used[within_window] <- TRUE
  }
  do.call(rbind, leads)
}

# ---------------------------------------------------------------------------
# Helper: match EUR and EAS loci within LOCUS_MATCH_KB
# Returns data.frame of shared loci with columns from both ancestries.
# ---------------------------------------------------------------------------
match_loci <- function(eur_leads, eas_leads, eur_ss, eas_ss,
                       match_kb = LOCUS_MATCH_KB, flank_kb = FLANK_KB,
                       eur_study, eas_study, trait_pair) {
  if (nrow(eur_leads) == 0 || nrow(eas_leads) == 0) return(NULL)
  shared <- vector("list", 0)
  for (i in seq_len(nrow(eur_leads))) {
    el <- eur_leads[i, ]
    # Find EAS leads on same chr within match_kb
    hits <- eas_leads[eas_leads$chr == el$chr &
                      abs(eas_leads$pos - el$pos) <= match_kb * 1000, ]
    if (nrow(hits) == 0) next
    # Use the nearest EAS lead if multiple match
    hits$dist <- abs(hits$pos - el$pos)
    ea_lead <- hits[which.min(hits$dist), ]

    # Merged window: encompass both lead SNP positions ±flank_kb
    win_start <- max(0, min(el$pos, ea_lead$pos) - flank_kb * 1000)
    win_end   <- max(el$pos, ea_lead$pos) + flank_kb * 1000

    # Find all EUR/EAS variants in the window for window-boundary definition
    eur_win <- eur_ss[CHR == el$chr & POS >= win_start & POS <= win_end]
    eas_win <- eas_ss[CHR == el$chr & POS >= win_start & POS <= win_end]
    if (nrow(eur_win) < 10 || nrow(eas_win) < 10) {
      warning(sprintf("Skipping shared locus chr%d:%d — too few variants (EUR:%d EAS:%d)",
                      el$chr, el$pos, nrow(eur_win), nrow(eas_win)))
      next
    }

    locus_id <- sprintf("locus_%s_chr%d_%d", trait_pair, el$chr, el$pos)
    shared <- c(shared, list(data.frame(
      locus_id       = locus_id,
      chr            = el$chr,
      window_start   = win_start,
      window_end     = win_end,
      eur_gwas       = eur_study,
      eas_gwas       = eas_study,
      eur_lead_snp   = sprintf("%d:%d:%s:%s", el$chr, el$pos,
                               eur_win[POS == el$pos]$A1[1],
                               eur_win[POS == el$pos]$A2[1]),
      eur_lead_pos   = el$pos,
      eur_lead_p     = el$pval,
      eas_lead_snp   = sprintf("%d:%d:%s:%s", ea_lead$chr, ea_lead$pos,
                               eas_win[POS == ea_lead$pos]$A1[1],
                               eas_win[POS == ea_lead$pos]$A2[1]),
      eas_lead_pos   = ea_lead$pos,
      eas_lead_p     = ea_lead$pval,
      trait_pair     = trait_pair,
      stringsAsFactors = FALSE
    )))
  }
  if (length(shared) == 0) return(NULL)
  do.call(rbind, shared)
}

# ---------------------------------------------------------------------------
# Main loop over trait pairs
# ---------------------------------------------------------------------------
all_shared <- vector("list", nrow(TRAIT_PAIRS))

for (k in seq_len(nrow(TRAIT_PAIRS))) {
  pair   <- TRAIT_PAIRS[k, ]
  cat("\n=== Trait pair:", pair$trait_pair, "(", pair$eur_study, "<->",
      pair$eas_study, ") ===\n")

  eur_ss <- tryCatch(
    load_sumstats(pair$eur_study, registry, FM_DIR),
    error = function(e) { message("  ERROR loading EUR: ", e$message); NULL }
  )
  eas_ss <- tryCatch(
    load_sumstats(pair$eas_study, registry, FM_DIR),
    error = function(e) { message("  ERROR loading EAS: ", e$message); NULL }
  )
  if (is.null(eur_ss) || is.null(eas_ss)) {
    cat("  Skipping pair due to load error\n")
    next
  }

  cat("  Clumping EUR GW-sig loci (p <", GW_PVAL, ") ...\n")
  eur_leads <- clump_loci(eur_ss)
  cat("    EUR lead loci:", nrow(eur_leads), "\n")

  cat("  Clumping EAS GW-sig loci ...\n")
  eas_leads <- clump_loci(eas_ss)
  cat("    EAS lead loci:", nrow(eas_leads), "\n")

  shared <- match_loci(
    eur_leads, eas_leads, eur_ss, eas_ss,
    eur_study  = pair$eur_study,
    eas_study  = pair$eas_study,
    trait_pair = pair$trait_pair
  )
  if (!is.null(shared)) {
    cat("    Shared loci:", nrow(shared), "\n")
    all_shared[[k]] <- shared
  } else {
    cat("    No shared loci found\n")
  }
}

result <- do.call(rbind, Filter(Negate(is.null), all_shared))

if (is.null(result) || nrow(result) == 0) {
  cat("\nNo shared loci identified across any trait pair.\n")
  # Write empty file so downstream scripts don't error
  result <- data.frame(
    locus_id = character(), chr = integer(),
    window_start = integer(), window_end = integer(),
    eur_gwas = character(), eas_gwas = character(),
    eur_lead_snp = character(), eur_lead_pos = integer(), eur_lead_p = numeric(),
    eas_lead_snp = character(), eas_lead_pos = integer(), eas_lead_p = numeric(),
    trait_pair = character()
  )
} else {
  cat("\n=== Summary ===\n")
  cat("Total shared loci:", nrow(result), "\n")
  print(table(result$trait_pair))
}

out_path <- file.path(OUTPUT_DIR, "shared_loci.csv")
write_csv(result, out_path)
cat("\nOutput written to:", out_path, "\n")
