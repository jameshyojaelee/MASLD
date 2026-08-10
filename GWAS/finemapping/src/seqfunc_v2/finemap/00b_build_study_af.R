#!/usr/bin/env Rscript

# Attach effect-allele frequency to each study's prepared summary statistics.
#
# Design rule: NEVER regenerate `data/sumstats/`.  Those files are read by the
# canonical COLOC pipeline and their paths are hashed into the frozen release
# manifest, so re-deriving them would breach the mutation firewall and risk
# silent drift (the formatters apply a maf >= 0.01 filter the moment an `af`
# column exists, which would change the variant set).
#
# Instead this LEFT-JOINS an `af` column onto the existing file and writes the
# result to `data/sumstats_af/<STUDY>_af.tsv`.  The seven core columns are
# therefore unchanged BY CONSTRUCTION, whatever the join does; AF coverage is
# just a reported statistic, and a study whose source has no frequency simply
# gets af = NA and is resolved downstream by its strand certificate instead.
#
# Sources fall into three shapes:
#   existing  -- the file already carries `af`; copied through unchanged
#   join      -- a source on the SAME build (hg19); joined on chr/pos/alleles
#   join_lift -- a source on hg38; lifted to hg19 with the project's own
#                hg38ToHg19 chain (the same chain format_mvp_for_coloc.R uses)
#                and then joined
#   none      -- no allele frequency exists in any on-disk source
#
# Usage:  Rscript 00b_build_study_af.R [study_name ...]   (default: all main)

source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT", unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "GWAS", "finemapping", "src", "seqfunc_v2", "finemap", "common.R"
))

OUT_DIR <- file.path(FM_ROOT, "data", "sumstats_af")
SOURCES <- file.path(FM_ROOT, "config", "gwas_af_sources.tsv")
CHAIN <- file.path(PROJECT_ROOT, "data", "broadaway_eqtl", "hg38ToHg19.over.chain")
OAHMED <- "/gpfs/commons/groups/sanjana_lab/oahmed/liver/gwas/250203_Batch_1/raw_sumstats"
MRD <- file.path(PROJECT_ROOT, "GWAS", "MR_Data")

# --- per-study source specification ---------------------------------------
# chr/pos/ea/oa/af name the columns IN THE SOURCE.  `ea` must be the allele the
# frequency refers to, and it must correspond to `allele1` in the reformatted
# file (verified for every entry below).
spec <- function(study, mode, path = NA_character_, chr = NA_character_, pos = NA_character_,
                 ea = NA_character_, oa = NA_character_, af = NA_character_, note = "") {
  data.table(study_name = study, mode = mode, source_path = path,
             src_chr = chr, src_pos = pos, src_ea = ea, src_oa = oa, src_af = af, note = note)
}

# The harmonised intermediate, NOT the raw `biomarkers-*.bgz` -- the raw file has
# a different schema (chr/pos/ref/alt plus per-population af_AFR, af_EUR, ...)
# and no pos_hg19.  The intermediate is what the finemapping formatter consumes.
panukbb <- function(pop, tr) spec(
  sprintf("PanUKBB_%s_%s", pop, tr), "join",
  file.path(MRD, "PanUKBB", sprintf("PanUKBB_%s_%s_harmonised_hg38.tsv.gz", pop, tr)),
  "chromosome", "pos_hg19", "effect_allele", "other_allele", "EAF",
  "source carries pos_hg19; no liftover")
bbj <- function(tr) spec(
  sprintf("BBJ_%s", tr), "join",
  file.path(MRD, "BBJ", sprintf("BBJ_%s_harmonised_hg38.tsv.gz", tr)),
  "chromosome", "pos_hg19", "effect_allele", "other_allele", "EAF",
  "source carries pos_hg19; no liftover")
# Binary disease endpoints are phecodes under raw/phecodes; the quantitative
# labs are inverse-normalised and live under raw/labs.  Both layouts share the
# SNP_ID/chrom/pos/ref/alt/ea/af prefix, and in both `af` is the frequency of
# `ea`, which the formatter maps to allele1.  Phecode digits are taken from
# GWAS/MR_Data/MVP/build_mvp_registry.R:28, the single mapping of record.
MVP_PHECODE <- c(NAFLD = "Phe_571_5", Cirrhosis = "Phe_571_51", ChronLiver = "Phe_571")
mvp <- function(phe, pop) spec(
  sprintf("MVP_%s_%s", phe, pop), "join_lift",
  if (phe %in% names(MVP_PHECODE))
    file.path(MRD, "MVP", "raw", "phecodes",
              sprintf("MVP_R4.1000G_AGR.%s.%s.GIA.dbGaP.txt.gz", MVP_PHECODE[[phe]], pop))
  else
    file.path(MRD, "MVP", "raw", "labs",
              sprintf("MVP_R4.1000G_AGR.%s_Mean_INT.%s.GIA.dbGaP.txt.gz", phe, pop)),
  "chrom", "pos", "ea", NA_character_, "af", "hg38 dbGaP; lifted with hg38ToHg19")
# The GWAS Catalog harmonised UKBB files DECLARE effect_allele_frequency but the
# column is 100% NA (verified over 200k rows in all three) -- the same phantom
# column as the Anstee harmonised file.  Treat as no AF rather than pretending.
# A future option is the PanUKBB raw biomarker file's af_EUR for the same cohort,
# but that is a cross-source positional join and AF is confirmatory here, not the
# primary resolver.
ukbb <- function(tr, gcst) spec(
  sprintf("UKBB_%s", tr), "none",
  note = sprintf("GCST%s declares effect_allele_frequency but it is 100%% NA", gcst))
oahmed <- function(study) spec(
  study, "join", file.path(OAHMED, paste0(study, ".tsv")),
  "chromosome", "base_pair_location", "effect_allele", "other_allele",
  "effect_allele_frequency", "batch-1 raw original; hg19 native, 1:1 with the preprocessed file")
finngen <- function(study, raw) spec(
  study, "join_lift", file.path(MRD, "FinnGen", raw),
  "#chrom", "pos", "alt", "ref", "af_alt", "FinnGen R12 hg38; lifted with hg38ToHg19")

SPECS <- rbindlist(list(
  spec("2023_36280732_NAFLD_deCode_EUR", "existing", note = "reformatted file already carries af"),
  spec("2023_36280732_NAFLD_Intermountain_EUR", "existing", note = "reformatted file already carries af"),
  spec("2023_36280732_NAFLD_UKBB_EUR", "existing", note = "reformatted file already carries af"),
  oahmed("2021_34128465_PDFF_EUR"), oahmed("2021_34841290_NAFLD_EUR"), oahmed("2021_34957434_PDFF_EUR"),
  spec("2019_31311600_NAFLD_EUR", "none", note = "GCST008468 has no frequency in any on-disk source"),
  spec("2020_32298765_NAFLD_EUR", "none", note = "harmonised effect_allele_frequency is 100% NA"),
  spec("2022_36402844_PDFF_EUR", "none", note = "z_value/n only; GCST90267352 copy also lacks AF"),
  finngen("FinnGen_NAFLD", "finngen_R12_NAFLD.gz"),
  finngen("FinnGen_NASH", "finngen_R12_CHIRHEP_NAS.gz"),
  bbj("ALT"), bbj("AST"), bbj("GGT"),
  ukbb("ALT", "90019492"), ukbb("AST", "90019497"), ukbb("GGT", "90019507"),
  panukbb("AFR", "ALT"), panukbb("AFR", "AST"), panukbb("AFR", "GGT"),
  panukbb("CSA", "ALT"), panukbb("CSA", "AST"), panukbb("CSA", "GGT"),
  mvp("NAFLD", "EUR"), mvp("NAFLD", "AFR"), mvp("NAFLD", "AMR"), mvp("NAFLD", "EAS"),
  mvp("ALT", "EUR"), mvp("ALT", "AFR"), mvp("ALT", "AMR"), mvp("ALT", "EAS"),
  mvp("AST", "EUR"), mvp("AST", "AFR"), mvp("AST", "AMR"), mvp("AST", "EAS"),
  # placement=="supp" (tier 3/4).  These are outside the finemap portfolio but
  # ARE in the 50-study registry that canonical COLOC reads, so without them the
  # palindrome frequency guard in src/06_susie_coloc.R is silently inactive for
  # 15 of 50 studies.  MVP_Cirrhosis has no EAS stratum in the registry.
  mvp("Cirrhosis", "EUR"), mvp("Cirrhosis", "AFR"), mvp("Cirrhosis", "AMR"),
  mvp("ChronLiver", "EUR"), mvp("ChronLiver", "AFR"), mvp("ChronLiver", "AMR"), mvp("ChronLiver", "EAS"),
  mvp("Albumin", "EUR"), mvp("Albumin", "AFR"), mvp("Albumin", "AMR"), mvp("Albumin", "EAS"),
  mvp("Platelet", "EUR"), mvp("Platelet", "AFR"), mvp("Platelet", "AMR"), mvp("Platelet", "EAS")
), fill = TRUE)

registry <- fread(file.path(FM_ROOT, "config", "gwas_registry.tsv"))
tiers <- fread(file.path(FM_ROOT, "config", "gwas_trait_tier.tsv"))
# Build AF for every registry study that has a source spec, not just the
# placement=="main" finemap portfolio.  The finemap branch consumes only the 35
# main studies, but canonical COLOC (src/06_susie_coloc.R) reads the full
# 50-study registry and uses `af` for its palindrome frequency guard; scoping
# this table to `main` left that guard inactive for the 15 supp studies without
# anything signalling it.
universe <- merge(registry[, .(study_name, sumstats_path, ancestry)],
                  tiers[, .(study_name, tier, placement)], by = "study_name", all.x = TRUE)
main <- universe[study_name %in% SPECS$study_name]
main[, reformatted := ifelse(grepl("^/", sumstats_path), sumstats_path, file.path(FM_ROOT, sumstats_path))]

FRAG_DIR <- file.path(FM_ROOT, "config", "af_sources")

args <- commandArgs(trailingOnly = TRUE)
# --merge collects the per-study fragments into the override table.  Each study
# writes only its own fragment so a 35-task array cannot race on the shared file.
if (length(args) == 1L && identical(args[[1L]], "--merge")) {
  frags <- list.files(FRAG_DIR, pattern = "[.]tsv$", full.names = TRUE)
  assert_true(length(frags) > 0L, "No AF fragments found in %s", FRAG_DIR)
  res <- rbindlist(lapply(frags, fread), fill = TRUE)
  setorder(res, study_name)
  atomic_fwrite(res, SOURCES)
  missing <- setdiff(main$study_name, res$study_name)
  cat(sprintf("Merged %d fragments -> %s\n", nrow(res), SOURCES))
  cat(sprintf("Studies with any AF: %d/%d | median coverage among them %.1f%%\n",
              sum(res$n_with_af > 0), nrow(res), 100 * median(res[n_with_af > 0]$af_coverage)))
  if (length(missing)) cat(sprintf("MISSING (not yet built): %s\n", paste(missing, collapse = ", ")))
  print(res[, .N, by = .(mode, has_af = n_with_af > 0)])
  quit(save = "no", status = if (length(missing)) 1L else 0L)
}

want <- if (length(args)) args else main$study_name
assert_true(all(want %in% main$study_name), "Unknown study: %s", paste(setdiff(want, main$study_name), collapse = ", "))
# Scoped to the finemap portfolio deliberately: `main` is now filtered BY spec
# membership, so setdiff(main, SPECS) is vacuously empty and would assert
# nothing.  The contract worth enforcing is that no placement=="main" study is
# left without an AF source.
missing_spec <- setdiff(tiers[placement == "main", study_name], SPECS$study_name)
assert_true(length(missing_spec) == 0L,
            "No AF source specified for: %s", paste(missing_spec, collapse = ", "))
unspecced <- setdiff(universe$study_name, SPECS$study_name)
if (length(unspecced))
  cat(sprintf("NOTE: %d registry studies have no AF spec and keep af=NA: %s\n",
              length(unspecced), paste(unspecced, collapse = ", ")))

lift_hg38_to_hg19 <- function(dt, chr_col, pos_col) {
  suppressPackageStartupMessages({ library(rtracklayer); library(GenomicRanges) })
  assert_true(file.exists(CHAIN), "Missing liftover chain: %s", CHAIN)
  chain <- import.chain(CHAIN)
  gr <- GRanges(seqnames = paste0("chr", dt[[chr_col]]),
                ranges = IRanges(start = dt[[pos_col]], width = 1L))
  mcols(gr)$row <- seq_len(nrow(dt))
  lifted <- unlist(liftOver(gr, chain))
  out <- data.table(row = mcols(lifted)$row,
                    chr19 = normalize_chr(as.character(seqnames(lifted))),
                    pos19 = start(lifted))
  # A position that maps to more than one hg19 coordinate is ambiguous; drop it
  # rather than pick arbitrarily.
  out <- out[!duplicated(row) & !(row %in% out[duplicated(row), row])]
  out
}

build_one <- function(study) {
  sp <- SPECS[study_name == study][1L]
  info <- main[study_name == study][1L]
  ref <- fread(info$reformatted, na.strings = c("", "NA", "NaN"))
  core <- c("chromosome", "position", "allele1", "allele2", "beta", "se", "pval")
  assert_true(all(core %in% names(ref)), "Reformatted schema mismatch for %s", study)
  ref[, `:=`(chromosome = normalize_chr(chromosome), position = as.integer(position),
             allele1 = toupper(as.character(allele1)), allele2 = toupper(as.character(allele2)))]

  if (identical(sp$mode, "existing")) {
    assert_true("af" %in% names(ref), "%s was declared 'existing' but has no af column", study)
    out <- ref[, c(core, "af"), with = FALSE]
  } else if (identical(sp$mode, "none")) {
    out <- ref[, c(core), with = FALSE][, af := NA_real_]
  } else {
    assert_true(file.exists(sp$source_path), "AF source missing for %s: %s", study, sp$source_path)
    keep <- c(sp$src_chr, sp$src_pos, sp$src_ea, sp$src_af)
    if (!is.na(sp$src_oa)) keep <- c(keep, sp$src_oa)
    raw <- fread(sp$source_path, select = keep, na.strings = c("", "NA", "NaN"), showProgress = FALSE)
    setnames(raw, sp$src_chr, "chromosome"); setnames(raw, sp$src_pos, "position")
    setnames(raw, sp$src_ea, "ea"); setnames(raw, sp$src_af, "src_af_value")
    if (!is.na(sp$src_oa)) setnames(raw, sp$src_oa, "oa") else raw[, oa := NA_character_]
    raw[, `:=`(chromosome = normalize_chr(chromosome), position = as.integer(position),
               ea = toupper(as.character(ea)), oa = toupper(as.character(oa)),
               src_af_value = suppressWarnings(as.numeric(src_af_value)))]
    raw <- raw[!is.na(chromosome) & !is.na(position) & is.finite(src_af_value)]
    # A source can declare a frequency column that is entirely NA (the GWAS
    # Catalog harmonised files do this).  Degrade to af = NA and record it rather
    # than crashing downstream on an empty GRanges.
    if (!nrow(raw)) {
      warning(sprintf("%s: AF source yielded no usable rows; emitting af = NA", study))
      out <- ref[, c(core), with = FALSE][, af := NA_real_]
      raw <- NULL
    }

    if (!is.null(raw) && identical(sp$mode, "join_lift")) {
      lift <- lift_hg38_to_hg19(raw, "chromosome", "position")
      raw <- raw[lift$row]
      raw[, `:=`(chromosome = lift$chr19, position = lift$pos19)]
      raw <- raw[!is.na(chromosome) & !is.na(position)]
    }
    # Join on the effect allele only.  `ea` corresponds to allele1 in every
    # reformatted file (the formatters all map effect_allele -> allele1), so the
    # frequency needs no flipping.  Matching allele1 alone tolerates sources
    # whose other-allele representation differs (indel normalisation).
    if (!is.null(raw)) {
      raw <- unique(raw, by = c("chromosome", "position", "ea"))
      out <- merge(ref[, c(core), with = FALSE], raw[, .(chromosome, position, ea, src_af_value)],
                   by.x = c("chromosome", "position", "allele1"),
                   by.y = c("chromosome", "position", "ea"),
                   all.x = TRUE, sort = FALSE)
      setnames(out, "src_af_value", "af")
      setcolorder(out, c(core, "af"))
    }
  }

  out[!is.finite(af) | af < 0 | af > 1, af := NA_real_]
  # --- verification A1: the seven core columns must be unchanged -----------
  a <- copy(ref[, c(core), with = FALSE]); b <- copy(out[, c(core), with = FALSE])
  setorderv(a, core); setorderv(b, core)
  assert_true(nrow(a) == nrow(b) && all.equal(a, b, check.attributes = FALSE) == TRUE,
              "AF augmentation altered the core columns for %s -- refusing to write", study)

  ensure_dirs(OUT_DIR)
  out_path <- file.path(OUT_DIR, paste0(study, "_af.tsv"))
  atomic_fwrite(out, out_path)
  n_af <- sum(!is.na(out$af))
  data.table(
    study_name = study, tier = info$tier, ancestry = info$ancestry, mode = sp$mode,
    af_sumstats_path = out_path, source_path = sp$source_path,
    n_variants = nrow(out), n_with_af = n_af,
    af_coverage = if (nrow(out)) n_af / nrow(out) else NA_real_,
    af_na_fraction = if (nrow(out)) 1 - n_af / nrow(out) else NA_real_,
    sha256 = sha256_file(out_path), note = sp$note
  )
}

ensure_dirs(FRAG_DIR)
for (s in want) {
  r <- build_one(s)
  atomic_fwrite(r, file.path(FRAG_DIR, paste0(s, ".tsv")))
  cat(sprintf("  %-40s %-10s %9d variants  AF %.1f%%\n",
              r$study_name, r$mode, r$n_variants, 100 * r$af_coverage))
}
cat(sprintf("\nWrote %d study AF tables and fragments under %s\n", length(want), FRAG_DIR))
cat("Run with --merge once all studies are built to produce the override table.\n")
