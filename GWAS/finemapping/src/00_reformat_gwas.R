#!/usr/bin/env Rscript
# 00_reformat_gwas.R
# Reformats our existing GWAS to the finemapping pipeline input format:
#   chromosome, position, allele1, allele2, beta, se, pval (hg19)
#
# Handles three formats:
#   1. GWAS Catalog harmonised (UKBB ALT/AST/GGT) — hg38, needs liftover
#   2. FinnGen R12 — hg38, needs liftover
#   3. BBJ — already has pos_hg19 column, just needs column rename
#
# Usage: Rscript 00_reformat_gwas.R
# Or:    sbatch 00_reformat_gwas.sh (for SLURM)

library(data.table)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
OUT_DIR <- file.path(FM_DIR, "data/sumstats")

CHAIN_FILE <- file.path(BASE_DIR, "data/broadaway_eqtl/hg38ToHg19.over.chain")
chain <- import.chain(CHAIN_FILE)
cat("Loaded liftover chain (hg38 -> hg19):", length(chain), "chains\n")

# ---------------------------------------------------------------------------
# Helper: liftover hg38 -> hg19
# ---------------------------------------------------------------------------
liftover_to_hg19 <- function(dt) {
  gr <- GRanges(seqnames = paste0("chr", dt$chromosome),
                ranges = IRanges(start = dt$position, width = 1))
  lifted <- liftOver(gr, chain)
  # Keep only 1:1 mappings
  keep <- lengths(lifted) == 1
  cat("  Liftover:", sum(keep), "/", length(keep), "variants mapped (1:1)\n")
  dt_out <- dt[keep, ]
  lifted_pos <- unlist(lifted[keep])
  dt_out$position <- start(lifted_pos)
  return(dt_out)
}

# ---------------------------------------------------------------------------
# Shared QC gate (review U01 / RR01 Winkler-2014). Catches corrupted/placeholder
# effect sizes — e.g. the deCODE/Intermountain NAFLD source had beta=0.018 sentinel
# for ~35% of variants and |beta| up to 102, which the prior `abs(beta)<Inf` no-op let
# through into coloc.abf. Applies: finite se>0 & beta!=0; bounded |beta| (binary log-OR
# <=5, quant <=beta_cap); sentinel detection (a single beta value repeated implausibly
# often => placeholder, dropped); optional INFO>=0.3 and MAF>=0.01 when those columns
# survive. Logs per-study rejection counts and warns if >50% dropped.
# ---------------------------------------------------------------------------
qc_filter <- function(dt, study, trait_type = "quant", beta_cap = 10,
                      sentinel_frac = 0.01, info_min = 0.3, maf_min = 0.01) {
  dt <- as.data.table(dt)
  n0 <- nrow(dt); rej <- integer(0)
  ok <- is.finite(dt$beta) & is.finite(dt$se) & dt$se > 0 & dt$beta != 0
  rej["nonfinite_se_beta"] <- sum(!ok); dt <- dt[ok]
  cap <- if (trait_type %in% c("binary", "cc")) min(beta_cap, 5) else beta_cap
  ok <- abs(dt$beta) <= cap
  rej["beta_exceeds_cap"] <- sum(!ok); dt <- dt[ok]
  if (nrow(dt) > 0) {
    top <- dt[, .N, by = beta][order(-N)][1]   # data.table grouping: fast on 30M rows
    if (!is.na(top$N) && (top$N / nrow(dt)) > sentinel_frac) {
      sv <- top$beta
      warning(sprintf("[%s] QC: beta=%g repeats in %.1f%% of rows — placeholder/sentinel; dropping.",
                      study, sv, 100 * top$N / nrow(dt)), call. = FALSE)
      ok <- dt$beta != sv; rej["sentinel_placeholder_beta"] <- sum(!ok); dt <- dt[ok]
    }
  }
  if ("info" %in% names(dt)) { ok <- is.na(dt$info) | dt$info >= info_min; rej["low_info"] <- sum(!ok); dt <- dt[ok] }
  if ("af" %in% names(dt))   { maf <- pmin(dt$af, 1 - dt$af); ok <- is.na(maf) | maf >= maf_min; rej["low_maf"] <- sum(!ok); dt <- dt[ok] }
  cat(sprintf("  [%s] QC gate: %d -> %d variants (%s)\n", study, n0, nrow(dt),
              paste(sprintf("%s=%d", names(rej), rej), collapse = ", ")))
  if (n0 > 0 && nrow(dt) < 0.5 * n0)
    warning(sprintf("[%s] QC dropped >50%% of variants — inspect the source file.", study), call. = FALSE)
  dt
}

# ---------------------------------------------------------------------------
# Format 1: GWAS Catalog harmonised (UKBB liver enzymes)
# Columns: chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, p_value
# ---------------------------------------------------------------------------
reformat_gwas_catalog <- function(input_path, output_name) {
  cat("\n=== Reformatting GWAS Catalog:", output_name, "===\n")
  dt <- fread(input_path)
  cat("  Loaded", nrow(dt), "variants\n")

  dt_fmt <- dt %>%
    select(chromosome, position = base_pair_location,
           allele1 = effect_allele, allele2 = other_allele,
           beta, se = standard_error, pval = p_value) %>%
    filter(!is.na(pval), nchar(allele1) >= 1, nchar(allele2) >= 1) %>%
    mutate(chromosome = as.integer(chromosome)) %>%
    filter(!is.na(chromosome), chromosome >= 1, chromosome <= 22)

  dt_fmt <- qc_filter(dt_fmt, output_name, trait_type = "quant")  # UKBB liver enzymes are quantitative

  dt_hg19 <- liftover_to_hg19(dt_fmt)

  out_path <- file.path(OUT_DIR, paste0(output_name, "_reformatted_hg19.tsv"))
  fwrite(dt_hg19, out_path, sep = "\t")
  cat("  Written to:", out_path, "(", nrow(dt_hg19), "variants)\n")
}

# ---------------------------------------------------------------------------
# Format 2: FinnGen R12
# Columns: #chrom, pos, ref, alt, rsids, nearest_genes, pval, mlogp, beta, sebeta, af_alt
# ---------------------------------------------------------------------------
reformat_finngen <- function(input_path, output_name) {
  cat("\n=== Reformatting FinnGen:", output_name, "===\n")
  dt <- fread(input_path)
  cat("  Loaded", nrow(dt), "variants\n")

  # FinnGen uses ref/alt where alt is effect allele
  # fread strips '#' from #chrom -> chrom
  chrom_col <- intersect(c("#chrom", "chrom", "X.chrom"), colnames(dt))[1]
  if (is.na(chrom_col)) stop("Cannot find chrom column in FinnGen file")
  setnames(dt, chrom_col, "chromosome")
  dt_fmt <- dt %>%
    select(chromosome, position = pos,
           allele1 = alt, allele2 = ref,  # alt = effect allele in FinnGen
           beta, se = sebeta, pval, af = af_alt) %>%   # retain AF for QC/concordance
    filter(!is.na(pval), nchar(allele1) >= 1, nchar(allele2) >= 1) %>%
    mutate(chromosome = as.integer(chromosome)) %>%
    filter(!is.na(chromosome), chromosome >= 1, chromosome <= 22)

  dt_fmt <- qc_filter(dt_fmt, output_name, trait_type = "binary")  # FinnGen NAFLD/NASH/HCC are case-control

  dt_hg19 <- liftover_to_hg19(dt_fmt)

  out_path <- file.path(OUT_DIR, paste0(output_name, "_reformatted_hg19.tsv"))
  fwrite(dt_hg19, out_path, sep = "\t")
  cat("  Written to:", out_path, "(", nrow(dt_hg19), "variants)\n")
}

# ---------------------------------------------------------------------------
# Format 3: BBJ (already has pos_hg19)
# Columns: chromosome, base_pair_location, pos_hg19, effect_allele, other_allele, beta, standard_error, p_value, EAF
# ---------------------------------------------------------------------------
reformat_bbj <- function(input_path, output_name) {
  cat("\n=== Reformatting BBJ:", output_name, "===\n")
  dt <- fread(input_path)
  cat("  Loaded", nrow(dt), "variants\n")

  dt_fmt <- dt %>%
    select(chromosome, position = pos_hg19,
           allele1 = effect_allele, allele2 = other_allele,
           beta, se = standard_error, pval = p_value, af = EAF) %>%   # retain AF for QC/concordance
    filter(!is.na(pval), !is.na(position), nchar(allele1) >= 1, nchar(allele2) >= 1) %>%
    mutate(chromosome = as.integer(chromosome)) %>%
    filter(!is.na(chromosome), chromosome >= 1, chromosome <= 22)

  dt_fmt <- qc_filter(dt_fmt, output_name, trait_type = "quant")  # BBJ liver enzymes are quantitative

  out_path <- file.path(OUT_DIR, paste0(output_name, "_reformatted_hg19.tsv"))
  fwrite(dt_fmt, out_path, sep = "\t")
  cat("  Written to:", out_path, "(", nrow(dt_fmt), "variants)\n")
}

# ---------------------------------------------------------------------------
# Format 4: deCODE / Sveinbjornsson 2022 NAFLD (PMID 36280732; deCode/Intermountain/UKBB cohorts).
# REPLACES the corrupted *_preprocessed.tsv (0.018-sentinel for ~35% of variants, |beta| up to ~102).
# Authoritative files: deCODE data-use form (decode.com/summarydata, "Multiomics study of NAFLD").
# CONFIRMED format (nafld2022readme.txt, case-control NAFL split by cohort):
#   Pval | Effect = ODDS RATIO (NO beta/SE column!) | marker | rsName | MAF_PC = minor-allele-freq in
#   PERCENT | Info | Chrom | Pos (build hg38) | Amin = minor/EFFECT allele | Amaj = major/OTHER allele.
# => beta = log(Effect); SE is back-computed from beta+pval; af = MAF_PC/100; hg38 -> liftover hg19.
# (The bundle's FinnGen file also has Standard_Error + Beta(log-OR), auto-detected if ever used; but we
#  use our own FinnGen R12 entry, so NAFL_FINNGEN is NOT downloaded.) effect_scale='or' is the default.
# ---------------------------------------------------------------------------
pick_col <- function(nms, patterns) {
  for (p in patterns) {
    hit <- grep(p, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}
reformat_decode <- function(input_path, output_name, trait_type = "binary", build = "hg38", effect_scale = "or") {
  cat("\n=== Reformatting deCODE/Sveinbjornsson:", output_name, "===\n")
  if (!file.exists(input_path)) {
    cat("  SKIP (staged): file not downloaded yet:", input_path, "\n")
    return(invisible(NULL))
  }
  dt <- fread(input_path); nms <- names(dt)
  cat("  Loaded", nrow(dt), "variants. Columns:", paste(nms, collapse = ", "), "\n")
  c_chr  <- pick_col(nms, c("^chrom$", "^chr$", "chromosome", "^#chrom$"))
  c_pos  <- pick_col(nms, c("^pos$", "^bp$", "position", "base_pair", "pos_b3[78]", "pos_hg3[78]"))
  c_ea   <- pick_col(nms, c("^amin$", "^effectall", "effect_?allele", "allele_minor", "^minor", "^alt$", "^a1$"))
  c_oa   <- pick_col(nms, c("^amaj$", "^canonref", "other_?allele", "allele_major", "^major", "^ref$", "^a2$"))
  c_beta <- pick_col(nms, c("^beta$", "log_?or", "logor"))                       # genuine log-OR (FinnGen file)
  c_eff  <- pick_col(nms, c("^effect$", "^or$", "odds_?ratio", "^oddsratio"))    # ODDS RATIO (deCODE case-control)
  c_se   <- pick_col(nms, c("standard_?error", "^se$", "stderr", "sebeta", "^se_"))
  c_pval <- pick_col(nms, c("^pval", "p_value", "^p$", "pvalue"))
  c_af   <- pick_col(nms, c("maf_?pc", "freq_?pc", "^eaf$", "effect_allele_freq", "^af$", "^maf$", "impmaf", "allele_freq", "freq"))
  c_info <- pick_col(nms, c("^info$", "imp_?info", "^r2$", "rsq", "imputation"))
  stopifnot("missing chrom" = !is.na(c_chr), "missing pos" = !is.na(c_pos),
            "missing effect allele" = !is.na(c_ea), "missing other allele" = !is.na(c_oa),
            "missing pval" = !is.na(c_pval),
            "missing beta AND effect/OR" = !(is.na(c_beta) & is.na(c_eff)))
  # beta: prefer a genuine log-OR column; else convert the OR (Effect) to log-OR.
  if (!is.na(c_beta)) {
    beta_vec <- as.numeric(dt[[c_beta]]); beta_src <- c_beta
  } else {
    eff <- as.numeric(dt[[c_eff]])
    if (tolower(effect_scale) == "or") { beta_vec <- log(eff); beta_src <- paste0("log(", c_eff, ")") }
    else { beta_vec <- eff; beta_src <- c_eff }
  }
  cat(sprintf("  Mapped: chr=%s pos=%s ea=%s oa=%s beta=%s se=%s pval=%s af=%s info=%s\n",
              c_chr, c_pos, c_ea, c_oa, beta_src, c_se, c_pval, c_af, c_info))
  out <- data.table(
    chromosome = dt[[c_chr]], position = dt[[c_pos]],
    allele1 = toupper(as.character(dt[[c_ea]])), allele2 = toupper(as.character(dt[[c_oa]])),
    beta = beta_vec, pval = as.numeric(dt[[c_pval]]))
  if (!is.na(c_se)) {
    out[, se := as.numeric(dt[[c_se]])]
  } else {
    warning(sprintf("[%s] no SE column — back-computing se=|beta|/qnorm(1-p/2) from beta+pval.",
                    output_name), call. = FALSE)
    z <- qnorm(1 - pmin(pmax(out$pval, 1e-300), 1) / 2)
    out[, se := abs(beta) / z]
  }
  if (!is.na(c_af)) {
    afv <- as.numeric(dt[[c_af]])
    if (grepl("pc|percent", c_af, ignore.case = TRUE) || max(afv, na.rm = TRUE) > 1.5) {
      afv <- afv / 100; cat("  AF/MAF appears to be in percent -> divided by 100\n")
    }
    out[, af := afv]   # deCODE MAF_PC is minor-allele freq = effect-allele(Amin) freq
  }
  if (!is.na(c_info)) out[, info := as.numeric(dt[[c_info]])]
  out[, chromosome := as.integer(sub("^chr", "", as.character(chromosome), ignore.case = TRUE))]
  out <- out[!is.na(chromosome) & chromosome >= 1 & chromosome <= 22 &
             nchar(allele1) >= 1 & nchar(allele2) >= 1 & !is.na(pval)]
  out <- qc_filter(out, output_name, trait_type = trait_type)  # sentinel detection guards residual placeholders
  if (tolower(build) == "hg38") out <- liftover_to_hg19(out)
  out_path <- file.path(OUT_DIR, paste0(output_name, "_reformatted_hg19.tsv"))
  fwrite(out, out_path, sep = "\t")
  cat("  Written to:", out_path, "(", nrow(out), "variants)\n")
}

# ---------------------------------------------------------------------------
# Run all reformatting
# ---------------------------------------------------------------------------
cat("============================================================\n")
cat("Reformatting GWAS summary statistics for finemapping pipeline\n")
cat("============================================================\n")

GWAS_DIR <- file.path(BASE_DIR, "GWAS/MR_Data")

# Group selector: REFORMAT_SET=all (default) | comma list of {gwas_catalog,finngen,bbj,decode}.
# Lets a targeted re-run (e.g. REFORMAT_SET=decode) reformat ONLY the new deCODE NAFL files
# without re-lifting + overwriting the UKBB/FinnGen/BBJ outputs that downstream COLOC consumes.
do_set  <- tolower(Sys.getenv("REFORMAT_SET", "all"))
run_grp <- function(g) do_set == "all" || g %in% trimws(strsplit(do_set, ",")[[1]])
cat("REFORMAT_SET =", do_set, "\n")

if (run_grp("gwas_catalog")) {  # UKBB liver enzymes (GWAS Catalog format, hg38)
  reformat_gwas_catalog(file.path(GWAS_DIR, "GCST90019492_UKBB_ALT_harmonised.tsv.gz"), "UKBB_ALT")
  reformat_gwas_catalog(file.path(GWAS_DIR, "GCST90019497_UKBB_AST_harmonised.tsv.gz"), "UKBB_AST")
  reformat_gwas_catalog(file.path(GWAS_DIR, "GCST90019507_UKBB_GGT_harmonised.tsv.gz"), "UKBB_GGT")
}
if (run_grp("finngen")) {  # FinnGen R12 (hg38)
  reformat_finngen(file.path(GWAS_DIR, "FinnGen/finngen_R12_NAFLD.gz"), "FinnGen_NAFLD")
  reformat_finngen(file.path(GWAS_DIR, "FinnGen/finngen_R12_CHIRHEP_NAS.gz"), "FinnGen_NASH")
  reformat_finngen(file.path(GWAS_DIR, "FinnGen/finngen_R12_C3_HEPATOCELLU_CARC_EXALLC.gz"), "FinnGen_HCC")
}
if (run_grp("bbj")) {  # BBJ (has pos_hg19, EAS ancestry)
  reformat_bbj(file.path(GWAS_DIR, "BBJ/BBJ_ALT_harmonised_hg38.tsv.gz"), "BBJ_ALT")
  reformat_bbj(file.path(GWAS_DIR, "BBJ/BBJ_AST_harmonised_hg38.tsv.gz"), "BBJ_AST")
  reformat_bbj(file.path(GWAS_DIR, "BBJ/BBJ_GGT_harmonised_hg38.tsv.gz"), "BBJ_GGT")
}
if (run_grp("decode")) {
  # deCODE / Sveinbjornsson 2022 NAFL (PMID 36280732) — CONFIRMED filenames (decode.com/summarydata
  # "Multiomics study of NAFLD" bundle). Effect=ODDS RATIO, no SE, Amin/Amaj, MAF_PC %, hg38 (handled by
  # reformat_decode). (NAFL_FINNGEN is NOT used — we have our own FinnGen R12 NAFLD entry.)
  DECODE_DIR <- file.path(GWAS_DIR, "Sveinbjornsson2022")
  reformat_decode(file.path(DECODE_DIR, "NAFL_deCODE_sumstat.txt"),        "2023_36280732_NAFLD_deCode_EUR",        "binary", "hg38")
  reformat_decode(file.path(DECODE_DIR, "NAFL_INTERMOUNTAIN_sumstat.txt"), "2023_36280732_NAFLD_Intermountain_EUR", "binary", "hg38")
  reformat_decode(file.path(DECODE_DIR, "NAFL_UKBB_sumstat.txt"),          "2023_36280732_NAFLD_UKBB_EUR",          "binary", "hg38")
}

cat("\n============================================================\n")
cat("All reformatting complete\n")
cat("============================================================\n")
