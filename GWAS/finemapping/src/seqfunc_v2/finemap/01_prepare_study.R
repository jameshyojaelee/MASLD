#!/usr/bin/env Rscript

source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT", unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "GWAS", "finemapping", "src", "seqfunc_v2", "finemap", "common.R"
))

args <- commandArgs(trailingOnly = TRUE)
assert_true(length(args) == 1L, "Usage: 01_prepare_study.R <study_index>")
study_index <- suppressWarnings(as.integer(args[[1L]]))
assert_true(!is.na(study_index), "study_index must be an integer")

manifest_path <- file.path(RUN_ROOT, "config", "study_manifest.tsv")
assert_true(file.exists(manifest_path), "Frozen study manifest is missing: %s", manifest_path)
manifest <- fread(manifest_path)
assert_true(study_index >= 1L && study_index <= nrow(manifest), "study_index out of range: %d", study_index)
study <- manifest[study_index]
assert_true(study$study_index == study_index, "Manifest index mismatch for task %d", study_index)
assert_true(identical(sha256_file(study$sumstats_abs), study$sumstats_sha256),
            "Summary-statistic hash changed for %s", study$study_name)

cat(sprintf("Preparing %s (%s, %s)\n", study$study_name, study$ancestry, study$trait_type))
ss <- fread(study$sumstats_abs, na.strings = c("", "NA", "NaN", "Inf", "-Inf"), showProgress = TRUE)
required <- c("chromosome", "position", "allele1", "allele2", "beta", "se", "pval")
assert_true(all(required %in% names(ss)), "Summary-statistic schema mismatch for %s", study$study_name)

# `af` is carried through as a mandatory column that MAY be NA: three studies in
# the portfolio have no allele frequency in any on-disk source, and those rely on
# the study strand certificate instead.  Sources are routed here by
# config/gwas_af_sources.tsv (built by 00b_build_study_af.R).
if (!"af" %in% names(ss)) ss[, af := NA_real_] else ss[, af := suppressWarnings(as.numeric(af))]
ss[!is.finite(af) | af < 0 | af > 1, af := NA_real_]

ss[, chromosome := normalize_chr(chromosome)]
ss[, position := suppressWarnings(as.integer(position))]
ss[, `:=`(
  allele1 = toupper(as.character(allele1)), allele2 = toupper(as.character(allele2)),
  beta = suppressWarnings(as.numeric(beta)), se = suppressWarnings(as.numeric(se)),
  pval = suppressWarnings(as.numeric(pval))
)]
ss <- ss[
  chromosome %between% c(1L, 22L) & !is.na(position) & position > 0L &
    is.finite(beta) & is.finite(se) & se > 0 & is.finite(pval) & pval >= 0 & pval <= 1 &
    nzchar(allele1) & nzchar(allele2)
]
ss[, variant_key := paste(chromosome, position, allele1, allele2, sep = ":")]
setorder(ss, variant_key, pval)
ss <- ss[!duplicated(variant_key)]

signals <- ss[pval < GWS_P, .(
  chromosome, position, allele1, allele2, beta, se, pval, variant_key, af
)]
setorder(signals, chromosome, position, pval)

fragment_dir <- file.path(RUN_ROOT, "work", "prep_fragments")
summary_dir <- file.path(RUN_ROOT, "work", "summary_stats", study$study_name)
ensure_dirs(fragment_dir, summary_dir)
study_status_path <- file.path(fragment_dir, sprintf("study_%02d_status.tsv", study_index))

if (nrow(signals) == 0L) {
  status <- data.table(
    study_index = study_index, study_name = study$study_name, prep_status = "no_signal",
    n_input_variants = nrow(ss), n_gws_variants = 0L, n_loci = 0L, message = "No p<5e-8 autosomal variant"
  )
  immutable_fwrite(status, study_status_path)
  immutable_fwrite(data.table(
    study_index = integer(), study_name = character(), locus_id = character()
  ), file.path(fragment_dir, sprintf("study_%02d_loci.tsv", study_index)))
  quit(save = "no", status = 0L)
}

blocks_path <- file.path(RUN_ROOT, "config", "ld_blocks", paste0(study$ancestry, ".tsv"))
blocks <- fread(blocks_path)
blocks <- blocks[, .(
  chromosome = as.integer(chr), block_start = as.integer(start), block_stop = as.integer(stop),
  ancestry, ld_panel_id, ld_panel_n, block_prefix, bim_exists, ld_exists
)]
study_ld_panel <- as.character(blocks$ld_panel_id[[1L]])
assert_true(uniqueN(blocks$ld_panel_id) == 1L,
            "%s maps to more than one LD panel", study$study_name)
# Strand certificate: derived from NON-palindromic SNVs against the panel AF
# sidecar, so it is independent of the palindrome resolution it later supports.
cert <- compute_strand_certificate(ss, study_ld_panel)
cat(sprintf("  strand certificate: %s (same=%d opp=%d, opp_fraction=%.2e)\n",
            cert$strand_certificate, cert$n_same, cert$n_opp, cert$opp_fraction))
immutable_fwrite(data.table(
  study_index = study_index, study_name = study$study_name,
  ld_panel_id = study_ld_panel, strand_certificate = cert$strand_certificate,
  n_same = cert$n_same, n_opp = cert$n_opp, opp_fraction = cert$opp_fraction,
  n_informative = cert$n_informative,
  af_present = sum(!is.na(ss$af)) > 0L,
  af_na_fraction = mean(is.na(ss$af))
), file.path(fragment_dir, sprintf("study_%02d_strand.tsv", study_index)))

blocks[, `:=`(query_start = block_start, query_stop = block_stop - 1L)]
signal_ranges <- signals[, .(chromosome, query_start = position, query_stop = position, signal_row = .I)]
setkey(blocks, chromosome, query_start, query_stop)
setkey(signal_ranges, chromosome, query_start, query_stop)
assigned <- foverlaps(signal_ranges, blocks, type = "within", nomatch = NA_integer_)
if (anyNA(assigned$block_start)) {
  missing_rows <- signals[assigned[is.na(block_start), signal_row]]
  atomic_fwrite(missing_rows, file.path(fragment_dir, sprintf("study_%02d_unassigned_signals.tsv", study_index)))
  stopf("%s has %d GWS variants outside the frozen LD blocks", study$study_name, nrow(missing_rows))
}
assert_true(!anyDuplicated(assigned$signal_row), "A GWS variant mapped to multiple LD blocks for %s", study$study_name)

assigned <- cbind(assigned, signals[assigned$signal_row, .(position, allele1, allele2, pval, variant_key)])
locus_defs <- assigned[, .(
  n_gws_variants = .N,
  top_position = position[which.min(pval)],
  top_allele1 = allele1[which.min(pval)],
  top_allele2 = allele2[which.min(pval)],
  top_p = min(pval),
  gws_variant_keys = paste(sort(unique(variant_key)), collapse = ";")
), by = .(chromosome, block_start, block_stop, ancestry, ld_panel_id, ld_panel_n, block_prefix, bim_exists, ld_exists)]
setorder(locus_defs, chromosome, block_start, block_stop)

count_lines <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  out <- system2("wc", c("-l", path), stdout = TRUE, stderr = TRUE)
  suppressWarnings(as.integer(strsplit(trimws(out[[1L]]), "[[:space:]]+")[[1L]][[1L]]))
}

rows <- vector("list", nrow(locus_defs))
for (i in seq_len(nrow(locus_defs))) {
  loc <- locus_defs[i]
  locus_id <- sprintf("%s__chr%d__%d_%d", study$study_name, loc$chromosome, loc$block_start, loc$block_stop)
  locus_ss <- ss[
    chromosome == loc$chromosome & position >= loc$block_start & position < loc$block_stop,
    .(chromosome, position, allele1, allele2, beta, se, pval, variant_key, af)
  ]
  setorder(locus_ss, position, allele1, allele2)
  ss_path <- file.path(summary_dir, paste0(locus_id, ".tsv.gz"))
  if (!file.exists(ss_path)) atomic_fwrite(locus_ss, ss_path)
  ss_hash <- sha256_file(ss_path)
  bim_path <- paste0(loc$block_prefix, ".bim")
  ld_path <- paste0(loc$block_prefix, ".ld")
  bim_exists <- file.exists(bim_path)
  ld_exists <- file.exists(ld_path)
  bim_data <- if (bim_exists) fread(bim_path, header = FALSE, select = c(1L, 4L, 5L, 6L)) else data.table()
  if (nrow(bim_data)) setnames(bim_data, c("ld_chr", "ld_position", "ld_allele1", "ld_allele2"))
  loc_signals <- signals[chromosome == loc$chromosome & position >= loc$block_start & position < loc$block_stop]
  represented <- if (nrow(bim_data)) loc_signals$position %in% bim_data$ld_position else rep(FALSE, nrow(loc_signals))
  unrepresented_keys <- loc_signals$variant_key[!represented]
  ld_min_position <- if (nrow(bim_data)) min(bim_data$ld_position, na.rm = TRUE) else NA_integer_
  ld_max_position <- if (nrow(bim_data)) max(bim_data$ld_position, na.rm = TRUE) else NA_integer_
  nominal_width <- as.integer(loc$block_stop - loc$block_start)
  covered_width <- if (nrow(bim_data)) as.integer(ld_max_position - ld_min_position + 1L) else NA_integer_
  coverage_fraction <- if (nrow(bim_data)) covered_width / nominal_width else NA_real_
  hash_stem <- file.path(
    RUN_ROOT, "work", "ld_hashes", loc$ld_panel_id,
    sprintf("chr%d_%d_%d", loc$chromosome, loc$block_start, loc$block_stop)
  )
  bim_sha <- if (bim_exists) cached_sha256(bim_path, paste0(hash_stem, ".bim.tsv")) else NA_character_
  ld_sha <- if (ld_exists) cached_sha256(ld_path, paste0(hash_stem, ".ld.tsv")) else NA_character_
  rows[[i]] <- data.table(
    study_index = study_index, study_name = study$study_name, trait = study$trait,
    tier = study$tier, ancestry = study$ancestry, trait_type = study$trait_type,
    N_tot = study$N_tot, N_cases = study$N_cases,
    locus_id = locus_id, chromosome = loc$chromosome,
    block_start = loc$block_start, block_stop = loc$block_stop,
    n_gws_variants = loc$n_gws_variants, top_position = loc$top_position,
    top_allele1 = loc$top_allele1, top_allele2 = loc$top_allele2,
    top_p = loc$top_p, gws_variant_keys = loc$gws_variant_keys,
    n_sumstats_variants = nrow(locus_ss), summary_stats_path = ss_path,
    summary_stats_sha256 = ss_hash, ld_panel_id = loc$ld_panel_id,
    ld_panel_n = loc$ld_panel_n, block_prefix = loc$block_prefix,
    bim_path = bim_path, ld_path = ld_path,
    bim_exists = bim_exists, ld_exists = ld_exists,
    bim_bytes = if (bim_exists) file.info(bim_path)$size else NA_real_,
    ld_bytes = if (ld_exists) file.info(ld_path)$size else NA_real_,
    bim_sha256 = bim_sha, ld_sha256 = ld_sha,
    n_ld_variants = if (bim_exists) nrow(bim_data) else NA_integer_,
    ld_min_position = ld_min_position, ld_max_position = ld_max_position,
    nominal_block_width = nominal_width, panel_covered_width = covered_width,
    panel_coverage_fraction = coverage_fraction,
    panel_truncated = if (nrow(bim_data)) coverage_fraction < 0.90 else NA,
    n_gws_positionally_represented = sum(represented),
    n_gws_positionally_unrepresented = sum(!represented),
    unrepresented_gws_variant_keys = paste(unrepresented_keys, collapse = ";")
  )
}
loci <- rbindlist(rows, fill = TRUE)
loci[, memory_tier := vapply(n_ld_variants, memory_tier, character(1))]
immutable_fwrite(loci, file.path(fragment_dir, sprintf("study_%02d_loci.tsv", study_index)))

status <- data.table(
  study_index = study_index, study_name = study$study_name, prep_status = "prepared",
  n_input_variants = nrow(ss), n_gws_variants = nrow(signals), n_loci = nrow(loci),
  message = sprintf("%d loci; %d missing BIM; %d missing LD", nrow(loci), sum(!loci$bim_exists), sum(!loci$ld_exists))
)
immutable_fwrite(status, study_status_path)
cat(sprintf("Prepared %s: %d GWS variants in %d unique LD blocks\n", study$study_name, nrow(signals), nrow(loci)))
