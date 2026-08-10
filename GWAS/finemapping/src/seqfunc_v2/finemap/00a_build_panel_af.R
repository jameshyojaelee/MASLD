#!/usr/bin/env Rscript

# Build per-population allele-frequency sidecars for the 1000 Genomes LD panels.
#
# Why this exists: strand-ambiguous (palindromic) variants cannot be oriented
# from letters alone, and NO frequency file ships with any LD panel in this
# repository.  PolyFun in particular has no genotypes at all, so EUR loci use
# 1000G EUR (n = 379) as a frequency PROXY; the two EUR panels were verified to
# agree on strand, which is the only property the proxy is asked to support.
#
# The join key is built exactly as `ld_key` is in 02_run_susie_locus.R: from BIM
# columns 5 and 6.  The variant ID is NEVER parsed -- in these files the ID
# encodes chr:pos:ref:alt while columns 5/6 are PLINK's A1/A2, and the two are
# inverted relative to each other.  Getting that backwards would silently
# invert every frequency and make palindrome resolution worse than useless, so
# an independent .bed allele counter re-derives a sample of frequencies from the
# raw genotype bytes and this script refuses to write unless they agree exactly.
#
# Usage:  Rscript 00a_build_panel_af.R [pop ...]      (default: all five)
# Env:    PANEL_AF_ROOT   output root (default data/ld_ref/panel_af)
#         PLINK2_BIN      plink2 binary (must self-report "PLINK v2")
#         PANEL_AF_VERIFY_N  variants re-derived per chromosome (default 25)

source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT", unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "GWAS", "finemapping", "src", "seqfunc_v2", "finemap", "common.R"
))

POPS <- c("eur", "afr", "amr", "eas", "sas")
args <- commandArgs(trailingOnly = TRUE)
pops <- if (length(args)) args else POPS
assert_true(all(pops %in% POPS), "Unknown population requested: %s", paste(setdiff(pops, POPS), collapse = ", "))

# NOTE: /nfs/sw/plink/plink-2.0a5.13/ contains a binary named `plink` that is in
# fact PLINK v2; there is no `plink2` in that directory.  Pin it explicitly and
# assert the reported version, because PLINK 1.9 writes a different --freq
# schema and a different A1 convention -- silently mis-taking one for the other
# is exactly the inversion this script exists to rule out.
PLINK2 <- Sys.getenv("PLINK2_BIN", unset = "/nfs/sw/plink/plink-2.0a5.13/plink")
VERIFY_N <- as.integer(Sys.getenv("PANEL_AF_VERIFY_N", unset = "25"))
LD_REF <- file.path(FM_ROOT, "data", "ld_ref")
OUT_ROOT <- PANEL_AF_ROOT

assert_true(file.exists(PLINK2), "plink2 binary not found: %s", PLINK2)
plink_version <- system2(PLINK2, "--version", stdout = TRUE, stderr = TRUE)[[1L]]
assert_true(grepl("PLINK v2", plink_version, fixed = TRUE),
            "PLINK2_BIN does not report PLINK v2 (got: %s)", plink_version)
cat(sprintf("plink2: %s -- %s\n", PLINK2, plink_version))

# ---------------------------------------------------------------------------
# Independent allele counter, straight from the .bed genotype bytes.
#
# PLINK 1 .bed is variant-major: a 3-byte magic header, then ceil(N/4) bytes per
# variant, 2 bits per sample, least-significant bits first.  Codes are
#   00 = homozygous A1 (BIM column 5)   01 = missing
#   10 = heterozygous                   11 = homozygous A2 (BIM column 6)
# so freq(A1) = (2*n00 + n10) / (2 * n_nonmissing).
#
# This deliberately shares no code with plink2.  If the two disagree, the BIM
# column convention has been misread somewhere and the sidecar is untrustworthy.
# ---------------------------------------------------------------------------
bed_allele1_freq <- function(bed_path, variant_index, n_samples) {
  bytes_per_variant <- ceiling(n_samples / 4)
  con <- file(bed_path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, "raw", n = 3L)
  assert_true(identical(as.integer(magic), c(0x6cL, 0x1bL, 0x01L)),
              "%s is not a variant-major PLINK 1 .bed", bed_path)
  vapply(variant_index, function(v) {
    seek(con, where = 3 + (as.numeric(v) - 1) * bytes_per_variant, origin = "start")
    raw_bytes <- readBin(con, "raw", n = bytes_per_variant)
    codes <- integer(0)
    ints <- as.integer(raw_bytes)
    for (shift in c(0L, 2L, 4L, 6L)) codes <- c(codes, bitwAnd(bitwShiftR(ints, shift), 3L))
    # Re-interleave to sample order, then drop the padding in the final byte.
    codes <- as.vector(t(matrix(codes, nrow = length(ints))))
    codes <- codes[seq_len(n_samples)]
    n00 <- sum(codes == 0L); n10 <- sum(codes == 2L); n11 <- sum(codes == 3L)
    denom <- 2 * (n00 + n10 + n11)
    if (denom == 0) return(NA_real_)
    (2 * n00 + n10) / denom
  }, numeric(1))
}

run_chromosome <- function(pop, chr) {
  prefix <- file.path(LD_REF, paste0("1kg_", pop), sprintf("chr%d_%s", chr, pop))
  bim_path <- paste0(prefix, ".bim")
  bed_path <- paste0(prefix, ".bed")
  fam_path <- paste0(prefix, ".fam")
  for (p in c(bim_path, bed_path, fam_path)) assert_true(file.exists(p), "Missing panel file: %s", p)

  bim <- fread(bim_path, header = FALSE, select = c(1L, 2L, 4L, 5L, 6L))
  setnames(bim, c("chromosome", "variant_id", "position", "bim_a1", "bim_a2"))
  bim[, `:=`(chromosome = normalize_chr(chromosome), position = as.integer(position),
             bim_a1 = toupper(as.character(bim_a1)), bim_a2 = toupper(as.character(bim_a2)))]
  n_samples <- length(readLines(fam_path, warn = FALSE))

  tmp_prefix <- file.path(tempdir(), sprintf("panelaf_%s_chr%d_%d", pop, chr, Sys.getpid()))
  on.exit(unlink(paste0(tmp_prefix, c(".afreq", ".log"))), add = TRUE)
  status <- system2(PLINK2,
                    c("--bfile", prefix, "--freq", "--out", tmp_prefix),
                    stdout = FALSE, stderr = FALSE)
  assert_true(identical(as.integer(status), 0L), "plink2 --freq failed for %s chr%d", pop, chr)
  afreq <- fread(paste0(tmp_prefix, ".afreq"))
  # plink2 loading a PLINK 1 fileset treats BIM column 5 as ALT, so ALT_FREQS is
  # the frequency of bim_a1 -- which is the orientation 02_run_susie_locus.R
  # expects.  The .bed verification below is what actually proves it.
  assert_true(all(c("ID", "ALT_FREQS", "OBS_CT") %in% names(afreq)),
              "Unexpected plink2 --freq schema for %s chr%d: %s", pop, chr, paste(names(afreq), collapse = ","))
  assert_true(nrow(afreq) == nrow(bim), "plink2 --freq row count differs from BIM for %s chr%d", pop, chr)
  assert_true(identical(as.character(afreq$ID), as.character(bim$variant_id)),
              "plink2 --freq is not in BIM order for %s chr%d", pop, chr)

  out <- data.table(
    chromosome = bim$chromosome, position = bim$position,
    bim_a1 = bim$bim_a1, bim_a2 = bim$bim_a2,
    af_a1 = as.numeric(afreq$ALT_FREQS), n_obs_alleles = as.integer(afreq$OBS_CT),
    panel_pop = pop
  )

  # ---- anti-inversion guard --------------------------------------------
  set.seed(42L + chr)
  informative <- which(is.finite(out$af_a1) & out$af_a1 > 0 & out$af_a1 < 1)
  probe <- if (length(informative) > VERIFY_N) sort(sample(informative, VERIFY_N)) else informative
  n_checked <- 0L; max_abs_diff <- 0
  if (length(probe)) {
    ours <- bed_allele1_freq(bed_path, probe, n_samples)
    max_abs_diff <- max(abs(ours - out$af_a1[probe]), na.rm = TRUE)
    n_checked <- sum(is.finite(ours))
    assert_true(max_abs_diff < 1e-6,
                paste("Independent .bed allele count disagrees with plink2 for %s chr%d",
                      "(max |diff| = %.3e over %d variants).  The BIM A1/A2 convention",
                      "has been misread; refusing to write a possibly inverted sidecar."),
                pop, chr, max_abs_diff, n_checked)
  }

  af_dir <- file.path(OUT_ROOT, pop)
  ensure_dirs(af_dir)
  af_path <- file.path(af_dir, sprintf("chr%d.af.tsv.gz", chr))
  atomic_fwrite(out, af_path)

  data.table(
    panel_pop = pop, chromosome = chr, n_variants = nrow(out), n_samples = n_samples,
    n_verified = n_checked, max_abs_diff = max_abs_diff,
    source_bim_sha256 = sha256_file(bim_path), source_fam_n = n_samples,
    af_path = af_path, af_sha256 = sha256_file(af_path)
  )
}

chrs <- Sys.getenv("PANEL_AF_CHRS", unset = "")
chrs <- if (nzchar(chrs)) as.integer(strsplit(chrs, ",", fixed = TRUE)[[1L]]) else 1:22
assert_true(all(chrs %in% 1:22), "PANEL_AF_CHRS must be autosomes")

rows <- list()
for (pop in pops) {
  for (chr in chrs) {
    r <- run_chromosome(pop, chr)
    rows[[length(rows) + 1L]] <- r
    cat(sprintf("  %s chr%-2d  %8d variants  n=%d  verified %d/%d (max diff %.2e)\n",
                pop, chr, r$n_variants, r$n_samples, r$n_verified, VERIFY_N, r$max_abs_diff))
  }
}
manifest <- rbindlist(rows)
ensure_dirs(OUT_ROOT)
atomic_fwrite(manifest, file.path(OUT_ROOT, "panel_af_manifest.tsv"))

atomic_json(list(
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  populations = pops,
  plink2_binary = PLINK2, plink2_version = plink_version,
  allele_orientation = "af_a1 is the frequency of BIM column 5 (PLINK A1), keyed on columns 5/6; the variant ID is never parsed",
  verification = "independent .bed allele counter re-derives af_a1 for a random sample per chromosome; exact agreement required",
  n_verified_per_chromosome = VERIFY_N,
  max_abs_diff_observed = max(manifest$max_abs_diff, na.rm = TRUE),
  total_variants = sum(manifest$n_variants),
  proxy_note = "PolyFun EUR ships no genotypes; EUR loci use the 1kg_eur sidecar as a frequency proxy and are flagged panel_af_is_proxy",
  canonical_outputs_mutated = FALSE, apply_only_firewall = TRUE
), file.path(OUT_ROOT, "panel_af_contract.json"))

cat(sprintf("\nWrote %d sidecars (%d variants total); worst |plink2 - bed| = %.3e\n",
            nrow(manifest), sum(manifest$n_variants), max(manifest$max_abs_diff, na.rm = TRUE)))
