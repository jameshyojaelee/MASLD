#!/usr/bin/env Rscript
# 00c_compare_overlapping_gwas.R
# ---------------------------------------------------------------------------
# Characterize OVERLAP and DIFFERENCES between the Sveinbjornsson-2022 NAFLD bundle
# (decode.com data-use form) and the GWAS we already use. Several of our sumstats share
# cohorts/samples with the bundle (Ghodsian is an EHR meta-analysis that INCLUDES deCODE +
# UKBB; our FinnGen R12 NAFLD is a newer release of the bundle's FinnGen NAFL; our PDFF is
# the same UKBB MRI cohort). For COLOC the question is two-fold:
#   (1) Do they measure the SAME biology?  -> effect-size correlation at GW-sig SNPs.
#   (2) Do they SHARE SAMPLES (so not independent evidence)? -> correlation of Z at NULL
#       SNPs (p>0.5). Under no sample overlap, null Z are ~uncorrelated; shared controls/
#       cases drive a positive null-Z correlation (the bivariate-LDSC intercept logic).
# Output: results/gwas_overlap/overlap_comparison.csv + a console report.
# Run on a compute node (16M-row files): sbatch wrapper 00c_compare.sbatch (bigmem/cpu, ~32G).
# ---------------------------------------------------------------------------
suppressMessages({ library(data.table) })
setDTthreads(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GW   <- file.path(BASE, "GWAS/MR_Data")           # data dir (rename to GWAS/data pending)
SV   <- file.path(GW, "Sveinbjornsson2022")
OUT  <- file.path(BASE, "GWAS/finemapping/results/gwas_overlap")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

COMPL <- c(A="T", T="A", C="G", G="C")
flip_strand <- function(a) { out <- COMPL[a]; out[is.na(out)] <- a[is.na(out)]; unname(out) }

# Standardized reader: returns DT(chr int, pos int, ea, oa, beta, se, pval, af) on autosomes.
# `spec` names the columns + an effect_scale ('or'|'beta'|'sd') and whether ea/oa need deriving.
read_sumstats <- function(path, spec, label, max_rows = Inf) {
  if (!file.exists(path)) { cat(sprintf("  [%s] MISSING: %s\n", label, path)); return(NULL) }
  sel <- unlist(spec$cols, use.names = TRUE)
  hdr <- names(fread(path, nrows = 0L))                # tolerate optional cols absent in some files
  present <- sel[sel %in% hdr]
  missing <- names(sel)[!(sel %in% hdr)]
  req <- c("pval", "beta", "chr", "pos", "ea", "oa")
  if (length(intersect(req, missing))) {
    cat(sprintf("  [%s] MISSING required column(s): %s -> skipping\n", label,
                paste(intersect(req, missing), collapse = ", "))); return(NULL) }
  if (length(missing)) cat(sprintf("  [%s] note: optional column(s) absent: %s\n", label, paste(missing, collapse = ", ")))
  dt <- fread(path, select = unname(present), nThread = getDTthreads(), nrows = max_rows)
  setnames(dt, unname(present), names(present))
  dt[, chr := suppressWarnings(as.integer(sub("^chr", "", as.character(chr), ignore.case = TRUE)))]
  dt[, pos := suppressWarnings(as.integer(pos))]
  dt[, pval := suppressWarnings(as.numeric(pval))]
  # effect -> beta (log-OR / SD / already-beta)
  if (spec$effect_scale == "or")      dt[, beta := log(suppressWarnings(as.numeric(beta)))]
  else                                dt[, beta := suppressWarnings(as.numeric(beta))]
  # other allele: derive from canonRef/canonAlt vs effect allele when needed
  if (isTRUE(spec$derive_oa)) dt[, oa := fifelse(toupper(ea) == toupper(canonAlt), toupper(canonRef), toupper(canonAlt))]
  dt[, `:=`(ea = toupper(as.character(ea)), oa = toupper(as.character(oa)))]
  if (!"se" %in% names(dt)) dt[, se := NA_real_]
  if (!"af" %in% names(dt)) dt[, af := NA_real_]
  dt <- dt[!is.na(chr) & chr >= 1 & chr <= 22 & !is.na(pos) & !is.na(beta) & is.finite(beta) &
             nchar(ea) >= 1 & nchar(oa) >= 1 & ea != oa]
  dt[, z := beta / se]
  cat(sprintf("  [%s] %d autosomal variants\n", label, nrow(dt)))
  dt[, .(chr, pos, ea, oa, beta, se, z, pval, af)]
}

# Harmonize two standardized DTs on chr:pos and align B's effect allele to A's.
# Returns matched DT with betaA/betaB/zA/zB (B aligned). Drops ambiguous palindromes by AF
# only when AF present; otherwise keeps identity-matched orientation.
harmonize_pair <- function(a, b) {
  setkey(a, chr, pos); setkey(b, chr, pos)
  m <- merge(a, b, by = c("chr", "pos"), suffixes = c("A", "B"), allow.cartesian = FALSE)
  if (!nrow(m)) return(m)
  same <- m$eaA == m$eaB & m$oaA == m$oaB
  flip <- m$eaA == m$oaB & m$oaA == m$eaB
  # strand-flip (palindrome-safe-ish): try complement match
  cflip <- (flip_strand(m$eaB) == m$eaA & flip_strand(m$oaB) == m$oaA)
  cswap <- (flip_strand(m$eaB) == m$oaA & flip_strand(m$oaB) == m$eaA)
  sign  <- rep(NA_real_, nrow(m))
  sign[same | cflip]  <-  1
  sign[flip | cswap]  <- -1
  m <- m[!is.na(sign)]; sign <- sign[!is.na(sign)]
  m[, `:=`(betaB = betaB * sign, zB = zB * sign)]
  m
}

report_pair <- function(a_spec, b_spec, name) {
  cat(sprintf("\n=== %s ===\n", name))
  a <- read_sumstats(a_spec$path, a_spec, a_spec$label)
  b <- read_sumstats(b_spec$path, b_spec, b_spec$label)
  if (is.null(a) || is.null(b)) { cat("  skipped (missing input)\n"); return(NULL) }
  m <- harmonize_pair(a, b)
  if (!nrow(m)) { cat("  no shared variants after harmonization\n"); return(NULL) }
  gw <- 5e-8
  sigEither <- m[pvalA < gw | pvalB < gw]
  nullBoth  <- m[pvalA > 0.5 & pvalB > 0.5]
  # null-Z correlation = sample-overlap proxy. ONLY valid when BOTH sides have a REAL SE column;
  # a back-computed SE makes z = sign(beta)*qnorm(1-p/2), a deterministic fn of p (not a true Z), so
  # the proxy is meaningless there (deCODE/Intermountain/UKBB NAFL have no SE -> reported N/A).
  both_real_se <- ("se" %in% names(a_spec$cols)) && ("se" %in% names(b_spec$cols))
  zr <- if (both_real_se && sum(is.finite(nullBoth$zA) & is.finite(nullBoth$zB)) > 1000)
          cor(nullBoth$zA, nullBoth$zB, use = "complete.obs") else NA_real_
  br_sig <- if (nrow(sigEither) > 10) cor(sigEither$betaA, sigEither$betaB, use = "complete.obs") else NA_real_
  conc   <- if (nrow(sigEither) > 0) mean(sign(sigEither$betaA) == sign(sigEither$betaB), na.rm = TRUE) else NA_real_
  res <- data.table(
    comparison = name, A = a_spec$label, B = b_spec$label,
    n_A = nrow(a), n_B = nrow(b), n_shared = nrow(m),
    n_gwsig_A = a[pval < gw, .N], n_gwsig_B = b[pval < gw, .N],
    n_sig_either = nrow(sigEither),
    r_beta_gwsig = round(br_sig, 4), dir_concordance_gwsig = round(conc, 4),
    r_z_null = round(zr, 4), n_null = nrow(nullBoth))
  print(res[, .(n_A, n_B, n_shared, n_gwsig_A, n_gwsig_B, r_beta_gwsig,
                dir_concordance_gwsig, r_z_null)])
  cat(sprintf("  INTERPRETATION: effect-r(GWsig)=%.3f conc=%.2f | null-Z-r=%s (sample-overlap proxy; ~0=independent)%s\n",
              br_sig, conc, ifelse(is.na(zr), "NA", sprintf("%.3f", zr)),
              if (!both_real_se) " [N/A: a side has back-computed SE]" else ""))
  res
}

# ---- column specs ---------------------------------------------------------
# Sveinbjornsson case-control (NAFL/Cirrhosis/HCC _{deCODE,INTERMOUNTAIN,UKBB}): Effect=OR, Amin=ea, Amaj=oa, no SE
sv_cc <- function(f, lab) list(path=file.path(SV,f), label=lab, effect_scale="or", derive_oa=FALSE,
  cols=c(pval="Pval", beta="Effect", chr="Chrom", pos="Pos", ea="Amin", oa="Amaj", af="MAF_PC"))
# Sveinbjornsson FinnGen (NAFL/Cirrhosis _FINNGEN): Beta=log-OR, effectAll=ea, derive oa from canonRef/canonAlt, SE present
sv_fg <- function(f, lab) list(path=file.path(SV,f), label=lab, effect_scale="beta", derive_oa=TRUE,
  cols=c(pval="Pval", beta="Beta", se="Standard_Error", chr="Chrom", pos="Pos", ea="effectAll",
         canonRef="canonRef", canonAlt="canonAlt", af="Freq_PC"))
# Sveinbjornsson PDFF UKBB: Effect in SD (quantitative), Amin=ea, Amaj=oa, no SE
sv_pdff <- function(f, lab) list(path=file.path(SV,f), label=lab, effect_scale="sd", derive_oa=FALSE,
  cols=c(pval="Pval", beta="Effect", chr="Chrom", pos="Pos", ea="Amin", oa="Amaj", af="MAF_PC"))
# our FinnGen R12 (tab, #chrom/pos/ref/alt, alt=effect allele): beta present, sebeta present
our_fg <- function(f, lab) list(path=file.path(GW,"FinnGen",f), label=lab, effect_scale="beta", derive_oa=FALSE,
  cols=c(pval="pval", beta="beta", se="sebeta", chr="#chrom", pos="pos", ea="alt", oa="ref", af="af_alt"))
# our GWAS-Catalog harmonised (Ghodsian / UKBB enzymes / Pazoki PDFF): beta present, SE present
our_gc <- function(f, lab) list(path=file.path(GW,f), label=lab, effect_scale="beta", derive_oa=FALSE,
  cols=c(pval="p_value", beta="beta", se="standard_error", chr="chromosome", pos="base_pair_location",
         ea="effect_allele", oa="other_allele", af="effect_allele_frequency"))

pairs <- list(
  list(report_pair, sv_fg("NAFL_FINNGEN_sumstats.txt","Sveinb_FinnGen_NAFL"),
                    our_fg("finngen_R12_NAFLD.gz","our_FinnGenR12_NAFLD"),
       "FinnGen NAFL: Sveinbjornsson-2022 release vs our R12"),
  list(report_pair, sv_pdff("Proton_density_fat_fraction_UKBB_sumstats.txt","Sveinb_PDFF_UKBB"),
                    our_gc("GCST90267352_PDFF_Pazoki2022.tsv.gz","our_Pazoki2022_PDFF"),
       "PDFF UKBB: Sveinbjornsson vs Pazoki-2022"),
  list(report_pair, sv_cc("NAFL_deCODE_sumstat.txt","Sveinb_deCODE_NAFL"),
                    our_gc("Ghodsian_2021_NAFLD_harmonised.tsv.gz","our_Ghodsian_NAFLD"),
       "NAFL deCODE (Sveinb) vs Ghodsian meta (Ghodsian INCLUDES deCODE)"),
  list(report_pair, sv_cc("NAFL_UKBB_sumstat.txt","Sveinb_UKBB_NAFL"),
                    our_gc("Ghodsian_2021_NAFLD_harmonised.tsv.gz","our_Ghodsian_NAFLD"),
       "NAFL UKBB (Sveinb) vs Ghodsian meta (Ghodsian INCLUDES UKBB)"),
  list(report_pair, sv_fg("Cirrhosis_FINNGEN_sumstats.txt","Sveinb_FinnGen_Cirrhosis"),
                    our_fg("finngen_R12_CHIRHEP_NAS.gz","our_FinnGenR12_Cirrhosis"),
       "Cirrhosis FinnGen: Sveinbjornsson vs our R12 CHIRHEP_NAS"),
  list(report_pair, sv_cc("NAFL_UKBB_sumstat.txt","Sveinb_UKBB_NAFL"),
                    our_gc("GCST90019492_UKBB_ALT_harmonised.tsv.gz","our_UKBB_ALT"),
       "Same UKBB cohort, different trait: NAFL diagnosis vs ALT enzyme (genetic corr)")
)

all_res <- rbindlist(lapply(pairs, function(p) tryCatch(p[[1]](p[[2]], p[[3]], p[[4]]),
                              error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })),
                     fill = TRUE)
if (nrow(all_res)) {
  fwrite(all_res, file.path(OUT, "overlap_comparison.csv"))
  cat("\nWritten:", file.path(OUT, "overlap_comparison.csv"), "\n")
  cat("\n=== SUMMARY ===\n"); print(all_res)
}
