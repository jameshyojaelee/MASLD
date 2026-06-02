#!/usr/bin/env Rscript
# 10b_run_mesusie.R
# Run MESuSiE multi-ancestry fine-mapping for one shared locus.
# Generalized to N-ancestry (default EUR+EAS, opt-in EUR+EAS+AFR+SAS via
# ANCESTRIES env var). Per-ancestry GWAS derived from trait_pair (ALT/AST/GGT):
#   EUR  → UKBB_<trait>            (n=343,850)
#   EAS  → BBJ_<trait>             (n=160,000)
#   AFR  → PanUKBB_AFR_<trait>     (n=6,636)
#   SAS  → PanUKBB_CSA_<trait>     (n=8,876)
# Arms with <MIN_SNPS in window are silently dropped; MESuSiE runs on the rest.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(MESuSiE)
})

# ---------------------------------------------------------------------------
# Paths and parameters
# ---------------------------------------------------------------------------

FM_DIR <- Sys.getenv("FM_DIR",
  unset = file.path(Sys.getenv("MASLD_PROJECT_ROOT",
    unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
    "GWAS/finemapping"))

source(file.path(FM_DIR, "src/finemapping_functions.R"))

SHARED_LOCI_FILE <- file.path(FM_DIR, "results/susiex/shared_loci.csv")
MESUSIE_RESULTS_SUFFIX <- Sys.getenv("MESUSIE_RESULTS_SUFFIX", unset = "")
OUTPUT_BASE      <- file.path(FM_DIR, paste0("results/mesusie", MESUSIE_RESULTS_SUFFIX))
MESUSIE_L        <- 10L
LD_REGULARIZE    <- 1e-3
MIN_SNPS         <- 50L

# ---------------------------------------------------------------------------
# Per-ancestry registry
# ---------------------------------------------------------------------------
ANCESTRY_REGISTRY <- list(
  EUR = list(template = "UKBB_{trait}",         n = 343850),
  EAS = list(template = "BBJ_{trait}",          n = 160000),
  AFR = list(template = "PanUKBB_AFR_{trait}",  n = 6636),
  SAS = list(template = "PanUKBB_CSA_{trait}",  n = 8876)
)

build_study_path <- function(anc, trait) {
  template <- ANCESTRY_REGISTRY[[anc]]$template
  study <- gsub("\\{trait\\}", trait, template)
  list(study = study,
       path  = file.path(FM_DIR, paste0("data/sumstats/", study, "_reformatted_hg19.tsv")))
}

get_active_ancestries <- function() {
  raw <- Sys.getenv("ANCESTRIES", unset = "EUR,EAS")
  ancs <- toupper(trimws(unlist(strsplit(raw, ","))))
  ancs <- ancs[nchar(ancs) > 0]
  bad <- setdiff(ancs, names(ANCESTRY_REGISTRY))
  if (length(bad) > 0) stop("Unsupported ancestries: ", paste(bad, collapse = ","))
  ancs
}

# ---------------------------------------------------------------------------
# Locus row
# ---------------------------------------------------------------------------

LOCUS_ROW <- as.integer(Sys.getenv("LOCUS_ROW", unset = ""))
if (is.na(LOCUS_ROW) || LOCUS_ROW == 0L) {
  LOCUS_ROW <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = ""))
}
if (is.na(LOCUS_ROW) || LOCUS_ROW < 1L) {
  stop("LOCUS_ROW must be set via env var LOCUS_ROW or SLURM_ARRAY_TASK_ID")
}

shared_loci <- fread(SHARED_LOCI_FILE)
if (LOCUS_ROW > nrow(shared_loci)) stop(sprintf("LOCUS_ROW %d > %d", LOCUS_ROW, nrow(shared_loci)))
locus      <- shared_loci[LOCUS_ROW]
locus_id   <- locus$locus_id
chrom      <- as.integer(locus$chr)
win_start  <- as.integer(locus$window_start)
win_end    <- as.integer(locus$window_end)
trait_pair <- locus$trait_pair

ancestries <- get_active_ancestries()

out_dir <- file.path(OUTPUT_BASE, trait_pair)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rds_path <- file.path(out_dir, paste0(locus_id, "_mesusie.rds"))
if (file.exists(rds_path)) {
  cat("Already completed:", rds_path, "— skipping\n"); quit(save = "no", status = 0)
}

cat(sprintf("\n%s\nLocus: %s  chr%d:%d-%d\nTrait: %s  Ancestries: %s\n",
            strrep("=", 60), locus_id, chrom, win_start, win_end,
            trait_pair, paste(ancestries, collapse = ",")))

# ---------------------------------------------------------------------------
# Load sumstats per ancestry
# ---------------------------------------------------------------------------

load_sumstats_window <- function(study, path, chr, start, end) {
  if (!file.exists(path)) return(NULL)
  ss <- fread(path)
  nms <- tolower(names(ss)); names(ss) <- nms
  if ("pos_hg19" %in% nms) {
    setnames(ss, "pos_hg19",       "POS",  skip_absent = TRUE)
    setnames(ss, "chromosome",     "CHR",  skip_absent = TRUE)
    setnames(ss, "effect_allele",  "A1",   skip_absent = TRUE)
    setnames(ss, "other_allele",   "A2",   skip_absent = TRUE)
    setnames(ss, "beta",           "BETA", skip_absent = TRUE)
    setnames(ss, "standard_error", "SE",   skip_absent = TRUE)
  } else if ("position" %in% nms) {
    setnames(ss, "chromosome", "CHR",  skip_absent = TRUE)
    setnames(ss, "position",   "POS",  skip_absent = TRUE)
    setnames(ss, "allele1",    "A1",   skip_absent = TRUE)
    setnames(ss, "allele2",    "A2",   skip_absent = TRUE)
    setnames(ss, "beta",       "BETA", skip_absent = TRUE)
    setnames(ss, "se",         "SE",   skip_absent = TRUE)
  } else {
    stop("Unrecognised sumstats format: ", study)
  }
  ss <- ss[, c("CHR", "POS", "A1", "A2", "BETA", "SE"), with = FALSE]
  ss$BETA <- as.numeric(ss$BETA); ss$SE <- as.numeric(ss$SE)
  ss$CHR  <- as.integer(ss$CHR);  ss$POS <- as.integer(ss$POS)
  ss <- ss[!is.na(POS) & !is.na(BETA) & !is.na(SE) & SE > 0]
  ss <- ss[CHR == chr & POS >= start & POS <= end]
  cat(sprintf("  [%s] %d variants in window\n", study, nrow(ss)))
  ss
}

ss_per_anc <- list(); study_per_anc <- list()
for (anc in ancestries) {
  bp <- build_study_path(anc, trait_pair)
  ss <- load_sumstats_window(bp$study, bp$path, chrom, win_start, win_end)
  if (is.null(ss) || nrow(ss) < MIN_SNPS) {
    cat(sprintf("  [%s] dropped (only %d variants)\n", anc,
                if (is.null(ss)) 0 else nrow(ss)))
    next
  }
  ss_per_anc[[anc]] <- ss
  study_per_anc[[anc]] <- bp$study
}
if (length(ss_per_anc) < 2) {
  cat(sprintf("  Fewer than 2 viable arms (%s) — skipping\n",
              paste(names(ss_per_anc), collapse = ","))); quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# Load LD per ancestry
# ---------------------------------------------------------------------------

load_ld_for_window <- function(chr, start, end, ancestry) {
  path_to_LD <- find_LD_block(chr, start, end, ancestry)
  if (is.null(path_to_LD)) return(NULL)
  if (length(path_to_LD) > 1) {
    # review A06#5: bdiag() below zeros inter-block LD for loci spanning >1 Berisa-Pickrell
    # block (~24% of loci). Blocks are ~independent by construction, but flag affected loci.
    warning(sprintf("MESuSiE locus chr%s:%s-%s spans %d LD blocks; bdiag() zeros cross-block LD.",
                    chr, start, end, length(path_to_LD)))
  }
  bim_all <- data.frame(); ld_all <- matrix(ncol = 0, nrow = 0)
  for (file in path_to_LD) {
    bim_f <- paste0(file, ".bim"); ld_f <- paste0(file, ".ld")
    if (!file.exists(bim_f) || !file.exists(ld_f)) next
    BIM <- fread(bim_f); LD <- fread(ld_f)
    bim_all <- rbind(bim_all, BIM); ld_all <- bdiag(ld_all, as.matrix(LD))
  }
  if (nrow(bim_all) == 0) return(NULL)
  colnames(bim_all) <- c("chr", "rsid", "dk", "pos", "alt", "ref")
  bim_all$SNP <- paste(bim_all$chr, bim_all$pos, bim_all$alt, bim_all$ref, sep = ":")
  ld_mat <- as.matrix(ld_all)
  colnames(ld_mat) <- rownames(ld_mat) <- bim_all$SNP
  keep <- !duplicated(bim_all$SNP)
  bim_all <- bim_all[keep, ]; ld_mat <- ld_mat[keep, keep]
  list(bim = bim_all, ld = ld_mat)
}

ld_per_anc <- list()
for (anc in names(ss_per_anc)) {
  cat(sprintf("  Loading %s LD...\n", anc))
  ld <- tryCatch(load_ld_for_window(chrom, win_start, win_end, anc),
                 error = function(e) { cat("    Error:", conditionMessage(e), "\n"); NULL })
  if (is.null(ld)) {
    cat(sprintf("  [%s] LD load failed — dropping arm\n", anc))
    ss_per_anc[[anc]] <- NULL; study_per_anc[[anc]] <- NULL
    next
  }
  ld_per_anc[[anc]] <- ld
}
if (length(ss_per_anc) < 2) {
  cat("  Fewer than 2 arms after LD load — skipping\n"); quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# Per-ancestry alignment of sumstats to LD bim (forward + flip)
# ---------------------------------------------------------------------------
align_ss_to_ld <- function(ss, ld_data) {
  bim <- ld_data$bim
  ss$SNP_fwd <- paste(ss$CHR, ss$POS, ss$A1, ss$A2, sep = ":")
  ss$SNP_rev <- paste(ss$CHR, ss$POS, ss$A2, ss$A1, sep = ":")
  fwd <- ss[SNP_fwd %in% bim$SNP]; fwd$SNP <- fwd$SNP_fwd
  rev <- ss[SNP_rev %in% bim$SNP & !(SNP_fwd %in% bim$SNP)]
  if (nrow(rev) > 0) { rev$BETA <- -rev$BETA; rev$SNP <- rev$SNP_rev }
  aligned <- rbind(fwd, rev, fill = TRUE)
  aligned <- aligned[!duplicated(SNP)]
  aligned <- aligned[SNP %in% bim$SNP]
  idx <- match(aligned$SNP, bim$SNP)
  aligned <- aligned[order(idx)]
  keep <- bim$SNP %in% aligned$SNP
  ld_sub <- ld_data$ld[keep, keep]; bim_sub <- bim[keep, ]
  ld_order <- match(aligned$SNP, bim_sub$SNP)
  ld_sub <- ld_sub[ld_order, ld_order]
  list(ss = aligned, ld = ld_sub)
}

aligned_per_anc <- list()
for (anc in names(ss_per_anc)) {
  aligned_per_anc[[anc]] <- align_ss_to_ld(ss_per_anc[[anc]], ld_per_anc[[anc]])
}

# ---------------------------------------------------------------------------
# Cross-ancestry intersection on chr:pos with allele reconciliation against EUR
# ---------------------------------------------------------------------------
# Pick "anchor" = the first ancestry in the active list (preserves caller order;
# usually EUR). All other arms are oriented to anchor's A1/A2.
anchor <- names(aligned_per_anc)[1]

for (anc in names(aligned_per_anc)) {
  aligned_per_anc[[anc]]$ss$KEY_POS <- paste(aligned_per_anc[[anc]]$ss$CHR,
                                             aligned_per_anc[[anc]]$ss$POS, sep = ":")
  k <- !duplicated(aligned_per_anc[[anc]]$ss$KEY_POS)
  aligned_per_anc[[anc]]$ss <- aligned_per_anc[[anc]]$ss[k]
  aligned_per_anc[[anc]]$ld <- aligned_per_anc[[anc]]$ld[k, k]
}

common_pos <- aligned_per_anc[[anchor]]$ss$KEY_POS
for (anc in setdiff(names(aligned_per_anc), anchor)) {
  common_pos <- intersect(common_pos, aligned_per_anc[[anc]]$ss$KEY_POS)
}

if (length(common_pos) < MIN_SNPS) {
  cat(sprintf("  Common variants %d < %d — skipping\n", length(common_pos), MIN_SNPS))
  quit(save = "no", status = 0)
}

# Anchor alleles
a_idx <- match(common_pos, aligned_per_anc[[anchor]]$ss$KEY_POS)
anchor_alleles <- paste(aligned_per_anc[[anchor]]$ss$A1[a_idx],
                        aligned_per_anc[[anchor]]$ss$A2[a_idx], sep = ":")

# For each non-anchor ancestry: keep variants where alleles match (forward or flip)
keep_mask <- rep(TRUE, length(common_pos))
flip_records <- list()
for (anc in setdiff(names(aligned_per_anc), anchor)) {
  idx <- match(common_pos, aligned_per_anc[[anc]]$ss$KEY_POS)
  fwd_alleles <- paste(aligned_per_anc[[anc]]$ss$A1[idx],
                       aligned_per_anc[[anc]]$ss$A2[idx], sep = ":")
  rev_alleles <- paste(aligned_per_anc[[anc]]$ss$A2[idx],
                       aligned_per_anc[[anc]]$ss$A1[idx], sep = ":")
  ori <- ifelse(anchor_alleles == fwd_alleles, 1L,
         ifelse(anchor_alleles == rev_alleles, -1L, NA_integer_))
  flip_records[[anc]] <- list(idx = idx, ori = ori)
  keep_mask <- keep_mask & !is.na(ori)
}
n_dropped <- sum(!keep_mask)
if (n_dropped > 0) cat(sprintf("  Allele reconciliation dropped %d variants\n", n_dropped))

common_pos <- common_pos[keep_mask]
if (length(common_pos) < MIN_SNPS) {
  cat(sprintf("  After allele reconciliation %d < %d — skipping\n",
              length(common_pos), MIN_SNPS))
  quit(save = "no", status = 0)
}

# Re-derive anchor-orientation indices and apply BETA flips on non-anchor arms
a_idx <- match(common_pos, aligned_per_anc[[anchor]]$ss$KEY_POS)
common_snps <- paste(aligned_per_anc[[anchor]]$ss$CHR[a_idx],
                     aligned_per_anc[[anchor]]$ss$POS[a_idx],
                     aligned_per_anc[[anchor]]$ss$A1[a_idx],
                     aligned_per_anc[[anchor]]$ss$A2[a_idx], sep = ":")

ss_final <- list(); ld_final <- list()
ss_final[[anchor]] <- aligned_per_anc[[anchor]]$ss[a_idx]
ld_final[[anchor]] <- aligned_per_anc[[anchor]]$ld[a_idx, a_idx]

for (anc in setdiff(names(aligned_per_anc), anchor)) {
  rec <- flip_records[[anc]]
  idx_kept <- rec$idx[keep_mask]
  ori_kept <- rec$ori[keep_mask]
  ss_anc <- aligned_per_anc[[anc]]$ss[idx_kept]
  flip_rows <- which(ori_kept == -1L)
  if (length(flip_rows) > 0) ss_anc$BETA[flip_rows] <- -ss_anc$BETA[flip_rows]
  ss_final[[anc]] <- ss_anc
  ld_final[[anc]] <- aligned_per_anc[[anc]]$ld[idx_kept, idx_kept]
}

cat(sprintf("  Common variants across %d arms: %d\n",
            length(ss_final), length(common_snps)))

# ---------------------------------------------------------------------------
# Build MESuSiE inputs
# ---------------------------------------------------------------------------

fix_ld_na <- function(R) {
  if (any(is.na(R))) { R[is.na(R)] <- 0; R <- (R + t(R)) / 2 }
  R
}

R_mat_list <- list(); summary_list <- list()
for (anc in names(ss_final)) {
  R <- as.matrix(ld_final[[anc]]) + LD_REGULARIZE * diag(length(common_snps))
  colnames(R) <- rownames(R) <- common_snps
  R_mat_list[[anc]] <- fix_ld_na(R)
  # MESuSiE's meSuSieData expects post-canonicalization column names
  # (organize_gwas renames lowercase → SNP/Beta/Se/Z/N). We bypass
  # organize_gwas (it's hard-wired to 2 ancestries) so we must produce
  # the canonical names ourselves.
  beta <- ss_final[[anc]]$BETA
  se   <- ss_final[[anc]]$SE
  summary_list[[anc]] <- data.frame(
    SNP  = common_snps,
    Beta = beta,
    Se   = se,
    Z    = beta / se,
    N    = ANCESTRY_REGISTRY[[anc]]$n,
    stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Run MESuSiE
# ---------------------------------------------------------------------------

cat(sprintf("  Running MESuSiE (L=%d, %d arms, %d variants)...\n",
            MESUSIE_L, length(summary_list), length(common_snps)))

# MESuSiE's organize_gwas/organize_ld helpers are hard-wired to 2 ancestries,
# but meSuSie_core accepts arbitrary-length R_mat_list / summary_stat_list.
# We've already aligned variants across all arms, so build the named lists
# directly and pass them to meSuSie_core.
fit <- tryCatch({
  # Validate column names of LD matrices match SNP order (mirror organize_ld checks)
  for (anc in names(summary_list)) {
    if (!identical(colnames(R_mat_list[[anc]]), summary_list[[anc]]$SNP)) {
      stop(sprintf("LD/SS SNP order mismatch for arm '%s'", anc))
    }
  }
  meSuSie_core(R_mat_list = R_mat_list, summary_stat_list = summary_list,
               L = MESUSIE_L, estimate_residual_variance = FALSE,
               max_iter = 200, cor_method = "min.abs.corr", cor_threshold = 0.5)
}, error = function(e) { cat("  MESuSiE error:", conditionMessage(e), "\n"); NULL })

if (is.null(fit)) {
  cat("  SKIP: MESuSiE failed\n"); quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------

pips <- fit$pip; names(pips) <- common_snps
cs_info <- fit$cs
n_cs <- if (!is.null(cs_info$cs)) length(cs_info$cs) else 0L
max_pip <- if (length(pips) > 0) max(pips, na.rm = TRUE) else NA_real_

output <- list(
  locus_id    = locus_id,
  chr         = chrom,
  window      = c(win_start, win_end),
  ancestries  = names(summary_list),
  studies     = unlist(study_per_anc[names(summary_list)]),
  trait_pair  = trait_pair,
  n_variants  = length(common_snps),
  snp_ids     = common_snps,
  pip         = pips,
  pip_config  = fit$pip_config,
  credible_sets = cs_info$cs,
  cs_category   = cs_info$cs_category,
  cs_purity     = cs_info$purity,
  n_cs        = n_cs,
  max_pip     = max_pip,
  converged   = TRUE
)
saveRDS(output, rds_path)
cat(sprintf("  Done: %d arms, %d CS, max PIP %.3f → %s\n",
            length(summary_list), n_cs, max_pip, rds_path))
