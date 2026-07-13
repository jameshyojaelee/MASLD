#!/usr/bin/env Rscript
# 09b_identify_shared_loci_mvp.R
# WITHIN-MVP multi-ancestry shared-locus identification for cross-ancestry joint
# fine-mapping (SuSiEx / MESuSiE). Distinct from 09_identify_shared_loci.R, which
# is HARDCODED to the enzyme UKBB(EUR) <-> BBJ(EAS) pairs — this script does NOT
# mutate that gate. Instead it enumerates MVP per-ancestry same-trait strata
# (MVP_<Trait>_<ANC>, ANC in {EUR,AFR,AMR,EAS}) and defines, per trait, loci
# where a genome-wide-significant (p < 5e-8) lead is shared across >= 2 MVP
# ancestries within LOCUS_MATCH_KB.
#
# Cross-ancestry fine-mapping is GWAS-only (no eQTL), so AMR is scientifically
# valid here (the EUR-eQTL caveat applies to COLOC only). AMR is included.
#
# Input:
#   config/gwas_registry.tsv                 — registry (study_name, sumstats_path, ancestry, N_tot, ...)
#   data/sumstats/MVP_<Trait>_<ANC>_reformatted_hg19.tsv  — hg19 sumstats (chromosome/position/allele1/allele2/beta/se/pval)
# Output:
#   results/susiex_mvp/shared_loci.csv        — MVP shared loci (schema below)
#   results/susiex_mvp/gwsig_cache/<study>_gwsig.tsv — cached GW-sig subsets (streamed via awk)
#
# shared_loci.csv schema (consumed by 10_run_susiex.py + 10b_run_mesusie.R):
#   locus_id, chr, window_start, window_end, trait_pair, cohort(="MVP"),
#   n_ancestries, shared_ancestries, anchor_ancestry, anchor_pos, anchor_p, anchor_snp
#
# Usage (compute node only — streams ~30GB of sumstats the first time):
#   MASLD_PROJECT_ROOT=/path Rscript 09b_identify_shared_loci_mvp.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

# ---------------------------------------------------------------------------
# Paths / parameters
# ---------------------------------------------------------------------------
FM_DIR <- Sys.getenv("FM_DIR",
  unset = file.path(Sys.getenv("MASLD_PROJECT_ROOT",
    unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
    "GWAS/finemapping"))

REGISTRY_FILE  <- file.path(FM_DIR, "config/gwas_registry.tsv")
OUTPUT_DIR     <- file.path(FM_DIR, "results/susiex_mvp")
CACHE_DIR      <- file.path(OUTPUT_DIR, "gwsig_cache")
GW_PVAL        <- 5e-8       # genome-wide significance threshold (matches 09)
CLUMP_KB       <- 500        # clumping window (each side, kb) — matches 09
LOCUS_MATCH_KB <- 500        # max distance to call cross-ancestry leads "shared" — matches 09
FLANK_KB       <- 100        # extra flanking added to merged window — matches 09
MVP_ANCESTRIES <- c("EUR", "AFR", "AMR", "EAS", "SAS")  # recognized MVP ancestry suffixes

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR,  recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load registry (skip comment lines)
# ---------------------------------------------------------------------------
load_registry <- function(path) {
  lines  <- readLines(path)
  active <- lines[!grepl("^\\s*#", lines) & nchar(trimws(lines)) > 0]
  fread(text = paste(active, collapse = "\n"), sep = "\t")
}
registry <- load_registry(REGISTRY_FILE)
cat("Registry loaded:", nrow(registry), "active studies\n")

mvp <- registry[grepl("^MVP_", study_name)]
if (nrow(mvp) == 0) stop("No MVP_* studies in registry: ", REGISTRY_FILE)

# Parse trait + ancestry from study_name = MVP_<Trait>_<ANC>
parse_trait_anc <- function(study) {
  stem <- sub("^MVP_", "", study)                 # <Trait>_<ANC>
  anc  <- sub(".*_", "", stem)                     # <ANC>
  trait<- sub(paste0("_", anc, "$"), "", stem)     # <Trait>
  list(trait = trait, anc = anc)
}
pa <- lapply(mvp$study_name, parse_trait_anc)
mvp[, trait := vapply(pa, `[[`, character(1), "trait")]
mvp[, anc   := vapply(pa, `[[`, character(1), "anc")]
mvp <- mvp[anc %in% MVP_ANCESTRIES]

cat("MVP strata:", nrow(mvp), "across traits:",
    paste(sort(unique(mvp$trait)), collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# Stream-extract genome-wide-significant variants via awk (low memory), cache.
# MVP sumstats columns: chromosome position allele1 allele2 beta se pval (col 7).
# ---------------------------------------------------------------------------
extract_gwsig <- function(study, sumstats_rel) {
  cache <- file.path(CACHE_DIR, paste0(study, "_gwsig.tsv"))
  if (file.exists(cache)) {
    dt <- fread(cache)
    return(if (nrow(dt) == 0) NULL else dt)
  }
  sumpath <- file.path(FM_DIR, sumstats_rel)
  if (!file.exists(sumpath)) {
    cat("    [", study, "] sumstats missing:", sumpath, "\n")
    fwrite(data.table(CHR = integer(), POS = integer(),
                      A1 = character(), A2 = character(), PVAL = numeric()), cache)
    return(NULL)
  }
  thr <- sprintf("%.0e", GW_PVAL)   # "5e-08"
  awk <- sprintf("awk -F'\\t' 'NR>1 && $7!=\"NA\" && ($7+0)<%s'", thr)
  cmd <- paste(awk, shQuote(sumpath))
  dt  <- tryCatch(fread(cmd = cmd, header = FALSE), error = function(e) {
    cat("    [", study, "] awk/fread error:", conditionMessage(e), "\n"); NULL
  })
  if (is.null(dt) || nrow(dt) == 0) {
    fwrite(data.table(CHR = integer(), POS = integer(),
                      A1 = character(), A2 = character(), PVAL = numeric()), cache)
    return(NULL)
  }
  # columns: chromosome position allele1 allele2 beta se pval
  setnames(dt, 1:7, c("CHR", "POS", "A1", "A2", "BETA", "SE", "PVAL"))
  dt <- dt[, .(CHR = suppressWarnings(as.integer(CHR)),
               POS = suppressWarnings(as.integer(POS)),
               A1  = as.character(A1), A2 = as.character(A2),
               PVAL= suppressWarnings(as.numeric(PVAL)))]
  dt <- dt[!is.na(CHR) & !is.na(POS) & !is.na(PVAL)]
  fwrite(dt, cache)
  if (nrow(dt) == 0) NULL else dt
}

# ---------------------------------------------------------------------------
# Clump GW-sig variants into lead SNPs (500kb sliding window; min-p lead) — mirrors 09
# ---------------------------------------------------------------------------
clump_leads <- function(dt) {
  empty <- data.table(chr = integer(), pos = integer(), pval = numeric(),
                      a1 = character(), a2 = character())
  if (is.null(dt) || nrow(dt) == 0) return(empty)
  sig  <- dt[order(PVAL)]
  used <- rep(FALSE, nrow(sig))
  leads <- vector("list", 0)
  for (i in seq_len(nrow(sig))) {
    if (used[i]) next
    ld <- sig[i]
    leads[[length(leads) + 1]] <- data.table(chr = ld$CHR, pos = ld$POS,
                                             pval = ld$PVAL, a1 = ld$A1, a2 = ld$A2)
    w <- sig$CHR == ld$CHR & abs(sig$POS - ld$POS) <= CLUMP_KB * 1000
    used[w] <- TRUE
  }
  rbindlist(leads)
}

# ---------------------------------------------------------------------------
# Per-stratum leads
# ---------------------------------------------------------------------------
leads_by_stratum <- list()
cat("\n=== Extracting GW-sig leads per MVP stratum ===\n")
for (i in seq_len(nrow(mvp))) {
  study <- mvp$study_name[i]
  cat(sprintf("  [%2d/%2d] %-22s ", i, nrow(mvp), study))
  gw <- extract_gwsig(study, mvp$sumstats_path[i])
  lds <- clump_leads(gw)
  leads_by_stratum[[study]] <- lds
  cat(sprintf("GW-sig=%d  leads=%d\n",
              if (is.null(gw)) 0L else nrow(gw), nrow(lds)))
}

# ---------------------------------------------------------------------------
# Per-trait cross-ancestry clustering: pooled greedy over leads (by p-value).
# A locus qualifies if >= 2 distinct MVP ancestries have a GW-sig lead within
# LOCUS_MATCH_KB of the anchor (mirrors the enzyme 09 "shared" definition,
# generalized from pairwise EUR<->EAS to N-way within-MVP).
# ---------------------------------------------------------------------------
shared_all <- list()

for (tr in sort(unique(mvp$trait))) {
  strata <- mvp[trait == tr]
  pooled <- vector("list", 0)
  for (j in seq_len(nrow(strata))) {
    st <- strata$study_name[j]; an <- strata$anc[j]
    lds <- leads_by_stratum[[st]]
    if (is.null(lds) || nrow(lds) == 0) next
    lds <- copy(lds); lds[, `:=`(anc = an, study = st)]
    pooled[[length(pooled) + 1]] <- lds
  }
  if (length(pooled) == 0) { cat(sprintf("\nTrait %-10s: no leads\n", tr)); next }
  pool <- rbindlist(pooled)
  setorder(pool, pval)
  n <- nrow(pool)
  used <- rep(FALSE, n)
  emitted <- 0L

  for (i in seq_len(n)) {
    if (used[i]) next
    anchor <- pool[i]
    near <- which(pool$chr == anchor$chr &
                  abs(pool$pos - anchor$pos) <= LOCUS_MATCH_KB * 1000 &
                  !used)
    ancs_here <- unique(pool$anc[near])
    if (length(ancs_here) >= 2) {
      contrib_pos <- pool$pos[near]
      win_start <- max(0L, min(contrib_pos) - FLANK_KB * 1000L)
      win_end   <- max(contrib_pos) + FLANK_KB * 1000L
      locus_id  <- sprintf("locus_MVP_%s_chr%d_%d", tr, anchor$chr, anchor$pos)
      shared_all[[length(shared_all) + 1]] <- data.table(
        locus_id          = locus_id,
        chr               = as.integer(anchor$chr),
        window_start      = as.integer(win_start),
        window_end        = as.integer(win_end),
        trait_pair        = tr,
        cohort            = "MVP",
        n_ancestries      = length(ancs_here),
        shared_ancestries = paste(sort(ancs_here), collapse = ";"),
        anchor_ancestry   = anchor$anc,
        anchor_pos        = as.integer(anchor$pos),
        anchor_p          = anchor$pval,
        anchor_snp        = sprintf("%d:%d:%s:%s", anchor$chr, anchor$pos,
                                    anchor$a1, anchor$a2)
      )
      emitted <- emitted + 1L
      # Consume all pooled leads within CLUMP_KB of anchor (same chr) into this locus
      mask <- pool$chr == anchor$chr &
              abs(pool$pos - anchor$pos) <= CLUMP_KB * 1000
      used[mask] <- TRUE
    } else {
      used[i] <- TRUE   # singleton ancestry — drop
    }
  }
  cat(sprintf("Trait %-10s: strata=%d  pooled_leads=%d  shared_loci=%d\n",
              tr, nrow(strata), n, emitted))
}

# ---------------------------------------------------------------------------
# Assemble + write
# ---------------------------------------------------------------------------
if (length(shared_all) == 0) {
  cat("\nNo MVP shared loci identified across any trait.\n")
  result <- data.table(
    locus_id = character(), chr = integer(),
    window_start = integer(), window_end = integer(),
    trait_pair = character(), cohort = character(),
    n_ancestries = integer(), shared_ancestries = character(),
    anchor_ancestry = character(), anchor_pos = integer(),
    anchor_p = numeric(), anchor_snp = character())
} else {
  result <- rbindlist(shared_all)
  setorder(result, trait_pair, chr, window_start)
  cat("\n=== Summary ===\n")
  cat("Total MVP shared loci:", nrow(result), "\n")
  cat("\nPer trait:\n"); print(table(result$trait_pair))
  cat("\nPer n_ancestries:\n"); print(table(result$n_ancestries))
  cat("\nPer shared-ancestry combination:\n"); print(sort(table(result$shared_ancestries), decreasing = TRUE))
}

out_path <- file.path(OUTPUT_DIR, "shared_loci.csv")
write_csv(result, out_path)
cat("\nOutput written to:", out_path, "\n")
