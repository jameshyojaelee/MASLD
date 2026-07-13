#!/usr/bin/env Rscript
# =============================================================================
# 63_broadaway_bench_prep.R
# DELIVERABLE 5 (Stage-2 data prep) — Broadaway concordance-benchmark TRUTH table
# -----------------------------------------------------------------------------
# Builds the ground-truth table that Borzoi / seq-function predictions will be
# scored against: credible-set variants that (a) COLOCALIZE (coloc PP.H4 > 0.5)
# and (b) have a MEASURED Broadaway liver cis-eQTL, carrying the measured eQTL
# effect size + DIRECTION (sign) with a fully-documented allele orientation.
#
# ADDITIVE ONLY. Reads existing substrate; writes a NEW file under
# GWAS/finemapping/results/seqfunc/. Edits no existing script.
#
# -----------------------------------------------------------------------------
# ALLELE ORIENTATION CONTRACT (the #1 sign-flip trap — read before scoring)
# -----------------------------------------------------------------------------
#  * Broadaway marginal files carry columns EA (effect allele), NEA (non-effect
#    allele) and Beta. Beta is oriented to EA: Beta>0 => the EA allele INCREASES
#    expression of the tested gene. (Confirmed against the exact merge in
#    06_susie_coloc.R:256-279, which sets eqtl_ea=EA, eqtl_nea=NEA and negates
#    Beta when the eQTL EA/NEA are reversed vs the GWAS a1/a2.)
#  * Each credible-set variant_id is chr:pos:ref:alt, with cs_member allele1=ref,
#    allele2=alt. We RE-ORIENT the Broadaway Beta to the variant's ALT allele so
#    that it matches Borzoi's ref->alt convention:
#         orientation == "match"  (EA==alt & NEA==ref)  ->  eqtl_beta =  Beta
#         orientation == "flip"   (EA==ref & NEA==alt)  ->  eqtl_beta = -Beta
#         orientation == "mismatch"(neither)            ->  DROPPED (alleles do
#                                                            not reconcile; indel
#                                                            representation diff.)
#    Output column `eqtl_beta` and `eqtl_sign` are therefore in ref->alt space:
#         eqtl_beta > 0  (eqtl_sign = +1)  =>  ALT allele INCREASES `gene` expr.
#         eqtl_beta < 0  (eqtl_sign = -1)  =>  ALT allele DECREASES `gene` expr.
#    `effect_allele` is set to `alt` for every retained row (post-reorientation).
#  * DOWNSTREAM CONTRACT: the Borzoi step MUST score the effect of substituting
#    THIS ref -> THIS alt (the ref/alt columns here), NOT genome-reference->other.
#    ref/alt here are the finemapping-harmonized alleles (allele1/allele2), which
#    reconcile with Broadaway EA/NEA at ~100% for non-palindromic SNVs (0 mismatch
#    observed on the chr22 pilot). Keep the same ref->alt polarity on both sides.
#  * STRAND-AMBIGUOUS (palindromic A/T or C/G) variants are FLAGGED
#    (strand_ambiguous=TRUE) but NOT dropped: letter-matching cannot detect a
#    strand flip for palindromes, so their sign may be unreliable — exclude them
#    (or LD-verify) before a strict sign-auROC. Raw Broadaway columns
#    (eqtl_ea, eqtl_nea, eqtl_beta_ea_raw) are passed through so any downstream
#    re-orientation is auditable.
#
# PP.H4 SOURCE NOTE: the task named lead_causal_annotation.csv for pp4_best, but
# that file's pp4_susie/pp4_abf/pp4_best columns are EMPTY (all NA). The populated
# per-variant coloc PP.H4 lives in cs_member_annotation.csv (coloc_pp4_best) — the
# same quantity — so pp4_best is sourced there. cs_member is already gated to
# coloc_pp4_best>0.5, so it IS the colocalizing credible-set variant substrate.
#
# SIGNAL vs CS-MEMBER granularity: one output row per Broadaway-tested
# colocalizing (variant x gene) pair. This is CS-member-level and therefore
# LD-pseudoreplicated (dozens of tagging variants per signal). `signal_id`
# = gene::locus and `is_signal_lead` (TRUE = strongest measured eQTL, i.e. min
# eqtl_p, within the signal) mark ONE independent representative per gene x
# fine-map locus. Run the concordance benchmark on is_signal_lead rows to avoid
# LD pseudoreplication; the full set is emitted for transparency.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CS_FILE   <- file.path(BASE, "RNA-seq/results/coloc_variant_classes/cs_member_annotation.csv")
EQTL_DIR  <- file.path(BASE, "data/broadaway_eqtl")
CHAIN     <- file.path(EQTL_DIR, "hg19ToHg38.over.chain")
OUT_DIR   <- file.path(BASE, "GWAS/finemapping/results/seqfunc")
OUT_TSV   <- file.path(OUT_DIR, "broadaway_benchmark_truth.tsv")
OUT_PROV  <- file.path(OUT_DIR, "broadaway_benchmark_truth.README.txt")
PP4_THRESH <- 0.5
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

log <- function(...) cat(sprintf(...), "\n")

# -----------------------------------------------------------------------------
# 1. Colocalizing credible-set variant substrate (cs_member, pp4>0.5)
# -----------------------------------------------------------------------------
log("--- Loading colocalizing credible-set variants (cs_member) ---")
cs <- fread(CS_FILE)
cs <- cs[!is.na(coloc_pp4_best) & coloc_pp4_best > PP4_THRESH]
log("  cs_member rows pp4>%.2f: %d", PP4_THRESH, nrow(cs))

# parse variant_id = chr:pos:ref:alt
vid_parts <- tstrsplit(cs$variant_id, ":", fixed = TRUE)
cs[, `:=`(vid_chr = vid_parts[[1]],
          vid_pos = suppressWarnings(as.integer(vid_parts[[2]])),
          ref     = vid_parts[[3]],
          alt     = vid_parts[[4]])]
# sanity: cs allele1/allele2 should equal ref/alt (variant_id convention)
n_allele_consistent <- sum(cs$allele1 == cs$ref & cs$allele2 == cs$alt, na.rm = TRUE)
log("  variant_id ref/alt == allele1/allele2 : %d / %d", n_allele_consistent, nrow(cs))

# explode coloc_genes (';' , '|' or ',' separated) into one row per (variant, gene)
split_genes <- function(s) {
  s <- trimws(s)
  if (is.na(s) || s == "" || s == "NA") return(character(0))
  parts <- unlist(strsplit(s, "[;,|]"))
  parts <- trimws(parts); parts[parts != ""]
}
cs_list <- cs[, {
  gl <- split_genes(coloc_genes)
  if (length(gl) == 0) NULL else .(gene = gl)
}, by = .(variant_id, vid_chr, vid_pos, ref, alt, allele1, allele2,
          coloc_pp4_best, study, locus, susie_cs, ancestry, trait)]
setnames(cs_list, "coloc_pp4_best", "pp4_best")
cs_list <- unique(cs_list)
log("  exploded colocalizing (variant x gene) target rows: %d", nrow(cs_list))
log("  unique variants: %d | unique genes: %d",
    uniqueN(cs_list$variant_id), uniqueN(cs_list$gene))

# -----------------------------------------------------------------------------
# 2. Join to Broadaway marginal eQTL per chromosome (both hg19)
# -----------------------------------------------------------------------------
palindrome <- function(a, b) {
  key <- paste(pmin(a, b), pmax(a, b))
  key %in% c(paste("A", "T"), paste("C", "G"))
}

autosomes <- sort(unique(cs_list$vid_chr))
autosomes <- autosomes[autosomes %in% as.character(1:22)]
skipped_chr <- setdiff(unique(cs_list$vid_chr), autosomes)
if (length(skipped_chr) > 0)
  log("  NOTE: %d target rows on non-autosomal chr (%s) have no Broadaway marginal file -> unmatched",
      sum(cs_list$vid_chr %in% skipped_chr), paste(skipped_chr, collapse = ","))

hits_all <- list()
n_mismatch_total <- 0L
for (chr in autosomes) {
  eqtl_file <- file.path(EQTL_DIR, sprintf("chr%s_marginal_summary_results.tsv", chr))
  if (!file.exists(eqtl_file)) { log("  chr%s: no Broadaway file, skip", chr); next }
  tgt <- cs_list[vid_chr == chr]
  tgt_pos <- unique(tgt$vid_pos)

  eq <- fread(eqtl_file,
              select = c("POS", "NEA", "EA", "Beta", "SE", "PVAL",
                         "GeneSymbol", "ENSG"))
  # Broadaway PVAL uses an "e-1" exponent style that fread may type as character;
  # force numeric so downstream ordering (is_signal_lead) is by NUMERIC p, not
  # lexical string order (silent lead-selection bug otherwise).
  eq[, `:=`(POS  = as.integer(POS),
            Beta = as.numeric(Beta),
            SE   = as.numeric(SE),
            PVAL = as.numeric(PVAL))]
  eq <- eq[POS %in% tgt_pos]
  if (nrow(eq) == 0) { log("  chr%s: 0 Broadaway rows at target positions", chr); next }

  # merge on (pos, gene symbol)
  m <- merge(tgt, eq,
             by.x = c("vid_pos", "gene"),
             by.y = c("POS", "GeneSymbol"),
             allow.cartesian = TRUE)
  if (nrow(m) == 0) { log("  chr%s: 0 (pos,gene) merges", chr); next }

  # orient Broadaway Beta (oriented to EA) into ref->alt (alt) space
  m[, orientation := fifelse(EA == alt & NEA == ref, "match",
                      fifelse(EA == ref & NEA == alt, "flip", "mismatch"))]
  n_mm <- sum(m$orientation == "mismatch")
  n_mismatch_total <- n_mismatch_total + n_mm
  m <- m[orientation != "mismatch"]
  if (nrow(m) == 0) { log("  chr%s: all merges allele-mismatch (%d)", chr, n_mm); next }

  m[, eqtl_beta_ea_raw := Beta]                                   # Broadaway native (EA-oriented)
  m[, eqtl_beta := fifelse(orientation == "flip", -Beta, Beta)]   # ref->alt (alt) space
  m[, eqtl_sign := fifelse(eqtl_beta >= 0, 1L, -1L)]
  m[, strand_ambiguous := palindrome(ref, alt)]
  # a (variant,gene) may hit >1 Broadaway row only via indel POS collisions;
  # keep the allele-reconciling row with the strongest eQTL
  setorder(m, variant_id, gene, PVAL)
  m <- m[, .SD[1], by = .(variant_id, gene)]
  hits_all[[chr]] <- m
  log("  chr%s: %d Broadaway-tested colocalizing (variant x gene) rows (%d mismatch dropped)",
      chr, nrow(m), n_mm)
}
truth <- rbindlist(hits_all, use.names = TRUE, fill = TRUE)
log("--- Merged: %d Broadaway-tested colocalizing (variant x gene) rows ---", nrow(truth))
log("    allele-mismatch dropped (all chr): %d", n_mismatch_total)

if (nrow(truth) == 0) stop("No overlap rows produced — investigate before scoring.")

# -----------------------------------------------------------------------------
# 3. LiftOver hg19 -> hg38 (rtracklayer; same chain/idiom as 55/55d)
# -----------------------------------------------------------------------------
log("--- LiftOver hg19 -> hg38 ---")
chain <- import.chain(CHAIN)
uniq_pos <- unique(truth[, .(vid_chr, vid_pos)])
gr19 <- GRanges(paste0("chr", uniq_pos$vid_chr),
                IRanges(uniq_pos$vid_pos, width = 1))
lifted <- liftOver(gr19, chain)
n_map  <- lengths(lifted)
uniq_pos[, pos_hg38 := NA_integer_]
ok <- which(n_map == 1L)
uniq_pos[ok, pos_hg38 := start(unlist(lifted[ok]))]
log("  lifted %d/%d unique positions (%d unmapped/multi)",
    length(ok), nrow(uniq_pos), sum(n_map != 1L))
truth <- merge(truth, uniq_pos, by = c("vid_chr", "vid_pos"), all.x = TRUE)

# -----------------------------------------------------------------------------
# 4. Independent-signal flag (gene x fine-map locus; lead = strongest eQTL)
# -----------------------------------------------------------------------------
truth[, signal_id := paste(gene, locus, sep = "::")]
setorder(truth, signal_id, PVAL)
truth[, is_signal_lead := seq_len(.N) == 1L, by = signal_id]

# -----------------------------------------------------------------------------
# 5. Assemble + write
# -----------------------------------------------------------------------------
out <- truth[, .(
  variant_id,
  chr           = vid_chr,
  pos_hg19      = vid_pos,
  pos_hg38,
  ref, alt,
  effect_allele = alt,          # eqtl_beta re-oriented to alt (ref->alt space)
  gene,
  ensembl       = ENSG,
  eqtl_beta,                     # ref->alt: >0 => alt increases gene expr
  eqtl_sign,                     # sign(eqtl_beta): +1 / -1
  eqtl_se       = SE,
  eqtl_p        = PVAL,
  pp4_best,
  # ---- provenance / audit passthrough (allele-orientation trail) ----
  eqtl_ea          = EA,
  eqtl_nea         = NEA,
  eqtl_beta_ea_raw,              # Broadaway native beta (oriented to EA)
  orientation,                   # match | flip
  strand_ambiguous,              # TRUE => palindromic, sign may be unreliable
  study, locus, susie_cs, ancestry, trait,
  signal_id, is_signal_lead
)]
setorder(out, chr, pos_hg19, gene)
fwrite(out, OUT_TSV, sep = "\t")
log("--- WROTE %s (%d rows) ---", OUT_TSV, nrow(out))

# -----------------------------------------------------------------------------
# 6. Summary + provenance sidecar
# -----------------------------------------------------------------------------
n_rows    <- nrow(out)
n_var     <- uniqueN(out$variant_id)
n_gene    <- uniqueN(out$gene)
n_signal  <- uniqueN(out$signal_id)
n_lead    <- sum(out$is_signal_lead)
n_palin   <- sum(out$strand_ambiguous)
n_flip    <- sum(out$orientation == "flip")
n_indel   <- sum(nchar(out$ref) != 1L | nchar(out$alt) != 1L)
n_hg38na  <- sum(is.na(out$pos_hg38))

log("")
log("=================== OVERLAP-SET SUMMARY ===================")
log("  rows (Broadaway-tested colocalizing variant x gene): %d", n_rows)
log("  unique variants                                     : %d", n_var)
log("  unique genes                                        : %d", n_gene)
log("  independent signals (gene x locus)                  : %d", n_signal)
log("  is_signal_lead rows (recommended benchmark subset)  : %d", n_lead)
log("  allele-flip rows (Beta negated)                     : %d", n_flip)
log("  strand-ambiguous (palindromic) rows                 : %d", n_palin)
log("  indel rows (ref/alt not single-base)                : %d", n_indel)
log("  rows with unmapped hg38                             : %d", n_hg38na)
log("===========================================================")

prov <- c(
  "Broadaway concordance-benchmark TRUTH table — provenance",
  "=========================================================",
  sprintf("Generated: %s", Sys.time()),
  sprintf("Script   : GWAS/finemapping/src/63_broadaway_bench_prep.R"),
  sprintf("Output   : %s", OUT_TSV),
  "",
  "WHAT THIS IS: credible-set variants that COLOCALIZE (coloc PP.H4 > 0.5) AND",
  "have a MEASURED Broadaway liver cis-eQTL, with the measured eQTL effect size",
  "and direction (sign) to score Borzoi/seq-function predictions against.",
  "",
  "SOURCES:",
  sprintf("  colocalizing CS variants : %s (coloc_pp4_best > %.2f)", CS_FILE, PP4_THRESH),
  sprintf("  measured eQTL            : %s/chr{1..22}_marginal_summary_results.tsv (hg19)", EQTL_DIR),
  sprintf("  liftover chain           : %s", CHAIN),
  "",
  "PP.H4 SOURCE: lead_causal_annotation.csv (named in the task) has EMPTY pp4",
  "columns; the populated per-variant coloc PP.H4 is coloc_pp4_best in",
  "cs_member_annotation.csv (same quantity) and is used as `pp4_best` here.",
  "",
  "ALLELE ORIENTATION (the sign-flip contract):",
  "  Broadaway Beta is oriented to EA (effect allele): Beta>0 => EA increases",
  "  expression. We re-orient to the variant ALT allele (ref->alt, Borzoi",
  "  convention): match(EA=alt,NEA=ref)->beta=Beta ; flip(EA=ref,NEA=alt)->",
  "  beta=-Beta ; mismatch->dropped. Output `eqtl_beta`/`eqtl_sign` are in",
  "  ref->alt space: eqtl_beta>0 (sign=+1) => ALT allele INCREASES gene expr.",
  "  `effect_allele`=alt for every row. ref/alt are finemapping-harmonized",
  "  alleles (variant_id chr:pos:ref:alt; allele1=ref, allele2=alt).",
  "  DOWNSTREAM: score Borzoi as substituting THIS ref -> THIS alt (same polarity)",
  "  or the sign-auROC will silently flip. Palindromic (A/T,C/G) variants are",
  "  flagged strand_ambiguous=TRUE (letter-matching can't detect a strand flip)",
  "  — exclude or LD-verify before a strict sign-auROC. Raw EA-oriented beta is",
  "  kept as eqtl_beta_ea_raw for audit.",
  "",
  "GRANULARITY: one row per Broadaway-tested colocalizing (variant x gene) pair",
  "(CS-member level = LD-pseudoreplicated). Use is_signal_lead==TRUE (one row per",
  "gene x fine-map locus, strongest measured eQTL) for the independent-signal",
  "benchmark to avoid LD pseudoreplication.",
  "",
  "COUNTS:",
  sprintf("  rows                         : %d", n_rows),
  sprintf("  unique variants              : %d", n_var),
  sprintf("  unique genes                 : %d", n_gene),
  sprintf("  independent signals          : %d", n_signal),
  sprintf("  is_signal_lead rows          : %d", n_lead),
  sprintf("  allele-flip rows             : %d", n_flip),
  sprintf("  strand-ambiguous rows        : %d", n_palin),
  sprintf("  indel rows                   : %d", n_indel),
  sprintf("  allele-mismatch dropped      : %d", n_mismatch_total),
  sprintf("  unmapped hg38                : %d", n_hg38na)
)
writeLines(prov, OUT_PROV)
log("--- WROTE provenance %s ---", OUT_PROV)
log("DONE.")
