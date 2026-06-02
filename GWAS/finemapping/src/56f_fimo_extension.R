#!/usr/bin/env Rscript
# 56f_fimo_extension.R
# A5 ATAC improvement: extend FIMO validation from top-500 cap to:
#   Tier 1 (canonical): ALL strong motifbreakR hits (effect == "strong")
#   Tier 2 (extended) : Non-strong hits with |alleleDiff| > 0.3 (alleleDiff-driven)
# Quantitative direction concordance per |alleleDiff| stratum.
#
# Does NOT modify Script 56. Reads motif_disruption_scores.csv (Script 56 output).
# Usage: Rscript 56f_fimo_extension.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(MotifDb)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE_DIR, "GWAS/finemapping")
ATAC_DIR <- file.path(BASE_DIR, "Analysis/ATAC/Human_Multiome")
OUT_DIR  <- file.path(FM_DIR, "results/gwas_atac")
TMP_ROOT <- Sys.getenv("TMPDIR", unset = file.path(OUT_DIR, "fimo_tmp_56f"))
dir.create(TMP_ROOT, showWarnings = FALSE, recursive = TRUE)
TMP_DIR  <- file.path(TMP_ROOT, paste0("fimo56f_", Sys.getpid()))
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("56f_fimo_extension.R\n")
cat("FIMO validation for ALL motifbreakR strong + alleleDiff-driven hits\n")
cat("TMP_DIR =", TMP_DIR, "\n")
cat("============================================================\n\n")

# ── 1. Load motif disruption scores + variant annotation ────────────────────
scores_file <- file.path(OUT_DIR, "motif_disruption_scores.csv")
ann_file    <- file.path(OUT_DIR, "gwas_atac_variant_annotation.csv")
stopifnot(file.exists(scores_file), file.exists(ann_file))

scores <- fread(scores_file)
ann    <- fread(ann_file)
cat("Loaded", nrow(scores), "motif disruption rows;",
    nrow(ann), "variant annotations\n")

# Recover REF/ALT per SNP_id from BSgenome (same logic as Script 56)
genome <- BSgenome.Hsapiens.UCSC.hg38

variants <- ann %>%
  distinct(chromosome, position, allele1, allele2, chr_hg38, pos_hg38, max_pip) %>%
  mutate(SNP_id = paste(chromosome, position, allele1, allele2, sep = ":"))

gr_v <- GRanges(seqnames = variants$chr_hg38,
                ranges = IRanges(start = variants$pos_hg38, width = 1))
ref_alleles <- as.character(getSeq(genome, gr_v))
variants$ref_genome <- ref_alleles
variants$REF <- ifelse(variants$allele1 == ref_alleles, variants$allele1,
                ifelse(variants$allele2 == ref_alleles, variants$allele2, NA_character_))
variants$ALT <- ifelse(variants$allele1 == ref_alleles, variants$allele2,
                ifelse(variants$allele2 == ref_alleles, variants$allele1, NA_character_))
variants <- variants %>% filter(!is.na(REF), nchar(REF) == 1, nchar(ALT) == 1)
cat("SNV variants with REF/ALT resolved:", nrow(variants), "\n")

# ── 2. Tier definitions ─────────────────────────────────────────────────────
tier1 <- scores %>% filter(effect == "strong")
tier2_pool <- scores %>% filter(effect != "strong", abs(alleleDiff) > 0.3)

# Tier 2 stratified sample: random 500 + top 500 by |alleleDiff|, dedupe
set.seed(42)
if (nrow(tier2_pool) <= 1000) {
  tier2 <- tier2_pool
} else {
  s1 <- tier2_pool %>% slice_sample(n = 500)
  s2 <- tier2_pool %>% arrange(desc(abs(alleleDiff))) %>% head(500)
  tier2 <- bind_rows(s1, s2) %>% distinct()
}
tier1$tier <- "tier1_strong"
tier2$tier <- "tier2_alleleDiff_driven"
work <- bind_rows(tier1, tier2)
cat("Tier 1 (strong) rows:        ", nrow(tier1), "\n")
cat("Tier 2 (alleleDiff-driven):  ", nrow(tier2), "\n")
cat("Total pairs to FIMO-validate:", nrow(work), "\n")

# ── 3. Build REF and ALT sequences for unique variants ──────────────────────
flank_size <- 25L
work_var_ids <- unique(work$SNP_id)
vmap <- variants %>% filter(SNP_id %in% work_var_ids) %>%
  distinct(SNP_id, .keep_all = TRUE)

# Some SNP_ids in scores may not appear in vmap (e.g., scores carries seqnames/start derived differently).
# Use seqnames + start from scores as fallback for those.
missing <- setdiff(work_var_ids, vmap$SNP_id)
if (length(missing) > 0) {
  cat("WARNING:", length(missing),
      "SNP_ids not resolvable from variant annotation; recovering via scores coords\n")
  fb <- work %>% filter(SNP_id %in% missing) %>%
    distinct(SNP_id, seqnames, start) %>%
    mutate(chr_hg38 = as.character(seqnames), pos_hg38 = as.integer(start))
  gr_fb <- GRanges(seqnames = fb$chr_hg38,
                   ranges = IRanges(start = fb$pos_hg38, width = 1))
  fb$ref_genome <- as.character(getSeq(genome, gr_fb))
  # Parse SNP_id "chr:pos:A1:A2"
  parts <- do.call(rbind, strsplit(fb$SNP_id, ":", fixed = TRUE))
  fb$allele1 <- parts[, 3]; fb$allele2 <- parts[, 4]
  fb$REF <- ifelse(fb$allele1 == fb$ref_genome, fb$allele1,
            ifelse(fb$allele2 == fb$ref_genome, fb$allele2, NA_character_))
  fb$ALT <- ifelse(fb$allele1 == fb$ref_genome, fb$allele2,
            ifelse(fb$allele2 == fb$ref_genome, fb$allele1, NA_character_))
  fb <- fb %>% filter(!is.na(REF), nchar(REF) == 1, nchar(ALT) == 1)
  vmap <- bind_rows(vmap %>% select(SNP_id, chr_hg38, pos_hg38, REF, ALT, max_pip),
                    fb %>% mutate(max_pip = NA_real_) %>%
                          select(SNP_id, chr_hg38, pos_hg38, REF, ALT, max_pip))
}
vmap <- vmap %>% filter(SNP_id %in% work_var_ids) %>% distinct(SNP_id, .keep_all = TRUE)
cat("Resolved", nrow(vmap), "unique variants for FIMO sequence extraction\n")

# FIMO truncates FASTA names at colons/whitespace, so use a colon-free surrogate id
# and map back via a lookup table.
vmap <- as.data.table(vmap)
vmap[, fasta_id := sprintf("snp%07d", seq_len(.N))]
seq_lookup <- vmap[, .(fasta_id, SNP_id)]

gr_seq <- GRanges(seqnames = vmap$chr_hg38,
                  ranges = IRanges(start = vmap$pos_hg38 - flank_size,
                                   end   = vmap$pos_hg38 + flank_size))
ref_seqs <- getSeq(genome, gr_seq)
names(ref_seqs) <- vmap$fasta_id
alt_seqs <- ref_seqs
pos_in_seq <- flank_size + 1L
alt_chars  <- as.character(ref_seqs)
for (i in seq_along(alt_chars)) {
  s <- alt_chars[i]
  substr(s, pos_in_seq, pos_in_seq) <- vmap$ALT[i]
  alt_chars[i] <- s
}
alt_seqs <- DNAStringSet(alt_chars)
names(alt_seqs) <- vmap$fasta_id

ref_fa <- file.path(TMP_DIR, "ref.fa")
alt_fa <- file.path(TMP_DIR, "alt.fa")
writeXStringSet(ref_seqs, ref_fa)
writeXStringSet(alt_seqs, alt_fa)
cat("Wrote REF/ALT FASTA (", length(ref_seqs), " sequences each)\n")

# ── 4. Build MEME motif file matching motifbreakR's selection ───────────────
cat("\n--- Building MEME motif file ---\n")
jaspar_motifs <- query(MotifDb, andStrings = c("hsapiens"),
                       orStrings = c("jaspar2022", "jaspar2024"))
if (length(jaspar_motifs) == 0) {
  jaspar_motifs <- query(MotifDb, "jaspar")
  jaspar_motifs <- jaspar_motifs[grep("Hsapiens", names(jaspar_motifs))]
}
if (length(jaspar_motifs) == 0) {
  jaspar_motifs <- query(MotifDb, "HOCOMOCOv11-core-A")
  if (length(jaspar_motifs) == 0) jaspar_motifs <- query(MotifDb, "HOCOMOCO")
}
cat("Loaded", length(jaspar_motifs), "motifs from MotifDb\n")

mdb_md <- values(jaspar_motifs)
mdb_geneSymbol <- toupper(as.character(mdb_md$geneSymbol))
mdb_providerName <- toupper(as.character(mdb_md$providerName))
mdb_providerId <- as.character(mdb_md$providerId)
mdb_full_names <- names(jaspar_motifs)

# tf_name in scores is what Script 56 set: geneSymbol (or providerName/Id fallback)
need_tfs <- unique(toupper(work$tf_name))
keep_idx <- which(mdb_geneSymbol %in% need_tfs | mdb_providerName %in% need_tfs)
hit_motifs <- jaspar_motifs[keep_idx]
cat("Matched", length(hit_motifs), "MotifDb entries for",
    length(need_tfs), "unique tf_names in input\n")

# Map: motif_id (we'll use full MotifDb name as MEME id) -> uppercased tf
motif_tf_map <- data.table(
  meme_id    = mdb_full_names[keep_idx],
  tf_upper   = mdb_geneSymbol[keep_idx],
  prov_name  = mdb_providerName[keep_idx],
  prov_id    = mdb_providerId[keep_idx]
)
# Fallback: if geneSymbol is NA/empty, use providerName
motif_tf_map[is.na(tf_upper) | tf_upper == "", tf_upper := prov_name]

# Write MEME file. Use shortened id (no spaces) for FIMO.
motif_tf_map[, fimo_id := paste0("M", sprintf("%05d", seq_len(.N)), "_",
                                  gsub("[^A-Za-z0-9_.-]", "_", tf_upper))]

meme_file <- file.path(TMP_DIR, "motifs.meme")
con <- file(meme_file, "w")
writeLines("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n", con)
writeLines("Background letter frequencies", con)
writeLines("A 0.25 C 0.25 G 0.25 T 0.25\n", con)
for (i in seq_along(hit_motifs)) {
  pwm <- as.matrix(hit_motifs[[i]])  # rows=A,C,G,T; cols=positions
  # MEME letter-probability matrix expects rows=positions, cols=A,C,G,T
  if (nrow(pwm) == 4 && all(rownames(pwm) %in% c("A", "C", "G", "T"))) {
    pwm <- t(pwm[c("A", "C", "G", "T"), , drop = FALSE])
  }
  # Normalize columns just in case
  pwm <- pmax(pwm, 0)
  rs <- rowSums(pwm)
  rs[rs == 0] <- 1
  pwm <- pwm / rs
  writeLines(paste("MOTIF", motif_tf_map$fimo_id[i]), con)
  writeLines(sprintf("letter-probability matrix: alength= 4 w= %d nsites= 20 E= 0",
                     nrow(pwm)), con)
  writeLines(apply(pwm, 1, function(r) paste(sprintf("%.6f", r), collapse = " ")), con)
  writeLines("", con)
}
close(con)
cat("Wrote MEME file with", length(hit_motifs), "motifs:", meme_file, "\n")

# ── 5. Run FIMO on REF and ALT ──────────────────────────────────────────────
fimo_bin <- Sys.which("fimo")
if (fimo_bin == "") fimo_bin <- "/nfs/sw/easybuild/software/meme/5.5.9/bin/fimo"
stopifnot(file.exists(fimo_bin))
cat("FIMO bin:", fimo_bin, "\n")

ref_out <- file.path(TMP_DIR, "fimo_ref.tsv")
alt_out <- file.path(TMP_DIR, "fimo_alt.tsv")

run_fimo <- function(meme, fa, out) {
  cmd <- sprintf("%s --thresh 1e-4 --text --max-stored-scores 1000000 %s %s > %s 2>/dev/null",
                 fimo_bin, meme, fa, out)
  cat("  ", cmd, "\n")
  t0 <- Sys.time()
  rc <- system(cmd)
  cat("   exit=", rc, " elapsed=", format(Sys.time() - t0), "\n")
  invisible(rc)
}
cat("\n--- Running FIMO on REF ---\n"); run_fimo(meme_file, ref_fa, ref_out)
cat("\n--- Running FIMO on ALT ---\n"); run_fimo(meme_file, alt_fa, alt_out)

read_fimo_tsv <- function(p) {
  if (!file.exists(p) || file.info(p)$size < 2) return(data.table())
  # FIMO --text writes a header line. Read and find header.
  lines <- readLines(p)
  # Skip leading comment / blank lines and parse header
  hdr_idx <- grep("^motif_id\\t|^motif_id ", lines)
  if (length(hdr_idx) == 0) {
    # try fread default
    dt <- tryCatch(fread(p, header = TRUE), error = function(e) data.table())
    return(dt)
  }
  body <- lines[hdr_idx[1]:length(lines)]
  body <- body[body != ""]
  dt <- fread(text = paste(body, collapse = "\n"), header = TRUE)
  dt
}

fimo_ref <- read_fimo_tsv(ref_out)
fimo_alt <- read_fimo_tsv(alt_out)
cat("FIMO REF hits:", nrow(fimo_ref), "  ALT hits:", nrow(fimo_alt), "\n")

# Standardise columns. FIMO --text typical cols:
# motif_id  motif_alt_id  sequence_name  start  stop  strand  score  p-value  q-value  matched_sequence
fix_cols <- function(dt) {
  if (nrow(dt) == 0) return(dt)
  cn <- colnames(dt)
  setnames(dt, cn[1], "motif_id")
  if ("sequence_name" %in% cn) {
    # ok
  } else if ("sequence name" %in% cn) {
    setnames(dt, "sequence name", "sequence_name")
  } else {
    setnames(dt, cn[3], "sequence_name")
  }
  if ("p-value" %in% colnames(dt)) setnames(dt, "p-value", "pvalue")
  if (!"score" %in% colnames(dt)) {
    # find numeric col 7
    setnames(dt, colnames(dt)[7], "score")
  }
  dt
}
fimo_ref <- fix_cols(fimo_ref)
fimo_alt <- fix_cols(fimo_alt)

# Best (max score) per (sequence_name, motif_id)
best_per <- function(dt, side) {
  if (nrow(dt) == 0) return(data.table(sequence_name = character(), motif_id = character()))
  dt[, score := as.numeric(score)]
  out <- dt[, .(score_best = max(score, na.rm = TRUE),
                n_hits = .N,
                pval_best = min(as.numeric(pvalue), na.rm = TRUE)),
            by = .(sequence_name, motif_id)]
  setnames(out, c("score_best", "n_hits", "pval_best"),
           paste0(side, "_", c("score", "nhits", "pval")))
  out
}
ref_best <- best_per(fimo_ref, "ref")
alt_best <- best_per(fimo_alt, "alt")

cmp <- merge(ref_best, alt_best, by = c("sequence_name", "motif_id"),
             all = TRUE)
cmp[is.na(ref_score), `:=`(ref_score = 0, ref_nhits = 0L, ref_pval = NA_real_)]
cmp[is.na(alt_score), `:=`(alt_score = 0, alt_nhits = 0L, alt_pval = NA_real_)]
cmp[, fimo_disrupted := ref_nhits > 0 & alt_nhits == 0]
cmp[, fimo_created   := ref_nhits == 0 & alt_nhits > 0]
cmp[, fimo_changed   := fimo_disrupted | fimo_created]
cmp[, score_delta    := ref_score - alt_score]
cat("FIMO compare rows:", nrow(cmp),
    " disrupted:", sum(cmp$fimo_disrupted),
    " created:",   sum(cmp$fimo_created), "\n")

# ── 6. Join FIMO compare onto Tier 1 / Tier 2 motifbreakR rows ──────────────
# Each motifbreakR row is (SNP_id, tf_name). FIMO motif_id = fimo_id with prefix.
# Need to map: for each (SNP_id, tf_name), get the union of fimo_id with matching tf.
cmp[, fimo_tf := toupper(sub("^M\\d+_", "", motif_id))]
# Map FIMO sequence_name (our surrogate fasta_id) back to SNP_id
cmp[, sequence_name := as.character(sequence_name)]
cmp <- merge(cmp, seq_lookup, by.x = "sequence_name", by.y = "fasta_id", all.x = TRUE)
# Build long pair join: per (SNP_id, tf_upper) -> aggregate FIMO best across all motifs for that TF.
agg_cmp <- cmp[, .(
  ref_score = max(ref_score, na.rm = TRUE),
  alt_score = max(alt_score, na.rm = TRUE),
  ref_nhits = max(ref_nhits, na.rm = TRUE),
  alt_nhits = max(alt_nhits, na.rm = TRUE),
  any_disrupted = any(fimo_disrupted),
  any_created   = any(fimo_created),
  any_changed   = any(fimo_changed),
  n_motifs_for_tf = .N
), by = .(SNP_id, tf_upper = fimo_tf)]
agg_cmp[, SNP_id := as.character(SNP_id)]
# FIX (review F110): orient score_delta as ALT-REF to match motifbreakR alleleDiff = scoreAlt-scoreRef.
# Previously ref-alt, so the sign(alleleDiff)==sign(score_delta) test below measured ANTI-agreement
# (~1.5% TRUE on disk) when the true motifbreakR<->FIMO direction agreement is ~98.5%.
agg_cmp[, score_delta := alt_score - ref_score]
agg_cmp[, fimo_validated := ref_nhits > 0 | alt_nhits > 0]

annotate_tier <- function(tdt) {
  if (nrow(tdt) == 0) return(tdt)
  tdt[, tf_upper := toupper(tf_name)]
  tdt[, SNP_id := as.character(SNP_id)]
  out <- merge(tdt, agg_cmp, by = c("SNP_id", "tf_upper"), all.x = TRUE)
  out[, fimo_dir_concord := fifelse(
    !is.na(score_delta) & fimo_validated,
    sign(alleleDiff) == sign(score_delta),
    NA
  )]
  out
}
tier1_t <- annotate_tier(as.data.table(tier1))
tier2_t <- annotate_tier(as.data.table(tier2))

fwrite(tier1_t, file.path(OUT_DIR, "fimo_comparison_full.csv"))
fwrite(tier2_t, file.path(OUT_DIR, "fimo_comparison_extended.csv"))
cat("Wrote fimo_comparison_full.csv  (", nrow(tier1_t), " rows)\n")
cat("Wrote fimo_comparison_extended.csv (", nrow(tier2_t), " rows)\n")

# ── 7. Stratified concordance summary ───────────────────────────────────────
all_t <- rbind(tier1_t, tier2_t, fill = TRUE)
bins <- c(0, 0.3, 0.5, 0.7, 1.0, 1.5, Inf)
labs <- c("0-0.3", "0.3-0.5", "0.5-0.7", "0.7-1.0", "1.0-1.5", ">1.5")
all_t[, ad_bin := cut(abs(alleleDiff), breaks = bins, labels = labs,
                       include.lowest = TRUE, right = FALSE)]

summarise_block <- function(dt, label) {
  s <- dt[, .(
    n_tested              = .N,
    n_fimo_validated      = sum(fimo_validated, na.rm = TRUE),
    n_disrupted_or_created = sum(any_changed,   na.rm = TRUE),
    n_dir_concord         = sum(fimo_dir_concord, na.rm = TRUE),
    n_dir_eval            = sum(!is.na(fimo_dir_concord))
  ), by = .(ad_bin)]
  s[, pct_validated := round(100 * n_fimo_validated / n_tested, 1)]
  s[, pct_changed   := round(100 * n_disrupted_or_created / n_tested, 1)]
  s[, pct_dir_concord := round(100 * n_dir_concord / pmax(n_dir_eval, 1), 1)]
  s[, group := label]
  setcolorder(s, c("group", "ad_bin"))
  s
}
summ <- rbind(
  summarise_block(tier1_t, "tier1_strong"),
  summarise_block(tier2_t, "tier2_alleleDiff_driven"),
  summarise_block(all_t,   "combined")
)
fwrite(summ, file.path(OUT_DIR, "fimo_concordance_by_stratum.csv"))
cat("\nStratified concordance summary:\n")
print(summ)

# ── 8. Headline + disease-regulon TF survival ───────────────────────────────
cat("\n=== Headline ===\n")
t1_validated <- sum(tier1_t$fimo_validated, na.rm = TRUE)
t1_changed   <- sum(tier1_t$any_changed,    na.rm = TRUE)
t1_concord_n <- sum(tier1_t$fimo_dir_concord, na.rm = TRUE)
t1_concord_d <- sum(!is.na(tier1_t$fimo_dir_concord))
cat(sprintf("Tier 1 (strong, n=%d): %d FIMO-validated (%.1f%%); %d disrupted-or-created (%.1f%%); direction-concordance %d/%d (%.1f%%)\n",
            nrow(tier1_t), t1_validated, 100*t1_validated/nrow(tier1_t),
            t1_changed,    100*t1_changed/nrow(tier1_t),
            t1_concord_n,  t1_concord_d,
            ifelse(t1_concord_d > 0, 100*t1_concord_n/t1_concord_d, NA_real_)))

# Disease-regulon TF survival (the 12 GWAS-disrupted disease-regulon TFs)
regulon_file <- file.path(ATAC_DIR, "scenic_plus/disease_regulons.csv")
if (file.exists(regulon_file)) {
  regulons <- fread(regulon_file)
  disease_tfs <- toupper(regulons$tf_name)
  dr_t1 <- tier1_t[toupper(tf_name) %in% disease_tfs]
  dr_summary <- dr_t1[, .(
    n_pairs              = .N,
    n_fimo_validated     = sum(fimo_validated, na.rm = TRUE),
    n_changed            = sum(any_changed,    na.rm = TRUE),
    n_dir_concord        = sum(fimo_dir_concord, na.rm = TRUE),
    n_dir_eval           = sum(!is.na(fimo_dir_concord))
  ), by = .(tf = toupper(tf_name))]
  dr_summary[, pct_validated := round(100 * n_fimo_validated / n_pairs, 1)]
  dr_summary[, pct_dir_concord := round(100 * n_dir_concord / pmax(n_dir_eval, 1), 1)]
  dr_summary <- dr_summary[order(-n_pairs)]
  fwrite(dr_summary, file.path(OUT_DIR, "fimo_disease_regulon_tf_survival.csv"))
  cat("\nDisease-regulon TF FIMO survival (Tier 1):\n")
  print(dr_summary)
  cat(sprintf("\nDisease-regulon TFs in Tier 1: %d unique; %d with >=1 FIMO-validated pair\n",
              nrow(dr_summary), sum(dr_summary$n_fimo_validated > 0)))
}

# ── 9. Cleanup ──────────────────────────────────────────────────────────────
unlink(TMP_DIR, recursive = TRUE)
cat("\nDone. Outputs in", OUT_DIR, "\n")
cat("============================================================\n")
