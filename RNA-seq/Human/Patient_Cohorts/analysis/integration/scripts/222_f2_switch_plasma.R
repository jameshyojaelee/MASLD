#!/usr/bin/env Rscript
# 222_f2_switch_plasma.R
# LEGACY SCRIPT — "F2 switch" framing retired 2026-05-09.
# Current thesis: multi-step cellular cascade across 4 CRN transitions (see paper_outline.md).
# Script name retained for provenance; outputs are audit trail only.
#
# Original purpose: F2 Metabolic-to-Inflammatory Switch Detection in Plasma Proteomics.
# Tests whether the tissue-level F2→F3 consecutive-contrast signature is detectable
# in circulating plasma proteins (Olink, 1,460 proteins × 218 subjects).
#
# Defines three stage-transition signatures from fibrosis consecutive contrasts:
#   Onset : genes UP in F1_vs_F0 OR F2_vs_F1  (early fibrosis onset)
#   Switch: genes UP in F3_vs_F2               (F2→F3 transition; legacy "switch" label)
#   Late  : genes UP in F4_vs_F3               (end-stage remodeling)
#
# Scores plasma subjects by mean NPX across signature proteins, then tests
# whether switch_score separates F0-2 vs F3 vs F4 subjects.
#
# Inputs:
#   fibrosis_consecutive_dream.csv   — 614 K rows, padj + logFC per contrast
#   dream_results_ashr.csv           — Ensembl → symbol mapping (34 K genes)
#   olink.qc.finished.mendeley.data.txt — 1,460 proteins × 218 subjects (NPX)
#   gse276114_disease_metadata.csv   — 177 liver samples with disease_group
#                                      (F0-2, F3, F4) and sample_number
#
# Outputs (→ results/multiprogram/):
#   plasma_switch_scores.csv         — per-subject scores + stage
#   plasma_switch_signatures.csv     — gene lists with tissue effect + Olink flag
#   plasma_switch_statistics.csv     — KW p-values, pairwise Wilcoxon, effect sizes
#
# SLURM: --partition=cpu --cpus-per-task=4 --mem=32G --time=48:00:00
# Env:   micromamba activate rnaseq

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(42)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
DS_DIR <- file.path(INT,  "results/disease_signatures")
RES    <- file.path(INT,  "results/integration")
OLINK_DIR <- file.path(BASE, "Analysis/Proteomics/data/olink_plasma")
PROT_DIR  <- file.path(BASE, "Analysis/Proteomics/results")
OUTDIR    <- file.path(INT,  "results/multiprogram")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 222: F2 Switch Detection in Plasma Proteomics ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ── STEP 1: Load consecutive contrast results ─────────────────────────────────
cat("=== STEP 1: Load fibrosis consecutive contrasts ===\n")

consec_file <- file.path(DS_DIR, "fibrosis_consecutive_dream.csv")
stopifnot(file.exists(consec_file))
consec <- fread(consec_file, na.strings = c("", "NA"))
cat("  Loaded:", nrow(consec), "rows x", ncol(consec), "cols\n")
cat("  Contrasts:", paste(consec[, .N, by = contrast][, paste0(contrast, "=", N)], collapse = ", "), "\n")

# Strip Ensembl version suffix (ENSG00000XXX.YY → ENSG00000XXX)
consec[, gene_id := sub("\\.\\d+$", "", gene)]

# Threshold: padj < 0.05, logFC > 0.2 (upregulated)
consec[, sig_up := !is.na(padj) & !is.na(logFC) &
                    as.numeric(padj) < 0.05 &
                    as.numeric(logFC) > 0.2]

cat("  Sig-up per contrast:\n")
print(consec[, .(n_sig_up = sum(sig_up, na.rm = TRUE)), by = contrast])

# ── STEP 2: Build Ensembl → symbol map ───────────────────────────────────────
cat("\n=== STEP 2: Build Ensembl → symbol map ===\n")

ashr_file <- file.path(RES, "dream_results_ashr.csv")
stopifnot(file.exists(ashr_file))
id_map <- fread(ashr_file, select = c("gene", "symbol"))
id_map[, gene_id := sub("\\.\\d+$", "", gene)]
id_map <- id_map[symbol != "" & !is.na(symbol)]
# Deduplicate: keep one symbol per stripped Ensembl ID
id_map <- unique(id_map, by = "gene_id")
cat("  ID map:", nrow(id_map), "Ensembl → symbol pairs\n")

# Merge symbols into consecutive contrast table
consec <- merge(consec, id_map[, .(gene_id, symbol)], by = "gene_id", all.x = TRUE)
cat("  Genes with symbol:", consec[!is.na(symbol), uniqueN(gene_id)],
    "of", consec[, uniqueN(gene_id)], "\n")

# ── STEP 3: Define onset / switch / late signatures ───────────────────────────
cat("\n=== STEP 3: Define tissue signatures ===\n")

# For each gene, take the row with the smallest padj per contrast to handle
# any per-gene duplicate rows (multiple probe/version collisions).
best <- consec[!is.na(symbol)][
  , .SD[which.min(as.numeric(padj))], by = .(gene_id, symbol, contrast)
]

# Onset: UP in F1_vs_F0 OR F2_vs_F1
onset_genes <- best[contrast %in% c("F1_vs_F0", "F2_vs_F1") & sig_up == TRUE,
                    unique(symbol)]

# Switch: UP in F3_vs_F2
switch_genes <- best[contrast == "F3_vs_F2" & sig_up == TRUE, unique(symbol)]

# Late: UP in F4_vs_F3
late_genes <- best[contrast == "F4_vs_F3" & sig_up == TRUE, unique(symbol)]

cat("  Onset genes (F0→F1 or F1→F2 UP):", length(onset_genes), "\n")
cat("  Switch genes (F2→F3 UP):", length(switch_genes), "\n")
cat("  Late genes (F3→F4 UP):", length(late_genes), "\n")

# Build annotation table for all three signature sets
build_sig_table <- function(sym, sig_label, best_dt) {
  # For each symbol, pick representative row (switch contrast first, else minimum padj)
  target_contrast <- switch(sig_label,
    onset  = c("F1_vs_F0", "F2_vs_F1"),
    switch = "F3_vs_F2",
    late   = "F4_vs_F3"
  )
  rows <- best_dt[symbol %in% sym & contrast %in% target_contrast]
  rows <- rows[, .SD[which.min(as.numeric(padj))], by = symbol]
  data.table(
    gene_symbol  = rows$symbol,
    signature    = sig_label,
    tissue_logFC = as.numeric(rows$logFC),
    tissue_padj  = as.numeric(rows$padj),
    contrast     = rows$contrast
  )
}

sig_table <- rbindlist(list(
  build_sig_table(onset_genes,  "onset",  best),
  build_sig_table(switch_genes, "switch", best),
  build_sig_table(late_genes,   "late",   best)
))
cat("  Signature table:", nrow(sig_table), "rows\n")

# ── STEP 4: Load Olink plasma data ────────────────────────────────────────────
cat("\n=== STEP 4: Load Olink plasma data ===\n")

olink_file <- file.path(OLINK_DIR, "olink.qc.finished.mendeley.data.txt")
stopifnot(file.exists(olink_file))
olink_raw <- fread(olink_file, header = TRUE, sep = "\t")
cat("  Olink raw:", nrow(olink_raw), "proteins x", ncol(olink_raw) - 1, "subjects\n")

olink_proteins <- olink_raw$Assay
olink_mat <- as.matrix(olink_raw[, -1, with = FALSE])
rownames(olink_mat) <- olink_proteins

# Transpose to subjects × proteins
npx <- t(olink_mat)
cat("  NPX matrix:", nrow(npx), "subjects x", ncol(npx), "proteins\n")

# Extract subject number from rownames ("Subject N" → N)
subject_nums <- as.integer(sub("^Subject\\s+", "", rownames(npx)))

# ── STEP 5: Load fibrosis stage metadata and match to Olink subjects ──────────
cat("\n=== STEP 5: Match subjects to fibrosis stage metadata ===\n")

# Search for metadata file; fall back to searching Analysis/Proteomics/ if needed
meta_file <- file.path(PROT_DIR, "gse276114_disease_metadata.csv")
if (!file.exists(meta_file)) {
  alt_files <- list.files(file.path(BASE, "Analysis/Proteomics"),
                          pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  fib_candidates <- alt_files[grep("disease|fibrosis|stage|meta", basename(alt_files),
                                   ignore.case = TRUE)]
  if (length(fib_candidates) > 0) {
    meta_file <- fib_candidates[1]
    cat("  Falling back to:", meta_file, "\n")
  } else {
    stop("No fibrosis metadata file found. Searched: ", PROT_DIR)
  }
}

meta <- fread(meta_file)
cat("  Metadata:", nrow(meta), "samples, columns:", paste(names(meta), collapse = ", "), "\n")
cat("  Disease groups:", paste(meta[, .N, by = disease_group][, paste0(disease_group, "=", N)], collapse = ", "), "\n")

# Map Olink "Subject N" → metadata via sample_number = N
label_dt <- data.table(
  subject_idx  = seq_len(nrow(npx)),
  subject_num  = subject_nums,
  subject_name = rownames(npx)
)
label_dt <- merge(label_dt, meta[, .(sample_number, disease_group, disease)],
                  by.x = "subject_num", by.y = "sample_number",
                  all.x = TRUE)

# Restrict to subjects with fibrosis staging metadata
label_dt_matched <- label_dt[!is.na(disease_group)]
cat("  Olink subjects with fibrosis metadata:", nrow(label_dt_matched), "\n")
cat("  Stage breakdown:", paste(label_dt_matched[, .N, by = disease_group][
                                  order(disease_group),
                                  paste0(disease_group, "=", N)],
                                collapse = ", "), "\n")

# Subset NPX to matched subjects
npx_matched <- npx[label_dt_matched$subject_idx, , drop = FALSE]

# ── STEP 6: Intersect signatures with Olink proteins ─────────────────────────
cat("\n=== STEP 6: Intersect tissue signatures with Olink ===\n")

olink_protein_names <- colnames(npx_matched)

n_onset_olink  <- sum(onset_genes  %in% olink_protein_names)
n_switch_olink <- sum(switch_genes %in% olink_protein_names)
n_late_olink   <- sum(late_genes   %in% olink_protein_names)

cat(sprintf("  Onset:  %d tissue genes → %d in Olink (%.1f%%)\n",
            length(onset_genes),  n_onset_olink,
            100 * n_onset_olink / max(length(onset_genes), 1)))
cat(sprintf("  Switch: %d tissue genes → %d in Olink (%.1f%%)\n",
            length(switch_genes), n_switch_olink,
            100 * n_switch_olink / max(length(switch_genes), 1)))
cat(sprintf("  Late:   %d tissue genes → %d in Olink (%.1f%%)\n",
            length(late_genes),   n_late_olink,
            100 * n_late_olink / max(length(late_genes), 1)))

# Annotate signature table with in_olink flag
sig_table[, in_olink := gene_symbol %in% olink_protein_names]

# ── STEP 7: Score plasma subjects ────────────────────────────────────────────
cat("\n=== STEP 7: Score plasma subjects ===\n")

score_subjects <- function(npx_mat, gene_set, set_name) {
  avail <- intersect(gene_set, colnames(npx_mat))
  if (length(avail) == 0) {
    cat(sprintf("  WARNING: no %s genes found in Olink — returning NA scores\n", set_name))
    return(rep(NA_real_, nrow(npx_mat)))
  }
  cat(sprintf("  %s: %d/%d genes available\n", set_name, length(avail), length(gene_set)))
  sub_mat <- npx_mat[, avail, drop = FALSE]
  # Impute remaining NA values with column median before scoring
  for (j in seq_len(ncol(sub_mat))) {
    na_idx <- is.na(sub_mat[, j])
    if (any(na_idx)) {
      sub_mat[na_idx, j] <- median(sub_mat[, j], na.rm = TRUE)
    }
  }
  rowMeans(sub_mat, na.rm = TRUE)
}

onset_score  <- score_subjects(npx_matched, onset_genes,  "onset")
switch_score <- score_subjects(npx_matched, switch_genes, "switch")
late_score   <- score_subjects(npx_matched, late_genes,   "late")

# switch_ratio: switch relative to onset (guards against global expression shifts)
switch_ratio <- switch_score / (onset_score + 0.01)

scores_dt <- data.table(
  subject_id     = label_dt_matched$subject_name,
  subject_num    = label_dt_matched$subject_num,
  fibrosis_stage = label_dt_matched$disease_group,
  disease        = label_dt_matched$disease,
  onset_score    = onset_score,
  switch_score   = switch_score,
  late_score     = late_score,
  switch_ratio   = switch_ratio
)

cat("  Score summary by stage:\n")
print(scores_dt[, .(
  n = .N,
  onset_median  = round(median(onset_score,  na.rm = TRUE), 3),
  switch_median = round(median(switch_score, na.rm = TRUE), 3),
  late_median   = round(median(late_score,   na.rm = TRUE), 3),
  ratio_median  = round(median(switch_ratio, na.rm = TRUE), 3)
), by = fibrosis_stage][order(fibrosis_stage)])

# ── STEP 8: Statistical tests ─────────────────────────────────────────────────
cat("\n=== STEP 8: Statistical tests ===\n")

stage_factor <- factor(scores_dt$fibrosis_stage,
                       levels = c("F0-2", "F3", "F4"))
scores_dt[, stage_f := stage_factor]

# Helper: Kruskal-Wallis + epsilon-squared effect size
run_kruskal <- function(score_vec, stage_vec, score_name) {
  kw <- kruskal.test(score_vec ~ stage_vec)
  # Epsilon-squared: KW H / (n - 1)
  n_total <- length(score_vec)
  eps2 <- kw$statistic / (n_total - 1)
  data.table(
    score       = score_name,
    test        = "Kruskal-Wallis",
    comparison  = "F0-2 vs F3 vs F4",
    statistic   = round(kw$statistic, 4),
    df          = kw$parameter,
    p_value     = kw$p.value,
    effect_size = round(eps2, 4),
    effect_name = "epsilon_squared"
  )
}

# Helper: pairwise Wilcoxon (F0-2 vs F3+) + r effect size
run_wilcox_pair <- function(score_vec, stage_vec, g1_labels, g2_labels, score_name) {
  g1 <- score_vec[stage_vec %in% g1_labels]
  g2 <- score_vec[stage_vec %in% g2_labels]
  g1 <- g1[!is.na(g1)]
  g2 <- g2[!is.na(g2)]
  wt <- wilcox.test(g1, g2, exact = FALSE)
  # r = Z / sqrt(n)
  z_val <- qnorm(wt$p.value / 2) * sign(median(g2, na.rm = TRUE) - median(g1, na.rm = TRUE))
  r_eff <- abs(z_val) / sqrt(length(g1) + length(g2))
  data.table(
    score       = score_name,
    test        = "Wilcoxon rank-sum",
    comparison  = paste0(paste(g1_labels, collapse = "/"), " vs ", paste(g2_labels, collapse = "/")),
    statistic   = round(wt$statistic, 4),
    df          = NA_real_,
    p_value     = wt$p.value,
    effect_size = round(r_eff, 4),
    effect_name = "r"
  )
}

stat_rows <- list()

for (sc_name in c("onset_score", "switch_score", "late_score", "switch_ratio")) {
  sc_vec <- scores_dt[[sc_name]]
  valid  <- !is.na(sc_vec) & !is.na(stage_factor)
  sv     <- sc_vec[valid]
  sf     <- stage_factor[valid]

  # Kruskal-Wallis across all three groups
  stat_rows[[length(stat_rows) + 1]] <- run_kruskal(sv, sf, sc_name)

  # Pairwise: F0-2 vs F3
  stat_rows[[length(stat_rows) + 1]] <- run_wilcox_pair(sv, sf, "F0-2", "F3", sc_name)

  # Pairwise: F0-2 vs F4
  stat_rows[[length(stat_rows) + 1]] <- run_wilcox_pair(sv, sf, "F0-2", "F4", sc_name)

  # Pairwise: F0-2 vs F3+
  stat_rows[[length(stat_rows) + 1]] <- run_wilcox_pair(sv, sf,
                                                         "F0-2", c("F3", "F4"), sc_name)

  # Pairwise: F3 vs F4
  stat_rows[[length(stat_rows) + 1]] <- run_wilcox_pair(sv, sf, "F3", "F4", sc_name)
}

stats_dt <- rbindlist(stat_rows)

# BH correction across all tests
stats_dt[, padj := p.adjust(p_value, method = "BH")]

# ── STEP 8b: Size-matched random null for switch panel (T0.12 / H5-inv HOLD) ──
# Added 2026-04-22 — addresses pre-registered HOLD H5-inv:
#   "Is the switch panel's KW p size-specific or generic to any 180-protein draw?"
# Constructs 3 null distributions (1000 draws each) over different universes:
#   full  — all 1,460 Olink proteins (is ANY random 180-protein panel significant?)
#   deg   — tissue-DEG-detected-in-Olink (are tissue DEGs special beyond DEG status?)
#   nondeg— non-tissue-DEG Olink proteins (are tissue DEGs needed at all?)
# CRITICAL: scoring function is byte-identical to score_subjects() used for the
# observed panel — same NA imputation (per-column median), same rowMeans, same
# subject ordering (npx_matched is the filtered matrix already).
cat("\n=== STEP 8b: Size-matched random null (T0.12 / H5-inv HOLD) ===\n")

set.seed(42)
N_DRAWS    <- 1000L
# Panel size = number of SWITCH proteins actually scored (Olink-detected subset).
# This is the panel whose KW p=2.06e-5 is the pre-registered claim.
PANEL_SIZE <- sum(switch_genes %in% colnames(npx_matched))
cat(sprintf("  Panel size (switch genes in Olink): %d\n", PANEL_SIZE))

# Universes (over Olink columns, not tissue genes)
olink_all    <- colnames(npx_matched)                          # 1,460 proteins
switch_olink <- intersect(switch_genes, olink_all)             # the observed panel
# Tissue DEGs = all genes significantly up in ANY consecutive fibrosis contrast
tissue_degs_all <- unique(c(onset_genes, switch_genes, late_genes))
deg_olink    <- intersect(tissue_degs_all, olink_all)
nondeg_olink <- setdiff(olink_all, tissue_degs_all)
# "Full universe" draws should NOT include the observed panel itself
full_pool    <- setdiff(olink_all, switch_olink)
# DEG pool: use ALL tissue DEGs in Olink (onset ∪ switch ∪ late). We do NOT
# exclude the switch panel because deg_olink \ switch_olink = 122 < 180 = k.
# This makes the DEG null strictly harder on the observed panel (it can
# re-sample panel members) — a conservative test of "are tissue DEGs special
# beyond DEG status?" Duplicates within a draw are not possible (sample()
# without replacement); but switch members CAN appear in a random DEG draw.
deg_pool     <- deg_olink

cat(sprintf("  Universe sizes: full=%d (excl. panel), deg=%d (tissue DEGs, incl. panel), nondeg=%d\n",
            length(full_pool), length(deg_pool), length(nondeg_olink)))

# Mirror score_subjects() exactly for null draws (same NA handling, same rowMeans)
score_draw <- function(npx_mat, gene_set) {
  sub_mat <- npx_mat[, gene_set, drop = FALSE]
  for (j in seq_len(ncol(sub_mat))) {
    na_idx <- is.na(sub_mat[, j])
    if (any(na_idx)) {
      sub_mat[na_idx, j] <- median(sub_mat[, j], na.rm = TRUE)
    }
  }
  rowMeans(sub_mat, na.rm = TRUE)
}

# Valid mask for KW — identical to the observed-panel branch above
obs_valid <- !is.na(switch_score) & !is.na(stage_factor)
sf_valid  <- stage_factor[obs_valid]
observed_kw <- unname(kruskal.test(switch_score[obs_valid] ~ sf_valid)$statistic)
cat(sprintf("  Observed KW H = %.4f\n", observed_kw))

null_kw_stats <- function(pool, k, n_draws) {
  if (length(pool) < k) {
    warning(sprintf("Pool size %d < panel size %d; skipping", length(pool), k))
    return(rep(NA_real_, n_draws))
  }
  vapply(seq_len(n_draws), function(i) {
    genes <- sample(pool, k)
    score <- score_draw(npx_matched, genes)
    sv    <- score[obs_valid]
    unname(kruskal.test(sv ~ sf_valid)$statistic)
  }, numeric(1))
}

cat(sprintf("  Running %d draws x 3 universes...\n", N_DRAWS))
t0 <- Sys.time()
null_full   <- null_kw_stats(full_pool,    PANEL_SIZE, N_DRAWS)
null_deg    <- null_kw_stats(deg_pool,     PANEL_SIZE, N_DRAWS)
null_nondeg <- null_kw_stats(nondeg_olink, PANEL_SIZE, N_DRAWS)
cat(sprintf("  Null draws complete (%.1f min)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

emp_p_full   <- mean(null_full   >= observed_kw, na.rm = TRUE)
emp_p_deg    <- mean(null_deg    >= observed_kw, na.rm = TRUE)
emp_p_nondeg <- mean(null_nondeg >= observed_kw, na.rm = TRUE)

cat(sprintf("  emp_p (full universe, 1,460 Olink): %.4f  [null median H = %.2f]\n",
            emp_p_full,   median(null_full,   na.rm = TRUE)))
cat(sprintf("  emp_p (tissue-DEG subset):          %.4f  [null median H = %.2f]\n",
            emp_p_deg,    median(null_deg,    na.rm = TRUE)))
cat(sprintf("  emp_p (non-DEG subset):             %.4f  [null median H = %.2f]\n",
            emp_p_nondeg, median(null_nondeg, na.rm = TRUE)))

null_out <- data.frame(
  observed_kw         = observed_kw,
  panel_size          = PANEL_SIZE,
  n_draws             = N_DRAWS,
  emp_p_full_universe = emp_p_full,
  emp_p_deg_subset    = emp_p_deg,
  emp_p_nondeg_subset = emp_p_nondeg,
  null_full_median    = median(null_full,   na.rm = TRUE),
  null_full_q95       = quantile(null_full,   0.95, na.rm = TRUE),
  null_deg_median     = median(null_deg,    na.rm = TRUE),
  null_deg_q95        = quantile(null_deg,    0.95, na.rm = TRUE),
  null_nondeg_median  = median(null_nondeg, na.rm = TRUE),
  null_nondeg_q95     = quantile(null_nondeg, 0.95, na.rm = TRUE),
  n_subjects_scored   = sum(obs_valid),
  universe_full_size  = length(full_pool),
  universe_deg_size   = length(deg_pool),
  universe_nondeg_size= length(nondeg_olink)
)
fwrite(null_out, file.path(OUTDIR, "plasma_switch_random_null.csv"))
cat("  Saved plasma_switch_random_null.csv\n")

saveRDS(list(
  null_full   = null_full,
  null_deg    = null_deg,
  null_nondeg = null_nondeg,
  observed_kw = observed_kw,
  panel_size  = PANEL_SIZE,
  n_draws     = N_DRAWS
), file.path(OUTDIR, "plasma_switch_random_null_draws.rds"))
cat("  Saved plasma_switch_random_null_draws.rds (full null vectors for figure)\n")
# ── end STEP 8b ───────────────────────────────────────────────────────────────

cat("  Key results (switch_score):\n")
print(stats_dt[score == "switch_score",
               .(comparison, p_value = signif(p_value, 3),
                 padj = signif(padj, 3), effect_size)])

cat("\n  Key results (switch_ratio):\n")
print(stats_dt[score == "switch_ratio",
               .(comparison, p_value = signif(p_value, 3),
                 padj = signif(padj, 3), effect_size)])

# ── STEP 9: Save outputs ──────────────────────────────────────────────────────
cat("\n=== STEP 9: Save outputs ===\n")

# Remove helper column before saving
scores_dt[, stage_f := NULL]

fwrite(scores_dt,
       file.path(OUTDIR, "plasma_switch_scores.csv"))
cat("  Saved plasma_switch_scores.csv:", nrow(scores_dt), "subjects\n")

fwrite(sig_table,
       file.path(OUTDIR, "plasma_switch_signatures.csv"))
cat("  Saved plasma_switch_signatures.csv:", nrow(sig_table), "genes\n")

fwrite(stats_dt,
       file.path(OUTDIR, "plasma_switch_statistics.csv"))
cat("  Saved plasma_switch_statistics.csv:", nrow(stats_dt), "tests\n")

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n=== Summary ===\n")
cat(sprintf("  Onset signature:  %d tissue genes, %d in Olink\n",
            length(onset_genes), n_onset_olink))
cat(sprintf("  Switch signature: %d tissue genes, %d in Olink\n",
            length(switch_genes), n_switch_olink))
cat(sprintf("  Late signature:   %d tissue genes, %d in Olink\n",
            length(late_genes), n_late_olink))
cat(sprintf("  Plasma subjects scored: %d (F0-2=%d, F3=%d, F4=%d)\n",
            nrow(scores_dt),
            sum(scores_dt$fibrosis_stage == "F0-2"),
            sum(scores_dt$fibrosis_stage == "F3"),
            sum(scores_dt$fibrosis_stage == "F4")))

kw_switch <- stats_dt[score == "switch_score" & test == "Kruskal-Wallis"]
cat(sprintf("  Switch score KW p = %.4g, eps2 = %.3f\n",
            kw_switch$p_value, kw_switch$effect_size))

wil_switch <- stats_dt[score == "switch_score" & comparison == "F0-2 vs F3/F4"]
cat(sprintf("  Switch score F0-2 vs F3+: W = %.0f, p = %.4g, r = %.3f\n",
            wil_switch$statistic, wil_switch$p_value, wil_switch$effect_size))

cat(sprintf("  [H5-inv] Observed KW H = %.3f vs random-draw null:\n", observed_kw))
cat(sprintf("           emp_p (full Olink universe) = %.4f (%d draws, panel_size=%d)\n",
            emp_p_full, N_DRAWS, PANEL_SIZE))
cat(sprintf("           emp_p (tissue-DEG subset)   = %.4f\n", emp_p_deg))
cat(sprintf("           emp_p (non-DEG subset)      = %.4f\n", emp_p_nondeg))

cat("\nFinished:", as.character(Sys.time()), "\n")
