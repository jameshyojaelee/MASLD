#!/usr/bin/env Rscript
# 56d_motif_quality_filter.R
# Track A3: JASPAR A/B-grade motif quality filter
#
# Grades each JASPAR motif used by Script 56 by experimental data type +
# cross-validates against HOCOMOCO consensus PWM. Produces motif quality
# tiers (A / B / C) and re-counts disease-regulon disruptions per tier.
#
# Output:
#   results/gwas_atac/motif_quality_tier.csv
#   results/gwas_atac/motif_tier_comparison.csv
#
# Grading:
#   A = ChIP-seq/ChIP-exo data AND HOCOMOCO Pearson cor >0.7
#   B = ChIP-seq/ChIP-exo data only OR (HOCOMOCO Pearson cor >0.7) OR
#       (SELEX-family AND information content >12 bits)
#   C = otherwise (low-confidence inferred / sparse / unsupported)
#
# Env: motifbreakr
#
# NEW script. Does NOT modify Script 56.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(MotifDb)
  library(TFBSTools)
  library(httr)
  library(jsonlite)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
OUT_DIR <- file.path(FM_DIR, "results/gwas_atac")
CACHE_DIR <- file.path(OUT_DIR, "jaspar_api_cache")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("56d_motif_quality_filter.R\n")
cat("JASPAR A/B/C-grade motif quality classifier\n")
cat("============================================================\n\n")

# Reproduce Script 56's motif query exactly
jaspar_motifs <- query(MotifDb, andStrings = c("Hsapiens"),
                       orStrings = c("jaspar2022", "jaspar2024"))
cat("JASPAR Hsapiens motifs:", length(jaspar_motifs), "\n")

# HOCOMOCO Hsapiens for consensus check
hocomoco_motifs <- query(MotifDb, andStrings = c("Hsapiens"),
                         orStrings = c("HOCOMOCO"))
cat("HOCOMOCO Hsapiens motifs:", length(hocomoco_motifs), "\n\n")

# ── 1. Build metadata table from MotifDb ────────────────────────────────────
jdf <- as.data.frame(values(jaspar_motifs))
jdf$motif_id <- names(jaspar_motifs)
jdf$matrix_id <- jdf$providerName  # MA0030.1 etc.

# ── 2. Compute IC per motif via TFBSTools ───────────────────────────────────
cat("Computing information content per PWM...\n")
ic_vec <- numeric(length(jaspar_motifs))
for (i in seq_along(jaspar_motifs)) {
  m <- as.matrix(jaspar_motifs[[i]])
  # MotifDb PWMs are normalized probabilities. Convert to PFM (synthetic counts)
  # then ICM via TFBSTools.
  counts <- round(m * 100)
  if (any(is.na(counts)) || ncol(counts) < 4) {
    ic_vec[i] <- NA_real_
    next
  }
  pfm <- PFMatrix(ID = jdf$matrix_id[i], name = jdf$geneSymbol[i],
                  profileMatrix = counts)
  icm <- toICM(pfm)
  ic_vec[i] <- sum(TFBSTools::totalIC(icm))
}
jdf$information_content <- ic_vec
cat("  IC range:", round(range(ic_vec, na.rm = TRUE), 2), "\n")
cat("  IC mean:",  round(mean(ic_vec,  na.rm = TRUE), 2), "\n")

# ── 3. JASPAR REST API: confirm experimentType per motif (with cache) ───────
# https://jaspar.elixir.no/api/v1/matrix/{matrix_id}/
# We use MotifDb's experimentType as the primary signal; API fetch is a
# sanity cross-check, capped + cached.
cat("\nJASPAR API cross-check (with cache, capped at 200 motifs)...\n")
api_base <- "https://jaspar.elixir.no/api/v1/matrix/"

fetch_jaspar_meta <- function(matrix_id) {
  cache_file <- file.path(CACHE_DIR, paste0(matrix_id, ".json"))
  if (file.exists(cache_file)) {
    return(tryCatch(fromJSON(cache_file), error = function(e) NULL))
  }
  url <- paste0(api_base, matrix_id, "/")
  resp <- tryCatch(GET(url, timeout(10)), error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) return(NULL)
  obj <- tryCatch(fromJSON(content(resp, as = "text", encoding = "UTF-8")),
                  error = function(e) NULL)
  if (!is.null(obj)) {
    writeLines(content(resp, as = "text", encoding = "UTF-8"), cache_file)
  }
  obj
}

# We hit API only for motifs whose TFs appear in the existing motifbreakR output
# (most relevant for tier comparison) — keeps API load minimal.
motif_score_file <- file.path(OUT_DIR, "motif_disruption_scores.csv")
if (file.exists(motif_score_file)) {
  disrupted <- fread(motif_score_file)
  used_tfs <- unique(disrupted$tf_name)
  use_api_idx <- which(toupper(jdf$geneSymbol) %in% toupper(used_tfs))
  cat("  Will query API for", length(use_api_idx), "motifs whose TF was disrupted\n")
} else {
  use_api_idx <- integer(0)
  cat("  motif_disruption_scores.csv not found — skipping API queries\n")
}

# Cap API hits to keep runtime bounded
use_api_idx <- head(use_api_idx, 300)

jdf$api_data_type <- NA_character_
jdf$api_source    <- NA_character_

for (i in use_api_idx) {
  mat_id <- jdf$matrix_id[i]
  obj <- fetch_jaspar_meta(mat_id)
  if (!is.null(obj)) {
    # JASPAR API field is "type" or "data_type" depending on version
    if (!is.null(obj$type))      jdf$api_data_type[i] <- as.character(obj$type)
    if (!is.null(obj$data_type)) jdf$api_data_type[i] <- as.character(obj$data_type)
    if (!is.null(obj$source))    jdf$api_source[i]    <- as.character(obj$source)
  }
}

n_api_hits <- sum(!is.na(jdf$api_data_type))
cat("  JASPAR API: successful fetches:", n_api_hits, "/", length(use_api_idx), "\n")

# ── 4. HOCOMOCO consensus check ─────────────────────────────────────────────
# For each TF in JASPAR, find the HOCOMOCO PWM for the same gene and compute
# PWMSimilarity (Pearson). High cor (>0.7) means JASPAR PWM is consistent with
# the HOCOMOCO database independent of source.
cat("\nHOCOMOCO consensus check (PWMSimilarity Pearson)...\n")
hdf <- as.data.frame(values(hocomoco_motifs))
hdf$motif_id <- names(hocomoco_motifs)

# Aggregate HOCOMOCO PWMs by geneSymbol — take the first PWM per gene
# (HOCOMOCO has multiple "M"/"P"/"S"/"B" subscores per TF; first is fine)
hocomoco_by_tf <- list()
for (gene in unique(hdf$geneSymbol)) {
  if (is.na(gene) || gene == "") next
  idx <- which(hdf$geneSymbol == gene)[1]
  pwm_mat <- as.matrix(hocomoco_motifs[[idx]])
  if (ncol(pwm_mat) < 4) next
  hocomoco_by_tf[[gene]] <- pwm_mat
}
cat("  HOCOMOCO PWMs available for", length(hocomoco_by_tf), "TFs\n")

jdf$hocomoco_cor    <- NA_real_
jdf$hocomoco_match  <- NA_character_

for (i in seq_along(jaspar_motifs)) {
  tf <- jdf$geneSymbol[i]
  if (is.na(tf) || tf == "") next
  # Try exact match
  h_mat <- hocomoco_by_tf[[tf]]
  if (is.null(h_mat)) {
    # Try uppercase / case-insensitive
    matches <- grep(paste0("^", tf, "$"), names(hocomoco_by_tf), ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) h_mat <- hocomoco_by_tf[[matches[1]]]
  }
  if (is.null(h_mat)) next
  j_mat <- as.matrix(jaspar_motifs[[i]])
  if (ncol(j_mat) < 4) next
  jpwm <- PWMatrix(ID = jdf$matrix_id[i], name = tf, profileMatrix = j_mat)
  hpwm <- PWMatrix(ID = "HOCOMOCO",        name = tf, profileMatrix = h_mat)
  cor_val <- tryCatch(PWMSimilarity(jpwm, hpwm, method = "Pearson"),
                      error = function(e) NA_real_)
  jdf$hocomoco_cor[i]   <- cor_val
  jdf$hocomoco_match[i] <- if (!is.null(matches <- jdf$geneSymbol[i])) tf else NA_character_
}

# Replace populated hocomoco_match values
jdf$hocomoco_match <- ifelse(!is.na(jdf$hocomoco_cor), jdf$geneSymbol, NA_character_)

n_hocomoco_match <- sum(!is.na(jdf$hocomoco_cor))
n_high_cor       <- sum(jdf$hocomoco_cor > 0.7, na.rm = TRUE)
cat("  JASPAR motifs with HOCOMOCO counterpart:", n_hocomoco_match, "\n")
cat("  Motifs with Pearson cor > 0.7:", n_high_cor, "\n")

# ── 5. Grade each motif ─────────────────────────────────────────────────────
cat("\nAssigning A/B/C grades...\n")
chip_types  <- c("ChIP-seq", "ChIP-exo")
selex_types <- c("HT-SELEX", "CAP-SELEX", "NCAP-SELEX", "SELEX", "SMiLE-seq",
                 "PBM", "bacterial 1-hybrid")
api_chip_keywords <- c("chip", "ChIP", "CHIP", "ChIP-seq", "ChIP-exo")

# Prefer experimentType from MotifDb metadata; fall back to API answer if any
jdf$effective_data_type <- jdf$experimentType
no_meta <- is.na(jdf$effective_data_type) | jdf$effective_data_type == ""
jdf$effective_data_type[no_meta] <- jdf$api_data_type[no_meta]

is_chip   <- jdf$effective_data_type %in% chip_types |
             grepl(paste(api_chip_keywords, collapse = "|"),
                   jdf$effective_data_type, ignore.case = TRUE)
is_selex  <- jdf$effective_data_type %in% selex_types |
             grepl("SELEX|PBM|SMiLE", jdf$effective_data_type, ignore.case = TRUE)
high_cor  <- !is.na(jdf$hocomoco_cor) & jdf$hocomoco_cor > 0.7
good_ic   <- !is.na(jdf$information_content) & jdf$information_content > 12

jdf$grade <- "C"
jdf$grade[is_chip] <- "B"
jdf$grade[is_selex & good_ic] <- "B"
jdf$grade[high_cor] <- "B"
jdf$grade[is_chip & high_cor] <- "A"

cat("Grade distribution:\n")
print(table(jdf$grade, useNA = "ifany"))

# ── 6. Write motif_quality_tier.csv ─────────────────────────────────────────
out_quality <- jdf %>%
  transmute(
    motif_id          = motif_id,
    matrix_id         = matrix_id,
    tf_name           = geneSymbol,
    data_source       = dataSource,
    jaspar_datatype   = effective_data_type,
    motifdb_exp_type  = experimentType,
    api_data_type     = api_data_type,
    jaspar_ic         = information_content,
    hocomoco_cor      = hocomoco_cor,
    hocomoco_match    = hocomoco_match,
    grade             = grade
  )

fwrite(out_quality, file.path(OUT_DIR, "motif_quality_tier.csv"))
cat("\nWrote motif_quality_tier.csv:", nrow(out_quality), "motifs\n")

# ── 7. Per-TF disruption-count comparison across tiers ──────────────────────
if (file.exists(motif_score_file)) {
  cat("\nBuilding motif_tier_comparison.csv...\n")
  disrupted <- fread(motif_score_file)
  # Disrupted file has tf_name (geneSymbol). Join via uppercase to handle case.
  qual_lookup <- out_quality %>%
    distinct(tf_name, grade, .keep_all = TRUE) %>%
    mutate(tf_upper = toupper(tf_name)) %>%
    select(tf_upper, grade)
  # If a TF has multiple motifs (different grades), we keep the BEST grade per TF.
  # Order: A > B > C
  grade_rank <- c(A = 1, B = 2, C = 3)
  qual_lookup <- out_quality %>%
    mutate(tf_upper = toupper(tf_name),
           rank = grade_rank[grade]) %>%
    group_by(tf_upper) %>%
    slice_min(rank, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(tf_upper, grade)

  disrupted_tf <- disrupted %>%
    mutate(tf_upper = toupper(tf_name)) %>%
    left_join(qual_lookup, by = "tf_upper") %>%
    mutate(grade = ifelse(is.na(grade), "C", grade))

  per_tf <- disrupted_tf %>%
    group_by(tf_name) %>%
    summarise(
      n_disruptions     = n(),
      n_strong          = sum(effect == "strong", na.rm = TRUE),
      in_disease_regulon = any(motif_in_disease_regulon, na.rm = TRUE),
      tf_best_grade     = dplyr::first(grade[order(grade_rank[grade])]),
      .groups = "drop"
    )

  per_tf$keep_all      <- TRUE                       # all
  per_tf$keep_AB       <- per_tf$tf_best_grade %in% c("A", "B")
  per_tf$keep_Aonly    <- per_tf$tf_best_grade == "A"

  fwrite(per_tf, file.path(OUT_DIR, "motif_tier_comparison.csv"))

  cat("Per-TF disruption tier counts:\n")
  cat("  All tiers (kept by Script 56):       ", sum(per_tf$keep_all),    "TFs\n")
  cat("  A/B-grade only:                       ", sum(per_tf$keep_AB),     "TFs\n")
  cat("  A-grade only (strict):                ", sum(per_tf$keep_Aonly),  "TFs\n")

  cat("\nDisease regulon TFs surviving each tier:\n")
  cat("  All tiers (current Script 56 output): ",
      sum(per_tf$in_disease_regulon),                                "TFs\n")
  cat("  A/B-grade only:                       ",
      sum(per_tf$in_disease_regulon & per_tf$keep_AB),               "TFs\n")
  cat("  A-grade only (strict):                ",
      sum(per_tf$in_disease_regulon & per_tf$keep_Aonly),            "TFs\n")
} else {
  cat("\nNote: motif_disruption_scores.csv not present; skipping comparison.\n")
}

cat("\n============================================================\n")
cat("Script 56d complete.\n")
cat("============================================================\n")
