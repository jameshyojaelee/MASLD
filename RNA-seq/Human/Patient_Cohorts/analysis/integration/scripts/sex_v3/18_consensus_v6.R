#!/usr/bin/env Rscript
# sex_v3/18_consensus_v6.R
# ---------------------------------------------------------------------------
# Cross-pillar consensus — joins v5 + P1 random-slope + P2 perm FDR +
# P4 bootstrap stability + P6 cohort Q + 8B magnitude LOCO into a single
# canonical per-gene table with class_v6_consensus.
#
# GRACEFUL DEGRADATION: every pillar input is optional. Missing pillars =>
# all derived cols NA + tier_evidence_incomplete = TRUE.
#
# Inputs (all under SEXV3):
#   interaction_classifier_v5.csv             (always required)
#   interaction_classifier_v5_random.csv      (P1)
#   perm_fdr_v6.csv                           (P2)
#   bootstrap_stability_v6.csv                (P4)
#   per_cohort_forest_data_v6.csv             (P6) -- per-gene Q/I2 summary
#   loco_concordance_per_class_v6.csv         (8B) -- per-class metric, joined by class
#   sensitivity_sweep_v6.csv                  (P7) -- canonical-cell counts only
#
# Tier rules:
#   Strong:    q_emp_per_gene<0.05 AND padj_int_5k_random<0.05 AND
#              stability_v6>=0.8 AND loco_rho>=0.5 AND I2<75
#   Moderate:  q_emp_per_gene<0.10 AND padj_int_5k_random<0.10 AND
#              stability_v6>=0.6 AND loco_rho>=0.3 AND I2<90
#   Suggestive: |beta_int|>0.5 AND evidence_tier_A_B_C == "C" AND stability_v6>=0.6
#   Power-limited: class_v5_interaction ends with "_underpowered"
#   Uncertain: else
# Inherits direction (F/M-biased) from v5 if available.
#
# Output: sex_deg_classification_v6.csv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
})

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SEXV3 <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "results/integration/sex_v3")

# Contrast routing (sex_v3_utils.R::contrast_paths) — overrides SEXV3 / IDIR
# when CONTRAST_NAME != "disease_vs_ctrl"; default preserves legacy layout.
if (!exists("contrast_paths")) source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT",
             "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3/sex_v3_utils.R"))
.cpaths <- contrast_paths()
SEXV3 <- .cpaths$sexv3
IDIR  <- .cpaths$idir
dir.create(IDIR, recursive = TRUE, showWarnings = FALSE)
UTILS <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "scripts/sex_v3/sex_v3_utils.R")
if (file.exists(UTILS)) source(UTILS)

OUT_CSV <- file.path(SEXV3, "sex_deg_classification_v6.csv")
LOG_TXT <- file.path(SEXV3, "logs/consensus_v6.log")
dir.create(dirname(LOG_TXT), showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  msg <- sprintf("[%s] %s",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_TXT, append = TRUE)
}

# ---------------------------------------------------------------------------
# 1. Base table: v5 classifier (disease_vs_ctrl canonical) OR the P1
# random-slope output (interaction_classifier_v5_random.csv) for non-canonical
# contrasts where the legacy v5 pipeline did not run.
# ---------------------------------------------------------------------------
v5_csv      <- file.path(SEXV3, "interaction_classifier_v5.csv")
v5rand_csv  <- file.path(SEXV3, "interaction_classifier_v5_random.csv")
if (file.exists(v5_csv)) {
  v5 <- fread(v5_csv)
  log_msg("v5 classifier: ", nrow(v5), " genes")
  base <- v5[, .(gene)]
  if ("gene_symbol" %in% names(v5)) base[, gene_symbol := v5$gene_symbol]
  if ("ensembl_base" %in% names(v5)) base[, ensembl_base := v5$ensembl_base]
  if ("chr" %in% names(v5)) base[, chr := v5$chr]
  if ("gene_biotype" %in% names(v5)) base[, gene_biotype := v5$gene_biotype]
  for (col in c("beta_F","beta_M","beta_int","se_int","padj_int_5k",
                "padj_int_full","lfsr_int_ashr","class_v5_interaction",
                "evidence_tier_A_B_C")) {
    if (col %in% names(v5)) base[, (col) := v5[[col]]]
  }
} else if (file.exists(v5rand_csv)) {
  log_msg("v5 classifier absent; using P1 random-slope output as base (",
          v5rand_csv, ")")
  v5r <- fread(v5rand_csv)
  log_msg("base from P1 random-slope: ", nrow(v5r), " genes")
  base <- v5r[, .(gene)]
  if ("gene_symbol" %in% names(v5r)) base[, gene_symbol := v5r$gene_symbol]
  if ("ensembl_base" %in% names(v5r)) base[, ensembl_base := v5r$ensembl_base]
  if ("chr" %in% names(v5r)) base[, chr := v5r$chr]
  if ("gene_biotype" %in% names(v5r)) base[, gene_biotype := v5r$gene_biotype]
  # Map random-slope output names to the names downstream code expects:
  rename_map <- c(beta_int_random      = "beta_int",
                  se_int_random        = "se_int",
                  padj_int_5k_random   = "padj_int_5k",
                  padj_int_full_random = "padj_int_full",
                  class_v5_random      = "class_v5_interaction",
                  beta_F_strat         = "beta_F",
                  beta_M_strat         = "beta_M")
  for (src in names(rename_map)) {
    if (src %in% names(v5r)) base[, (rename_map[[src]]) := v5r[[src]]]
  }
  # Missing in P1: lfsr_int_ashr, evidence_tier_A_B_C (set NA)
  if (!"lfsr_int_ashr" %in% names(base))
    base[, lfsr_int_ashr := NA_real_]
  if (!"evidence_tier_A_B_C" %in% names(base))
    base[, evidence_tier_A_B_C := NA_character_]
} else {
  stop("Both v5 classifier and P1 random-slope output are missing — ",
       "cannot build consensus base. Looked for:\n  ", v5_csv,
       "\n  ", v5rand_csv)
}

# Direction
infer_direction <- function(cls, bF, bM) {
  if (is.na(cls)) return(NA_character_)
  if (grepl("Female", cls, fixed = TRUE)) return("F_biased")
  if (grepl("Male",   cls, fixed = TRUE)) return("M_biased")
  if (grepl("Divergent", cls, fixed = TRUE)) return("Divergent")
  if (grepl("Concordant", cls, fixed = TRUE)) return("Concordant")
  NA_character_
}
base[, direction_v5 := mapply(infer_direction, class_v5_interaction, beta_F, beta_M)]

# ---------------------------------------------------------------------------
# 2. Optional pillars — load with graceful degradation
# ---------------------------------------------------------------------------
missing_pillars <- character(0)

load_optional <- function(path, label) {
  if (!file.exists(path)) {
    log_msg("MISSING pillar input: ", label, "  (", path, ")")
    missing_pillars <<- c(missing_pillars, label)
    return(NULL)
  }
  dt <- tryCatch(fread(path), error = function(e) {
    log_msg("ERROR loading ", label, ": ", conditionMessage(e))
    missing_pillars <<- c(missing_pillars, label)
    NULL
  })
  if (!is.null(dt)) log_msg("Loaded ", label, ": ", nrow(dt), " rows")
  dt
}

# Per meta-review B0-9: load P6's *aggregated* cochranQ_v6.csv (per-gene Q,
# df, p_Q, I²) rather than per_cohort_forest_data_v6.csv (per-gene-per-cohort
# long-form). The forest file has no Q/p_Q/I² columns, so loading it silently
# zeroed the cohort-heterogeneity gate across all 34,453 genes.
p1 <- load_optional(file.path(SEXV3, "interaction_classifier_v5_random.csv"), "P1_random_slope")
p2 <- load_optional(file.path(SEXV3, "perm_fdr_v6.csv"),                      "P2_perm_fdr")
p4 <- load_optional(file.path(SEXV3, "bootstrap_stability_v6.csv"),           "P4_bootstrap_dream")
p6 <- load_optional(file.path(SEXV3, "cochranQ_v6.csv"),                      "P6_cochranQ")
p8b<- load_optional(file.path(SEXV3, "loco_concordance_per_class_v6.csv"),    "P8B_loco_magnitude")
p7 <- load_optional(file.path(SEXV3, "sensitivity_sweep_v6.csv"),             "P7_threshold_sweep")

# ---------------------------------------------------------------------------
# 3. Join optional cols
# ---------------------------------------------------------------------------
add_cols <- function(base, dt, key = "gene", cols_map) {
  for (out_col in names(cols_map)) {
    base[, (out_col) := NA_real_]
  }
  if (is.null(dt) || nrow(dt) == 0) return(base)
  if (!(key %in% names(dt))) return(base)
  for (out_col in names(cols_map)) {
    src <- cols_map[[out_col]]
    if (!(src %in% names(dt))) next
    sub <- dt[, c(key, src), with = FALSE]
    setnames(sub, src, "_val_")
    base[sub, on = setNames(key, key), (out_col) := i._val_]
  }
  base
}

# P1: random slope coefficients
base <- add_cols(base, p1, "gene", list(
  beta_int_random       = "beta_int_random",
  se_int_random         = "se_int_random",
  padj_int_5k_random    = "padj_int_5k_random",
  padj_int_full_random  = "padj_int_full_random",
  re_variant_used       = "re_variant_used"
))
if (!is.null(p1) && "class_v5_random" %in% names(p1)) {
  base[, class_v5_random := NA_character_]
  sub <- p1[, .(gene, class_v5_random)]
  base[sub, on = "gene", class_v5_random := i.class_v5_random]
}

# P2: permutation FDR
base <- add_cols(base, p2, "gene", list(
  q_emp_per_gene = "q_emp_per_gene",
  q_emp_pooled   = "q_emp_pooled",
  p_emp_per_gene = "p_emp_per_gene"
))

# P4: bootstrap stability — `13b_boot_aggregate.R` writes `stability_modal`
# (the rate at which each gene's modal class was recovered across 50 boots).
# Map it onto the consensus's `stability_v6` slot. Fallback ladder retained.
base <- add_cols(base, p4, "gene", list(
  stability_v6 = "stability_modal"
))
if (is.null(p4) || !("stability_modal" %in% names(p4))) {
  if (!is.null(p4)) {
    alt <- intersect(c("stability_v6","stability","fraction_same_class","n_same_class"),
                     names(p4))[1]
    if (!is.na(alt)) {
      base[, stability_v6 := NA_real_]
      sub <- p4[, c("gene", alt), with = FALSE]
      setnames(sub, alt, "_val_")
      base[sub, on = "gene", stability_v6 := i._val_]
    }
  }
}

# P6: cohort heterogeneity — `15b_cochranQ_aggregate.R` emits `p_Q` (not
# `Q_pval`), `Q_padj_BH`, `I2_pct`. The earlier `Q_pval` lookup missed,
# zeroing the I² gate across the genome.
base <- add_cols(base, p6, "gene", list(
  cochranQ_pval    = "p_Q",
  cochranQ_padj_BH = "Q_padj_BH",
  I2_pct           = "I2_pct"
))
if (!is.null(p6) && "cohort_heterogeneous" %in% names(p6)) {
  base[, cohort_heterogeneous := FALSE]
  sub <- p6[, .(gene, cohort_heterogeneous)]
  base[sub, on = "gene", cohort_heterogeneous := i.cohort_heterogeneous]
}

# 8B: magnitude LOCO — join by class
if (!is.null(p8b) && "class" %in% names(p8b)) {
  # Pool across folds per class
  pool <- p8b[, .(
    loco_rho = weighted.mean(spearman_rho, n_genes, na.rm = TRUE),
    loco_sign = weighted.mean(sign_concord, n_genes, na.rm = TRUE),
    loco_magnitude = weighted.mean(magnitude_concord, n_genes, na.rm = TRUE),
    loco_holdout_padj = weighted.mean(holdout_padj_lt_0.20, n_genes, na.rm = TRUE)
  ), by = class]
  base[, loco_rho := NA_real_]
  base[, loco_sign := NA_real_]
  base[, loco_magnitude := NA_real_]
  base[, loco_holdout_padj := NA_real_]
  base[pool, on = c("class_v5_interaction" = "class"), `:=`(
    loco_rho = i.loco_rho,
    loco_sign = i.loco_sign,
    loco_magnitude = i.loco_magnitude,
    loco_holdout_padj = i.loco_holdout_padj
  )]
} else {
  base[, loco_rho := NA_real_]
  base[, loco_sign := NA_real_]
  base[, loco_magnitude := NA_real_]
  base[, loco_holdout_padj := NA_real_]
}

# ---------------------------------------------------------------------------
# 4. Compute class_v6_consensus per gene
# ---------------------------------------------------------------------------
classify <- function(q_pg, padj_rand, stab, rho, i2, cls_v5, evid_tier, beta_int) {
  # Power-limited carve-out first
  if (!is.na(cls_v5) && grepl("_underpowered$", cls_v5)) return("Power_limited")
  # Use defensive NA-handling — missing input => fail gate
  pass_strong <- (!is.na(q_pg)     && q_pg     < 0.05) &&
                 (!is.na(padj_rand)&& padj_rand< 0.05) &&
                 (!is.na(stab)     && stab    >= 0.8 ) &&
                 (!is.na(rho)      && rho     >= 0.5 ) &&
                 (!is.na(i2)       && i2      <  75  )
  if (pass_strong) return("Strong")
  pass_moderate <- (!is.na(q_pg)     && q_pg     < 0.10) &&
                   (!is.na(padj_rand)&& padj_rand< 0.10) &&
                   (!is.na(stab)     && stab    >= 0.6 ) &&
                   (!is.na(rho)      && rho     >= 0.3 ) &&
                   (!is.na(i2)       && i2      <  90  )
  if (pass_moderate) return("Moderate")
  pass_suggestive <- (!is.na(beta_int)  && abs(beta_int) > 0.5) &&
                     (!is.na(evid_tier) && evid_tier == "C")    &&
                     (!is.na(stab)      && stab >= 0.6)
  if (pass_suggestive) return("Suggestive")
  "Uncertain"
}

# Diagnostic emit (meta-review B0-9): track per-gene which pillars actually
# carried non-NA values so the consensus quality can be inspected gene-by-gene
# instead of relying on a single run-wide tier_evidence_incomplete scalar.
base[, n_pillars_used := (!is.na(q_emp_per_gene))      +
                         (!is.na(padj_int_5k_random))  +
                         (!is.na(stability_v6))        +
                         (!is.na(loco_rho))            +
                         (!is.na(I2_pct))]

# Defensive length assertion — mapply silently recycles length-1 inputs which
# would broadcast NA across every gene if a column collapsed to a singleton.
stopifnot(length(base$q_emp_per_gene)       == nrow(base),
          length(base$padj_int_5k_random)   == nrow(base),
          length(base$stability_v6)         == nrow(base),
          length(base$loco_rho)             == nrow(base),
          length(base$I2_pct)               == nrow(base),
          length(base$class_v5_interaction) == nrow(base),
          length(base$evidence_tier_A_B_C)  == nrow(base),
          length(base$beta_int)             == nrow(base))

class_v6 <- mapply(classify,
                   base$q_emp_per_gene,
                   base$padj_int_5k_random,
                   base$stability_v6,
                   base$loco_rho,
                   base$I2_pct,
                   base$class_v5_interaction,
                   base$evidence_tier_A_B_C,
                   base$beta_int)
stopifnot(length(class_v6) == nrow(base))
base[, class_v6_consensus := class_v6]

# Per-gene quality label: distinguishes a real cross-pillar calibrated null
# (Uncertain with n_pillars_used >= 3) from an artifact (Uncertain with
# n_pillars_used <= 1 — most gates were NA so all genes default to Uncertain).
base[, consensus_quality := fifelse(
  class_v6_consensus %in% c("Strong","Moderate","Suggestive"), class_v6_consensus,
  fifelse(n_pillars_used >= 3, "calibrated_null",
  fifelse(n_pillars_used >= 1, "partial_evidence", "incomplete"))
)]

# Combine with direction for canonical
base[, class_v6_combined := fifelse(
  class_v6_consensus %in% c("Strong","Moderate","Suggestive"),
  paste0(class_v6_consensus, "_", direction_v5),
  class_v6_consensus
)]

# Tier-evidence-incomplete flag (run-wide; kept for backward compat)
base[, tier_evidence_incomplete := length(missing_pillars) > 0]
base[, missing_pillars := paste(missing_pillars, collapse = ";")]

# ---------------------------------------------------------------------------
# 5. Tally + write
# ---------------------------------------------------------------------------
log_msg("class_v6_consensus tally:")
print(base[, .N, by = class_v6_consensus][order(-N)])
log_msg("class_v6_combined tally:")
print(base[, .N, by = class_v6_combined][order(-N)])

if (exists("write_atomic_csv")) {
  write_atomic_csv(base, OUT_CSV)
} else {
  fwrite(base, OUT_CSV)
}
log_msg("Wrote: ", OUT_CSV)
log_msg("Missing pillars: ",
        if (length(missing_pillars) == 0) "NONE" else paste(missing_pillars, collapse = ", "))

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
