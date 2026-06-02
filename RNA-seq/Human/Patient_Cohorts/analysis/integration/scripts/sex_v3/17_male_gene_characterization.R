#!/usr/bin/env Rscript
# sex_v3/17_male_gene_characterization.R
# ---------------------------------------------------------------------------
# Pillar 8A — Characterize Male-biased genes from v5 classifier.
#
# For each gene with class_v5_interaction %in% {Male_biased, Male_biased_F_underpowered}
# in interaction_classifier_v5.csv, gather:
#   - gene symbol, ensembl_id, chr, biotype
#   - β_F, β_M, padj_int, lfsr_F_strat, lfsr_M_strat, power_M_at_F, evidence_tier
#   - Cross-references against external sb-gene tables (graceful if missing):
#       * Oliva2020 GTEx sb-genes liver
#       * Mostafavi2023 sb-eQTL liver
#       * Khramtsova2019 sb-genes
#       * ENCODE AR HepG2 ChIP targets
#       * Tukiainen2017 XCI status
#   - Manual stub biological notes (header-curated)
#   - Pre-registered suspicious flag rule:
#       chrY gene OR (no Mostafavi/AR support AND 1-cohort-dominated)
#         => suspicious_flag = "chrY_pseudo_artifact" / "single_cohort_artifact"
#
# Output: male_biased_characterization.csv
# Runtime: <30 min on login or io 1-core.
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
EXT   <- file.path(BASE, "data/external")
UTILS <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "scripts/sex_v3/sex_v3_utils.R")
if (file.exists(UTILS)) source(UTILS)

OUT_CSV <- file.path(SEXV3, "male_biased_characterization.csv")
LOG_TXT <- file.path(SEXV3, "logs/male_gene_characterization.log")
dir.create(dirname(LOG_TXT), showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  msg <- sprintf("[%s] %s",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_TXT, append = TRUE)
}

# ---------------------------------------------------------------------------
# 1. Load v5 classifier and filter to Male-biased genes
# ---------------------------------------------------------------------------
v5_csv <- file.path(SEXV3, "interaction_classifier_v5.csv")
# interaction_classifier_v5.csv is a legacy disease_vs_ctrl artifact (v5
# pipeline). For non-canonical contrasts (mash_vs_masl etc.) it is absent.
# Gracefully degrade: write the empty stub and exit so the downstream
# aggregator + consensus can proceed.
if (!file.exists(v5_csv)) {
  log_msg("v5 classifier not present (",
          basename(v5_csv),
          "); contrast lacks a v5 sex-class call. Writing empty stub.")
  m_dt <- data.table()
} else {
  v5 <- fread(v5_csv)
  log_msg("Loaded v5 classifier: ", nrow(v5), " rows")
  male_classes <- c("Male_biased", "Male_biased_F_underpowered")
  m_dt <- v5[class_v5_interaction %in% male_classes]
  log_msg("Male-biased genes (any class): ", nrow(m_dt))
}
if (nrow(m_dt) == 0) {
  log_msg("WARNING: zero Male-biased genes — writing empty stub")
  fwrite(data.table(
    gene_symbol = character(0),
    ensembl_id = character(0),
    chr = character(0),
    biotype = character(0),
    class_v5_interaction = character(0),
    beta_F = numeric(0),
    beta_M = numeric(0),
    padj_int_5k = numeric(0),
    padj_int_full = numeric(0),
    lfsr_F_strat = numeric(0),
    lfsr_M_strat = numeric(0),
    power_M_at_F = numeric(0),
    evidence_tier_A_B_C = character(0),
    oliva2020_sb = logical(0),
    mostafavi2023_sb_eqtl = logical(0),
    khramtsova2019_sb = logical(0),
    ar_hepg2_target = logical(0),
    xci_status = character(0),
    suspicious_flag = character(0),
    bio_note = character(0)
  ), OUT_CSV)
  log_msg("Wrote empty stub: ", OUT_CSV)
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# 2. External cross-reference tables — graceful degradation
# ---------------------------------------------------------------------------
load_tsv_safe <- function(path, key_cols = "gene_symbol") {
  if (!file.exists(path)) {
    log_msg("MISSING external file (will write NA): ", path)
    return(NULL)
  }
  dt <- tryCatch(fread(path), error = function(e) {
    log_msg("ERROR reading ", path, ": ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  log_msg("Loaded external table ", basename(path), ": ", nrow(dt), " rows")
  dt
}

oliva_path     <- file.path(EXT, "oliva2020_gtex_sbgenes_liver.tsv")
mostafavi_path <- file.path(EXT, "mostafavi2023_sb_eqtl_liver.tsv")
khramtsova_path<- file.path(EXT, "khramtsova2019_sbgenes.tsv")
ar_path        <- file.path(EXT, "encode_liver/AR_HepG2_targets.tsv")
xci_path       <- file.path(EXT, "xci_status/tukiainen2017_xci.tsv")

oliva     <- load_tsv_safe(oliva_path)
mostafavi <- load_tsv_safe(mostafavi_path)
khramtsova<- load_tsv_safe(khramtsova_path)
ar        <- load_tsv_safe(ar_path)
xci       <- load_tsv_safe(xci_path)

annotate_in_set <- function(symbol, ref_dt, sym_col_candidates = c("gene_symbol","gene","symbol","SYMBOL")) {
  if (is.null(ref_dt)) return(NA)
  sc <- intersect(sym_col_candidates, names(ref_dt))[1]
  if (is.na(sc)) return(NA)
  symbol %in% unique(ref_dt[[sc]])
}

xci_lookup <- function(symbol) {
  if (is.null(xci)) return(NA_character_)
  sc <- intersect(c("gene_symbol","gene","symbol","Gene"), names(xci))[1]
  cc <- intersect(c("xci_status","status","Combined_XCI_status","XCI"), names(xci))[1]
  if (is.na(sc) || is.na(cc)) return(NA_character_)
  hit <- xci[get(sc) == symbol]
  if (nrow(hit) == 0) return(NA_character_)
  as.character(hit[[cc]][1])
}

# ---------------------------------------------------------------------------
# 3. Manual biology stubs — gene symbol keyed; user/curator to expand later.
# ---------------------------------------------------------------------------
bio_stub_lookup <- function(symbol, chr) {
  notes <- c(
    POU6F2 = "POU-class homeodomain TF, retina/neural enriched; sparse liver expression. Male-biased here may reflect rare-cell or contamination signal.",
    CCL8   = "CC chemokine (MCP-2); monocyte/T-cell recruitment. Liver expression context-dependent; testosterone-modulated leukocyte recruitment plausible.",
    ASS1   = "Argininosuccinate synthase 1, urea cycle; hepatocyte-essential. Sex-dimorphic urea-cycle flux reported in rodent liver. F_underpowered call.",
    DDX3Y  = "chrY DEAD-box helicase; male-exclusive expression. Should NOT appear as sex-modulated DEG — chrY-mapping artifact if present.",
    XIST   = "X-inactive specific transcript; constitutive F-only ncRNA. Should not be flagged as disease-modulated."
  )
  if (!is.null(symbol) && length(symbol) == 1 && symbol %in% names(notes)) {
    return(notes[[symbol]])
  }
  # Default: chromosome-keyed cue
  if (!is.na(chr) && chr == "chrY") {
    return("chrY-encoded gene; male-exclusive expression. Differential-expression call between sexes is artifactual — expect missing F counts.")
  }
  if (!is.na(chr) && chr == "chrX") {
    return("chrX gene; check XCI escape status. If escape-variable, sex-modulated DE plausible; if subject to XCI, F downregulation expected at baseline.")
  }
  paste0("Autosomal gene; mechanism unknown. Suggest cross-check against Oliva2020 sb-genes + Mostafavi2023 sb-eQTL + AR ChIP and report cohort-distribution forest.")
}

# ---------------------------------------------------------------------------
# 4. Build output rows
# ---------------------------------------------------------------------------
out_rows <- lapply(seq_len(nrow(m_dt)), function(i) {
  r <- m_dt[i]
  sym <- if ("gene_symbol" %in% names(r)) r$gene_symbol else NA_character_
  if (is.na(sym) || sym == "") sym <- r$ensembl_base
  chr <- if ("chr" %in% names(r)) r$chr else NA_character_
  bt  <- if ("gene_biotype" %in% names(r)) r$gene_biotype else NA_character_

  in_oliva     <- annotate_in_set(sym, oliva)
  in_mostafavi <- annotate_in_set(sym, mostafavi)
  in_khram     <- annotate_in_set(sym, khramtsova)
  in_ar        <- annotate_in_set(sym, ar)
  xci_st       <- xci_lookup(sym)

  # Pre-registered suspicious flag
  flag <- NA_character_
  if (!is.na(chr) && chr == "chrY") {
    flag <- "chrY_pseudo_artifact"
  } else if (isFALSE(in_mostafavi) && isFALSE(in_ar)) {
    # Cohort-dominated proxy not available without per-cohort fits;
    # leave a soft tag the curator can confirm against forest plot.
    flag <- "weak_external_support"
  }

  data.table(
    gene_symbol           = sym,
    ensembl_id            = r$ensembl_base,
    chr                   = chr,
    biotype               = bt,
    class_v5_interaction  = r$class_v5_interaction,
    beta_F                = if ("beta_F" %in% names(r)) r$beta_F else NA_real_,
    beta_M                = if ("beta_M" %in% names(r)) r$beta_M else NA_real_,
    padj_int_5k           = if ("padj_int_5k" %in% names(r)) r$padj_int_5k else NA_real_,
    padj_int_full         = if ("padj_int_full" %in% names(r)) r$padj_int_full else NA_real_,
    lfsr_F_strat          = if ("lfsr_F" %in% names(r)) r$lfsr_F else NA_real_,
    lfsr_M_strat          = if ("lfsr_M" %in% names(r)) r$lfsr_M else NA_real_,
    power_M_at_F          = if ("power_M_at_F" %in% names(r)) r$power_M_at_F else NA_real_,
    evidence_tier_A_B_C   = if ("evidence_tier_A_B_C" %in% names(r)) r$evidence_tier_A_B_C else NA_character_,
    oliva2020_sb          = in_oliva,
    mostafavi2023_sb_eqtl = in_mostafavi,
    khramtsova2019_sb     = in_khram,
    ar_hepg2_target       = in_ar,
    xci_status            = xci_st,
    suspicious_flag       = flag,
    bio_note              = bio_stub_lookup(sym, chr)
  )
})

out_dt <- rbindlist(out_rows, fill = TRUE)
log_msg("Built characterization table: ", nrow(out_dt), " rows, ", ncol(out_dt), " cols")

# Use atomic writer if loaded, else fallback
if (exists("write_atomic_csv")) {
  write_atomic_csv(out_dt, OUT_CSV)
} else {
  fwrite(out_dt, OUT_CSV)
}
log_msg("Wrote: ", OUT_CSV)

# Print summary table for log
print(out_dt[, .(gene_symbol, chr, class_v5_interaction,
                 beta_F = round(beta_F, 3), beta_M = round(beta_M, 3),
                 oliva2020_sb, mostafavi2023_sb_eqtl,
                 ar_hepg2_target, suspicious_flag)])

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
