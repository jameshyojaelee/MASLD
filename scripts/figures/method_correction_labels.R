# =============================================================================
# method_correction_labels.R
# -----------------------------------------------------------------------------
# SHARED human-readable label lookups for the multimethod / degx-factorial
# validation figures. Sourced by the figS_multimethod_* / figS_degx_* /
# figS_degcount_by_scheme / fig1_loo_cv_panels scripts so the cryptic engine,
# correction, and surrogate-variable codes (C0..C11, C6g/s/r, kbe/kna/kleek,
# engine short-codes) are NEVER shown raw on a figure axis / legend / title.
#
# WHAT MAPS TO WHAT
#   * engine_label(x)  : DE-engine short-code -> readable name (dream, DESeq2
#     (Wald), edgeR (QLF), limma-voom (QW), metafor (RE), ...).
#   * correction_label(x): batch-correction code -> the ACTUAL covariate scheme
#     it fits, derived verbatim from the degx harness correction registry
#     (~/degx/R/corrections.R, make_correction()). e.g. C2 = "~dataset+sex+group"
#     (the canonical model), C6g/C6s/C6r = RUVg/RUVs/RUVr on top of cohort, etc.
#   * ksv_label(x)     : SV/W count code -> readable (be="num.sv (BE)",
#     leek="num.sv (Leek)", na="-", k5="k=5").
#   * cell_label(id)   : a full cell_id "engine__correction__k" (e.g.
#     "limma_voom_qw__C2__kna") -> "limma-voom (QW) - ~dataset+sex+group".
#   * engine_family(x) / engine_family_palette : categorical colour axis by
#     engine family (the colour the figures use; never one colour per raw cell).
#
# SOURCE OF TRUTH for the correction covariate schemes:
#   /gpfs/commons/home/jameslee/degx/R/corrections.R  ::  make_correction()
#   (read 2026-06-10). The free-text covariate strings below are the literal
#   `fixed_terms` / `random_formula` each correction id builds, expressed in the
#   canonical "~dataset+sex+group" notation used throughout the manuscript.
#
# All labels ASCII; "~" reads as the model formula tilde, "+" joins covariates.
# =============================================================================

# ---------------------------------------------------------------------------
# ENGINE labels (DE-engine short-code -> readable). Includes both the harness
# factorial short-codes (deseq2_wald, edger_qlf, ...) AND the bare aliases used
# by the legacy resample tables (deseq2, edger, metafor).
# ---------------------------------------------------------------------------
ENGINE_LABELS <- c(
  # --- harness factorial engine ids (run_methods_real.R) ---
  dream             = "dream",
  deseq2_wald       = "DESeq2 (Wald)",
  deseq2_lrt        = "DESeq2 (LRT)",
  edger_qlf         = "edgeR (QLF)",
  edger_qlf_robust  = "edgeR (QLF-robust)",
  edger_lrt         = "edgeR (LRT)",
  limma_voom        = "limma-voom",
  limma_voom_qw     = "limma-voom (QW)",
  limma_trend       = "limma-trend",
  metafor_re        = "metafor (RE)",
  # --- bare aliases used by the legacy multimethod_validation resample data ---
  deseq2            = "DESeq2",
  edger             = "edgeR",
  limma             = "limma-voom",
  metafor           = "metafor (RE)",
  # --- R-exact / sim-resample extras seen in figS_degx_multimethod.R ---
  deseq2_apeglm     = "DESeq2 (apeglm)",
  deseq2_ashr       = "DESeq2 (ashr)",
  edger_exact       = "edgeR (exact)",
  metafor_deseq2_re = "metafor DESeq2 (RE)",
  metafor_deseq2_fe = "metafor DESeq2 (FE)",
  metafor_deseq2_hk = "metafor DESeq2 (HK)",
  metafor_voom_re   = "metafor voom (RE)",
  metafor_voom_fe   = "metafor voom (FE)",
  metafor_voom_hk   = "metafor voom (HK)",
  metafor_fe        = "metafor (FE)",
  metafor_hk        = "metafor (HK)",
  combatseq_deseq2  = "ComBat-seq + DESeq2",
  sva_limma         = "sva + limma",
  # --- NB-sim resample family (figS_degx_multimethod.R method_labels) ---
  voom              = "voom",
  voom_tmm          = "voom (TMM)",
  voom_robust       = "voom (robust)",
  voom_ashr         = "voom (ashr)",
  nb_glm            = "NB-GLM (Wald)",
  nb_glm_lrt        = "NB-GLM (LRT)",
  nb_glm_tmm        = "NB-GLM (TMM)",
  nb_glm_ashr       = "NB-GLM (ashr)",
  dream_tmm         = "dream (TMM)",
  dream_ashr        = "dream (ashr)",
  metafor_deseq2    = "metafor (DESeq2)"
)

#' Readable label for a DE-engine short-code (vectorised). Unknown codes pass
#' through unchanged so a future engine is at worst shown as its raw id.
engine_label <- function(x) {
  x <- as.character(x)
  out <- ENGINE_LABELS[x]
  ifelse(is.na(out), x, unname(out))
}

# ---------------------------------------------------------------------------
# CORRECTION labels (batch-correction code -> the actual covariate scheme).
# Derived from ~/degx/R/corrections.R::make_correction(). Two flavours:
#   CORRECTION_LABELS       : the literal model the correction fits (for in-cell
#                             / detailed contexts, e.g. cell_label()).
#   CORRECTION_LABELS_SHORT : a compact "code + 1-3 word gloss" for crowded axes
#                             (boxplot x-axis with 12+ corrections).
# ---------------------------------------------------------------------------
CORRECTION_LABELS <- c(
  C0   = "~group",                              # intercept + group only
  C1   = "~dataset+group",                      # dataset fixed effect
  C2   = "~dataset+sex+group",                  # CANONICAL
  C3   = "~dataset+sex+age+group",              # listwise-drops NA-age samples
  C4   = "~sex+SVs+group (SVs replace dataset)",
  C5   = "~dataset+sex+SVs+group (SVA on top)",
  C6g  = "~dataset+sex+W+group (RUVg)",
  C6s  = "~dataset+sex+W+group (RUVs)",
  C6r  = "~dataset+sex+W+group (RUVr)",
  C6hk = "~sex+SVs+group (SVA, housekeeping controls)",
  C7   = "ComBat-seq -> ~sex+group",
  C8   = "~sex+group+(1|dataset) (dream random intercept)",
  C9   = "~sex+SVs+group+(1|dataset) (dream + SVs)",
  C10  = "~sex+group+(1+group|dataset) (dream random slope)",
  C11  = "residualize dataset+sex+SVs -> test group (PCA-style)"
)

CORRECTION_LABELS_SHORT <- c(
  C0   = "C0: ~group",
  C1   = "C1: +dataset",
  C2   = "C2: +dataset+sex",
  C3   = "C3: +dataset+sex+age",
  C4   = "C4: SVA (replace cohort)",
  C5   = "C5: SVA (on top)",
  C6g  = "C6g: RUVg",
  C6s  = "C6s: RUVs",
  C6r  = "C6r: RUVr",
  C6hk = "C6hk: SVA (housekeeping)",
  C7   = "C7: ComBat-seq",
  C8   = "C8: dream (1|dataset)",
  C9   = "C9: dream + SVs",
  C10  = "C10: dream (1+group|dataset)",
  C11  = "C11: residualize-then-test"
)

#' Readable label for a batch-correction code (vectorised). `short=TRUE` returns
#' the compact "code: gloss" form for crowded axes; default returns the full
#' covariate scheme. Unknown codes pass through unchanged.
correction_label <- function(x, short = FALSE) {
  x <- as.character(x)
  tab <- if (short) CORRECTION_LABELS_SHORT else CORRECTION_LABELS
  out <- tab[x]
  ifelse(is.na(out), x, unname(out))
}

# ---------------------------------------------------------------------------
# k_sv labels (surrogate-variable / RUV-factor count code -> readable).
#   be   -> num.sv (BE)     (data-driven, Buja-Eyuboglu; the primary k)
#   leek -> num.sv (Leek)   (data-driven, Leek 2011; winner-eligible)
#   na   -> "-"             (correction is k-free)
#   2/5/...-> k=N           (fixed-k sensitivity arm)
# Accepts either the raw code ("kbe","k5","kna") or the stripped form
# ("be","5","na").
# ---------------------------------------------------------------------------
ksv_label <- function(x) {
  x <- sub("^k", "", as.character(x))          # kbe -> be, k5 -> 5, kna -> na
  vapply(x, function(v) {
    if (is.na(v) || v == "" || v == "na") return("—")  # em dash
    if (v == "be")   return("num.sv (BE)")
    if (v == "leek") return("num.sv (Leek)")
    if (grepl("^[0-9]+$", v)) return(paste0("k=", v))
    v
  }, character(1), USE.NAMES = FALSE)
}

# ---------------------------------------------------------------------------
# cell_label — a full factorial cell_id "engine__correction__k" rendered as
# "<engine readable> - <correction covariate scheme>". The k_sv suffix is
# appended only when it is an informative SV/W count (not 'na'/'-').
#   "limma_voom_qw__C2__kna"  -> "limma-voom (QW) - ~dataset+sex+group"
#   "dream__C9__kbe"          -> "dream - ~sex+SVs+group+(1|dataset) (dream + SVs) [num.sv (BE)]"
# `short=TRUE` uses the compact correction gloss.
# ---------------------------------------------------------------------------
cell_label <- function(id, short = FALSE) {
  id <- as.character(id)
  vapply(id, function(one) {
    if (is.na(one) || one == "") return(one)
    p <- strsplit(one, "__", fixed = TRUE)[[1]]
    eng  <- if (length(p) >= 1) engine_label(p[1]) else one
    corr <- if (length(p) >= 2) correction_label(p[2], short = short) else NA_character_
    k    <- if (length(p) >= 3) sub("^k", "", p[3]) else NA_character_
    lab <- eng
    if (!is.na(corr)) lab <- paste0(lab, " · ", corr)   # middle dot
    if (!is.na(k) && k != "" && k != "na")
      lab <- paste0(lab, " [", ksv_label(k), "]")
    lab
  }, character(1), USE.NAMES = FALSE)
}

# ---------------------------------------------------------------------------
# Engine FAMILY (categorical colour axis) + palette. Used wherever the figures
# colour by method family rather than by raw cell id. Colours match the
# guideline palette already used across the degx figures.
# ---------------------------------------------------------------------------
ENGINE_FAMILY_LEVELS <- c("DESeq2", "edgeR", "limma", "dream", "metafor",
                          "voom", "NB-GLM", "batch-corr", "other")
ENGINE_FAMILY_PALETTE <- c(
  DESeq2       = "#4aa2c2",
  edgeR        = "#518dc9",
  limma        = "#9b75d6",
  dream        = "#40b499",
  metafor      = "#e37faf",
  voom         = "#4aa2c2",
  "NB-GLM"     = "#C9265E",
  "batch-corr" = "#f7bf87",
  other        = "#BDBDBD"
)

#' Engine FAMILY for a method/engine short-code (vectorised), for the
#' categorical colour axis. Covers harness engines, bare aliases, R-exact,
#' and the NB-sim resample family.
engine_family <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^deseq2",  x), "DESeq2",
  ifelse(grepl("^edger",   x), "edgeR",
  ifelse(grepl("^limma",   x), "limma",
  ifelse(grepl("^dream",   x), "dream",
  ifelse(grepl("^metafor", x), "metafor",
  ifelse(grepl("^voom",    x), "voom",
  ifelse(grepl("^nb_glm",  x), "NB-GLM",
  ifelse(grepl("combatseq|^sva", x), "batch-corr", "other"))))))))
}
