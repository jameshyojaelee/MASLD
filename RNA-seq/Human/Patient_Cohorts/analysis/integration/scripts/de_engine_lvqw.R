#!/usr/bin/env Rscript
# ===========================================================================
# de_engine_lvqw.R — the SINGLE vetted limma-voom-quality-weighted (LVQW) DEG
# engine. LIBRARY-ONLY: this file defines functions and has NO side effects
# when sourced (no I/O, no model fits, no global mutation), so any analysis
# script can `source()` it safely.
#
# Why this exists
# ---------------
# All "other" bulk DEG analyses in this repo (per-study, strata, successive-stage
# transitions, etc.) should call ONE vetted engine instead of re-implementing
# voom/eBayes/ashr by hand. This engine is copied verbatim from the already-vetted
# canonical producer 05h_limma_voom_qw_canonical.R (engine = lines ~40-54):
#     voomWithQualityWeights(dge, design)  ->  eBayes(lmFit(v, design))
#     topTable(..., sort.by="none");  SE = abs(logFC / t)
#     ashr::ash(logFC, SE, mixcompdist="normal")  ->  shrunk_logFC + lfsr
#
# One engine, many designs
# -------------------------
# The DESIGN is the caller's responsibility. This file does NOT decide model
# formulas, contrasts, or which coefficient to test — it only fits whatever
# design/coef/contrast it is handed and returns a fixed two-tier schema. Callers
# build their own designs (helper build_design_guarded() is provided for the
# common case) and MUST verify their target coef survived design construction:
#     stopifnot(coef %in% colnames(design))
#
# ashr mixcompdist
# ----------------
# This engine standardizes on mixcompdist="normal", matching the 05h canonical.
# NOTE: the legacy dream-era ashr calls used mixcompdist="halfuniform"; we do
# NOT reproduce that here. "normal" is the canonical going forward.
#
# Output schema (matches dream_results_ashr.csv consumers)
# --------------------------------------------------------
# A data.table with EXACTLY these columns, in this order:
#   gene, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr, AveExpr
# topTable's native column names are preserved: effect = logFC, average
# expression = AveExpr, moderated t = t, raw p = P.Value (adj.P.Val -> padj).
# When do_ashr=FALSE, shrunk_logFC and lfsr are NA_real_.
# ===========================================================================
suppressMessages({
  library(edgeR)
  library(limma)
  library(ashr)
  library(data.table)
})

# ---------------------------------------------------------------------------
# add_ashr() — adaptive shrinkage on a fitted result table.
# Factored out of 05h lines 48-52 EXACTLY: shrink only finite logFC/SE with
# SE>0; write PosteriorMean -> shrunk_logFC and lfsr, leaving non-shrinkable
# rows as NA_real_. Mutates `dt` by reference (data.table) and returns it.
# ---------------------------------------------------------------------------
add_ashr <- function(dt, mixcompdist = "normal") {
  ok <- is.finite(dt$logFC) & is.finite(dt$SE) & dt$SE > 0
  a  <- ashr::ash(dt$logFC[ok], dt$SE[ok], mixcompdist = mixcompdist)
  dt[, shrunk_logFC := NA_real_][ok, shrunk_logFC := a$result$PosteriorMean]
  dt[, lfsr := NA_real_][ok, lfsr := a$result$lfsr]
  dt
}

# ---------------------------------------------------------------------------
# .lvqw_finalize() — internal: assemble topTable output into the fixed
# two-tier schema (and run ashr if requested). Shared by fit_lvqw() and
# fit_lvqw_contrast() so the schema is defined in exactly one place.
# ---------------------------------------------------------------------------
.lvqw_finalize <- function(res, do_ashr, mixcompdist) {
  res$gene <- rownames(res)
  dt <- as.data.table(res)
  setnames(dt, "adj.P.Val", "padj", skip_absent = TRUE)
  dt[, SE := abs(logFC / t)]
  if (do_ashr) {
    dt <- add_ashr(dt, mixcompdist)
  } else {
    dt[, shrunk_logFC := NA_real_]
    dt[, lfsr := NA_real_]
  }
  dt[, .(gene, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr, AveExpr)]
}

# ---------------------------------------------------------------------------
# fit_lvqw() — fit the LVQW engine for a single named coefficient.
#   dge          : edgeR DGEList (normalized; lib sizes set by caller)
#   design       : model.matrix the caller built
#   coef         : name of the coefficient column in `design` to test
#   do_ashr      : run ashr adaptive shrinkage (default TRUE)
#   mixcompdist  : ashr mixture component dist (default "normal", canonical)
#   weights      : TRUE -> voomWithQualityWeights; FALSE -> plain voom
# Returns the fixed two-tier schema data.table (see header).
# ---------------------------------------------------------------------------
fit_lvqw <- function(dge, design, coef, do_ashr = TRUE,
                     mixcompdist = "normal", weights = TRUE) {
  stopifnot(coef %in% colnames(design))
  v   <- if (weights) voomWithQualityWeights(dge, design) else voom(dge, design)
  fit <- eBayes(lmFit(v, design))
  res <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  .lvqw_finalize(res, do_ashr, mixcompdist)
}

# ---------------------------------------------------------------------------
# fit_lvqw_contrast() — fit the LVQW engine for a custom contrast vector.
# Supports successive-stage transition contrasts (e.g. F1->F2 vs F0->F1),
# replacing the dream makeContrastsDream pattern.
#   contrast_vec : numeric contrast over the columns of `design` (length =
#                  ncol(design)), or a matrix with one contrast column.
# Pipeline: voomWithQualityWeights -> lmFit -> contrasts.fit -> eBayes ->
#           topTable(coef = 1). Same fixed two-tier output schema.
# ---------------------------------------------------------------------------
fit_lvqw_contrast <- function(dge, design, contrast_vec, do_ashr = TRUE,
                              mixcompdist = "normal", weights = TRUE) {
  v    <- if (weights) voomWithQualityWeights(dge, design) else voom(dge, design)
  fit  <- lmFit(v, design)
  fit2 <- eBayes(contrasts.fit(fit, contrasts = contrast_vec))
  res  <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
  .lvqw_finalize(res, do_ashr, mixcompdist)
}

# ---------------------------------------------------------------------------
# build_design_guarded() — build a rank-safe model.matrix from `info`.
#
# In per-study / strata subsets a factor term (e.g. `dataset`) can collapse to
# a single observed level, which makes ~ term contribute no columns or makes
# the design rank-deficient. This helper:
#   1. droplevels() every factor named in `rhs_terms`,
#   2. DROPS any factor term that has <2 observed levels in the subset,
#   3. builds model.matrix(~ <surviving terms>),
#   4. checks qr(design)$rank < ncol(design); if rank-deficient, identifies
#      aliased columns via limma::nonEstimable() (with a caret/QR-pivot
#      fallback if nonEstimable is unavailable) and drops them (with a warning).
#
# Args:
#   info          : data.frame of sample-level covariates (rows = samples).
#   rhs_terms     : character vector of term names in `info` to put on the RHS
#                   of the formula (order preserved).
#   response_drop : optional character vector of terms to forcibly exclude
#                   (e.g. a response/outcome you do not want as a covariate).
#
# Returns list(design, formula_used, dropped) where:
#   design       = the (possibly column-reduced) model.matrix,
#   formula_used = the formula string actually built,
#   dropped      = character vector of dropped TERMS and/or aliased COLUMNS.
#
# CALLER CONTRACT: this helper may drop a term whose single level made it
# unidentifiable. The caller MUST verify its target coefficient survived:
#     stopifnot(coef %in% colnames(out$design))
# ---------------------------------------------------------------------------
build_design_guarded <- function(info, rhs_terms, response_drop = NULL) {
  stopifnot(is.data.frame(info), is.character(rhs_terms))
  dropped <- character(0)

  # forcibly exclude any response/outcome terms the caller flagged
  if (!is.null(response_drop)) {
    hit <- intersect(rhs_terms, response_drop)
    if (length(hit)) {
      dropped <- c(dropped, hit)
      rhs_terms <- setdiff(rhs_terms, response_drop)
    }
  }

  # work on a local copy; droplevels factors and drop single-level factor terms
  keep_terms <- character(0)
  for (trm in rhs_terms) {
    if (!trm %in% names(info)) {
      warning(sprintf("build_design_guarded: term '%s' not found in info; dropping", trm))
      dropped <- c(dropped, trm)
      next
    }
    x <- info[[trm]]
    if (is.factor(x) || is.character(x)) {
      x <- droplevels(as.factor(x))
      info[[trm]] <- x
      if (nlevels(x) < 2L) {
        warning(sprintf("build_design_guarded: factor term '%s' has <2 observed levels; dropping", trm))
        dropped <- c(dropped, trm)
        next
      }
    }
    keep_terms <- c(keep_terms, trm)
  }

  formula_used <- if (length(keep_terms)) {
    paste("~", paste(keep_terms, collapse = " + "))
  } else {
    "~ 1"
  }
  design <- model.matrix(as.formula(formula_used), data = info)

  # ---- rank check + aliased-column drop -----------------------------------
  if (qr(design)$rank < ncol(design)) {
    ne <- .find_aliased_cols(design)
    if (length(ne)) {
      warning(sprintf("build_design_guarded: design rank-deficient; dropping aliased column(s): %s",
                      paste(ne, collapse = ", ")))
      dropped <- c(dropped, ne)
      design  <- design[, setdiff(colnames(design), ne), drop = FALSE]
    } else {
      warning("build_design_guarded: design rank-deficient but no aliased columns identified")
    }
  }

  list(design = design, formula_used = formula_used, dropped = dropped)
}

# ---------------------------------------------------------------------------
# .find_aliased_cols() — internal: return names of columns to drop to make a
# rank-deficient design full-rank. Primary path = limma::nonEstimable()
# (present in this env, limma 3.62.2). Fallbacks (in order): caret::findLinearCombos
# (if installed), then a manual QR-pivot column drop. Always returns the
# intercept-preserving choice that nonEstimable/QR pivots imply.
# ---------------------------------------------------------------------------
.find_aliased_cols <- function(design) {
  # 1) limma::nonEstimable (preferred; available in limma 3.62.2)
  if (exists("nonEstimable", where = asNamespace("limma"), inherits = FALSE)) {
    ne <- limma::nonEstimable(design)
    if (!is.null(ne) && length(ne)) return(ne)
    if (is.null(ne)) return(character(0))
  }
  # 2) caret::findLinearCombos fallback (only if caret is installed)
  if (requireNamespace("caret", quietly = TRUE)) {
    lc <- caret::findLinearCombos(design)
    if (length(lc$remove)) return(colnames(design)[lc$remove])
    return(character(0))
  }
  # 3) manual QR-pivot fallback: columns beyond the numerical rank, in pivot order
  qrd <- qr(design)
  r   <- qrd$rank
  if (r < ncol(design)) return(colnames(design)[qrd$pivot[(r + 1L):ncol(design)]])
  character(0)
}
