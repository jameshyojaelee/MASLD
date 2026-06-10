#!/usr/bin/env Rscript
# =============================================================================
# build_ranked_table.R
# -----------------------------------------------------------------------------
# SELECTION scorer for the factorial DEG-method benchmark. Applies the FROZEN
# calibration-gated composite rule from preregistration.yaml to pick the winning
# (engine x correction x n_sv) cell, with PROGRAMMATIC anti-circularity asserts.
#
# This is a PURE CONSUMER (login-safe, CSV math only). It does NOT fit any model.
# It reads:
#   * manifest_disease_vs_control*.csv  (cell identity + status; GLOB + dedupe)
#   * metrics_calibration_disease_vs_control.csv  (GATE; from --permute_null)
#   * metrics_loco.csv                            (RANK_R1; from the LOCO scorer)
#   * metrics_external.csv                        (RANK_R2; from the truth driver)
#   * metrics_report.csv                          (REPORT; from the truth driver)
#
# The FROZEN RULE (preregistration.yaml step1_gate / step2_rank / step3_report /
# truth_assignment is authoritative):
#   STEP 1 GATE  -- a cell PASSES iff permutation empirical type-I in [0.04,0.06].
#                   >0.06 -> FILTERED (inflated, removed from winner race).
#                   <0.04 -> over-conservative: EXCLUDED from race, RETAINED in
#                   reporting. Single pre-committed fallback: if ZERO cells pass,
#                   relax ONCE to [0.035,0.065], set gate_band_relaxed=TRUE,
#                   emit a DEVIATION note.
#   ELIGIBILITY  -- only k_sv in {be,na} can WIN. Fixed-k {2,5,10,20} are
#                   sensitivity-only: scored + reported but eligible=FALSE.
#   STEP 2 RANK  -- over gate-passing AND eligible cells:
#                   R1 = loco_repro_scalar (held-out LOCO reproducibility).
#                   R2 = mean z-score of {ot_enrichment_or_up (UP-directional),
#                        mouse_precision, coloc_or}.
#                   Authoritative combine = rank-of-ranks within survivors:
#                     composite = 0.5*rank(R1) +
#                                 0.5*mean(rank(ot_up),rank(mouse),rank(coloc));
#                     lowest composite wins; ties -> better R1 -> larger n_deg.
#                   Sensitivity variant = z-score average (reported, NOT
#                   authoritative).
#   STEP 3 REPORT -- positive-control P/R/F1, genetic-control capture,
#                   Govaere/Hoang concordance, heterogeneity. JOINED to the
#                   output but NEVER used in gate or rank.
#
# ANTI-CIRCULARITY (structural, asserted with stopifnot):
#   * gate_function's inputs are RESTRICTED to calibration columns; it cannot
#     even see loco/external/report columns.
#   * composite_function's inputs are RESTRICTED to loco + the 3 external
#     columns; it cannot see calibration or any REPORT column.
#   These restrictions are enforced as code (formals + a column-name firewall),
#   then re-asserted post-hoc against the column sets actually referenced.
#
# Env: micromamba run -n rnaseq Rscript build_ranked_table.R [--dir <degx_dir>]
# DO NOT COMMIT outputs (repo no-auto-commit rule).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# 0. paths + CLI
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEGX_DIR_DEFAULT <- file.path(BASE, "RNA-seq/results/degx_factorial")

CONTRAST <- "disease_vs_control"

# Gate band (FROZEN -- figS_degx_multimethod.R CAL_LO/CAL_HI).
GATE_LO <- 0.04
GATE_HI <- 0.06
# Single pre-committed fallback band.
GATE_LO_RELAX <- 0.035
GATE_HI_RELAX <- 0.065

# Winner-eligible k_sv values (data-driven k only; fixed-k are sensitivity-only).
ELIGIBLE_KSV <- c("be", "na")

# ----- FROZEN column-role partition (the anti-circularity firewall) ---------
# These vectors are the authoritative statement of which columns each step may
# touch. The gate/composite functions are restricted to their slice; an assert
# at the end re-verifies disjointness against truth_assignment.
GATE_COLS     <- c("perm_typeI_mean", "perm_typeI_mcse",
                   "null_pi0_mean", "lambda_null_mean", "n_perm")
RANK_R1_COLS  <- c("loco_repro_scalar")                  # primary R1 scalar
RANK_R1_SUB   <- c("loco_jaccard", "loco_lfc_spearman",
                   "loco_dir_concord", "loco_auc")       # constituents (R1 fallback)
RANK_R2_COLS  <- c("ot_enrichment_or_up", "mouse_precision", "coloc_or")
REPORT_COLS   <- c("pc_precision", "pc_recall", "pc_f1", "genetic_ctrl_captured",
                   "govaere_dir_pct", "govaere_lfc_r", "govaere_lfc_rho",
                   "hoang_lfc_r", "hoang_lfc_rho")

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  out <- list(dir = DEGX_DIR_DEFAULT)
  i <- 1
  while (i <= length(a)) {
    key <- sub("^--", "", a[i])
    if (key %in% names(out)) { out[[key]] <- a[i + 1]; i <- i + 2 } else { i <- i + 1 }
  }
  out
}
args <- parse_args()
DEGX_DIR <- args$dir
CELLS_DIR <- file.path(DEGX_DIR, "cells")
stopifnot(dir.exists(DEGX_DIR))

cat("=== build_ranked_table.R (SELECTION scorer) ===\n")
cat("dir      :", DEGX_DIR, "\n")
cat("contrast :", CONTRAST, "\n")
cat("gate band:", sprintf("[%.3f, %.3f]", GATE_LO, GATE_HI), "\n\n")

read_csv_safe <- function(path, what) {
  if (!file.exists(path)) {
    cat(sprintf("[input] %-22s MISSING -> %s\n", what, basename(path)))
    return(NULL)
  }
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) {
    cat(sprintf("[input] %-22s EMPTY/UNREADABLE -> %s\n", what, basename(path)))
    return(NULL)
  }
  cat(sprintf("[input] %-22s OK (%d rows) -> %s\n", what, nrow(dt), basename(path)))
  dt
}

# ---------------------------------------------------------------------------
# 1. MERGE manifests (glob concurrent shards -> concat -> dedupe by cell_id)
# ---------------------------------------------------------------------------
manifest_glob <- list.files(DEGX_DIR,
                            pattern = "^manifest_disease_vs_control.*\\.csv$",
                            full.names = TRUE)
if (length(manifest_glob) == 0) {
  stop("No manifest_disease_vs_control*.csv found in ", DEGX_DIR,
       " -- cannot enumerate cells. Run the engine factorial first.")
}
cat(sprintf("[manifest] %d shard manifest(s) globbed\n", length(manifest_glob)))
mf_list <- lapply(manifest_glob, function(f) {
  d <- fread(f)
  d[, .src_manifest := basename(f)]
  # Guard 1: capture each manifest file's mtime so dedupe can prefer the
  # MOST-RECENTLY-MODIFIED manifest. After a re-run (e.g. the RUVseq fix
  # 2026-06-07), a cell can appear in BOTH a stale shard manifest (broken
  # ~sex+W+group design, ncov=6) and a newer fix shard manifest
  # (~dataset+sex+W+group, ncov=10). The newer manifest's mtime wins.
  d[, ._mtime := file.info(f)$mtime]
  d
})
manifest <- rbindlist(mf_list, use.names = TRUE, fill = TRUE)

# ---- deterministic dedupe (belt + suspenders) -----------------------------
# Detect which cell_ids appear in >1 manifest row BEFORE collapsing, so we can
# WARN about exactly what got superseded.
dup_ids <- manifest[, .N, by = cell_id][N > 1L, cell_id]

# Tiebreak helpers (used only when mtimes are equal/missing):
#   - RUV-fix marker: message contains "ON TOP of cohort" (the fixed design).
#   - then higher n_covariates (the fixed C6 design adds the cohort terms).
manifest[, ._is_ok      := as.integer(!is.na(status) & status == "ok")]
manifest[, ._ruvfix     := as.integer(!is.na(message) &
                                        grepl("ON TOP of cohort", message, fixed = TRUE))]
manifest[, ._ncov       := suppressWarnings(as.numeric(n_covariates))]
manifest[is.na(._ncov), ._ncov := -Inf]
manifest[, ._mt         := suppressWarnings(as.numeric(._mtime))]
manifest[is.na(._mt), ._mt := -Inf]

# Priority order (all DESC so the FIRST row per cell_id is the keeper):
#   1. status=="ok"            -- never let a failed row supersede an ok one
#   2. most-recent manifest    -- Guard 1 (newer re-run wins automatically)
#   3. RUV-fix marker          -- Guard 2a (explicit fixed-design preference)
#   4. higher n_covariates     -- Guard 2b (fixed C6 design has more covariates)
setorder(manifest, cell_id, -._is_ok, -._mt, -._ruvfix, -._ncov)
manifest <- manifest[!duplicated(cell_id)]

# Guard 3: printed WARNING listing every superseded cell_id.
if (length(dup_ids)) {
  cat(sprintf("[manifest] WARNING: deduped %d cell_id(s) with multiple manifest rows; kept most-recent.\n",
              length(dup_ids)))
  kept <- manifest[cell_id %in% dup_ids,
                   .(cell_id, kept_src = .src_manifest, kept_ncov = n_covariates,
                     kept_mtime = ._mtime)]
  print(kept[order(cell_id)], row.names = FALSE)
}

manifest[, c(".src_manifest", "._mtime", "._is_ok", "._ruvfix",
             "._ncov", "._mt") := NULL]
cat(sprintf("[manifest] %d unique cells after dedupe (%d ok)\n",
            nrow(manifest), sum(manifest$status == "ok", na.rm = TRUE)))

# Score only status == "ok" cells.
ok <- manifest[status == "ok"]
if (nrow(ok) == 0) stop("No status=='ok' cells in merged manifest.")

# Identity columns carried into outputs.
id_cols <- intersect(c("cell_id", "engine", "correction_id", "k_sv",
                       "k_sv_source", "k_be", "k_leek", "n_genes",
                       "n_sig_padj05", "n_sig_lfc05", "n_covariates",
                       "n_samples", "seconds"),
                     names(ok))
base <- ok[, ..id_cols]

# ---------------------------------------------------------------------------
# 2. LOAD metric tables (graceful degradation -- any may be absent)
# ---------------------------------------------------------------------------
# Calibration is written PER-CELL by the DEGperm array: each task emits one
# single-row file metrics_calibration_<contrast>__<cell_id>.csv (this replaced a
# single shared metrics_calibration_<contrast>.csv, which clobbered under
# concurrent array tasks). GLOB + rbind + dedupe-by-cell_id, exactly like the
# manifest shards above. If NO per-cell file matches, fall through to the
# existing graceful-degradation path (calib = NULL -> gate_status="pending").
load_calibration_glob <- function() {
  pat  <- sprintf("^metrics_calibration_%s__.*\\.csv$", CONTRAST)
  cfs  <- list.files(DEGX_DIR, pattern = pat, full.names = TRUE)
  if (length(cfs) == 0) {
    cat(sprintf("[input] %-22s MISSING -> no metrics_calibration_%s__*.csv\n",
                "calibration (GATE)", CONTRAST))
    return(NULL)
  }
  parts <- lapply(cfs, function(f) tryCatch(fread(f), error = function(e) NULL))
  parts <- Filter(function(d) !is.null(d) && nrow(d) > 0, parts)
  if (length(parts) == 0) {
    cat(sprintf("[input] %-22s EMPTY/UNREADABLE -> %d per-cell file(s)\n",
                "calibration (GATE)", length(cfs)))
    return(NULL)
  }
  cb <- rbindlist(parts, use.names = TRUE, fill = TRUE)
  stopifnot("cell_id" %in% names(cb))
  cb <- cb[!duplicated(cell_id)]   # one row per cell; dedupe defensively
  cat(sprintf("[input] %-22s OK (%d per-cell file(s) -> %d unique cells) -> metrics_calibration_%s__*.csv\n",
              "calibration (GATE)", length(cfs), nrow(cb), CONTRAST))
  cb
}
calib <- load_calibration_glob()
loco  <- read_csv_safe(file.path(DEGX_DIR, "metrics_loco.csv"), "loco (RANK_R1)")
ext   <- read_csv_safe(file.path(DEGX_DIR, "metrics_external.csv"), "external (RANK_R2)")
rep_  <- read_csv_safe(file.path(DEGX_DIR, "metrics_report.csv"), "report (REPORT)")
cat("\n")

have_calib <- !is.null(calib)
have_loco  <- !is.null(loco)
have_ext   <- !is.null(ext)

# ---------------------------------------------------------------------------
# 3. ELIGIBILITY (k_sv in {be,na}); fixed-k are sensitivity-only.
# ---------------------------------------------------------------------------
base[, eligible := k_sv %in% ELIGIBLE_KSV]

# ===========================================================================
# 4. GATE FUNCTION  -- inputs RESTRICTED to calibration columns ONLY.
#    Anti-circularity: this function's formal arguments are exactly the gate
#    columns; it has no access to loco/external/report data. It returns the
#    pass/status decision plus the relaxation flag.
# ===========================================================================
gate_function <- function(perm_typeI_mean) {
  # FROZEN: PASS iff CAL_LO <= type-I <= CAL_HI. Below uses ONLY perm_typeI_mean.
  # First pass at the frozen band; the single fallback is decided by the caller
  # AFTER seeing whether zero cells passed (pre-committed, logged as DEVIATION).
  n <- length(perm_typeI_mean)
  res <- data.table(
    gate_pass   = rep(NA, n),
    gate_status = rep("pending", n),
    perm_typeI_mean = perm_typeI_mean
  )
  ok_idx <- is.finite(perm_typeI_mean)
  pass   <- ok_idx & perm_typeI_mean >= GATE_LO & perm_typeI_mean <= GATE_HI
  infl   <- ok_idx & perm_typeI_mean > GATE_HI
  cons   <- ok_idx & perm_typeI_mean < GATE_LO
  res[pass, c("gate_pass", "gate_status") := list(TRUE,  "pass")]
  res[infl, c("gate_pass", "gate_status") := list(FALSE, "inflated")]
  res[cons, c("gate_pass", "gate_status") := list(FALSE, "conservative")]
  res
}
# Hard restriction: the gate function may reference ONLY GATE_COLS.
GATE_FN_INPUTS <- names(formals(gate_function))   # -> "perm_typeI_mean"
stopifnot(all(GATE_FN_INPUTS %in% GATE_COLS))

# Apply the gate (or mark pending if calibration absent).
gate_band_relaxed <- FALSE
deviation_notes <- character(0)
if (have_calib) {
  stopifnot("cell_id" %in% names(calib), "perm_typeI_mean" %in% names(calib))
  cal <- merge(base[, .(cell_id)], calib[, .(cell_id, perm_typeI_mean)],
               by = "cell_id", all.x = TRUE, sort = FALSE)
  g <- gate_function(cal$perm_typeI_mean)
  g[, cell_id := cal$cell_id]

  # Pre-committed single fallback: if ZERO cells pass the frozen band, relax once.
  n_pass <- sum(g$gate_pass %in% TRUE)
  if (n_pass == 0L && any(is.finite(g$perm_typeI_mean))) {
    gate_band_relaxed <- TRUE
    # Relaxed band is still a GATE-only decision (uses only perm_typeI_mean).
    ptm <- g$perm_typeI_mean
    ok_idx <- is.finite(ptm)
    pass   <- ok_idx & ptm >= GATE_LO_RELAX & ptm <= GATE_HI_RELAX
    infl   <- ok_idx & ptm > GATE_HI_RELAX
    cons   <- ok_idx & ptm < GATE_LO_RELAX
    g[, gate_pass := NA]; g[, gate_status := "pending"]
    g[pass, c("gate_pass", "gate_status") := list(TRUE,  "pass")]
    g[infl, c("gate_pass", "gate_status") := list(FALSE, "inflated")]
    g[cons, c("gate_pass", "gate_status") := list(FALSE, "conservative")]
    msg <- sprintf("DEVIATION: gate band relaxed to [%.3f,%.3f] (pre-committed single fallback).",
                   GATE_LO_RELAX, GATE_HI_RELAX)
    deviation_notes <- c(deviation_notes, msg)
    cat("\n!! ", msg, "\n", sep = "")
  }
  base <- merge(base, g[, .(cell_id, perm_typeI_mean, gate_pass, gate_status)],
                by = "cell_id", all.x = TRUE, sort = FALSE)
  # Carry the GATE secondary diagnostics (mcse / pi0 / lambda / n_perm) for the
  # REPORTING long table only. These are NOT used by gate_function -- the gate
  # decision above consumed perm_typeI_mean exclusively.
  diag_cols <- intersect(c("perm_typeI_mcse", "null_pi0_mean",
                           "lambda_null_mean", "n_perm"), names(calib))
  if (length(diag_cols)) {
    base <- merge(base, calib[, c("cell_id", diag_cols), with = FALSE],
                  by = "cell_id", all.x = TRUE, sort = FALSE)
  }
} else {
  base[, perm_typeI_mean := NA_real_]
  base[, gate_pass := NA]
  base[, gate_status := "pending"]
}
for (cc in c("perm_typeI_mcse", "null_pi0_mean", "lambda_null_mean", "n_perm"))
  if (!cc %in% names(base)) base[, (cc) := NA_real_]
base[, gate_band_relaxed := gate_band_relaxed]

# ===========================================================================
# 5. R1 / R2 metric assembly (join loco + external; NO calibration, NO report).
# ===========================================================================
# R1 scalar (held-out LOCO reproducibility). If loco_repro_scalar absent, derive
# it from the four constituents = mean(jaccard,(spearman+1)/2,dir_concord/100,auc).
if (have_loco) {
  stopifnot("cell_id" %in% names(loco))
  lc <- copy(loco)
  if (!"loco_repro_scalar" %in% names(lc)) {
    need <- RANK_R1_SUB
    if (all(need %in% names(lc))) {
      lc[, loco_repro_scalar := rowMeans(cbind(
        loco_jaccard,
        (loco_lfc_spearman + 1) / 2,
        loco_dir_concord / 100,
        loco_auc), na.rm = TRUE)]
      cat("[loco] loco_repro_scalar absent -> derived from 4 constituents.\n")
    } else {
      lc[, loco_repro_scalar := NA_real_]
      cat("[loco] loco_repro_scalar absent AND constituents incomplete -> NA.\n")
    }
  }
  base <- merge(base, lc[, .(cell_id, loco_repro_scalar)],
                by = "cell_id", all.x = TRUE, sort = FALSE)
} else {
  base[, loco_repro_scalar := NA_real_]
}

# R2 external constituents.
if (have_ext) {
  stopifnot("cell_id" %in% names(ext))
  ext_keep <- intersect(c("cell_id", "n_deg", RANK_R2_COLS), names(ext))
  base <- merge(base, ext[, ..ext_keep], by = "cell_id", all.x = TRUE, sort = FALSE)
}
for (cc in c("n_deg", RANK_R2_COLS)) if (!cc %in% names(base)) base[, (cc) := NA_real_]

# ===========================================================================
# 6. COMPOSITE FUNCTION  -- inputs RESTRICTED to R1 + the 3 R2 columns ONLY.
#    Anti-circularity: formals are exactly {loco_repro_scalar,
#    ot_enrichment_or_up, mouse_precision, coloc_or, n_deg(tiebreak)}.
#    It has NO access to calibration columns or any REPORT column.
# ===========================================================================
composite_function <- function(loco_repro_scalar,
                               ot_enrichment_or_up,
                               mouse_precision,
                               coloc_or,
                               n_deg) {
  # Operates ONLY on the survivor subset passed in (the caller restricts to
  # gate-pass & eligible cells so ranks are computed WITHIN survivors).
  n <- length(loco_repro_scalar)
  rank_best <- function(x) {
    # rank ascending so 1 = best for "higher is better" metrics; NA -> worst.
    r <- rank(-x, na.last = TRUE, ties.method = "average")
    r
  }
  R1_rank <- rank_best(loco_repro_scalar)
  r_ot    <- rank_best(ot_enrichment_or_up)
  r_mouse <- rank_best(mouse_precision)
  r_coloc <- rank_best(coloc_or)
  R2_rank <- rowMeans(cbind(r_ot, r_mouse, r_coloc))
  composite_ranksum <- 0.5 * R1_rank + 0.5 * R2_rank

  # Sensitivity variant (NOT authoritative): z-score average.
  zsc <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }
  z_R1 <- zsc(loco_repro_scalar)
  z_R2 <- rowMeans(cbind(zsc(ot_enrichment_or_up), zsc(mouse_precision), zsc(coloc_or)),
                   na.rm = TRUE)
  composite_zavg <- 0.5 * z_R1 + 0.5 * z_R2     # HIGHER = better (variant)

  data.table(R1_rank = R1_rank, R2_rank = R2_rank,
             composite_ranksum = composite_ranksum,
             composite_zavg = composite_zavg)
}
# Hard restriction: composite function may reference ONLY R1 + R2 columns
# (plus n_deg, a descriptive tiebreak, which is neither calibration nor report).
COMPOSITE_FN_INPUTS <- names(formals(composite_function))
ALLOWED_COMPOSITE_INPUTS <- c(RANK_R1_COLS, RANK_R2_COLS, "n_deg")
stopifnot(all(COMPOSITE_FN_INPUTS %in% ALLOWED_COMPOSITE_INPUTS))

# ---------------------------------------------------------------------------
# 7. Apply ranking over SURVIVORS = gate-pass AND eligible.
#    Requires both calibration (gate) and loco+external (rank) to be present.
# ---------------------------------------------------------------------------
base[, c("R1_rank", "R2_rank", "composite_ranksum", "composite_zavg",
         "winner_rank", "is_winner") :=
       list(NA_real_, NA_real_, NA_real_, NA_real_, NA_integer_, FALSE)]

rank_done <- FALSE
can_rank <- have_calib && have_loco && have_ext
if (can_rank) {
  surv <- base[gate_pass %in% TRUE & eligible == TRUE]
  if (nrow(surv) >= 1L) {
    comp <- composite_function(
      loco_repro_scalar   = surv$loco_repro_scalar,
      ot_enrichment_or_up = surv$ot_enrichment_or_up,
      mouse_precision     = surv$mouse_precision,
      coloc_or            = surv$coloc_or,
      n_deg               = surv$n_deg)
    surv[, names(comp) := comp]
    # Winner ordering: lowest composite_ranksum; ties -> better (lower) R1_rank,
    # then larger n_deg.
    setorder(surv, composite_ranksum, R1_rank, -n_deg)
    surv[, winner_rank := seq_len(.N)]
    surv[, is_winner := winner_rank == 1L]
    # write back by cell_id
    rk <- surv[, .(cell_id, R1_rank, R2_rank, composite_ranksum, composite_zavg,
                   winner_rank, is_winner)]
    base[rk, on = "cell_id",
         c("R1_rank", "R2_rank", "composite_ranksum", "composite_zavg",
           "winner_rank", "is_winner") :=
           list(i.R1_rank, i.R2_rank, i.composite_ranksum, i.composite_zavg,
                i.winner_rank, i.is_winner)]
    rank_done <- TRUE
  } else {
    cat("[rank] No gate-pass & eligible survivors -> no winner selected.\n")
  }
} else {
  miss <- c(if (!have_calib) "calibration", if (!have_loco) "loco",
            if (!have_ext) "external")
  cat(sprintf("[rank] PENDING -- missing input(s): %s. Gate/eligibility applied; rank skipped.\n",
              paste(miss, collapse = ", ")))
}

# ===========================================================================
# 8. ANTI-CIRCULARITY ASSERTS (hard stopifnot, printed).
#    Verify structurally that:
#      (a) the gate consumed ONLY calibration columns;
#      (b) the composite consumed ONLY loco + the 3 external columns (+n_deg);
#      (c) NO REPORT column ever entered gate OR composite;
#      (d) NO rank/loco/external column ever entered the gate.
# ===========================================================================
cat("\n==================== ANTI-CIRCULARITY ====================\n")

# (a) gate inputs subset of GATE_COLS, disjoint from rank/report.
assert_gate_clean <- all(GATE_FN_INPUTS %in% GATE_COLS) &&
  length(intersect(GATE_FN_INPUTS, c(RANK_R1_COLS, RANK_R1_SUB, RANK_R2_COLS, REPORT_COLS))) == 0L
cat(sprintf("  gate inputs               = {%s}\n", paste(GATE_FN_INPUTS, collapse = ", ")))
stopifnot("gate references a non-calibration column" = assert_gate_clean)

# (b) composite inputs subset of {R1,R2,n_deg}; disjoint from gate + report.
assert_comp_clean <- all(COMPOSITE_FN_INPUTS %in% ALLOWED_COMPOSITE_INPUTS) &&
  length(intersect(COMPOSITE_FN_INPUTS, GATE_COLS)) == 0L &&
  length(intersect(COMPOSITE_FN_INPUTS, REPORT_COLS)) == 0L
cat(sprintf("  composite inputs          = {%s}\n", paste(COMPOSITE_FN_INPUTS, collapse = ", ")))
stopifnot("composite references a calibration column" =
            length(intersect(COMPOSITE_FN_INPUTS, GATE_COLS)) == 0L)
stopifnot("composite references a REPORT column" =
            length(intersect(COMPOSITE_FN_INPUTS, REPORT_COLS)) == 0L)
stopifnot("composite references a disallowed column" = assert_comp_clean)

# (c) no REPORT column in either decision-maker.
stopifnot("REPORT column leaked into gate" =
            length(intersect(GATE_FN_INPUTS, REPORT_COLS)) == 0L)
stopifnot("REPORT column leaked into composite" =
            length(intersect(COMPOSITE_FN_INPUTS, REPORT_COLS)) == 0L)

# (d) no rank/loco/external column in the gate.
stopifnot("rank/external column leaked into gate" =
            length(intersect(GATE_FN_INPUTS,
                             c(RANK_R1_COLS, RANK_R1_SUB, RANK_R2_COLS))) == 0L)

# (e) the partition is mutually exclusive (truth_assignment: exactly one role).
all_role_cols <- list(GATE = GATE_COLS,
                      RANK_R1 = c(RANK_R1_COLS, RANK_R1_SUB),
                      RANK_R2 = RANK_R2_COLS,
                      REPORT = REPORT_COLS)
flat <- unlist(all_role_cols, use.names = FALSE)
stopifnot("a column is assigned to two roles" = !any(duplicated(flat)))

cat("ANTI-CIRCULARITY: PASS\n")
cat("==========================================================\n")

# ---------------------------------------------------------------------------
# 9. JOIN the REPORT columns (AFTER selection; never referenced above).
# ---------------------------------------------------------------------------
if (!is.null(rep_)) {
  rep_keep <- intersect(c("cell_id", REPORT_COLS), names(rep_))
  base <- merge(base, rep_[, ..rep_keep], by = "cell_id", all.x = TRUE, sort = FALSE)
}
for (cc in REPORT_COLS) if (!cc %in% names(base)) base[, (cc) := NA_real_]

# ---------------------------------------------------------------------------
# 10. WRITE degmethod_ranked_winners.csv  (one row per ok cell).
#     Sorted by winner_rank (NA -> last).
# ---------------------------------------------------------------------------
out_cols <- c(
  id_cols,
  "perm_typeI_mean", "gate_pass", "gate_status", "gate_band_relaxed", "eligible",
  "loco_repro_scalar", "ot_enrichment_or_up", "mouse_precision", "coloc_or",
  "n_deg",
  "R1_rank", "R2_rank", "composite_ranksum", "composite_zavg",
  "winner_rank", "is_winner",
  REPORT_COLS)
out_cols <- intersect(out_cols, names(base))   # tolerate any absent id col
winners <- base[, ..out_cols]
setorder(winners, winner_rank, na.last = TRUE)

winners_path <- file.path(DEGX_DIR, "degmethod_ranked_winners.csv")
fwrite(winners, winners_path)

# ---------------------------------------------------------------------------
# 11. WRITE degmethod_scoring_long.csv  (cell_id x metric x value x family x role)
#     Roles match prereg truth_assignment.
# ---------------------------------------------------------------------------
make_long <- function(dt, cols, family, role) {
  cols <- intersect(cols, names(dt))
  if (length(cols) == 0) return(NULL)
  d2 <- copy(dt[, c("cell_id", cols), with = FALSE])
  for (col in cols) set(d2, j = col, value = as.numeric(d2[[col]]))
  m <- melt(d2, id.vars = "cell_id", variable.name = "metric",
            value.name = "value", variable.factor = FALSE)
  m[, family := family]; m[, role := role]
  m
}
long_parts <- list(
  make_long(base, GATE_COLS,                     "calibration", "GATE"),
  make_long(base, c("loco_repro_scalar"),        "loco",        "RANK_R1"),
  make_long(base, RANK_R2_COLS,                  "external",    "RANK_R2"),
  make_long(base, REPORT_COLS,                   "report",      "REPORT"),
  # descriptive context (never gate/rank): identity counts + composite outputs.
  make_long(base, c("n_deg", "n_sig_padj05", "n_sig_lfc05", "n_covariates",
                    "n_samples", "R1_rank", "R2_rank", "composite_ranksum",
                    "composite_zavg", "winner_rank"),
            "descriptive", "DESCRIPTIVE")
)
scoring_long <- rbindlist(Filter(Negate(is.null), long_parts), use.names = TRUE)
setorder(scoring_long, cell_id, role, family, metric)
long_path <- file.path(DEGX_DIR, "degmethod_scoring_long.csv")
fwrite(scoring_long, long_path)

# ---------------------------------------------------------------------------
# 12. SUMMARY + WINNER line.
# ---------------------------------------------------------------------------
n_cells     <- nrow(base)
n_gatepass  <- sum(base$gate_pass %in% TRUE)
n_eligible  <- sum(base$eligible == TRUE)
n_filtered  <- sum(base$gate_status == "inflated")
n_cons      <- sum(base$gate_status == "conservative")
n_pending   <- sum(base$gate_status == "pending")

cat("\n=== wrote ===\n")
cat("  ranked winners :", winners_path, sprintf("(%d cells)\n", nrow(winners)))
cat("  scoring long   :", long_path, sprintf("(%d rows)\n", nrow(scoring_long)))

cat("\n=== SUMMARY ===\n")
cat(sprintf("  cells (ok)          : %d\n", n_cells))
cat(sprintf("  gate-pass           : %d\n", n_gatepass))
cat(sprintf("  gate-inflated (FILT): %d\n", n_filtered))
cat(sprintf("  gate-conservative   : %d\n", n_cons))
cat(sprintf("  gate-pending        : %d\n", n_pending))
cat(sprintf("  eligible (k in {be,na}): %d\n", n_eligible))
cat(sprintf("  gate-pass & eligible: %d\n", sum(base$gate_pass %in% TRUE & base$eligible == TRUE)))
cat(sprintf("  gate_band_relaxed   : %s\n", gate_band_relaxed))
if (length(deviation_notes)) for (m in deviation_notes) cat("  ", m, "\n", sep = "")

cat("\n=== WINNER ===\n")
if (rank_done && any(base$is_winner %in% TRUE)) {
  w <- winners[is_winner == TRUE][1]
  cat(sprintf("  cell_id : %s\n", w$cell_id))
  cat("  full metric line:\n")
  print(t(w), quote = FALSE)
} else {
  cat("  No winner selected (rank ", if (!can_rank) "PENDING (missing inputs)" else "produced no eligible survivor", ").\n", sep = "")
}
cat("\nDone.\n")
