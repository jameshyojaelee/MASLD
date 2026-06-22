#!/usr/bin/env bash
# CI gate for the 2026-06-13 mega-review A7 fix (prediction & F-stage leakage).
# FAILS (exit 1) if any ACTIVE script reintroduces one of the three confirmed
# leakage regressions that A7 remediated:
#
#   (a) HARD-CODED "reproduces the canonical" verdict in the NMF compare scripts
#       (94c_compare_refit_vs_canonical.R). The honest verdict MUST be COMPUTED
#       from ari_val / per-program Jaccard. A bare declarative "the re-fit
#       reproduces the canonical decomposition" literal that is NOT derived from
#       those metrics is a regression (k=6 is seed-unstable: ARI~0.30, agree
#       0.504, min Jaccard 0.06). Lines that build the verdict inside a sprintf
#       from ari_val / min_jacc / agree (the data-driven branch) are NOT flagged.
#
#   (b) FULL-COHORT div_* feature screen upstream of the LOCO-CV loop in the
#       bifurcation / prognosis predictor (117_bifurcation_divergence.py). The
#       leakage-free panel is recomputed per outer fold
#       (divergence_genes_per_fold.csv). A full-cohort detect_divergence_genes()
#       result wired into the *_per_fold / feature_matrix predictor input is a
#       regression. The retained DESCRIPTIVE full-cohort ranking is allowed and
#       is recognised by its leakage-fix annotation on a nearby line.
#
#   (c) ANY active use of F_stage_augmented (OR its remediated alias
#       F_stage_augmented_clean) as the F-stage axis / stage column / donor-
#       inclusion gate in WS3 scripts 346-349. The augmented column is LEAKED
#       (scVI QWK 0.74-0.76 -> jackknife 0.286 / held-out Andrews 0.0); the
#       cascade routes the axis through F_stage_inferred (343b, bootstrap-gated).
#       Three leak shapes are flagged: (i) a BARE assignment of
#       F_stage_augmented[_clean] to a stage-axis/column sink variable;
#       (ii) the COHORT_TAG switch indirection `augmented = "F_stage_augmented
#       [_clean]"` (348b's switch arm — the augmented arm still maps the leaked
#       column even when the default is documented); (iii) a pandas
#       `dropna(subset=["F_stage_augmented[_clean]"])` donor-inclusion filter
#       (348_pseudobulk gates the kept donor set on the leaked axis).
#       NOTE: F_stage_augmented_clean is the column actually consumed
#       downstream, so the previous `([^_]|$)` boundary that spared `_clean`
#       was a FALSE-GREEN — it is REMOVED. Provenance references (covar lists,
#       comments, the `_clean := F_stage_augmented else NA` fallback mapping,
#       the `donor_meta["..._clean"] = donor_meta["F_stage_augmented"]` column
#       set-up, .notna()/.astype()/print echoes, bootstrap-gate warnings) are
#       NOT flagged; the three leak shapes above ARE.
#
# Escape hatch: an intentional reference must be annotated on the SAME LINE with
# the sentinel `# A7-OK-provenance` (mirrors the `# C2-OK-sensitivity` sentinel
# used by check_no_dream_cols.sh).
#
# Usage: bash scripts/ci/check_no_leakage.sh   (run from repo root)
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 2

# Shared excludes: archives, docs, backups, build dirs, and this gate itself.
COMMON_EXCLUDE='docs/audit/|docs/archive/|_backup|_premigration|\.bak:|/\.next/|/out/|/node_modules/|/dist/'
SENTINEL='A7-OK-provenance'

fail=0

emit_fail () {  # $1 = human label, $2 = hits blob
  local label="$1" hits="$2"
  local n
  n=$(printf '%s\n' "$hits" | grep -c . )
  echo "FAIL [$label]: $n offending line(s):"
  printf '%s\n' "$hits"
  echo
  fail=1
}

# ---------------------------------------------------------------------------
# (a) Hard-coded "reproduces the canonical" NMF verdict.
#     Flag the declarative verdict literal; spare lines that compute it from the
#     real metrics (ari_val / ARI_REPRO / min_jacc / agree / reproduces var).
# ---------------------------------------------------------------------------
PATTERN_A='re-?fit[[:space:]]+(REPRODUCES|reproduces)[[:space:]]+(the[[:space:]]+)?canonical|reproduces[[:space:]]+the[[:space:]]+canonical[[:space:]]+decomposition'
hits_a=$(grep -rEn \
  --include='9[0-9][a-z]_*compare*.R' --include='9[0-9][a-z]_*refit*.R' \
  --include='*nmf*compare*.R' --include='*compare_refit*.R' \
  --exclude-dir=archive --exclude-dir=.claude --exclude-dir=worktrees --exclude-dir=.git \
  --exclude-dir=.next --exclude-dir=out --exclude-dir=node_modules --exclude-dir=dist \
  --exclude='check_no_leakage.sh' \
  "$PATTERN_A" . 2>/dev/null \
  | grep -vE "$COMMON_EXCLUDE" \
  | grep -v "$SENTINEL" \
  | grep -vE 'ari_val|ARI_REPRO|min_jacc|JACC_REPRO|\bagree\b|if \(reproduces\)|reproduces <-|ARI=%|Jaccard=%|%\.[0-9]f' \
  || true)
[ -n "$hits_a" ] && emit_fail "94c hard-coded reproduces verdict" "$hits_a"

# ---------------------------------------------------------------------------
# (b) Full-cohort div_* feature screen wired into the CV predictor input.
#     Flag a full-cohort detect_divergence_genes()/div_df result that is written
#     to or assigned into the predictor panel (*_per_fold / feature_matrix /
#     div_panel) WITHOUT the leakage-fix annotation. The retained descriptive
#     full-cohort line carries a DESCRIPTIVE / leakage-fix marker and is spared.
# ---------------------------------------------------------------------------
PATTERN_B='(divergence_genes_per_fold|feature_matrix|div_panel|div_features)[^=]*=[[:space:]]*detect_divergence_genes\(|detect_divergence_genes\([^)]*\)[[:space:]]*#[[:space:]]*full.?cohort'
hits_b=$(grep -rEn \
  --include='11[0-9]_*.py' --include='1[5-9][0-9]_*.py' --include='1[5-9][0-9]_*.R' \
  --include='*bifurcation*divergence*.py' \
  --exclude-dir=archive --exclude-dir=.claude --exclude-dir=worktrees --exclude-dir=.git \
  --exclude-dir=.next --exclude-dir=out --exclude-dir=node_modules --exclude-dir=dist \
  --exclude='check_no_leakage.sh' \
  "$PATTERN_B" . 2>/dev/null \
  | grep -vE "$COMMON_EXCLUDE" \
  | grep -v "$SENTINEL" \
  | grep -viE 'DESCRIPTIVE|leakage.?free|leakage fix|MUST NOT be used|training cohorts only' \
  || true)
[ -n "$hits_b" ] && emit_fail "117 full-cohort div_* screen upstream of CV" "$hits_b"

# ---------------------------------------------------------------------------
# (c) Active use of F_stage_augmented / F_stage_augmented_clean as the stage
#     axis, COHORT_TAG switch arm, or donor-inclusion gate in WS3 346-349.
#     Three leak shapes (alternation branches):
#       1. <stage_sink> (<-|:=|=) "F_stage_augmented[_clean]"   (bare axis assign)
#       2. augmented = "F_stage_augmented[_clean]"              (switch arm — 348b)
#       3. dropna(subset=["F_stage_augmented[_clean]"])         (pandas donor gate)
#     `_clean` is INCLUDED (it is the column actually consumed downstream; the
#     old `([^_]|$)` boundary that spared it was the false-green this gate fixes).
#     Spared by the comment-skip (#-lines), the SENTINEL, and by NOT matching:
#     the `_clean := F_stage_augmented else NA` fallback, the
#     `donor_meta["..._clean"] = donor_meta["F_stage_augmented"]` col set-up,
#     and .notna()/.astype()/print echoes.
# ---------------------------------------------------------------------------
PATTERN_C='(fstage_col|f_stage_col|F_STAGE_COL|stage_col|stage_axis|axis_term|F_STAGE_VAR|F_STAGE_AXIS|F_stage_numeric)[[:space:]]*(<-|:=|=)[[:space:]]*["'"'"']?F_stage_augmented(_clean)?["'"'"']?|(^|[^_[:alnum:]])augmented[[:space:]]*=[[:space:]]*["'"'"']F_stage_augmented(_clean)?["'"'"']|dropna\([[:space:]]*subset[[:space:]]*=[[:space:]]*\[?["'"'"']F_stage_augmented(_clean)?["'"'"']'
hits_c=$(grep -rEn \
  --include='34[6-9]*.R' --include='34[6-9]*.py' \
  --exclude-dir=archive --exclude-dir=.claude --exclude-dir=worktrees --exclude-dir=.git \
  --exclude-dir=.next --exclude-dir=out --exclude-dir=node_modules --exclude-dir=dist \
  --exclude='check_no_leakage.sh' \
  "$PATTERN_C" . 2>/dev/null \
  | grep -vE "$COMMON_EXCLUDE" \
  | grep -v "$SENTINEL" \
  | grep -vE '^[^:]+:[0-9]+:[[:space:]]*#' \
  || true)
[ -n "$hits_c" ] && emit_fail "WS3 346-349 active F_stage_augmented[_clean] axis/switch/dropna" "$hits_c"

# ---------------------------------------------------------------------------
if [ "$fail" -ne 0 ]; then
  echo "check_no_leakage: FAIL — one or more A7 leakage regressions detected."
  echo "If a flagged line is a deliberate provenance/descriptive reference,"
  echo "annotate it on the same line with the sentinel: # $SENTINEL"
  exit 1
fi
echo "PASS: no NMF hard-coded verdict, no full-cohort div_* screen upstream of CV, no active F_stage_augmented axis in WS3 346-349."
exit 0
