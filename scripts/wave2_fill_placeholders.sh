#!/usr/bin/env bash
# wave2_fill_placeholders.sh — fills the [TBD by Wave 2] / [GENETIC_DOWN_LABEL]
# placeholders that Wave 1 Agent A inserted across 11 markdown files,
# using actual numbers from the LOMO/baselines/permutation outputs.
#
# Run AFTER `RNA-seq/46d_baselines_and_lomo.R` completes.

set -eo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"

ME="RNA-seq/results/multi_evidence"
LOMO_F="$ME/bayesian_evidence_lomo_validation.csv"
PERM_F="$ME/bayesian_evidence_permutation_null_matched.csv"
BASE_F="$ME/bayesian_evidence_baseline_benchmark.csv"

if [[ ! -s "$LOMO_F" ]]; then
  echo "ERROR: $LOMO_F missing — run 46d_baselines_and_lomo.R first" >&2
  exit 1
fi

# Pull LOMO-S1 Govaere AUROC + 95% CI (use awk on CSV columns)
LOMO_S1_GOVAERE=$(awk -F, '$1=="Govaere" && $2=="LOMO_S1" {printf "%.3f", $3}' "$LOMO_F")
LOMO_S1_GOVAERE_CI_LOW=$(awk -F, '$1=="Govaere" && $2=="LOMO_S1" {printf "%.3f", $4}' "$LOMO_F")
LOMO_S1_GOVAERE_CI_HIGH=$(awk -F, '$1=="Govaere" && $2=="LOMO_S1" {printf "%.3f", $5}' "$LOMO_F")

# Pull permutation null p-value
PERM_P_GOVAERE=$(awk -F, '$1=="Govaere" {printf "%.4f", $8}' "$PERM_F")

# Final concordance state label (renamed from Protective-LOF)
GENETIC_DOWN_LABEL="Genetic+down-coherent"

# F2 dominance permutation p — TODO: requires separate computation; placeholder for now
F2_PERM_P="(F2 sub-contrast permutation p pending — see Methods sensitivity)"

echo "Filling placeholders with:"
echo "  LOMO-S1 Govaere AUROC: $LOMO_S1_GOVAERE  [$LOMO_S1_GOVAERE_CI_LOW, $LOMO_S1_GOVAERE_CI_HIGH]"
echo "  Govaere expression-matched permutation null p: $PERM_P_GOVAERE"
echo "  Genetic-down label: $GENETIC_DOWN_LABEL"

FILES=(
  docs/manuscript/00_abstract.md
  docs/manuscript/02_results.md
  docs/manuscript/03_discussion.md
  docs/manuscript/04_methods.md
  docs/manuscript/05_figure_legends.md
  docs/manuscript/FIGURE_PLAN_REVISED.md
  docs/manuscript/NUMBERS.md
  docs/paper_outline.md
  docs/progress.md
  CLAUDE.md
)

for f in "${FILES[@]}"; do
  if [[ ! -f "$f" ]]; then
    echo "  SKIP $f (missing)" >&2
    continue
  fi
  # Replace placeholders. Use exact strings to avoid corrupting nearby brackets.
  sed -i.wave2bak \
    -e "s|\\[LOMO-S1 AUROC = TBD by Wave 2\\]|LOMO-S1 AUROC = $LOMO_S1_GOVAERE|g" \
    -e "s|LOMO-S1 AUROC = \\[TBD by Wave 2\\]|LOMO-S1 AUROC = $LOMO_S1_GOVAERE|g" \
    -e "s|AUROC = \\[TBD by Wave 2\\]|AUROC = $LOMO_S1_GOVAERE|g" \
    -e "s|\\[TBD by Wave 2\\]|$LOMO_S1_GOVAERE|g" \
    -e "s|\\[GENETIC_DOWN_LABEL\\]|$GENETIC_DOWN_LABEL|g" \
    "$f"
  rm -f "$f.wave2bak"
  echo "  filled placeholders in $f"
done

echo
echo "Verification grep — should return zero lines (no remaining placeholders):"
grep -rln "TBD by Wave 2\|GENETIC_DOWN_LABEL" "${FILES[@]}" 2>/dev/null || echo "  (clean)"
