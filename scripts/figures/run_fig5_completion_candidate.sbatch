#!/usr/bin/env bash
#SBATCH --job-name=fig5-completion
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=scripts/figures/logs/fig5_completion_%j.out
#SBATCH --error=scripts/figures/logs/fig5_completion_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
CANDIDATE_ROOT="$PROJECT_ROOT/figures/candidates/fig5-completion-2026-08-24-v6"
CANONICAL_ROOT="$PROJECT_ROOT/figures/main/fig5_molecular_context"
R_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

test ! -e "$CANDIDATE_ROOT"
mkdir -p "$CANDIDATE_ROOT/figure5/panels/data" "$CANDIDATE_ROOT/figure5/data" "$CANDIDATE_ROOT/provenance"

export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
export FIGURE_CANDIDATE_ROOT="$CANDIDATE_ROOT"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "$PROJECT_ROOT"

"$R_BIN" scripts/figures/fig5a_input_firewall.R
"$R_BIN" scripts/figures/composite_mrna_protein.R
"$R_BIN" scripts/figures/fig5d_variant_consistent_accessibility.R

for panel in \
  fig5b_protein_triage.pdf \
  fig5e_multimodal_program_summary.pdf \
  fig5f_spatial_program_calibration.pdf
do
  cp "$CANONICAL_ROOT/panels/$panel" "$CANDIDATE_ROOT/figure5/panels/$panel"
done
for sidecar in \
  fig5b_protein_triage_summary.tsv \
  fig5e_multimodal_program_summary_selection.tsv \
  fig5f_spatial_program_calibration.tsv \
  fig5f_spatial_histology_overlay.tsv.gz \
  fig5f_spatial_histology_overlay_metadata.tsv
do
  cp "$CANONICAL_ROOT/panels/data/$sidecar" "$CANDIDATE_ROOT/figure5/panels/data/$sidecar"
done

"$PY_BIN" scripts/figures/assemble_fig5.py --candidate-root "$CANDIDATE_ROOT"
"$PY_BIN" scripts/figures/validate_fig5_completion.py --candidate-root "$CANDIDATE_ROOT"
FIG5_SESSION_INFO="$CANDIDATE_ROOT/provenance/sessionInfo.txt" \
  "$R_BIN" -e 'writeLines(capture.output(sessionInfo()), Sys.getenv("FIG5_SESSION_INFO"))'
"$PY_BIN" -m pip freeze > "$CANDIDATE_ROOT/provenance/pip_freeze.txt"
python3 scripts/manuscript/validate_resource_scope.py > "$CANDIDATE_ROOT/provenance/resource_scope_validation.txt"
touch "$CANDIDATE_ROOT/VALIDATED"

printf 'VALIDATED\t%s\n' "$CANDIDATE_ROOT"
