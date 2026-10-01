#!/usr/bin/env bash
#SBATCH --job-name=fig5-completion
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=scripts/figures/logs/fig5_completion_%j.out
#SBATCH --error=scripts/figures/logs/fig5_completion_%j.err

set -euo pipefail

# Required: CANDIDATE_ROOT (new directory, e.g. figures/candidates/fig5-completion-<date>-v7)
# and BULK_DEG_FILE (the bulk deg_results.csv that Figure 5C mRNA values come from).
# Panel 5C proteomics are rebuilt into $CANDIDATE_ROOT/proteomics from BULK_DEG_FILE.
: "${CANDIDATE_ROOT:?CANDIDATE_ROOT is required}"
: "${BULK_DEG_FILE:?BULK_DEG_FILE is required}"
PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
case "$CANDIDATE_ROOT" in /*) ;; *) CANDIDATE_ROOT="$PROJECT_ROOT/$CANDIDATE_ROOT" ;; esac
case "$BULK_DEG_FILE" in /*) ;; *) BULK_DEG_FILE="$PROJECT_ROOT/$BULK_DEG_FILE" ;; esac
test -s "$BULK_DEG_FILE"
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

PROTEOMICS_OUT_DIR="$CANDIDATE_ROOT/proteomics" BULK_DEG_FILE="$BULK_DEG_FILE" \
  "$R_BIN" Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R
export FIG5C_PROTEOMICS_DIR="$CANDIDATE_ROOT/proteomics"
sha256sum "$BULK_DEG_FILE" > "$CANDIDATE_ROOT/provenance/bulk_deg_file.sha256"

"$R_BIN" scripts/figures/fig5a_input_firewall.R
"$R_BIN" scripts/figures/composite_mrna_protein.R
"$R_BIN" scripts/figures/fig5d_variant_consistent_accessibility.R
if [[ -n "${FIG5B_BULK_DEG:-}" ]]; then
  "$R_BIN" scripts/figures/fig5b_protein_triage.R
else
  cp "$CANONICAL_ROOT/panels/fig5b_protein_triage.pdf" "$CANDIDATE_ROOT/figure5/panels/"
  cp "$CANONICAL_ROOT/panels/data/fig5b_protein_triage_summary.tsv" "$CANDIDATE_ROOT/figure5/panels/data/"
fi

for panel in \
  fig5e_multimodal_program_summary.pdf \
  fig5f_spatial_program_calibration.pdf
do
  cp "$CANONICAL_ROOT/panels/$panel" "$CANDIDATE_ROOT/figure5/panels/$panel"
done
for sidecar in \
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
