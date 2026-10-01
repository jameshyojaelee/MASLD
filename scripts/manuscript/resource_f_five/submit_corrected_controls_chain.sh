#!/usr/bin/env bash
# Runbook: corrected GSE130970 control chain (review item 1, Phase 1 steps 2-6).
#
#   1 fcorr_refit   refit_corrected_controls.R fit        -> candidates/bulk-corrected-controls-<TS>
#   2 fcorr_repro   refit_corrected_controls.R reproduce  -> candidates/bulk-corrected-controls-<TS>-reproduction
#   3 f0_arms       build_pi_stage_f0_arms.R              -> figures/candidates/f0-reference-arms-corrected-controls-<TS>
#   4 bulk_corr     stage extensions, NAS grid, Fig 3D, cascade
#                                                         -> figures/candidates/bulk-release-corrected-controls-<TS>
#   5 hac_bulk_corr continuum 10_bulk_transcript_deg.R    -> .../molecular_layers/hac-molecular-layers-corrected-controls-<TS>
#
# Each job starts only after the previous one succeeds (afterok). Job 4 needs the
# F0 arms from job 3 (the release sbatch checks f0_arm_results.tsv.gz on entry).
# Release packaging (build_bulk_release.py) is OFF: the retained Figure 3G and
# lncRNA inputs were built on the old fit and are rebuilt in Phase 3.
#
# Dry run by default. Submit only after approval:
#   bash scripts/manuscript/resource_f_five/submit_corrected_controls_chain.sh --submit
set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
TS="$(date -u +%Y%m%dT%H%M%SZ)"
REFIT_ROOT="${PROJECT_ROOT}/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-${TS}"
REPRO_ROOT="${REFIT_ROOT}-reproduction"
POOLED_ARM_ROOT="${REPRO_ROOT}/primary_source_controls"
F0_CANDIDATE="${PROJECT_ROOT}/figures/candidates/f0-reference-arms-corrected-controls-${TS}"
BULK_RELEASE_ID="corrected-controls-${TS}"
BULK_CANDIDATE="${PROJECT_ROOT}/figures/candidates/bulk-release-${BULK_RELEASE_ID}"
HAC_ROOT="${PROJECT_ROOT}/RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected-controls-${TS}"

for path in "${REFIT_ROOT}" "${REPRO_ROOT}" "${F0_CANDIDATE}" "${BULK_CANDIDATE}" "${HAC_ROOT}"; do
  if [[ -e "${path}" ]]; then echo "Refusing existing output: ${path}" >&2; exit 1; fi
done

cd "${PROJECT_ROOT}"
submit=0
[[ "${1:-}" == "--submit" ]] && submit=1
# Numeric job names only; JOB_INDEX_FILE holds the last index used.
next_job_name() {
  local f="${JOB_INDEX_FILE:?JOB_INDEX_FILE is required}" i
  i=$(( $(cat "$f") + 1 )); echo "$i" > "$f"; printf '%04d' "$i"
}
run() {
  if (( submit )); then sbatch --parsable "$@"; else echo "DRYRUN"; echo "sbatch --parsable $*" >&2; fi
}

j1=$(run --job-name="$(next_job_name)" \
  --export=ALL,REFIT_MODE=fit,REFIT_OUTPUT_ROOT="${REFIT_ROOT}" \
  scripts/manuscript/resource_f_five/run_refit_corrected_controls.sbatch)
j2=$(run --job-name="$(next_job_name)" --dependency="afterok:${j1}" \
  --export=ALL,REFIT_MODE=reproduce,REFIT_OUTPUT_ROOT="${REPRO_ROOT}",REFIT_SOURCE_ROOT="${REFIT_ROOT}" \
  scripts/manuscript/resource_f_five/run_refit_corrected_controls.sbatch)
j3=$(run --job-name="$(next_job_name)" --dependency="afterok:${j2}" \
  --export=ALL,STAGE_RELEASE_ROOT="${POOLED_ARM_ROOT}",FIGURE_CANDIDATE_ROOT="${F0_CANDIDATE}",RSCRIPT="${RNASEQ_ENV}/bin/Rscript" \
  scripts/figures/run_pi_stage_f0_arms.sbatch)
j4=$(run --job-name="$(next_job_name)" --dependency="afterok:${j3}" \
  --output="${PROJECT_ROOT}/logs/bulk_corr_%j.out" --error="${PROJECT_ROOT}/logs/bulk_corr_%j.err" \
  --export=ALL,BULK_RELEASE_ID="${BULK_RELEASE_ID}",POOLED_ROOT="${POOLED_ARM_ROOT}",F0_ROOT="${F0_CANDIDATE}/analysis/f0_arms",BULK_RELEASE_PACKAGE=0 \
  scripts/manuscript/bulk_release/run_bulk_release_20260909.sbatch)
j5=$(run --job-name="$(next_job_name)" --dependency="afterok:${j4}" \
  --export=ALL,HAC_ML_OUT_ROOT="${HAC_ROOT}",POOLED_ARM_ROOT="${POOLED_ARM_ROOT}",STAGE_EXTENSION_ROOT="${BULK_CANDIDATE}/analysis/stage_extensions" \
  scripts/analysis/histology_anchored_continuum/molecular_layers/run_bulk_corrected_controls.sbatch)

printf 'refit\t%s\t%s\n' "${j1}" "${REFIT_ROOT}"
printf 'reproduce\t%s\t%s\n' "${j2}" "${REPRO_ROOT}"
printf 'f0_arms\t%s\t%s\n' "${j3}" "${F0_CANDIDATE}"
printf 'bulk_figures\t%s\t%s\n' "${j4}" "${BULK_CANDIDATE}"
printf 'continuum\t%s\t%s\n' "${j5}" "${HAC_ROOT}"
(( submit )) || echo "Dry run only. Re-run with --submit after approval." >&2
