#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers"
RESULT_PARENT="${PROJECT_ROOT}/RNA-seq/results/histology_anchored_continuum/molecular_layers"

mkdir -p "${RESULT_PARENT}"
if [[ -z "${HAC_ML_OUT_ROOT:-}" ]]; then
  timestamp="$(date -u +%Y%m%dT%H%M%SZ)"
  export HAC_ML_OUT_ROOT="${RESULT_PARENT}/hac-molecular-layers-${timestamp}"
fi
if [[ -e "${HAC_ML_OUT_ROOT}" ]]; then
  echo "Refusing to overwrite existing candidate: ${HAC_ML_OUT_ROOT}" >&2
  exit 1
fi
mkdir -p "${HAC_ML_OUT_ROOT}"

echo "Candidate: ${HAC_ML_OUT_ROOT}"
echo "Concurrent workstreams: 52 CPUs, 480 GB; serial freeze/integration: 4 CPUs, 32 GB each"

freeze_job="$(sbatch --parsable --export=ALL,HAC_ML_OUT_ROOT "${SCRIPT_DIR}/run_preflight.sbatch")"
bulk_job="$(sbatch --parsable --dependency="afterok:${freeze_job}" --export=ALL,HAC_ML_OUT_ROOT "${SCRIPT_DIR}/run_bulk.sbatch")"
pathway_job="$(sbatch --parsable --dependency="afterok:${freeze_job}" --export=ALL,HAC_ML_OUT_ROOT "${SCRIPT_DIR}/run_pathway_tf.sbatch")"
program_job="$(sbatch --parsable --dependency="afterok:${freeze_job}" --export=ALL,HAC_ML_OUT_ROOT "${SCRIPT_DIR}/run_programs.sbatch")"
paired_job="$(sbatch --parsable --dependency="afterok:${freeze_job}" --export=ALL,HAC_ML_OUT_ROOT "${SCRIPT_DIR}/run_paired.sbatch")"
integration_job="$(sbatch --parsable --dependency="afterok:${bulk_job}:${pathway_job}:${program_job}:${paired_job}" --export=ALL,HAC_ML_OUT_ROOT "${SCRIPT_DIR}/run_integrate.sbatch")"

printf 'freeze\t%s\n' "${freeze_job}"
printf 'bulk\t%s\n' "${bulk_job}"
printf 'pathway_tf\t%s\n' "${pathway_job}"
printf 'programs\t%s\n' "${program_job}"
printf 'paired\t%s\n' "${paired_job}"
printf 'integrate\t%s\n' "${integration_job}"
printf 'candidate\t%s\n' "${HAC_ML_OUT_ROOT}"
