#!/bin/bash
#SBATCH --job-name=count_r8_f10
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=1-10
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/count_r8_f10_%a_%A.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/count_r8_f10_%a_%A.err

set -euo pipefail
task_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
manifest="${task_root}/scripts/analysis/bulk_refit_followup/library_manifest_gse126848_f10.tsv"
export ANNOTATION_ROOT="${task_root}/figures/candidates/resource-corrections-20261005-v1/annotation_resolution_r8"
COUNT_SAMPLE_ROOT="${task_root}/figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_gse126848_f10_20261118/samples"
mkdir -p "${COUNT_SAMPLE_ROOT}" "${task_root}/logs"

line=$(awk -F'\t' -v n="${SLURM_ARRAY_TASK_ID}" 'NR==n+1{print $1"\t"$3"\t"$4}' "${manifest}")
IFS=$'\t' read -r SAMPLE_ID SOURCE_BAM layout <<<"${line}"
export SAMPLE_ID SOURCE_BAM LIBRARY_LAYOUT="${layout^^}" COUNT_OUTPUT_ROOT="${COUNT_SAMPLE_ROOT}/${SAMPLE_ID}"
if [[ -s "${COUNT_OUTPUT_ROOT}/count_completion.json" && -s "${COUNT_OUTPUT_ROOT}/gene_counts.tsv" ]]; then
  echo "already counted: ${SAMPLE_ID}"
  exit 0
fi

unset SLURM_TMPDIR
SCRATCH_BASE="/scratch"
if [[ ! -d "${SCRATCH_BASE}" || ! -w "${SCRATCH_BASE}" ]]; then
  SCRATCH_BASE="/tmp"
fi
export COUNT_SCRATCH_ROOT=$(mktemp -d "${SCRATCH_BASE}/count-f10-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}-XXXXXX")
trap 'rm -rf -- "${COUNT_SCRATCH_ROOT}"' EXIT
bash "${task_root}/scripts/analysis/alt_contig_recount/run_verified_count.sbatch"
