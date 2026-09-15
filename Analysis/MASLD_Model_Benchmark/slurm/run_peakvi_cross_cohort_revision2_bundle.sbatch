#!/usr/bin/env bash
#SBATCH --job-name=model-training-333
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE="${ROOT}/slurm/run_peakvi_cross_cohort_bundle.sbatch"
SOURCE_SHA256=b55f8a56da7e3fa76849cdbfb3bb3cc6ac4f4eae6f48dfb627135e2c14b19e4a
DERIVED_SHA256=b94976d2630bde9326389825529898d40003cdb6affb0a82a48b9a0179948002
OLD_CONFIG=peakvi_gse296875_to_gse281367_fit_predict_20260825.json
NEW_CONFIG=peakvi_gse296875_to_gse281367_fit_predict_revision2_20260825.json
DERIVED=$(mktemp "${TMPDIR:-/tmp}/peakvi-cross-cohort-revision2.XXXXXX")

cleanup() {
  rm -f -- "${DERIVED}"
}
trap cleanup EXIT
printf '%s  %s\n' "${SOURCE_SHA256}" "${SOURCE}" | sha256sum --check --strict
if [[ $(grep -oF "${OLD_CONFIG}" "${SOURCE}" | wc -l) -ne 1 ]]; then
  printf 'PeakVI source-wrapper config occurrence differs\n' >&2
  exit 2
fi
sed "s/${OLD_CONFIG}/${NEW_CONFIG}/g" "${SOURCE}" > "${DERIVED}"
printf '%s  %s\n' "${DERIVED_SHA256}" "${DERIVED}" | sha256sum --check --strict
bash "${DERIVED}"
