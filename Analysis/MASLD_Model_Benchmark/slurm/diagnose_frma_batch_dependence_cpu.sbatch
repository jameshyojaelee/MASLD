#!/usr/bin/env bash
#SBATCH --job-name=model-check-087
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Diagnostic for the model-check-085 separation failure. Not a published output file
# and not an input to anything: it locates which fRMA stage carries the batch
# dependence so the requirements question can be put to the owner with evidence.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null
export R_LIBS_USER=/nonexistent/masld_bench_r_user_library

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/microarray-r-bioc-21104442
RUNTIME_SHA256=8115e0966804dfc4b7bfac3f95016896d099d82fc51f3428fab1ec34182a7d07
RAW_AUDIT=${ROOT}/executions/model-data-065-21079931
RAW_AUDIT_SHA256=8b30a2265c421f2917ed886d99cfaa08ee8d874abd40e68c8c72d48c726781d1
CEL_TAR=${RAW_AUDIT}/raw/GSE49541_RAW.tar
CEL_MANIFEST=${RAW_AUDIT}/GSE49541_cel_manifest.tsv
CEL_MANIFEST_SHA256=46204986a7bbefbcde169411d7fb5b6bded4b6babd560e49fc3e47e065c579f8
PLATFORM_AXIS=${ROOT}/executions/model-data-067-21080193-platform/GPL570_selected_annotation.tsv
PLATFORM_AXIS_SHA256=4bbd91765d1a80b12439e8c39b2aca9218a4135a80fc6a50f6f4eed276d61b69

TARGET=${ROOT}/executions/model-check-087-${SLURM_JOB_ID}-diagnostic
STAGE=${TARGET}.staging
CEL_WORK=${TMPDIR:-/tmp}/masld_gse49541_${SLURM_JOB_ID}

cleanup() { rm -rf "${CEL_WORK}"; }
trap cleanup EXIT
for path in "${TARGET}" "${STAGE}"; do
  if [[ -e "${path}" ]]; then
    printf 'Refusing to overwrite diagnostic: %s\n' "${path}" >&2
    exit 1
  fi
done
install -d -m 0750 "${STAGE}"
install -d -m 0700 "${CEL_WORK}"

printf '%s  %s\n' \
  "${RUNTIME_SHA256}"      "${RUNTIME}/ARTIFACTS.json" \
  "${RAW_AUDIT_SHA256}"    "${RAW_AUDIT}/ARTIFACTS.json" \
  "${CEL_MANIFEST_SHA256}" "${CEL_MANIFEST}" \
  "${PLATFORM_AXIS_SHA256}" "${PLATFORM_AXIS}" \
  | sha256sum --check --strict

cd "${ROOT}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
python3 scripts/stage_gse49541_cel_arrays.py \
  --tar "${CEL_TAR}" --manifest "${CEL_MANIFEST}" \
  --output "${CEL_WORK}" --receipt "${STAGE}/cel_staging_receipt.json" > /dev/null

module purge
module load R/4.5.2-mba
R_LIBS="${RUNTIME}/library" Rscript --vanilla \
  "${ROOT}/scripts/diagnose_frma_batch_dependence.R" \
  "${CEL_WORK}" "${CEL_MANIFEST}" "${PLATFORM_AXIS}" "${STAGE}" \
  > "${STAGE}/diagnostic_stdout.txt" 2> "${STAGE}/diagnostic_stderr.txt"

rm -rf "${CEL_WORK}"
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'fRMA batch-dependence diagnostic: %s\n' "${TARGET}"
