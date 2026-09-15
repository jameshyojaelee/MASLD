#!/usr/bin/env bash
#SBATCH --job-name=model-check-092
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Confirm the repaired matrix writer on the real GSE49541 values preserved in
# the failed stage of job 21109357, rather than re-running 47 minutes blind.
# Diagnostic; publishes no matrix and is not an input to anything.

set -euo pipefail
umask 027
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null
export R_LIBS_USER=/nonexistent/masld_bench_r_user_library

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/microarray-r-bioc-21104442
RUNTIME_SHA256=8115e0966804dfc4b7bfac3f95016896d099d82fc51f3428fab1ec34182a7d07
LOSSY=${ROOT}/executions/model-data-086-21109357-summarization.failed/summary/gse49541_gpl570_platform_feature_matrix.tsv

TARGET=${ROOT}/executions/model-check-092-${SLURM_JOB_ID}-roundtrip
STAGE=${TARGET}.staging
WORK=${TMPDIR:-/tmp}/masld_roundtrip_${SLURM_JOB_ID}

cleanup() { rm -rf "${WORK}"; }
trap cleanup EXIT
for path in "${TARGET}" "${STAGE}"; do
  if [[ -e "${path}" ]]; then
    printf 'Refusing to overwrite: %s\n' "${path}" >&2
    exit 1
  fi
done
install -d -m 0750 "${STAGE}"
install -d -m 0700 "${WORK}"

printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
if [[ ! -f "${LOSSY}" ]]; then
  printf 'Preserved failed-stage matrix is absent: %s\n' "${LOSSY}" >&2
  exit 1
fi

cd "${ROOT}"
module purge
module load R/4.5.2-mba
R_LIBS="${RUNTIME}/library" Rscript --vanilla \
  "${ROOT}/scripts/verify_matrix_text_roundtrip.R" "${LOSSY}" "${WORK}" \
  > "${STAGE}/roundtrip_stdout.txt" 2> "${STAGE}/roundtrip_stderr.txt"
cp "${WORK}/matrix_text_roundtrip.json" "${STAGE}/matrix_text_roundtrip.json"

rm -rf "${WORK}"
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Matrix text round-trip check: %s\n' "${TARGET}"
cat "${TARGET}/roundtrip_stdout.txt"
