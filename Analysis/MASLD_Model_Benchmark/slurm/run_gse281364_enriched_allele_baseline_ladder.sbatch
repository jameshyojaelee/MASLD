#!/usr/bin/env bash
#SBATCH --job-name=mpra-enriched-baseline
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
# scikit-learn pinned to 1.9.0, the version recorded in the frozen campaign's
# validation/pip-freeze.txt; numpy 2.2.6 and scipy 1.15.3 come from the campaign
# environment underneath and already match.
PYTHON=${ROOT}/executions/environments/mpra-enriched-baseline-sklearn190-20260905T192000Z/env/bin/python
PRESPEC_DIR=${ROOT}/executions/mpra-enriched-baseline-prespec-20260905T192000Z
PRESPEC=${PRESPEC_DIR}/PRESPEC.md
PRESPEC_SHA256=b670655d9fd3bafec3d11c3129c25f71f677a7805c10b3a37fb1e097bca4139c
TARGET=${ROOT}/executions/mpra-enriched-baseline-20260905T192000Z-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite enriched baseline lane: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH=${ROOT}:${ROOT}/src
cd "${ROOT}"

SOURCE_FILES=(
  scripts/fit_gse281364_enriched_allele_baseline_ladder.py
  scripts/fit_gse281364_enformer_sei_heads.py
  scripts/evaluate_gse281364_dna_language_seeded_heads.py
  scripts/evaluate_gse281364_dna_lm_common_lane.py
  slurm/run_gse281364_enriched_allele_baseline_ladder.sbatch
)
sha256sum "${SOURCE_FILES[@]}" > "${STAGE}/source.start.sha256"
sha256sum "${PRESPEC}" > "${STAGE}/prespec.sha256"

"${PYTHON}" "${ROOT}/scripts/fit_gse281364_enriched_allele_baseline_ladder.py" \
  --root "${ROOT}" \
  --prespec "${PRESPEC}" \
  --prespec-sha256 "${PRESPEC_SHA256}" \
  --output "${STAGE}/ladder" \
  > "${STAGE}/validation/ladder.stdout.json" \
  2> "${STAGE}/validation/ladder.stderr.txt"

cp "${PRESPEC}" "${STAGE}/PRESPEC.md"
"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
"${PYTHON}" - <<'PY' > "${STAGE}/validation/seeds.json"
import json
print(json.dumps({
    "fixed_seeds": [1103, 2909, 4721, 6673, 8111],
    "inner_bootstrap_salt_offset": 10000,
    "outer_bootstrap_salt_offset": 20000,
    "training_mean_bootstrap_salt_offset": 30000,
    "bootstrap_seed": 20260825,
    "bootstrap_resamples": 10000,
    "bootstrap_unit": "long_range_block_id",
}, indent=2, sort_keys=True))
PY
sha256sum "${SOURCE_FILES[@]}" > "${STAGE}/source.terminal.sha256"
cmp --silent "${STAGE}/source.start.sha256" "${STAGE}/source.terminal.sha256"
sha256sum --check --strict "${STAGE}/source.terminal.sha256" \
  > "${STAGE}/validation/terminal_source_recheck.stdout.txt"

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Completed enriched allele baseline ladder: %s\n' "${TARGET}"
