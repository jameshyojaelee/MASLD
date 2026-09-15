#!/usr/bin/env bash
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4

: "${ENDPOINT:?}" ; : "${VARIANCE_FRACTION:?}" ; : "${CRIT_NOMINAL:?}" ; : "${CRIT_BH:?}"
: "${REPLICATES:=200}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON="${PROJECT}/.agent/scanpy_env/bin/python"
PSEUDOBULK="${ROOT}/executions/gse296875-donor-lineage-pseudobulk-21107249"
PSEUDOBULK_SHA256=e381f60328645d8e0f2887976021094be35f3d5f8d703a78efe9f5ea7fadb1b7
TARGET="${ROOT}/executions/gse296875-realised-power-${ENDPOINT}-vf${VARIANCE_FRACTION}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite realised-power artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
mkdir -p "${STAGE}/validation"

cat > "${STAGE}/validation/bound_fixtures.sha256" <<SUMS
6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140  executions/gse296875-phenotype-endpoints-20260825/ARTIFACTS.json
e6bc9161acaab508c298c49036af1ad10a720949178a14cba42282933195d416  executions/gse296875-donor-folds-20260825/ARTIFACTS.json
474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58  executions/gse296875-phenotype-campaign-spec-21107179/ARTIFACTS.json
${PSEUDOBULK_SHA256}  ${PSEUDOBULK}/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/validation/bound_fixtures.sha256"

"${PYTHON}" scripts/simulate_gse296875_realised_power.py \
  --root "${ROOT}" \
  --pseudobulk "${PSEUDOBULK}" \
  --pseudobulk-sha256 "${PSEUDOBULK_SHA256}" \
  --endpoint "${ENDPOINT}" \
  --variance-fraction "${VARIANCE_FRACTION}" \
  --replicates "${REPLICATES}" \
  --critical-nominal "${CRIT_NOMINAL}" \
  --critical-bh "${CRIT_BH}" \
  --output "${STAGE}/power" \
  > "${STAGE}/validation/simulation_stdout.jsonl"

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true
sha256sum scripts/simulate_gse296875_realised_power.py \
          slurm/simulate_gse296875_realised_power_cpu.sbatch > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "power" / "audit.json").read_text())
if audit["observed_model_performance_read"] or audit["score_artifact_opened"]:
    raise SystemExit("a power simulation may not read an observed result")
if audit["effect_is_planted_not_observed"] is not True:
    raise SystemExit("the effect must be planted, not observed")
freeze_tree(root, {
    "artifact_class": "gse296875_realised_power_simulation",
    "dataset_id": "gse296875",
    "endpoint": audit["endpoint"],
    "planted_variance_fraction": audit["planted_variance_fraction"],
    "effect_is_planted_not_observed": True,
    "observed_model_performance_read": False,
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Realised power: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
