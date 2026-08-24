#!/usr/bin/env bash
#SBATCH --job-name=masld-stabmap-stability
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

: "${STABMAP_EVAL4817:?STABMAP_EVAL4817 is required}"
: "${STABMAP_EVAL5839:?STABMAP_EVAL5839 is required}"
: "${STABMAP_EVAL6857:?STABMAP_EVAL6857 is required}"
: "${STABMAP_EVAL4817_SHA:?STABMAP_EVAL4817_SHA is required}"
: "${STABMAP_EVAL5839_SHA:?STABMAP_EVAL5839_SHA is required}"
: "${STABMAP_EVAL6857_SHA:?STABMAP_EVAL6857_SHA is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/stabmap-retrieval-seed-stability-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed StabMap stability summary: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite StabMap stability summary: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "${STABMAP_EVAL4817_SHA} ${STABMAP_EVAL4817}/ARTIFACTS.json" \
  "${STABMAP_EVAL5839_SHA} ${STABMAP_EVAL5839}/ARTIFACTS.json" \
  "${STABMAP_EVAL6857_SHA} ${STABMAP_EVAL6857}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${STABMAP_EVAL4817}" "${STABMAP_EVAL5839}" "${STABMAP_EVAL6857}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_summarize_stabmap_seed_stability \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}/scripts" python \
  "${ROOT}/scripts/summarize_stabmap_seed_stability.py" \
  --seed-evaluation "4817=${STABMAP_EVAL4817}" \
  --seed-evaluation "5839=${STABMAP_EVAL5839}" \
  --seed-evaluation "6857=${STABMAP_EVAL6857}" \
  --output "${STAGE}/summary" \
  > "${STAGE}/validation/summary.stdout.json" \
  2> "${STAGE}/validation/summary.stderr.txt"
sha256sum \
  "${ROOT}/scripts/summarize_stabmap_seed_stability.py" \
  "${ROOT}/tests/unit/test_summarize_stabmap_seed_stability.py" \
  "${ROOT}/slurm/stabmap_summarize_seed_stability.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - \
  "${STAGE}" "${STABMAP_EVAL4817_SHA}" \
  "${STABMAP_EVAL5839_SHA}" "${STABMAP_EVAL6857_SHA}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
evaluation_sha = {
    "4817": sys.argv[2],
    "5839": sys.argv[3],
    "6857": sys.argv[4],
}
summary = json.loads((root / "summary/stability.json").read_text())
if (
    summary["status"] != "pass"
    or summary["model_id"] != "stabmap"
    or summary["seed_count"] != 3
    or set(summary["seeds"]) != {"4817", "5839", "6857"}
    or summary["hidden_pair_map_available_to_models"]
    or summary["cells_used_as_biological_replicates"]
    or summary["test_outcomes_read"]
    or summary["rna_conditioned_atac_prediction_claim"]
    or summary["sealed_inference_eligible"]
    or summary["champion_claim_allowed"]
):
    raise SystemExit("StabMap seed stability contract differs")
freeze_tree(root, {
    "artifact_class": "stabmap_retrieval_seed_stability",
    "model_id": "stabmap",
    "dataset_id": "gse296875",
    "evaluation_artifacts_sha256": evaluation_sha,
    "seeds_are_biological_replicates": False,
    "hidden_pair_map_available_to_models": False,
    "test_outcomes_read": False,
    "sealed_inference_eligible": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'StabMap seed stability: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
