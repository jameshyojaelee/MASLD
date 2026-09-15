#!/usr/bin/env bash
#SBATCH --job-name=model-check-629
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PLAN="${ROOT}/OVERALL_PLAN.md"
TARGET="${ROOT}/executions/model-check-629-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite LS-GKM plan validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
cd "${REPO}"
git diff --check -- Analysis/MASLD_Model_Benchmark/OVERALL_PLAN.md
cd "${ROOT}"
PYTHONPATH="${ROOT}/src" python - "${ROOT}" "${STAGE}" <<'PY'
from hashlib import sha256
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
stage = Path(sys.argv[2]).resolve(strict=True)
plan = root / "OVERALL_PLAN.md"
text = plan.read_text(encoding="utf-8")
required = (
    "model-cpu-train-614-21099919",
    "model-cpu-train-618-21099921",
    "model-data-623-21100461",
    "model-check-624-21100481",
    "model-check-625-21100532",
    "model-eval-626-21100570",
    "model-check-627-21100620",
    "The locked exposed-development evaluation did not promote either LS-GKM\nreadout.",
    "zero biological donors and four experimental replicates per",
    "not an external evaluation, shortlist,",
    "share the same 25 fits and remain one model family",
)
for token in required:
    if token not in text:
        raise SystemExit(f"required LS-GKM plan statement absent: {token}")
for stale in (
    "are running as five split bundles",
    "completion and outcome-blind aggregation of the active independent LS-GKM redesign",
    "That independent redesign is now active and remains outcome-blind",
):
    if stale in text:
        raise SystemExit(f"stale LS-GKM plan statement remains: {stale}")
authorities = {
    "executions/model-data-623-21100461/ARTIFACTS.json": "2d74167f7670bbcad7ec6e178f228b9f7efc80b308f827041d4df7df18357d1a",
    "executions/model-check-624-21100481/ARTIFACTS.json": "78c77e7d617593090bc4618c7c7276ec75475a9984357534edef1a9b88e4f87e",
    "executions/model-check-625-21100532/ARTIFACTS.json": "de04d91ea4b790e3de8f6b8bab2509d2b0ca3e9311f901f32982250f8e512d1f",
    "executions/model-eval-626-21100570/ARTIFACTS.json": "3a2db575d98c1fb22954c9f993ff3ba6390204437f5d285d0fda5ff5aa54276e",
    "executions/model-check-627-21100620/ARTIFACTS.json": "583cdca59e568adc05e0937c3853250cc232df72fa43f09286212c7270a40407",
}
for relative, expected in authorities.items():
    observed = sha256((root / relative).read_bytes()).hexdigest()
    if observed != expected or expected not in text:
        raise SystemExit(f"LS-GKM plan authority differs: {relative}")
receipt = {
    "schema_version": "masld-bench-lsgkm-overall-plan-update-validation-v1",
    "status": "pass_lsgkm_exposed_development_plan_update_validation",
    "plan_sha256": sha256(plan.read_bytes()).hexdigest(),
    "required_statements_checked": len(required),
    "stale_statements_rejected": 3,
    "artifact_authorities_checked": len(authorities),
    "development_only": True,
    "external_evaluation": False,
    "champion_claim": False,
}
(stage / "validation_receipt.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
freeze_tree(stage, {
    "artifact_class": "lsgkm_overall_plan_update_validation",
    **receipt,
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'LS-GKM plan validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
