#!/usr/bin/env bash
#SBATCH --job-name=model-check-224
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
CONFIG=${ROOT}/config/gse281364_open_sequence_head_campaign.toml
PRIOR_VALIDATION=${ROOT}/executions/model-check-223-21090508
PRIOR_VALIDATION_SHA=edaa354661b2d8c564d040472ac7f4e9d87099f04292842c00d2bd31b262bcf3
TARGET=${ROOT}/executions/model-check-224-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE281364 v2 binding validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH=${ROOT}:${ROOT}/src

printf '%s  %s\n' "${PRIOR_VALIDATION_SHA}" "${PRIOR_VALIDATION}/ARTIFACTS.json" \
  | sha256sum --check --strict
"${PYTHON}" - "${PRIOR_VALIDATION}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]).resolve(strict=True))
PY

bash -n slurm/run_gse281364_open_sequence_head_campaign_v2_cpu.sbatch
bash -n slurm/validate_gse281364_open_sequence_head_campaign_v2_cpu.sbatch
"${PYTHON}" -m py_compile \
  scripts/freeze_gse281364_open_sequence_taskspec.py \
  scripts/fit_gse281364_open_sequence_heads.py \
  scripts/evaluate_gse281364_open_sequence_heads.py \
  scripts/verify_gse281364_open_sequence_production_binding.py \
  tests/unit/test_gse281364_open_sequence_head_campaign.py \
  tests/unit/test_verify_gse281364_open_sequence_production_binding.py

"${PYTHON}" -m unittest -v \
  tests.unit.test_gse281364_open_sequence_head_campaign \
  tests.unit.test_verify_gse281364_open_sequence_production_binding \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" scripts/freeze_gse281364_open_sequence_taskspec.py \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/task_spec" \
  > "${STAGE}/validation/taskspec.stdout.json" \
  2> "${STAGE}/validation/taskspec.stderr.txt"

"${PYTHON}" - "${STAGE}/task_spec" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(root, {
    "artifact_class": "gse281364_open_sequence_taskspec",
    "dataset_id": "gse281364",
    "task_id": "variant_to_regulation",
    "candidate_count": 12,
    "biological_donors": 0,
    "experimental_replicates_per_context": 4,
    "outcomes_read": False,
    "model_fit": False,
    "status": "pass_prespecified_taskspec",
})
verify_frozen_tree(root)
PY

sha256sum \
  config/gse281364_open_sequence_head_campaign.toml \
  scripts/freeze_gse281364_open_sequence_taskspec.py \
  scripts/fit_gse281364_open_sequence_heads.py \
  scripts/evaluate_gse281364_open_sequence_heads.py \
  scripts/verify_gse281364_open_sequence_production_binding.py \
  tests/unit/test_gse281364_open_sequence_head_campaign.py \
  tests/unit/test_verify_gse281364_open_sequence_production_binding.py \
  slurm/run_gse281364_open_sequence_head_campaign_v2_cpu.sbatch \
  slurm/validate_gse281364_open_sequence_head_campaign_v2_cpu.sbatch \
  > "${STAGE}/source.sha256"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
"${PYTHON}" - "${STAGE}" "${PRIOR_VALIDATION}" "${PRIOR_VALIDATION_SHA}" <<'PY'
from hashlib import sha256
import json
from pathlib import Path
import sys

from scripts.verify_gse281364_open_sequence_production_binding import (
    CONFIG,
    MISSING_BASELINES,
    PRODUCTION_WRAPPER,
)


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


stage = Path(sys.argv[1]).resolve(strict=True)
root = Path.cwd().resolve(strict=True)
receipt = {
    "schema_version": "masld-bench-gse281364-open-sequence-validation-v2",
    "status": "passed_bound_production_wrapper",
    "campaign_scope": "partial_common_head_campaign",
    "candidate_count": 12,
    "unit_tests": 11,
    "mandatory_baselines_complete": False,
    "mandatory_baselines_missing": list(MISSING_BASELINES),
    "shortlist_blocked": True,
    "finalist_claim_blocked": True,
    "complementarity_blocked": True,
    "conditional_trigger_blocked": True,
    "outcomes_read": False,
    "model_fit": False,
    "predictions_generated": False,
    "production_wrapper_submitted": False,
    "production_wrapper_sha256": digest(root / PRODUCTION_WRAPPER),
    "config_sha256": digest(root / CONFIG),
    "source_manifest_sha256": digest(stage / "source.sha256"),
    "validated_source_count": 9,
    "prior_validation_root": Path(sys.argv[2]).resolve(strict=True).as_posix(),
    "prior_validation_artifacts_sha256": sys.argv[3],
}
(stage / "validation/receipt.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(receipt, sort_keys=True))
PY

"${PYTHON}" - "${STAGE}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads(
    (stage / "validation/receipt.json").read_text(encoding="utf-8")
)
if (
    receipt["status"] != "passed_bound_production_wrapper"
    or receipt["unit_tests"] != 11
    or receipt["validated_source_count"] != 9
    or receipt["mandatory_baselines_complete"]
    or len(receipt["mandatory_baselines_missing"]) != 5
    or not receipt["shortlist_blocked"]
    or not receipt["finalist_claim_blocked"]
    or not receipt["complementarity_blocked"]
    or not receipt["conditional_trigger_blocked"]
    or receipt["outcomes_read"]
    or receipt["model_fit"]
    or receipt["predictions_generated"]
    or receipt["production_wrapper_submitted"]
):
    raise SystemExit("GSE281364 v2 binding-validation receipt differs")
freeze_tree(stage, {
    "artifact_class": "gse281364_open_sequence_head_campaign_validation_v2",
    "status": "passed_bound_production_wrapper",
    "campaign_scope": "partial_common_head_campaign",
    "candidate_count": 12,
    "unit_tests": 11,
    "mandatory_baselines_complete": False,
    "mandatory_baselines_missing": receipt["mandatory_baselines_missing"],
    "shortlist_blocked": True,
    "finalist_claim_blocked": True,
    "complementarity_blocked": True,
    "conditional_trigger_blocked": True,
    "outcomes_read": False,
    "model_fit": False,
    "predictions_generated": False,
    "production_wrapper_submitted": False,
    "production_wrapper_sha256": receipt["production_wrapper_sha256"],
    "config_sha256": receipt["config_sha256"],
    "source_manifest_sha256": receipt["source_manifest_sha256"],
    "prior_validation_root": receipt["prior_validation_root"],
    "prior_validation_artifacts_sha256": receipt["prior_validation_artifacts_sha256"],
})
verify_frozen_tree(stage)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 open-sequence binding validation v2: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
