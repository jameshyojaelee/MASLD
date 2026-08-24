#!/usr/bin/env bash
#SBATCH --job-name=masld-ablation-specs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONTEXT="${ROOT}/executions/context-conditioning-precondition-21064508"
CONTEXT_ARTIFACTS_SHA256=3577ad32ac4b34252b78a4a8cd489375dcfaba3d50f3093f870eaedceabaccb5
TARGET="${ROOT}/executions/conditional-ablation-preconditions-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite conditional ablation preconditions: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${CONTEXT_ARTIFACTS_SHA256}" \
  "${CONTEXT}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="src:." python scripts/build_conditional_ablation_preconditions.py \
  --campaign config/campaigns/v1_conditional_model.toml \
  --context-fixture "${CONTEXT}" \
  --output "${STAGE}/conditional_ablation_preconditions.json" \
  > "${STAGE}/validation/build.stdout.txt" \
  2> "${STAGE}/validation/build.stderr.txt"
PYTHONPATH="src:." python - "${STAGE}/conditional_ablation_preconditions.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.hashing import canonical_sha256

payload = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
identity = dict(payload)
claimed = identity.pop("precondition_spec_id")
if canonical_sha256(identity) != claimed:
    raise SystemExit("conditional ablation identity does not rederive")
expected = [
    "sequence_only",
    "trans_only",
    "permuted_context",
    "reverse_complement",
    "matched_shuffle",
]
if payload["ablation_names"] != expected or set(payload["ablations"]) != set(expected):
    raise SystemExit("conditional ablation roster differs")
for name in expected:
    spec = payload["ablations"][name]
    if (
        spec["executable"]
        or spec["model_training_allowed"]
        or spec["sealed_outcomes_allowed"]
        or spec["independent_champion_identity"]
        or not spec["future_realization_requires_new_frozen_parent_bound_spec"]
    ):
        raise SystemExit(f"conditional ablation firewall differs: {name}")
if (
    payload["architecture_built"]
    or payload["model_training_started"]
    or payload["conditional_model_authorized"]
    or payload["sealed_outcomes_exposed"]
    or payload["conditional_campaign"]["status"] != "blocked_trigger"
    or payload["conditional_campaign"]["submit_enabled"]
):
    raise SystemExit("conditional campaign firewall differs")
print(f"conditional ablation preconditions verified: {claimed}")
PY
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/build_conditional_ablation_preconditions.py" \
  "${ROOT}/src/masld_bench/context_preconditions.py" \
  "${ROOT}/config/campaigns/v1_conditional_model.toml" \
  "${ROOT}/slurm/build_conditional_ablation_preconditions.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1])
payload = json.loads(
    (root / "conditional_ablation_preconditions.json").read_text(encoding="utf-8")
)
freeze_tree(
    root,
    metadata={
        "artifact_class": "conditional_ablation_preconditions",
        "status": "precondition_only_parent_unresolved",
        "precondition_spec_id": payload["precondition_spec_id"],
        "ablation_count": 5,
        "architecture_built": False,
        "model_training_started": False,
        "conditional_model_authorized": False,
        "sealed_outcomes_exposed": False,
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Conditional ablation preconditions: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
