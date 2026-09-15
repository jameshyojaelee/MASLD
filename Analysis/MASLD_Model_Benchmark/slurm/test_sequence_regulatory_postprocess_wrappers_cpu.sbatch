#!/usr/bin/env bash
#SBATCH --job-name=model-check
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-check-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-check-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SEQUENCE=${ROOT}/slurm/evaluate_sequence_task_native_five_seed_rectangle_cpu.sbatch
SCBASSET=${ROOT}/slurm/evaluate_scbasset_five_seed_rectangle_cpu.sbatch
GATE=${ROOT}/slurm/run_sequence_regulatory_completion_and_launch_io.sbatch
TARGET=${ROOT}/executions/sequence-regulatory-postprocess-wrapper-tests-${SLURM_JOB_ID}
STAGE=${TARGET}.staging

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${TARGET}.failed" ]]; then
    mv "${STAGE}" "${TARGET}.failed"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite wrapper acceptance: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

cd "${ROOT}"
printf '%s  %s\n' \
  c134bd7796af7b28d75e4a9cfbe4f75a045413d14e86695a6309e538f6482611 "${SEQUENCE}" \
  41d1c755497d105b979f7cca7349a0af936ef11e1ae5488dd3127d73f186392c "${SCBASSET}" \
  5db37a347e4afca769f1e8e92f2c05be87993653f4ef9d7bf70b4b7c262a5710 "${GATE}" \
  | sha256sum --check --strict
bash -n "${SEQUENCE}" "${SCBASSET}" "${GATE}"

module load python/3.11.5-GCCcore-13.2.0
python - "${SEQUENCE}" "${SCBASSET}" "${GATE}" \
  > "${STAGE}/validation/contracts.stdout.txt" <<'PY'
from pathlib import Path
import sys

sequence, scbasset, gate = [Path(value).read_text(encoding="utf-8") for value in sys.argv[1:]]
for label, source in (("sequence", sequence), ("scbasset", scbasset)):
    required = (
        "#SBATCH --job-name=model-eval",
        "#SBATCH --partition=cpu",
        "#SBATCH --account=nslab",
        "#SBATCH --qos=nslab",
        ': "${COMPLETION_ROOT:?',
        ': "${COMPLETION_SHA:?',
        "ready_both_rectangles_locked",
        "partial_rectangle_used",
        "cross_task_family_ranking_performed",
        "sealed_data_read",
        "champion_claim_allowed",
    )
    if any(value not in source for value in required) or "#SBATCH --array" in source:
        raise SystemExit(f"{label} evaluator wrapper contract differs")
if (
    "#SBATCH --job-name=model-audit" not in gate
    or "#SBATCH --partition=io" not in gate
    or "#SBATCH --array" in gate
    or gate.index('[[ "${STATUS}" != ready_both_rectangles_locked ]]')
    > gate.index("SEQUENCE_JOB=$(sbatch")
    or gate.count("COMPLETION_ROOT=") < 2
    or gate.count("COMPLETION_SHA=") < 2
    or "completion_artifacts_sha256" not in gate
    or "cross_task_family_ranking_performed" not in gate
):
    raise SystemExit("production completion-and-launch wrapper contract differs")
print("postprocessing wrappers pass static exact-coverage contracts")
PY
python --version > "${STAGE}/validation/python-version.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum "${SEQUENCE}" "${SCBASSET}" "${GATE}" "$0" > "${STAGE}/source.sha256"

PYTHONPATH=${ROOT}/src:${ROOT} python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(freeze_tree(root, {
    "artifact_class": "sequence_regulatory_postprocess_wrapper_acceptance",
    "status": "pass",
    "wrappers": 3,
    "shell_syntax_passed": True,
    "exact_coverage_guard_present": True,
    "arrays_used": False,
    "outcomes_read": False,
    "prediction_values_read": False,
    "sealed_data_read": False,
    "partial_ranking_performed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Postprocess wrapper acceptance: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
