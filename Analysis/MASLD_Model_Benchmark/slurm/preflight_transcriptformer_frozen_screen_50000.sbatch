#!/usr/bin/env bash
#SBATCH --job-name=model-check-032
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
FIXTURE=${ROOT}/executions/model-data-051-21077382
FIXTURE_SHA256=ba40ce9555bb6d7bc278be2a80015b993f9f3650fa06f59505ad279c88cc117a
CHECKPOINT=${ROOT}/executions/transcriptformer-materialized-21066178/bundle/tf_sapiens
PYTHON=${RUNTIME}/env/bin/python
TARGET=${ROOT}/executions/model-check-032-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite model preflight artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FIXTURE_SHA256}" "${FIXTURE}/ARTIFACTS.json" \
  | sha256sum --check --strict

cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
"${PYTHON}" -m unittest \
  tests.unit.test_build_transcriptformer_unlabeled_fixture \
  tests.unit.test_transcriptformer_extract_unlabeled_embeddings \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" - "${FIXTURE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

value = verify_frozen_tree(Path(sys.argv[1]))
if value["metadata"].get("development_rows") != 50_000:
    raise SystemExit("50k input fixture differs")
PY
"${PYTHON}" scripts/preflight_transcriptformer_cpu.py \
  --checkpoint-root "${CHECKPOINT}" \
  --fixture-root "${FIXTURE}/fixture" \
  --output "${STAGE}/cpu_preflight.json" \
  --expected-rows 50000 \
  > "${STAGE}/validation/preflight_stdout.json" \
  2> "${STAGE}/validation/preflight_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/preflight_transcriptformer_cpu.py \
  scripts/transcriptformer_extract_unlabeled_embeddings.py \
  tests/unit/test_transcriptformer_extract_unlabeled_embeddings.py \
  slurm/preflight_transcriptformer_frozen_screen_50000.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${FIXTURE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "cpu_preflight.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_cpu_preflight"
    or receipt.get("development_rows") != 50_000
    or receipt.get("model_parameter_count") != 429_870_637
    or receipt.get("checkpoint_strict_load") is not True
    or receipt.get("evaluation_label_columns_read")
    or receipt.get("sealed_outcomes_read")
    or receipt.get("downstream_head_fit")
):
    raise SystemExit("50k CPU preflight receipt differs")
freeze_tree(root, {
    "artifact_class": "cell_foundation_frozen_screen_cpu_preflight",
    "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
    "development_rows": 50_000,
    "checkpoint_strict_load": True,
    "evaluation_labels_read": False,
    "model_id": "transcriptformer_tf_sapiens",
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
