#!/usr/bin/env bash
#SBATCH --job-name=model-check-045
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
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
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/geneformer-v2-316m-compat-runtime-21077474
RUNTIME_SHA256=c034ee71fb2c67f04f11e91b41ea8c9bf6a7b9df3a4bb72378e5b2696058c307
RUNTIME_ARCHIVE_SHA256=a4fdc1f1d0b67c4753354c5a04ccc05c84a3bf7442d8655bc8ed6493a2487775
RUNTIME_VALIDATION=${ROOT}/executions/geneformer-v2-316m-preflight-validation-21077937
RUNTIME_VALIDATION_SHA256=7e936a9ac9c39e90bcb538f5b0562198f62ab20a4b4c3987c53d6c8d233cc063
FIXTURE=${ROOT}/executions/geneformer-v2-316m-frozen-screen-fixture-21078339
FIXTURE_SHA256=bb94045a46c7231d50a17ffd921a146ad6858b39bcfcda516da61230fcfb9bf6
BUNDLE=${ROOT}/executions/geneformer-acquisition-21014827/bundles/geneformer_v2_316m--11f311a5078fc4c8
BUNDLE_SHA256=b269cef53fa92cc7e7ab9c21a2b85de69904cdfaaf5881e9676c57dc2ae02ae0
TARGET=${ROOT}/executions/geneformer-v2-316m-frozen-screen-preflight-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
WORK=$(mktemp -d "${TMPDIR:-/tmp}/geneformer-50k-preflight.${SLURM_JOB_ID}.XXXXXX")

finish() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  if [[ -d "${WORK}" ]]; then
    rm -rf -- "${WORK}"
  fi
  exit "${exit_code}"
}
trap finish EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Geneformer 50k preflight.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${WORK}/runtime"
cd "${ROOT}"
printf '%s  %s\n' \
  "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  "${RUNTIME_ARCHIVE_SHA256}" "${RUNTIME}/payload/geneformer-decima-compat-env.tar.gz" \
  "${RUNTIME_VALIDATION_SHA256}" "${RUNTIME_VALIDATION}/ARTIFACTS.json" \
  "${FIXTURE_SHA256}" "${FIXTURE}/ARTIFACTS.json" \
  "${BUNDLE_SHA256}" "${BUNDLE}/ARTIFACTS.json" \
  > "${STAGE}/validation/input.sha256"
sha256sum --check --strict "${STAGE}/validation/input.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RUNTIME}" "${RUNTIME_VALIDATION}" "${FIXTURE}" "${BUNDLE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY
module purge
tar -xzf "${RUNTIME}/payload/geneformer-decima-compat-env.tar.gz" -C "${WORK}/runtime"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${WORK}/runtime/bin/python" -m unittest -v \
  tests.unit.test_geneformer_v2_316m_frozen_screen_50000 \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${WORK}/runtime/bin/python" - "${FIXTURE}/fixture" \
  "${BUNDLE}/upstream/Geneformer-V2-316M" "${STAGE}/cpu_preflight.json" <<'PY'
from pathlib import Path
import json
import sys

from scripts.extract_geneformer_v2_316m_frozen_screen_50000 import load_fixture
from scripts.preflight_geneformer_v2_316m_common import load_prior_compatible_model

fixture = load_fixture(Path(sys.argv[1]))
model, loading = load_prior_compatible_model(Path(sys.argv[2]))
parameters = sum(value.numel() for value in model.parameters())
receipt = {
    "schema_version": "masld-bench-geneformer-v2-316m-frozen-screen-cpu-preflight-v1",
    "status": "pass_cpu_preflight",
    "model_id": "geneformer_v2_316m",
    "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
    "rows": len(fixture["rows"]),
    "embedding_dimension": 1152,
    "model_parameter_count_including_unused_pooler": parameters,
    "prior_compatible_restore": loading,
    "one_batch_l40s_preflight_required": True,
    "one_batch_l40s_preflight_passed": False,
    "full_extraction_queued": False,
    "full_extraction_submitted": False,
    "evaluation_labels_read": False,
    "histology_read": False,
    "sealed_outcomes_read": False,
}
Path(sys.argv[3]).write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
PY
"${WORK}/runtime/bin/python" -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  scripts/preflight_geneformer_v2_316m_common.py \
  scripts/build_geneformer_v2_316m_frozen_screen_50000.py \
  scripts/extract_geneformer_v2_316m_frozen_screen_50000.py \
  tests/unit/test_geneformer_v2_316m_frozen_screen_50000.py \
  slurm/run_geneformer_v2_316m_frozen_screen_50000.sbatch \
  slurm/preflight_geneformer_v2_316m_frozen_screen_50000_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${RUNTIME_VALIDATION}/ARTIFACTS.json" \
  "${FIXTURE}/ARTIFACTS.json" \
  "${BUNDLE}/ARTIFACTS.json" \
  > "${STAGE}/validation/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((stage / "cpu_preflight.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_cpu_preflight"
    or receipt.get("rows") != 50_000
    or receipt.get("embedding_dimension") != 1152
    or receipt.get("one_batch_l40s_preflight_required") is not True
    or receipt.get("one_batch_l40s_preflight_passed") is not False
    or receipt.get("full_extraction_queued") is not False
    or receipt.get("full_extraction_submitted") is not False
    or receipt.get("evaluation_labels_read") is not False
    or receipt.get("histology_read") is not False
    or receipt.get("sealed_outcomes_read") is not False
):
    raise SystemExit("Geneformer 50k CPU preflight receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "geneformer_v2_316m_frozen_screen_cpu_preflight",
        "model_id": "geneformer_v2_316m",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "rows": 50_000,
        "one_batch_l40s_preflight_required": True,
        "one_batch_l40s_preflight_passed": False,
        "full_extraction_queued": False,
        "full_extraction_submitted": False,
        "evaluation_labels_read": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY

trap - EXIT
rm -rf -- "${WORK}"
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
