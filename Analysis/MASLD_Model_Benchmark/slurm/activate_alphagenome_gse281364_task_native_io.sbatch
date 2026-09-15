#!/usr/bin/env bash
#SBATCH --job-name=model-check-302
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export MPLBACKEND=Agg
export WANDB_MODE=offline
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome/bin/python
CONFIG="${ROOT}/config/alphagenome_gse281364_task_native_activation.json"
TARGET="${ROOT}/executions/model-check-302-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
export PYTHONPATH="${ROOT}:${ROOT}/src"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite task-native activation: %s\n' "${TARGET}" >&2
  exit 1
fi

install -d -m 0750 "${STAGE}/sources" "${STAGE}/validation"

curl --fail --silent --show-error --location --retry 4 --retry-delay 2 \
  --connect-timeout 30 --max-time 300 --proto '=https' --tlsv1.2 \
  https://api.github.com/repos/google-deepmind/alphagenome/git/ref/tags/v0.8.0 \
  --output "${STAGE}/sources/github_tag_v0.8.0.json"
curl --fail --silent --show-error --location --retry 4 --retry-delay 2 \
  --connect-timeout 30 --max-time 300 --proto '=https' --tlsv1.2 \
  https://pypi.org/pypi/alphagenome/0.8.0/json \
  --output "${STAGE}/sources/pypi_0.8.0.json"
curl --fail --silent --show-error --location --retry 4 --retry-delay 2 \
  --connect-timeout 30 --max-time 300 --proto '=https' --tlsv1.2 \
  https://huggingface.co/api/models/google/alphagenome-all-folds/revision/a8f293a76ee73d5b57f3bf2ae146510589fcf187 \
  --output "${STAGE}/sources/hf_checkpoint_revision.json"

curl --fail --silent --show-error --location --retry 4 --retry-delay 2 \
  --connect-timeout 30 --max-time 300 --proto '=https' --tlsv1.2 \
  --dump-header "${STAGE}/sources/api_terms.headers.txt" \
  https://deepmind.google.com/science/alphagenome/terms \
  --output "${STAGE}/sources/api_terms.html"
curl --fail --silent --show-error --location --retry 4 --retry-delay 2 \
  --connect-timeout 30 --max-time 300 --proto '=https' --tlsv1.2 \
  --dump-header "${STAGE}/sources/model_terms.headers.txt" \
  https://deepmind.google.com/science/alphagenome/model-terms \
  --output "${STAGE}/sources/model_terms.html"

curl --fail --silent --show-error --location --retry 4 --retry-delay 2 \
  --connect-timeout 30 --max-time 300 --proto '=https' --tlsv1.2 \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_report.txt \
  --output "${STAGE}/sources/GRCh38.p13_assembly_report.txt"
curl --fail --silent --show-error --location --retry 4 --retry-delay 2 \
  --connect-timeout 30 --max-time 300 --proto '=https' --tlsv1.2 \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_assembly_report.txt \
  --output "${STAGE}/sources/GRCh38.p14_assembly_report.txt"

"${PYTHON}" -m unittest -v \
  tests.unit.test_activate_alphagenome_gse281364_task_native \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/activate_alphagenome_gse281364_task_native.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --source-dir "${STAGE}/sources" \
  --output "${STAGE}/activation" \
  > "${STAGE}/validation/activation.stdout.json" \
  2> "${STAGE}/validation/activation.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
curl --version > "${STAGE}/validation/curl.version.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/activate_alphagenome_gse281364_task_native.py" \
  "${ROOT}/tests/unit/test_activate_alphagenome_gse281364_task_native.py" \
  "${ROOT}/slurm/activate_alphagenome_gse281364_task_native_io.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${TARGET}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace


stage = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2])
receipt = json.loads((stage / "activation/receipt.json").read_text(encoding="utf-8"))
if (
    not receipt["status"].startswith("pass_outcome_blind_task_native_activation_")
    or receipt["logical_model_version"] != "ALL_FOLDS"
    or receipt["request_manifest"]["requests"] != 1033
    or receipt["request_manifest"]["outer_locus_sequence_groups"] != 1033
    or receipt["request_manifest"]["input_length_bp"] != 524288
    or not receipt["reference_crosswalk"]["all_used_primary_sequence_accessions_equal"]
    or receipt["api_metadata_called"]
    or receipt["api_prediction_called"]
    or receipt["checkpoint_bytes_downloaded"]
    or receipt["checkpoint_deserialized"]
    or receipt["model_forward_executed"]
    or receipt["observed_outcomes_loaded"]
    or receipt["sealed_assets_loaded"]
    or receipt["metrics_calculated"]
    or receipt["open_champion_eligible"]
):
    raise SystemExit("AlphaGenome task-native activation receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "restricted_sequence_task_native_activation",
        "model_id": "alphagenome",
        "dataset_id": "gse281364",
        "logical_model_version": "ALL_FOLDS",
        "requests": 1033,
        "input_length_bp": 524288,
        "api_credential_present": receipt["credential_presence"]["api_credential_present"],
        "huggingface_credential_present": receipt["credential_presence"]["huggingface_credential_present"],
        "api_prediction_called": False,
        "checkpoint_bytes_downloaded": False,
        "model_forward_executed": False,
        "outcomes_read": False,
        "sealed_assets_read": False,
        "open_champion_eligible": False,
        "status": receipt["status"],
    },
)
publish_directory_noreplace(stage, target)
PY

trap - EXIT
printf 'Task-native activation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
