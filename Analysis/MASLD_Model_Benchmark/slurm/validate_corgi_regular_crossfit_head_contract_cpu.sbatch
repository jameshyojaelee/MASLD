#!/usr/bin/env bash
#SBATCH --job-name=model-check-256
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONTRACT="${ROOT}/config/evaluation/corgi_regular_crossfit_profile_head_smoke.toml"
TARGET="${ROOT}/executions/model-check-256-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi head-contract preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_freeze_corgi_regular_crossfit_head_contract \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="src:." python scripts/freeze_corgi_regular_crossfit_head_contract.py \
  --project-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/contract" \
  > "${STAGE}/validation/contract.stdout.json"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  config/evaluation/corgi_regular_crossfit_profile_head_smoke.toml \
  scripts/freeze_corgi_regular_crossfit_head_contract.py \
  tests/unit/test_freeze_corgi_regular_crossfit_head_contract.py \
  slurm/validate_corgi_regular_crossfit_head_contract_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "contract/contract_receipt.json").read_text(encoding="utf-8"))
freeze_tree(
    root,
    {
        "artifact_class": "corgi_regular_crossfit_profile_head_contract_preflight",
        "model_id": "corgi_regular",
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "biological_unit": "donor",
        "genomic_block_unit": "chromosome",
        "folds": 5,
        "tests": 6,
        "source_signal_arrays_read": False,
        "development_outcomes_read": False,
        "histology_or_disease_labels_read": False,
        "test_or_sealed_features_or_outcomes_read": False,
        "model_fit_performed": False,
        "evaluation_performed": False,
        "champion_or_external_claim_allowed": False,
        "global_census_modified": False,
        "contract_sha256": receipt["contract_sha256"],
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi cross-fit profile-head contract preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
