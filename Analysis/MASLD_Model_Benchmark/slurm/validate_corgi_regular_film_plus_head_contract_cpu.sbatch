#!/usr/bin/env bash
#SBATCH --job-name=model-check-260
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
CONTRACT="${ROOT}/config/evaluation/corgi_regular_film_plus_head_smoke.toml"
TARGET="${ROOT}/executions/model-check-260-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi FiLM-plus-head contract preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
for pair in \
  "615f3b37d094a6d6031c7efd1fd68d8f1f37e0446727aab2bf5abd1fa4ff24b9 config/evaluation/corgi_regular_film_plus_head_smoke.toml" \
  "61fc4adc976ef6f4c8cf5f3b35d349ae756a3a94a7d933ff5f058adae44c047b scripts/freeze_corgi_regular_film_plus_head_contract.py" \
  "56db3e870af4f738021bf1b8dbac3d10af15bd8e37c36d2abeef4bdb847c3d01 tests/unit/test_freeze_corgi_regular_film_plus_head_contract.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_freeze_corgi_regular_film_plus_head_contract \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="src:." python scripts/freeze_corgi_regular_film_plus_head_contract.py \
  --project-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/contract" \
  > "${STAGE}/validation/contract.stdout.json"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  config/evaluation/corgi_regular_film_plus_head_smoke.toml \
  scripts/freeze_corgi_regular_film_plus_head_contract.py \
  tests/unit/test_freeze_corgi_regular_film_plus_head_contract.py \
  slurm/validate_corgi_regular_film_plus_head_contract_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "contract/contract_receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_prospective_contract_execution_blocked"
    or receipt["crossfit_folds"] != 5
    or receipt["trainable_parameter_numel"] != 5869229
    or receipt["development_outcomes_read"]
    or receipt["histology_or_disease_labels_read"]
    or receipt["test_or_sealed_features_or_outcomes_read"]
    or receipt["training_execution_authorized"]
    or receipt["prediction_execution_authorized"]
    or receipt["evaluation_execution_authorized"]
    or receipt["global_census_modified"]
    or receipt["global_promotion_gate_modified"]
):
    raise SystemExit("Corgi FiLM-plus-head contract receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "corgi_regular_film_plus_head_prospective_contract",
        "model_id": "corgi_regular",
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "lineage_id": "hepatocyte",
        "adaptation_rung": "film_plus_head",
        "biological_unit": "donor",
        "genomic_block_unit": "chromosome",
        "folds": 5,
        "tests": 6,
        "trainable_parameter_numel": receipt["trainable_parameter_numel"],
        "development_outcomes_read": False,
        "histology_or_disease_labels_read": False,
        "test_or_sealed_features_or_outcomes_read": False,
        "training_execution_authorized": False,
        "prediction_execution_authorized": False,
        "evaluation_execution_authorized": False,
        "champion_or_external_claim_allowed": False,
        "global_census_modified": False,
        "global_promotion_gate_modified": False,
        "contract_sha256": receipt["contract_sha256"],
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi FiLM-plus-head prospective contract preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
