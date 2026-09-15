#!/usr/bin/env bash
#SBATCH --job-name=model-check-259
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
TARGET="${ROOT}/executions/model-check-259-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi negative-disposition audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
for pair in \
  "bff0cb47b7df5311859b38c0f0375db611c54b8d1a7bc7fe4f5896b21c54d9d1 config/artifacts/models/corgi/new_assay_head_development_disposition_20260825.json" \
  "4fcac8d09ecff9ac82aec359aa2b06c8f2a1044440832832b5f7d7aee1ea658a scripts/audit_corgi_regular_new_head_disposition.py" \
  "f43f3308cf1b8582a6cbf15d38848792fb8cb1d223248ced5ece7cbed1cba716 tests/unit/test_audit_corgi_regular_new_head_disposition.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_corgi_regular_new_head_disposition \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="src:." python scripts/audit_corgi_regular_new_head_disposition.py \
  > "${STAGE}/validation/audit.stdout.json"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  config/artifacts/models/corgi/new_assay_head_development_disposition_20260825.json \
  scripts/audit_corgi_regular_new_head_disposition.py \
  tests/unit/test_audit_corgi_regular_new_head_disposition.py \
  slurm/validate_corgi_regular_new_head_disposition_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/audit.stdout.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_frozen_negative_development_result"
    or receipt["selected_alpha_by_fold"] != [0.0] * 5
    or receipt["numeric_residual_is_biological_gain"]
    or receipt["incremental_corgi_profile_contribution"]
    or receipt["film_plus_head_execution_authorized"]
    or receipt["sealed_or_test_data_opened_by_audit"]
):
    raise SystemExit("Corgi negative-disposition receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "corgi_regular_new_assay_head_negative_development_disposition",
        "model_id": "corgi_regular",
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "lineage_id": "hepatocyte",
        "adaptation_rung": "new_assay_head",
        "biological_unit": "donor",
        "genomic_block_unit": "chromosome",
        "folds": 5,
        "tests": 3,
        "development_result": "negative_no_incremental_corgi_profile_contribution",
        "numeric_residual_is_biological_gain": False,
        "histology_or_disease_labels_read": False,
        "test_or_sealed_features_or_outcomes_read": False,
        "promotion_gate_evaluated": False,
        "global_census_modified": False,
        "global_promotion_gate_modified": False,
        "film_plus_head_execution_authorized": False,
        "authority_sha256": receipt["authority_sha256"],
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi negative new-assay-head disposition audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
