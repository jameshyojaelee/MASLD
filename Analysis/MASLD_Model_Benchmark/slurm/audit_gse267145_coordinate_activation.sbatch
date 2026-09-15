#!/usr/bin/env bash
#SBATCH --job-name=gse267145-coord-task
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SUPPLEMENT="${ROOT}/executions/gse267145-supplement3-21070402"
METADATA="${ROOT}/executions/geo-multimodal-source-audit-21064310"
MATRIX="${ROOT}/executions/gse267145-matrix-audit-21064333"
JOIN="${ROOT}/executions/gse267145-authoritative-join-21064930"
MEASUREMENT="${ROOT}/executions/gse267145-rna-measurement-21065097"
REFERENCE="${ROOT}/executions/gse267145-reference-crosswalk-21065154"
QC_RIGHTS="${ROOT}/executions/gse267145-qc-rights-21065336"
TASK_SPEC="${ROOT}/config/evaluation/paired_bulk_rna_h3k27ac_task.toml"
PROMOTION_GATE="${ROOT}/config/evaluation/paired_bulk_rna_h3k27ac_promotion_gate.json"
TARGET="${ROOT}/executions/gse267145-coordinate-activation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE267145 coordinate activation audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_gse267145_coordinate_activation
PYTHONPATH="src" python - "${SUPPLEMENT}" "${METADATA}" "${MATRIX}" "${JOIN}" "${MEASUREMENT}" "${REFERENCE}" "${QC_RIGHTS}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
print("input frozen-tree verification passed")
PY
PYTHONPATH="src:." python scripts/audit_gse267145_coordinate_activation.py \
  --supplement "${SUPPLEMENT}" \
  --supplement-artifacts-sha256 6bebedcc10f92840d62e6337cc2f96c38e8cc873c89d85f15e4c85ade70f7dce \
  --metadata "${METADATA}" \
  --metadata-artifacts-sha256 bac433336705b74dcb2a6ef6b6594af10fef362769d564adfea9788fe546e88e \
  --matrix "${MATRIX}" \
  --matrix-artifacts-sha256 4e0d17f392530f4e6f0b54346d5d9333e450b7606d9698fbf7253a4f2a0e0417 \
  --join "${JOIN}" \
  --join-artifacts-sha256 e6c539c5fb29cd357dc126073779948c73c27fd219ff7a4d80b3d20a6cddb580 \
  --measurement "${MEASUREMENT}" \
  --measurement-artifacts-sha256 38531e463837cafd2c8e321411048c7b8b493876c3c964c0f50790cac8806f2f \
  --reference "${REFERENCE}" \
  --reference-artifacts-sha256 da491302ea738a92017978518239bc3884a7635d580e01e3bd28cd4f55b630c2 \
  --qc-rights "${QC_RIGHTS}" \
  --qc-rights-artifacts-sha256 ac7c92fd01b2b4240b380e4f6a58d98f1a6f853f8d0ad98d7167d7af34baaa73 \
  --task-spec "${TASK_SPEC}" \
  --task-spec-sha256 8bf05fe693c787ffdd8fe11443e019e668ecd03318c708e390785edd837b1a28 \
  --promotion-gate "${PROMOTION_GATE}" \
  --promotion-gate-sha256 f26a2048346cb60f810407971c41bbc22fe862444ee5438c90da954eb44b21a1 \
  --output "${STAGE}"
python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
sha256sum \
  scripts/audit_gse267145_coordinate_activation.py \
  tests/unit/test_audit_gse267145_coordinate_activation.py \
  slurm/audit_gse267145_coordinate_activation.sbatch \
  config/evaluation/paired_bulk_rna_h3k27ac_task.toml \
  config/evaluation/paired_bulk_rna_h3k27ac_promotion_gate.json \
  > "${STAGE}/source.sha256"
PYTHONPATH="src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
coordinate = json.loads((root / "coordinate_audit.json").read_text(encoding="utf-8"))
contract = json.loads((root / "activation_contract.json").read_text(encoding="utf-8"))
if coordinate["coordinate_semantics_resolved"] or coordinate["sequence_extraction_allowed"]:
    raise SystemExit("coordinate audit did not fail closed")
if not contract["task_specific_model_fitting_allowed"]:
    raise SystemExit("coordinate-independent task was not activated")
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_coordinate_and_task_activation",
        "cohort_family_id": "gse267145_znf469_human_liver",
        "task_id": "paired_bulk_rna_h3k27ac",
        "participants": 99,
        "pairing": "same_sample_different_aliquot",
        "coordinate_semantics_resolved": False,
        "sequence_extraction_allowed": False,
        "task_specific_coordinate_independent_training_allowed": True,
        "source_outcomes_allowed_for_model_selection": False,
        "single_cell_methods_allowed": False,
        "champion_eligible": False,
        "derivative_weight_release_review_required": True,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE267145 coordinate/task activation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
