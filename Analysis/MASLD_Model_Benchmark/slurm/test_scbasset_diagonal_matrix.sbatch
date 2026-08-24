#!/usr/bin/env bash
#SBATCH --job-name=masld-scbasset-matrix-gate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/scbasset-diagonal-matrix-gate-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"
shell_scripts=(
  "${ROOT}/slurm/scbasset_build_donor0_genomic0_inputs.sbatch"
  "${ROOT}/slurm/scbasset_train_donor0_genomic0_seed11.sbatch"
  "${ROOT}/slurm/scbasset_predict_valid_donor0_genomic0_seed11.sbatch"
)
python_scripts=(
  "${ROOT}/scripts/scbasset_build_fold_inputs.py"
  "${ROOT}/scripts/scbasset_train_fold.py"
  "${ROOT}/scripts/scbasset_predict_fold.py"
  "${ROOT}/src/masld_bench/adapters/scbasset_sequence_profile.py"
)
for script in "${shell_scripts[@]}"; do bash -n "${script}"; done
printf '%s  %s\n' 02f76f4f7af62063a606943e2b1bcf11ca0c8d8a63f29879ed31bfc1fbc30dc0 \
  "${ROOT}/executions/scbasset-donor0-genomic0-inputs-21064455/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 06db5a37eed92ebb6fd95b47f10cb998121cc0a641cab42686c27d47660e9e00 \
  "${ROOT}/executions/scbasset-donor0-genomic0-seed11-21064694/ARTIFACTS.json" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${python_scripts[@]}" \
  > "${STAGE}/validation/gate.json" <<'PY'
from pathlib import Path
import ast
import json
import sys

paths = [Path(value).resolve(strict=True) for value in sys.argv[1:]]
for path in paths:
    ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
train = paths[1].read_text(encoding="utf-8")
predict = paths[2].read_text(encoding="utf-8")
for marker in ("--split-id", "--donor-test-fold", "--donor-valid-fold"):
    if marker not in train or marker not in predict:
        raise SystemExit(f"missing dynamic fold marker: {marker}")
print(json.dumps({
    "schema_version": "masld-bench-scbasset-diagonal-matrix-gate-v1",
    "status": "pass", "outer_folds": list(range(5)),
    "screening_seeds": [20260824, 20260825, 20260826],
    "models_per_fold": 3, "prediction_roles": ["valid"],
    "fold0_legacy_fixture_receipts_verified": True,
    "shell_syntax_passed": True, "python_ast_passed": True,
}, sort_keys=True, indent=2))
PY
sha256sum "${shell_scripts[@]}" "${python_scripts[@]}" \
  "${ROOT}/slurm/test_scbasset_diagonal_matrix.sbatch" > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/gate.json").read_text())
freeze_tree(root, {"artifact_class": "scbasset_diagonal_matrix_gate", **receipt})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
