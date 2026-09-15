#!/usr/bin/env bash
#SBATCH --job-name=g256-view-check
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET=${ROOT}/executions/gse256398-unlabeled-smoke-validation-${SLURM_JOB_ID}

if [[ -e "${TARGET}" ]]; then
  printf 'Refusing to overwrite GSE256398 view validation: %s\n' "${TARGET}" >&2
  exit 1
fi
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH=${ROOT}/src:${ROOT}
mkdir "${TARGET}"
python -m unittest -v \
  tests.unit.test_build_gse256398_qc_sensitivity \
  tests.unit.test_build_gse256398_unlabeled_smoke \
  tests.contract.test_cross_cohort_plan \
  > "${TARGET}/tests.stdout.txt" \
  2> "${TARGET}/tests.stderr.txt"
python - "${ROOT}" > "${TARGET}/validation.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree
from masld_bench.contracts import DatasetManifest
from masld_bench.registry import _validate_dataset_view, load_dataset_view_contract

root = Path(sys.argv[1]).resolve(strict=True)
config = root / "config"
parent_path = config / "datasets/gse256398.toml"
view_path = config / "dataset_views/gse256398_rna_unlabeled_smoke_1040_v1.toml"
parent = DatasetManifest.load_toml(parent_path)
view = load_dataset_view_contract(view_path)
_validate_dataset_view(
    view=view,
    parent=parent,
    parent_path=parent_path,
    config_root=config,
)
manifest = verify_frozen_tree(
    root / "executions/gse256398-unlabeled-smoke-21069823"
)
if (
    not parent.admission_blocking
    or parent.activation.ready
    or view.allowed_waves != ("smoke",)
    or view.row_count != 1040
    or view.biological_unit_count != 26
    or view.class_counts != {"unlabeled": 1040}
    or manifest["metadata"]["compatibility_smoke_only"] is not True
):
    raise SystemExit("GSE256398 parent/view fail-closed contract differs")
print(json.dumps({
    "dataset_id": parent.dataset_id,
    "full_activation_ready": parent.activation.ready,
    "parent_admission_blocking": parent.admission_blocking,
    "view_id": view.view_id,
    "allowed_waves": list(view.allowed_waves),
    "rows": view.row_count,
    "donors": view.biological_unit_count,
    "recursive_artifacts_verified": True,
}, sort_keys=True))
PY
python - "${TARGET}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1])
validation = json.loads((root / "validation.json").read_text(encoding="utf-8"))
freeze_tree(root, {
    "artifact_class": "gse256398_unlabeled_smoke_validation",
    "dataset_id": "gse256398",
    "view_id": validation["view_id"],
    "full_activation_ready": validation["full_activation_ready"],
    "smoke_view_ready": True,
})
verify_frozen_tree(root)
PY
sha256sum "${TARGET}/ARTIFACTS.json"
