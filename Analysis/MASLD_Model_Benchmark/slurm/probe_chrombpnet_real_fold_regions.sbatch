#!/usr/bin/env bash
#SBATCH --job-name=masld-real-fold-regions-probe
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
TARGET="${ROOT}/executions/chrombpnet-real-fold-regions-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite real fold-regions probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
  "${SPLIT}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="src:." python - "${STAGE}/training.narrowPeak" <<'PY'
from pathlib import Path
import sys

Path(sys.argv[1]).write_text(
    "chr3\t999990\t1000010\tfixture_train_peak\t100\t.\t10\t5\t5\t10\n",
    encoding="utf-8",
)
PY
PYTHONPATH="src:." python scripts/assemble_chrombpnet_fold_regions.py \
  --training-peaks "${STAGE}/training.narrowPeak" \
  --split-contract "${SPLIT}" \
  --split-id donor0_genomic0 \
  --output-peaks "${STAGE}/regions.bed" \
  --output-fold "${STAGE}/fold.json" \
  --output-summary "${STAGE}/summary.json"
PYTHONPATH="src:." python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "summary.json").read_text())
if summary["role_counts"] != {"train": 1, "valid": 16000, "test": 16000}:
    raise SystemExit(f"real role census differs: {summary['role_counts']!r}")
if summary["assembled_regions"] != 32001:
    raise SystemExit("real assembled-region total differs")
if summary["held_atac_used_for_region_selection"]:
    raise SystemExit("held ATAC was used for region selection")
with (root / "regions.bed").open(encoding="utf-8") as handle:
    region_rows = sum(1 for _line in handle)
if region_rows != 32001:
    raise SystemExit("real region-file row count differs")
print(json.dumps(summary, sort_keys=True))
PY
python -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/assemble_chrombpnet_fold_regions.py" \
  "${ROOT}/slurm/probe_chrombpnet_real_fold_regions.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "chrombpnet_real_fold_regions_probe",
            "split_id": "donor0_genomic0",
            "split_artifacts_sha256": "00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f",
            "assembled_regions": 32001,
            "held_atac_used_for_region_selection": False,
            "status": "passed",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet real fold-regions probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
