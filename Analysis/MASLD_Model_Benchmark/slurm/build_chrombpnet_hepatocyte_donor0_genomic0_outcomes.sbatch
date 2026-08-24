#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-hep-d0g0-outcomes
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
DONOR_BIGWIGS="${ROOT}/executions/gse296875-deduplicated-tn5-bigwigs-21063992"
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
TARGET="${ROOT}/executions/chrombpnet-hepatocyte-donor0-genomic0-outcomes-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite ChromBPNet outcome artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f \
  "${DONOR_BIGWIGS}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
  "${SPLIT}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${DONOR_BIGWIGS}" "${SPLIT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("input frozen-tree verification passed")
PY

module purge
module load singularity/4.2.2
cd "${ROOT}"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python scripts/build_chrombpnet_development_outcomes.py \
    --donor-bigwigs "${DONOR_BIGWIGS}" \
    --donor-bigwigs-artifacts-sha256 \
      01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f \
    --split-contract "${SPLIT}" \
    --split-artifacts-sha256 \
      00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
    --split-id donor0_genomic0 \
    --lineage hepatocyte \
    --output "${STAGE}/outcomes" \
    --workers 8 \
    > "${STAGE}/validation/build.stdout.txt" \
    2> "${STAGE}/validation/build.stderr.txt"

singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

import h5py
import numpy as np

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "outcomes" / "summary.json").read_text())
expected = {
    "valid": {"donors": 9, "windows": 16000, "donor_fold": 1, "genomic_fold": 1},
    "test": {"donors": 11, "windows": 16000, "donor_fold": 0, "genomic_fold": 0},
}
for role, contract in expected.items():
    for field, value in contract.items():
        if summary["roles"][role][field] != value:
            raise SystemExit(f"{role} {field} differs")
if (
    summary["model_inference_input_eligible"]
    or summary["benchmark_metrics_calculated"]
    or not summary["observed_atac_outcomes"]
):
    raise SystemExit("evaluator-only outcome firewall differs")
with h5py.File(root / "outcomes" / "outcomes.h5", "r") as handle:
    for role, contract in expected.items():
        counts = handle[role]["counts"]
        if counts.shape != (contract["donors"], contract["windows"]):
            raise SystemExit(f"{role} outcome matrix shape differs")
        if np.any(np.asarray(counts[:, :100]) < 0):
            raise SystemExit(f"{role} outcome matrix contains negative counts")
validation = {
    "schema_version": "masld-bench-chrombpnet-development-outcomes-probe-v1",
    "status": "pass",
    "split_id": "donor0_genomic0",
    "lineage_id": "hepatocyte",
    "valid_donors": 9,
    "test_donors": 11,
    "windows_per_role": 16000,
    "observed_atac_outcomes": True,
    "model_inference_input_eligible": False,
    "benchmark_metrics_calculated": False,
}
(root / "validation" / "probe.json").write_text(
    json.dumps(validation, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(validation, sort_keys=True))
PY
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/build_chrombpnet_development_outcomes.py" \
  "${ROOT}/slurm/build_chrombpnet_hepatocyte_donor0_genomic0_outcomes.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "chrombpnet_development_outcomes",
            "dataset_id": "gse296875",
            "split_id": "donor0_genomic0",
            "lineage_id": "hepatocyte",
            "biological_unit": "donor",
            "donor_bigwigs_artifacts_sha256": "01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f",
            "split_artifacts_sha256": "00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f",
            "observed_atac_outcomes": True,
            "model_inference_input_eligible": False,
            "benchmark_metrics_calculated": False,
            "status": "passed",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet hepatocyte donor0/genomic0 outcomes: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
