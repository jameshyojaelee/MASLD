#!/usr/bin/env bash
#SBATCH --job-name=masld-gse296875-train-pb
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=5
#SBATCH --mem=96G
#SBATCH --time=48:00:00
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
OVERLAY="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
FRAGMENTS="${ROOT}/executions/gse296875-donor-lineage-fragments-21063889"
FRAGMENTS_ARTIFACTS_SHA256=ca8f157669d5a6953765525b62ab593db1d4a352e99a73ed1b6db353aab35ef4
REFERENCE="${ROOT}/executions/sequence-reference-contract-21063670"
CHROM_SIZES="${REFERENCE}/GRCh38.p14.sequence_model.chrom.sizes"
TARGET="${ROOT}/executions/gse296875-training-fold-pseudobulk-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite training pseudobulk: %s\n' "${TARGET}" >&2
  exit 1
fi
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${OVERLAY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${FRAGMENTS_ARTIFACTS_SHA256}" "${FRAGMENTS}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' 7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d \
  "${REFERENCE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 1901777f65ef618569eabb62967c5ace55a63b4bef8795af7a8838fc7042858b \
  "${CHROM_SIZES}" | sha256sum --check --strict
cd "${ROOT}"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${FRAGMENTS}" "${REFERENCE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
verify_frozen_tree(Path(sys.argv[2]))
print("input frozen-tree verification passed")
PY

module purge
module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${OVERLAY}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python scripts/build_training_fold_pseudobulk.py \
    --fragments "${FRAGMENTS}" \
    --fragments-artifacts-sha256 "${FRAGMENTS_ARTIFACTS_SHA256}" \
    --chrom-sizes "${CHROM_SIZES}" \
    --chrom-sizes-sha256 1901777f65ef618569eabb62967c5ace55a63b4bef8795af7a8838fc7042858b \
    --workers 5 \
    --output "${STAGE}"
install -d -m 0750 "${STAGE}/validation"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${OVERLAY}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/build_training_fold_pseudobulk.py" \
  "${ROOT}/slurm/build_gse296875_training_fold_pseudobulk.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${FRAGMENTS}" "${REFERENCE}" \
  "${STAGE}/validation/input-verification.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

payload = {
    "fragments": verify_frozen_tree(Path(sys.argv[1])),
    "reference": verify_frozen_tree(Path(sys.argv[2])),
}
Path(sys.argv[3]).write_text(
    json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
PY
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "gse296875_training_fold_pseudobulk",
            "dataset_id": "gse296875",
            "biological_outer_unit": "donor",
            "outputs": 25,
            "workers": 5,
            "fragments_artifacts_sha256": "ca8f157669d5a6953765525b62ab593db1d4a352e99a73ed1b6db353aab35ef4",
            "reference_artifacts_sha256": "7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d",
            "algorithm": "training_fold_heap_merge_and_sweep_line_exact_Tn5",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'GSE296875 training-fold pseudobulk: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
