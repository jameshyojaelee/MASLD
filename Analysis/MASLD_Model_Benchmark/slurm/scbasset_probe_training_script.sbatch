#!/usr/bin/env bash
#SBATCH --job-name=masld-scbasset-train-probe
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
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
ACQUISITION="${ROOT}/executions/upstream-acquisitions/scbasset-aed3a6f713091fd988196b297e9c06e092ff1d22-21064184"
SOURCE="${ACQUISITION}/source/scBasset-aed3a6f713091fd988196b297e9c06e092ff1d22"
TARGET="${ROOT}/executions/scbasset-training-script-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite scBasset training probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/inputs" "${STAGE}/validation"

printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 9b0b418cea67bc788d7657913965c5e11ed62e36c35a8104c381c9007f6dffc5 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict

module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python - "${STAGE}/inputs" <<'PY'
from pathlib import Path
import csv
import json
import sys

import h5py
import numpy as np
from scipy import sparse

root = Path(sys.argv[1])
rng = np.random.default_rng(20260824)
for role, rows in (("train", 8), ("valid", 4)):
    with h5py.File(root / f"{role}_seqs.h5", "x") as handle:
        handle.create_dataset(
            "X", data=rng.integers(0, 4, size=(rows, 1344), dtype=np.int8)
        )
    values = rng.integers(0, 2, size=(rows, 15), dtype=np.int8)
    values[:, 0] = np.arange(rows) % 2
    sparse.save_npz(root / f"m_{role}.npz", sparse.csr_matrix(values))
with (root / "training_cells.tsv").open("x", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=("cell_id", "donor_id", "lineage", "outer_fold"),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    lineages = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
    for fold in (2, 3, 4):
        for lineage in lineages:
            writer.writerow(
                {
                    "cell_id": f"fold{fold}_{lineage}",
                    "donor_id": f"donor_{fold}",
                    "lineage": lineage,
                    "outer_fold": fold,
                }
            )
(root / "summary.json").write_text(
    json.dumps(
        {
            "status": "pass",
            "donor_train_folds": [2, 3, 4],
            "donor_valid_fold": 1,
            "donor_test_fold": 0,
            "held_donor_atac_read": False,
            "genomic_test_atac_parsed_only_to_exclude": True,
            "genomic_test_atac_used_for_matrix": False,
            "genomic_test_atac_exported": False,
            "model_training_input_contains_genomic_test_atac": False,
            "intermediate_lineage_matrices_removed": True,
            "missing_evidence_encoded_as_zero": False,
        },
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
PY

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python "${ROOT}/scripts/scbasset_train_fold.py" \
    --inputs "${STAGE}/inputs" \
    --source "${SOURCE}" \
    --output "${STAGE}/model" \
    --input-artifacts-sha256 aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa \
    --source-artifacts-sha256 9b0b418cea67bc788d7657913965c5e11ed62e36c35a8104c381c9007f6dffc5 \
    --seed 11 \
    --batch-size 4 \
    --learning-rate 0.01 \
    --max-epochs 2 \
    --early-stopping-patience 1 \
    > "${STAGE}/validation/training.stdout.txt"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "model" / "model_receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["n_training_cells"] != 15
    or receipt["public_tutorial_weights_used"]
    or receipt["held_donor_atac_used"]
    or receipt["genomic_test_atac_used"]
):
    raise SystemExit("scBasset training probe receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "scbasset_training_script_probe",
        "model_id": "scbasset",
        "source_revision": "aed3a6f713091fd988196b297e9c06e092ff1d22",
        "device_class": "NVIDIA_L40S_sm89",
        "epochs": 2,
        "public_tutorial_weights_used": False,
        "held_donor_atac_used": False,
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scBasset training script probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
