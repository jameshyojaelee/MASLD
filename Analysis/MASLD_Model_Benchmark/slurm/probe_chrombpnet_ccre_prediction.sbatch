#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-ccre-predict
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONHASHSEED=20260824

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/chrombpnet-v1.0.1-21063022"
SOURCE="${ACQUISITION}/source/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c"
MODEL_ROOT="${ROOT}/executions/chrombpnet-direct-training-probe-21064025"
FULL_MODEL="${MODEL_ROOT}/model/chrombpnet.h5"
NOBIAS_MODEL="${MODEL_ROOT}/model/chrombpnet_nobias.h5"
FOLD_ROOT="${ROOT}/executions/chrombpnet-real-fold-regions-probe-21063996"
FOLD="${FOLD_ROOT}/fold.json"
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
WINDOWS="${SPLIT}/ccre_evaluation_windows.tsv"
FASTA_ROOT="${ROOT}/executions/sequence-fasta-21063808"
FASTA="${FASTA_ROOT}/GRCh38.p14.sequence_model.fa"
TARGET="${ROOT}/executions/chrombpnet-ccre-prediction-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
PYTHONPATH_VALUE="${ROOT}:${CORE}/site-packages:${SOURCE}"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite cCRE prediction probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/predictions" "${STAGE}/validation"
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 04bd59c825874f2d2d7dbf38b741633ccacad4735e95997a7477288e9c9656c4 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 2ed8603984659578eda44b327f5aee478686bc2995e45dd13a26e833e4c08896 \
  "${MODEL_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 29c7e577a78600a43422dca8e29c69f622c0ecf32557886b3955844169c8cfda \
  "${FULL_MODEL}" | sha256sum --check --strict
printf '%s  %s\n' 9551a4ccf3ec3df6350cc273cedae92c77c5269d5feaec9df2ee9d53fa3864f3 \
  "${NOBIAS_MODEL}" | sha256sum --check --strict
printf '%s  %s\n' 4a49d56c1fd4084392fbbff0b85a1a4b65e9120f63f65dfa2e69423b8a01df91 \
  "${FOLD_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' c23fb8dd996a6fbc3f15d5f5b61c01a7db2da5d63fbf9fc65a1bf0dd1340ebee \
  "${FOLD}" | sha256sum --check --strict
printf '%s  %s\n' 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
  "${SPLIT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 3dbefbb961c040236431773fe3e2e32527622e250e9fc86133ed117dc4420c4c \
  "${WINDOWS}" | sha256sum --check --strict
printf '%s  %s\n' 7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 \
  "${FASTA_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' ed20dceeb315cffa2bc3de050d846eceefc552451fef741f0dd95b0149d10fce \
  "${FASTA}" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${MODEL_ROOT}" "${FOLD_ROOT}" "${SPLIT}" "${FASTA_ROOT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("input frozen-tree verification passed")
PY

module purge
module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
singularity --version > "${STAGE}/validation/singularity.version.txt"
cd "${ROOT}"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env PYTHONHASHSEED=20260824 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python scripts/predict_chrombpnet_ccre_windows.py \
    --model "${FULL_MODEL}" \
    --fasta "${FASTA}" \
    --windows "${WINDOWS}" \
    --fold "${FOLD}" \
    --output "${STAGE}/predictions/full_model" \
    --batch-size 512 \
    --seed 20260824 \
    --count-transform log1p_absolute \
    > "${STAGE}/validation/full-model-prediction.stdout.txt" \
    2> "${STAGE}/validation/full-model-prediction.stderr.txt"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env PYTHONHASHSEED=20260824 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python scripts/predict_chrombpnet_ccre_windows.py \
    --model "${NOBIAS_MODEL}" \
    --fasta "${FASTA}" \
    --windows "${WINDOWS}" \
    --fold "${FOLD}" \
    --output "${STAGE}/predictions/nobias_model" \
    --batch-size 512 \
    --seed 20260824 \
    --count-transform log_component \
    > "${STAGE}/validation/nobias-model-prediction.stdout.txt" \
    2> "${STAGE}/validation/nobias-model-prediction.stderr.txt"

singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

import h5py
import numpy as np

root = Path(sys.argv[1]).resolve(strict=True)
transforms = {"full_model": "log1p_absolute", "nobias_model": "log_component"}
for variant, transform in transforms.items():
    prediction = root / "predictions" / variant
    summary = json.loads((prediction / "summary.json").read_text())
    if (
        summary["windows"] != 32000
        or summary["role_counts"] != {"valid": 16000, "test": 16000}
        or summary["regional_mass_transform"] != transform
        or summary["observed_atac_input_exposed"]
        or summary["benchmark_metrics_calculated"]
    ):
        raise SystemExit(f"{variant} cCRE prediction summary differs")
    with (prediction / "regional_counts.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    masses = [float(row["strand_averaged_mass"]) for row in rows]
    if len(rows) != 32000 or any(value < 0 for value in masses) or max(masses) <= 0:
        raise SystemExit(f"{variant} cCRE regional-count output differs")
    with h5py.File(prediction / "profile_probabilities.h5", "r") as handle:
        values = handle["profile_probability"]
        if values.shape != (32000, 1000):
            raise SystemExit(f"{variant} cCRE profile tensor shape differs")
        sample = np.asarray(values[[0, 15999, 16000, 31999]], dtype=np.float64)
    if not np.isfinite(sample).all() or not np.allclose(
        sample.sum(axis=1), 1.0, rtol=0, atol=1e-6
    ):
        raise SystemExit(f"{variant} cCRE profile tensor values differ")
validation = {
    "schema_version": "masld-bench-chrombpnet-ccre-prediction-probe-v1",
    "status": "pass",
    "windows": 32000,
    "valid_windows": 16000,
    "test_windows": 16000,
    "reverse_complement_tta": True,
    "model_variants": ["full_model", "nobias_model"],
    "primary_biological_score": "nobias_model",
    "primary_regional_mass_transform": "exp(log_sequence_component)",
    "assay_qc_score": "full_model",
    "assay_qc_regional_mass_transform": "max(expm1(log1p_total),0)",
    "observed_atac_input_exposed": False,
    "benchmark_metrics_calculated": False,
    "source_model_role": "restricted_synthetic_smoke_comparator",
}
(root / "validation" / "probe.json").write_text(
    json.dumps(validation, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(validation, sort_keys=True))
PY
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/predict_chrombpnet_ccre_windows.py" \
  "${ROOT}/slurm/probe_chrombpnet_ccre_prediction.sbatch" \
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
            "artifact_class": "chrombpnet_ccre_prediction_probe",
            "model_id": "chrombpnet",
            "source_model_artifacts_sha256": "2ed8603984659578eda44b327f5aee478686bc2995e45dd13a26e833e4c08896",
            "split_id": "donor0_genomic0",
            "windows": 32000,
            "reverse_complement_tta": True,
            "model_variants": ["full_model", "nobias_model"],
            "primary_biological_score": "nobias_model",
            "primary_regional_mass_transform": "exp(log_sequence_component)",
            "assay_qc_score": "full_model",
            "assay_qc_regional_mass_transform": "max(expm1(log1p_total),0)",
            "observed_atac_input_exposed": False,
            "benchmark_metrics_calculated": False,
            "source_model_role": "restricted_synthetic_smoke_comparator",
            "status": "passed",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet cCRE prediction probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
