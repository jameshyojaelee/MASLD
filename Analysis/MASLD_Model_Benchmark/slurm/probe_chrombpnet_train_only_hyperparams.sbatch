#!/usr/bin/env bash
#SBATCH --job-name=masld-train-only-hp-probe
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

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/chrombpnet-v1.0.1-21063022"
SOURCE="${ACQUISITION}/source/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c"
BIAS_ROOT="${ROOT}/executions/chrombpnet-bias-models-zenodo7443683-21063957"
BIAS_MODEL="${BIAS_ROOT}/models/bias_models/ATAC/scATAC_dermal_fibroblast.h5"
FOLD_PROBE="${ROOT}/executions/chrombpnet-real-fold-regions-probe-21063996"
TARGET="${ROOT}/executions/chrombpnet-train-only-hyperparams-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
PYTHONPATH_VALUE="${ROOT}:${CORE}/site-packages:${SOURCE}"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite train-only hyperparameter probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/fixture" "${STAGE}/validation"
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 04bd59c825874f2d2d7dbf38b741633ccacad4735e95997a7477288e9c9656c4 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 8ecc40f452328df965be823d14f9ce76e2777407dba62cc7826685be635bb6e5 \
  "${BIAS_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 3fcadc521716156405650b1f439bce81c074193c226b4883471b328f1c7b27be \
  "${BIAS_MODEL}" | sha256sum --check --strict
printf '%s  %s\n' 4a49d56c1fd4084392fbbff0b85a1a4b65e9120f63f65dfa2e69423b8a01df91 \
  "${FOLD_PROBE}/ARTIFACTS.json" | sha256sum --check --strict
cp "${FOLD_PROBE}/fold.json" "${STAGE}/fixture/actual.fold.json"

module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
singularity --version > "${STAGE}/validation/singularity.version.txt"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python - "${STAGE}/fixture" <<'PY'
from pathlib import Path
import json
import sys

import pyBigWig
import pyfaidx

from scripts.prepare_chrombpnet_train_only_hyperparams import PRIMARY_CONTIGS

root = Path(sys.argv[1]).resolve(strict=True)
contig_length = 12000
sequence = "ACGT" * (contig_length // 4)
fasta = root / "fixture.fa"
with fasta.open("x", encoding="utf-8") as handle:
    for contig in PRIMARY_CONTIGS:
        handle.write(f">{contig}\n{sequence}\n")
reference = pyfaidx.Fasta(str(fasta), rebuild=True)
reference.close()
bigwig_path = root / "fixture.tn5.bw"
bigwig = pyBigWig.open(str(bigwig_path), "w")
bigwig.addHeader([(contig, contig_length) for contig in PRIMARY_CONTIGS])
for index, contig in enumerate(PRIMARY_CONTIGS):
    starts = [2000, 4500, 7000, 9500]
    ends = [3000, 5500, 8000, 10500]
    values = [
        1.0 + index * 0.05,
        0.2 + index * 0.01,
        2.0 + index * 0.05,
        0.3 + index * 0.01,
    ]
    bigwig.addEntries([contig] * 4, starts, ends=ends, values=values)
bigwig.close()
with (root / "peaks.bed").open("x", encoding="utf-8") as peaks:
    for contig in PRIMARY_CONTIGS:
        peaks.write(f"{contig}\t2490\t2510\tpeak_a\t100\t.\t10\t5\t5\t10\n")
        peaks.write(f"{contig}\t7490\t7510\tpeak_b\t100\t.\t10\t5\t5\t10\n")
with (root / "nonpeaks.bed").open("x", encoding="utf-8") as nonpeaks:
    for contig in PRIMARY_CONTIGS:
        nonpeaks.write(f"{contig}\t4990\t5010\tnonpeak_a\t0\t.\t0\t-1\t-1\t10\n")
        nonpeaks.write(f"{contig}\t9990\t10010\tnonpeak_b\t0\t.\t0\t-1\t-1\t10\n")
(root / "fixture.json").write_text(
    json.dumps(
        {
            "contigs": len(PRIMARY_CONTIGS),
            "contig_length": contig_length,
            "peaks": 48,
            "nonpeaks": 48,
            "signal_intervals_per_contig": 4,
        },
        sort_keys=True,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)
PY

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python scripts/prepare_chrombpnet_train_only_hyperparams.py \
    --genome "${STAGE}/fixture/fixture.fa" \
    --bigwig "${STAGE}/fixture/fixture.tn5.bw" \
    --peaks "${STAGE}/fixture/peaks.bed" \
    --nonpeaks "${STAGE}/fixture/nonpeaks.bed" \
    --actual-fold "${STAGE}/fixture/actual.fold.json" \
    --bias-model "${BIAS_MODEL}" \
    --output "${STAGE}/prepared" \
    --seed 20260824 \
    --negative-sampling-ratio 0.1 \
    --outlier-threshold 0.9999 \
    --input-length 2114 \
    --output-length 1000 \
    --max-jitter 500 \
    --filters 8 \
    --dilation-layers 2 \
    > "${STAGE}/validation/preparation.stdout.txt"

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

import numpy as np
import tensorflow as tf

from chrombpnet.training.utils.losses import multinomial_nll

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "prepared" / "summary.json").read_text())
expected_input = {"train": 28, "valid": 10, "test": 10}
if summary["input_peak_role_census"] != expected_input:
    raise SystemExit(f"synthetic peak role census differs: {summary!r}")
if summary["input_nonpeak_role_census"] != expected_input:
    raise SystemExit("synthetic nonpeak role census differs")
if summary["filtered_peak_role_census"]["valid"] != 10 or summary[
    "filtered_peak_role_census"
]["test"] != 10:
    raise SystemExit("held peak rows changed")
if summary["filtered_nonpeak_role_census"]["valid"] != 10 or summary[
    "filtered_nonpeak_role_census"
]["test"] != 10:
    raise SystemExit("held nonpeak rows changed")
if (
    summary["validation_outcomes_used_for_hyperparameter_fit"]
    or summary["test_outcomes_used_for_hyperparameter_fit"]
    or not summary["final_parameter_paths_are_artifact_relative"]
):
    raise SystemExit("train-only or artifact-relative contract failed")
parameters = summary["final_model_parameters"]
if parameters["chr_fold_path"] != "actual_fold.json" or parameters[
    "bias_model_path"
] != "train_only_bias_model_scaled.h5":
    raise SystemExit("final portable parameter bindings differ")
scaled = root / "prepared" / "train_only_bias_model_scaled.h5"
model = tf.keras.models.load_model(
    scaled,
    custom_objects={"multinomial_nll": multinomial_nll, "tf": tf},
    compile=False,
)
rng = np.random.default_rng(20260824)
sequence = np.eye(4, dtype=np.float32)[rng.integers(0, 4, size=(2, 2114))]
predictions = model.predict(sequence, verbose=0)
if len(predictions) != 2 or any(not np.isfinite(value).all() for value in predictions):
    raise SystemExit("scaled bias checkpoint forward pass differs")
validation = {
    "schema_version": "masld-bench-chrombpnet-train-only-hyperparams-probe-v1",
    "status": "pass",
    "input_role_census": expected_input,
    "filtered_peak_role_census": summary["filtered_peak_role_census"],
    "filtered_nonpeak_role_census": summary["filtered_nonpeak_role_census"],
    "validation_outcomes_used": False,
    "test_outcomes_used": False,
    "portable_parameter_bindings": True,
    "scaled_bias_forward_passed": True,
    "source_bias_role": "restricted_smoke_comparator",
}
(root / "validation" / "probe.json").write_text(
    json.dumps(validation, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(validation, sort_keys=True))
PY
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/prepare_chrombpnet_train_only_hyperparams.py" \
  "${ROOT}/slurm/probe_chrombpnet_train_only_hyperparams.sbatch" \
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
            "artifact_class": "chrombpnet_train_only_hyperparams_probe",
            "source_revision": "eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c",
            "source_bias_sha256": "3fcadc521716156405650b1f439bce81c074193c226b4883471b328f1c7b27be",
            "source_bias_role": "restricted_smoke_comparator",
            "hyperparameter_fit_outcome_roles": ["train"],
            "portable_parameter_bindings": True,
            "status": "passed",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet train-only hyperparameter probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
