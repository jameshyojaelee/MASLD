#!/usr/bin/env bash
#SBATCH --job-name=masld-bpnet-train-probe
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
ACQUISITION="${ROOT}/executions/upstream-acquisitions/bpnet-v2.0.0-21064186"
SOURCE="${ACQUISITION}/source/bpnet-f4593eedca51741f25b8c5cc0c0d647faf8c8a3c"
RUNTIME_PROBE="${ROOT}/executions/bpnet-native-runtime-probe-21064198"
TARGET="${ROOT}/executions/bpnet-training-inference-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
PYTHONPATH_VALUE="${CORE}/site-packages:${SOURCE}:${ROOT}/scripts"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite BPNet training probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/fixture" "${STAGE}/validation"

printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 8062574bfaa9067b0aca82ffa3eeae87188140ee5f5842399171ad4535cf88f6 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' cd25afab907693c048eb71b0764db81984c5c664d0b267999682dc10cb3b05e4 \
  "${RUNTIME_PROBE}/ARTIFACTS.json" | sha256sum --check --strict

module load singularity/4.2.2
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python - "${STAGE}/fixture" <<'PY'
from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import pyBigWig


root = Path(sys.argv[1]).resolve(strict=True)
contigs = [f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"]
length = 20_000
(root / "fixture.chrom.sizes").write_text(
    "".join(f"{contig}\t{length}\n" for contig in contigs), encoding="utf-8"
)
sequence = ("ACGT" * (length // 4 + 1))[:length]
(root / "fixture.fa").write_text(
    "".join(f">{contig}\n{sequence}\n" for contig in contigs), encoding="utf-8"
)

train_centers = [2_000 + 700 * index for index in range(16)]
valid_centers = [2_000 + 900 * index for index in range(8)]
test_centers = [3_000, 6_000]
with (root / "regions.narrowPeak").open("x", encoding="utf-8") as handle:
    for role, contig, centers in (
        ("train", "chr1", train_centers),
        ("valid", "chr21", valid_centers),
        ("test", "chr22", test_centers),
    ):
        for index, center in enumerate(centers):
            start = center - 100
            handle.write(
                f"{contig}\t{start}\t{center + 100}\t{role}_{index}\t0\t.\t0\t-1\t-1\t100\n"
            )
with (root / "nonpeaks.bed").open("x", encoding="utf-8") as handle:
    for contig, centers in (("chr1", train_centers), ("chr21", valid_centers)):
        for index, center in enumerate(centers):
            start = center + 250
            handle.write(f"{contig}\t{start}\t{start + 200}\tnegative_{index}\n")

bigwig = pyBigWig.open(str(root / "fixture.tn5.bw"), "w")
bigwig.addHeader([(contig, length) for contig in contigs], maxZooms=0)
for contig, centers in (("chr1", train_centers), ("chr21", valid_centers)):
    starts = []
    ends = []
    values = []
    for center in centers:
        for position in range(center - 50, center + 51):
            starts.append(position)
            ends.append(position + 1)
            values.append(float(1 + (position - center) % 5))
    bigwig.addEntries([contig] * len(starts), starts=starts, ends=ends, values=values)
bigwig.close()

fold = {
    "train": [f"chr{index}" for index in range(1, 21)],
    "valid": ["chr21"],
    "test": ["chr22", "chrX", "chrY"],
}
(root / "fold.json").write_text(
    json.dumps(fold, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)

fields = (
    "contig", "output_start", "output_end", "window_id", "genomic_fold",
    "window_class", "ccre_class", "ccre_id", "ccre_start", "ccre_end",
    "input_start", "input_end", "selection_hash",
)
with (root / "held_windows.tsv").open("x", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for role, contig, centers, genomic_fold in (
        ("valid", "chr21", [4_000, 8_000], "1"),
        ("test", "chr22", [4_000, 8_000], "0"),
    ):
        for index, center in enumerate(centers):
            writer.writerow({
                "contig": contig,
                "output_start": center - 500,
                "output_end": center + 500,
                "window_id": f"{role}_{index}",
                "genomic_fold": genomic_fold,
                "window_class": "fixed_ccre",
                "ccre_class": "PLS",
                "ccre_id": f"fixture_{role}_{index}",
                "ccre_start": center - 100,
                "ccre_end": center + 100,
                "input_start": center - 1057,
                "input_end": center + 1057,
                "selection_hash": f"fixture_hash_{role}_{index}",
            })
PY

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env PYTHONHASHSEED=20260824 \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python "${ROOT}/scripts/bpnet_run_training.py" \
    --genome "${STAGE}/fixture/fixture.fa" \
    --chrom-sizes "${STAGE}/fixture/fixture.chrom.sizes" \
    --bigwig "${STAGE}/fixture/fixture.tn5.bw" \
    --regions "${STAGE}/fixture/regions.narrowPeak" \
    --nonpeaks "${STAGE}/fixture/nonpeaks.bed" \
    --fold "${STAGE}/fixture/fold.json" \
    --output "${STAGE}/model" \
    --epochs 2 \
    --batch-size 4 \
    --threads 1 \
    --max-jitter 8 \
    --negative-ratio 0.5 \
    --early-stop 1 \
    --seed 20260824 \
    > "${STAGE}/validation/training.stdout.txt" \
    2> "${STAGE}/validation/training.stderr.txt"

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  "${IMAGE}" python "${ROOT}/scripts/bpnet_predict_ccre_windows.py" \
    --model "${STAGE}/model/bpnet.inference.h5" \
    --fasta "${STAGE}/fixture/fixture.fa" \
    --windows "${STAGE}/fixture/held_windows.tsv" \
    --fold "${STAGE}/fixture/fold.json" \
    --output "${STAGE}/predictions" \
    --batch-size 2 \
    --expected-windows 4 \
    --seed 20260824 \
    > "${STAGE}/validation/inference.stdout.txt" \
    2> "${STAGE}/validation/inference.stderr.txt"

singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${CORE}/site-packages:${ROOT}/scripts" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

import h5py
import numpy as np

root = Path(sys.argv[1]).resolve(strict=True)
training = json.loads((root / "model" / "training_receipt.json").read_text())
prediction = json.loads((root / "predictions" / "summary.json").read_text())
if training["architecture"]["parameter_count"] != 109124:
    raise SystemExit("BPNet parameter count differs")
if training["test_contigs_exposed"] or training["evaluator_outcomes_exposed"]:
    raise SystemExit("BPNet training leakage flag differs")
if prediction["role_counts"] != {"valid": 2, "test": 2}:
    raise SystemExit("BPNet prediction role census differs")
with h5py.File(root / "predictions" / "profile_probabilities.h5") as handle:
    values = handle["profile_probability"][:]
if values.shape != (4, 1000) or not np.isfinite(values).all():
    raise SystemExit("BPNet profile predictions differ")
np.testing.assert_allclose(values.sum(axis=1), 1.0, rtol=0, atol=1e-6)
with (root / "predictions" / "regional_counts.tsv").open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if len(rows) != 4 or any(float(row["strand_averaged_mass"]) < 0 for row in rows):
    raise SystemExit("BPNet count predictions differ")
(root / "validation" / "probe.json").write_text(
    json.dumps({
        "schema_version": "masld-bench-bpnet-training-inference-probe-v1",
        "status": "pass",
        "source_revision": "f4593eedca51741f25b8c5cc0c0d647faf8c8a3c",
        "seed": 20260824,
        "epochs_completed": training["completed_epochs"],
        "parameter_count": training["architecture"]["parameter_count"],
        "bias_or_control_input": False,
        "test_outcomes_used_for_training": False,
        "inference_observed_atac_exposed": False,
    }, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
PY

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/bpnet_contract.py" \
  "${ROOT}/scripts/bpnet_run_training.py" \
  "${ROOT}/scripts/bpnet_predict_ccre_windows.py" \
  "${ROOT}/slurm/bpnet_probe_training_inference.sbatch" \
  > "${STAGE}/validation/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(freeze_tree(root, {
    "artifact_class": "bpnet_training_inference_probe",
    "model_id": "bpnet",
    "source_revision": "f4593eedca51741f25b8c5cc0c0d647faf8c8a3c",
    "runtime_probe_artifacts_sha256": "cd25afab907693c048eb71b0764db81984c5c664d0b267999682dc10cb3b05e4",
    "seed": 20260824,
    "test_outcomes_used_for_training": False,
    "status": "passed",
}))
PY
mv "${STAGE}" "${TARGET}"
printf 'BPNet training/inference probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
