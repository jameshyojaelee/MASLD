#!/usr/bin/env bash
#SBATCH --job-name=model-check-337
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
SOURCE=${ROOT}/executions/model-training-500-21080428
HEADS=${ROOT}/executions/model-training-21082720
SPLIT=${ROOT}/executions/model-data-052-21077595/split
TARGET=${ROOT}/executions/model-check-337-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
trap 'status=$?; if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi; exit ${status}' EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then exit 1; fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="${ROOT}:${ROOT}/src"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${RUNTIME}/env/bin/python" - <<'PY' "${SOURCE}" "${HEADS}" "${SPLIT}" "${STAGE}"
import csv, json, sys
from pathlib import Path
import numpy as np
from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from scripts.predict_gse296875_transcriptformer_outer_ensemble import (
    ROSTER, SPLIT_SHA256, _head_logits, _softmax,
)
from masld_bench.hashing import sha256_file
source, heads, split, output = map(Path, sys.argv[1:])
expected = {
    source: "fc4601988491b086d5e89f0ff35f2f73d9c27cf74846a6ff7895efc130c825a5",
    heads: "f0cefd02290c2b2135c3578784076257e2a63213a99e670bf848a47f18a1e27e",
    split: SPLIT_SHA256,
}
for root, digest in expected.items():
    if sha256_file(root / "ARTIFACTS.json") != digest: raise SystemExit("input hash differs")
    verify_frozen_tree(root)
with np.load(source / "embeddings/common_embeddings.npz", allow_pickle=False) as values:
    embeddings = np.asarray(values["embeddings"], dtype=np.float32)
    row_ids = values["row_ids"].astype(str)
    activation = np.asarray(values["outer_folds"], dtype=np.int8)
with (split / "row_outer_folds.tsv").open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
outer = np.asarray([int(row["outer_fold"]) for row in rows], dtype=np.int8)
if [row["row_id"] for row in rows] != row_ids.tolist() or int(np.sum(activation != outer)) != 42584:
    raise SystemExit("study split alignment differs")
held = np.flatnonzero(outer == 0)
errors = {}
for head_id in ("linear", "two_layer_mlp"):
    shard = heads / "shards" / f"{head_id}__seed1103__fold0"
    receipt = json.loads((shard / "prediction_receipt.json").read_text())
    reproduced = _softmax(
        _head_logits(shard / "states/final_outer_training_refit.npz", embeddings[held], head_id),
        float(receipt["temperature"]),
    )
    with (shard / "predictions.tsv").open() as handle:
        observed_rows = list(csv.DictReader(handle, delimiter="\t"))
    if [row["row_id"] for row in observed_rows] != row_ids[held].tolist():
        raise SystemExit("source held rows differ")
    observed = np.asarray([[float(row[f"probability::{label}"]) for label in ROSTER] for row in observed_rows])
    error = float(np.max(np.abs(reproduced - observed)))
    if error > 2e-6: raise SystemExit(f"{head_id} reproduction error {error}")
    errors[head_id] = error
write_json_exclusive(output / "receipt.json", {
    "status": "passed", "outer_fold": 0, "screen_seed": 1103,
    "held_rows": len(held), "activation_vs_study_fold_mismatch_rows": 42584,
    "maximum_absolute_errors": errors, "evaluator_labels_read": False,
    "sealed_outcomes_read": False,
})
freeze_tree(output, {
    "artifact_class": "gse296875_transcriptformer_source_reproduction_preflight",
    "heads": ["linear", "two_layer_mlp"], "outer_fold": 0, "screen_seed": 1103,
    "held_rows": len(held), "evaluator_labels_read": False,
    "sealed_outcomes_read": False, "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
