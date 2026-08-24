#!/usr/bin/env bash
#SBATCH --job-name=masld-seq-screen-gate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
INPUT_ROSTER="${ROOT}/config/campaigns/sequence_native_diagonal_input_matrix_20260824.tsv"
MODEL_ROSTER="${ROOT}/config/campaigns/sequence_native_hepatocyte_screen_matrix_20260824.tsv"
TARGET="${ROOT}/executions/sequence-hepatocyte-screen-gate-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"

scripts=(
  "${ROOT}/slurm/build_chrombpnet_hepatocyte_donor0_genomic0_inputs.sbatch"
  "${ROOT}/scripts/sequence_control_production_driver.sh"
  "${ROOT}/slurm/sequence_control_train_cnn_hepatocyte_donor0_genomic0.sbatch"
  "${ROOT}/slurm/sequence_control_train_transformer_hepatocyte_donor0_genomic0.sbatch"
  "${ROOT}/slurm/bpnet_train_predict_hepatocyte_donor0_genomic0.sbatch"
  "${ROOT}/slurm/train_chrombpnet_hepatocyte_donor0_genomic0_seed20260824.sbatch"
)
for script in "${scripts[@]}"; do bash -n "${script}"; done

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${INPUT_ROSTER}" "${MODEL_ROSTER}" \
  > "${STAGE}/validation/gate.json" <<'PY'
from pathlib import Path
import csv
import json
import sys

input_roster, model_roster = map(Path, sys.argv[1:])
with input_roster.open(encoding="utf-8", newline="") as handle:
    inputs = list(csv.DictReader(handle, delimiter="\t"))
with model_roster.open(encoding="utf-8", newline="") as handle:
    models = list(csv.DictReader(handle, delimiter="\t"))
input_keys = {(row["lineage_id"], int(row["outer_fold"])) for row in inputs}
if len(inputs) != 24 or len(input_keys) != 24:
    raise SystemExit("diagonal input roster differs")
allowed_models = {
    "sequence_cnn_control", "sequence_transformer_control", "bpnet",
    "chrombpnet", "scbasset",
}
allowed_seeds = {20260824, 20260825, 20260826}
model_keys = {
    (row["model_id"], row["lineage_id"], int(row["outer_fold"]), int(row["seed"]))
    for row in models
}
if len(models) != 75 or len(model_keys) != 75:
    raise SystemExit("hepatocyte model roster cardinality differs")
for row in models:
    fold = int(row["outer_fold"])
    seed = int(row["seed"])
    if (
        row["model_id"] not in allowed_models
        or row["lineage_id"] != "hepatocyte"
        or fold not in range(5)
        or seed not in allowed_seeds
        or row["split_id"] != f"donor{fold}_genomic{fold}"
    ):
        raise SystemExit("hepatocyte model roster identity differs")
existing = [row for row in models if row["disposition"].startswith("completed_existing_")]
submitted = [row for row in models if row["disposition"] == "submit"]
if len(existing) != 4 or len(submitted) != 71:
    raise SystemExit("hepatocyte model roster disposition differs")
print(json.dumps({
    "schema_version": "masld-bench-sequence-hepatocyte-screen-gate-v1",
    "status": "pass",
    "input_rows": len(inputs),
    "model_rows": len(models),
    "existing_rows": len(existing),
    "submission_rows": len(submitted),
    "models": sorted(allowed_models),
    "folds": list(range(5)),
    "seeds": sorted(allowed_seeds),
    "frozen_upstream_payloads_not_rehashed": True,
}, sort_keys=True, indent=2))
PY
sha256sum "${INPUT_ROSTER}" "${MODEL_ROSTER}" "${scripts[@]}" \
  "${ROOT}/slurm/test_sequence_hepatocyte_screen_matrix.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/gate.json").read_text())
if receipt["status"] != "pass" or receipt["submission_rows"] != 71:
    raise SystemExit("screen gate receipt differs")
freeze_tree(root, {"artifact_class": "sequence_hepatocyte_screen_gate", **receipt})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
