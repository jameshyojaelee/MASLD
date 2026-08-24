#!/usr/bin/env bash
#SBATCH --job-name=masld-seq-fold0-retry-gate
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
INPUT="${ROOT}/executions/chrombpnet-hepatocyte-donor0-genomic0-inputs-21064244"
TARGET="${ROOT}/executions/sequence-fold0-legacy-input-retry-gate-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"
scripts=(
  "${ROOT}/slurm/sequence_control_train_cnn_hepatocyte_donor0_genomic0.sbatch"
  "${ROOT}/slurm/sequence_control_train_transformer_hepatocyte_donor0_genomic0.sbatch"
  "${ROOT}/slurm/bpnet_train_predict_hepatocyte_donor0_genomic0.sbatch"
  "${ROOT}/slurm/train_chrombpnet_hepatocyte_donor0_genomic0_seed20260824.sbatch"
)
for script in "${scripts[@]}"; do
  bash -n "${script}"
  rg -q 'hepatocyte:0:21064244' "${script}"
  rg -q 'chrombpnet-hepatocyte-donor0-genomic0-inputs-21064244' "${script}"
done
printf '%s  %s\n' fc26ec558ea51f74212bea120dab59053e5f2856494407bb0a717e9eaecc0208 \
  "${INPUT}/ARTIFACTS.json" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${INPUT}" > "${STAGE}/validation/gate.json" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(root)
print(json.dumps({
    "schema_version": "masld-bench-sequence-fold0-legacy-input-retry-gate-v1",
    "status": "pass", "lineage_id": "hepatocyte", "split_id": "donor0_genomic0",
    "legacy_path_token": "donor0-genomic0", "input_job_id": 21064244,
    "failed_attempts": [21067067, 21067068, 21067069, 21067070],
    "failure_stage": "pretraining_input_path_resolution",
}, sort_keys=True, indent=2))
PY
sha256sum "${scripts[@]}" "${ROOT}/slurm/test_sequence_fold0_legacy_input_retry.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/gate.json").read_text())
freeze_tree(root, {"artifact_class": "sequence_fold0_legacy_input_retry_gate", **receipt})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
