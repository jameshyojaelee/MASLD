#!/usr/bin/env bash
#SBATCH --job-name=model-check-303
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
runtime_artifact="${root}/executions/environments/scprint2-native-torch2.8.0-v2-21073597"
runtime="${runtime_artifact}/source/.venv"
source="${runtime_artifact}/source"
checkpoint="${root}/executions/scprint2-acquisition-21066044/18hebyht-final-small.ckpt"
target="${root}/executions/scprint2-native-runtime-preflight-v2-${SLURM_JOB_ID}"
stage="${target}.staging"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite scPRINT-2 preflight v2: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"

(
  cd "${runtime_artifact}"
  sha256sum -c SHA256SUMS
)
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

"${runtime}/bin/python" "${root}/scripts/preflight_scprint2_native_runtime.py" \
  --checkpoint "${checkpoint}" \
  --source "${source}" \
  --output "${stage}/runtime_preflight.json"

source_sha256="$(sha256sum "${root}/scripts/preflight_scprint2_native_runtime.py" | cut -d' ' -f1)"
report_sha256="$(sha256sum "${stage}/runtime_preflight.json" | cut -d' ' -f1)"
status="$("${runtime}/bin/python" -c \
  'import json,sys; print(json.load(open(sys.argv[1]))["status"])' \
  "${stage}/runtime_preflight.json")"
python3 - "${stage}/ARTIFACTS.json" "${source_sha256}" "${report_sha256}" "${status}" <<'PY'
import json
from pathlib import Path
import sys

output, source_sha, report_sha, status = sys.argv[1:]
payload = {
    "outcomes_or_labels_read": False,
    "preflight_source_sha256": source_sha,
    "runtime_artifact_sha256": "4474a0d069d09a9a805e7d6e33936af95eeab9fb72d6e3c131e7e7733350842d",
    "runtime_preflight_sha256": report_sha,
    "schema_version": "masld-bench-scprint2-native-runtime-preflight-artifact-v2",
    "status": status,
}
Path(output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
(
  cd "${stage}"
  sha256sum ARTIFACTS.json runtime_preflight.json > SHA256SUMS
)
chmod -R a-w "${stage}"
mv "${stage}" "${target}"
sha256sum "${target}/ARTIFACTS.json"
