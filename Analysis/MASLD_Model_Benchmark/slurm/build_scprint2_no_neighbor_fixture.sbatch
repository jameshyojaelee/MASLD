#!/usr/bin/env bash
#SBATCH --job-name=model-data-302
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
runtime="${root}/executions/environments/scprint2-native-torch2.8.0-v2-21073597/source/.venv"
fixture="${root}/executions/transcriptformer-activation-preflight-21070971/fixture"
target="${root}/executions/scprint2-no-neighbor-fixture-${SLURM_JOB_ID}"
stage="${target}.staging"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite scPRINT-2 fixture: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export LD_LIBRARY_PATH="/gpfs/commons/home/jameslee/micromamba/envs/scprint2x/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

"${runtime}/bin/python" "${root}/scripts/build_scprint2_no_neighbor_fixture.py" \
  --source "${fixture}/unlabeled_counts.h5ad" \
  --rows "${fixture}/row_contract.tsv" \
  --checkpoint "${root}/executions/scprint2-acquisition-21066044/18hebyht-final-small.ckpt" \
  --output "${stage}/fixture"

source_sha256="$(sha256sum "${root}/scripts/build_scprint2_no_neighbor_fixture.py" | cut -d' ' -f1)"
fixture_manifest_sha256="$(sha256sum "${stage}/fixture/ARTIFACTS.json" | cut -d' ' -f1)"
python3 - "${stage}/ARTIFACTS.json" "${source_sha256}" "${fixture_manifest_sha256}" <<'PY'
import json
from pathlib import Path
import sys

output, source_sha, fixture_sha = sys.argv[1:]
payload = {
    "evaluation_labels_read": False,
    "fixture_manifest_sha256": fixture_sha,
    "fixture_source_sha256": source_sha,
    "neighbor_tensors_allowed": False,
    "schema_version": "masld-bench-scprint2-no-neighbor-fixture-artifact-v1",
    "sealed_outcomes_read": False,
    "status": "pass",
}
Path(output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
(
  cd "${stage}"
  sha256sum ARTIFACTS.json fixture/ARTIFACTS.json fixture/aligned_counts.h5ad fixture/fixture_receipt.json > SHA256SUMS
)
chmod -R a-w "${stage}"
mv "${stage}" "${target}"
sha256sum "${target}/ARTIFACTS.json"
