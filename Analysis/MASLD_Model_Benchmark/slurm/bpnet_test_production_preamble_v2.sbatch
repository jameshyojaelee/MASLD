#!/usr/bin/env bash
#SBATCH --job-name=masld-bpnet-preamble-v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=executions/masld-bpnet-preamble-v2-%j.out
#SBATCH --error=executions/masld-bpnet-preamble-v2-%j.err

set -euo pipefail
umask 027

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ACQUISITION="${ROOT}/executions/upstream-acquisitions/bpnet-v2.0.0-21064186"
TARGET="${ROOT}/executions/bpnet-production-preamble-v2-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite BPNet preamble regression: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

bash -n "${ROOT}/slurm/bpnet_train_predict_hepatocyte_donor0_genomic0.sbatch"
module load python/3.11.5-GCCcore-13.2.0
python "${ROOT}/scripts/bpnet_verify_acquisition.py" \
  --acquisition "${ACQUISITION}" \
  --output "${STAGE}/acquisition-verification.json"
python - "${STAGE}/acquisition-verification.json" <<'PY'
from pathlib import Path
import json
import sys

receipt = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
if receipt != {
    "status": "pass",
    "artifact_count": 301,
    "manifest_sha256": "8062574bfaa9067b0aca82ffa3eeae87188140ee5f5842399171ad4535cf88f6",
}:
    raise SystemExit("BPNet acquisition regression receipt differs")
PY

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

freeze_tree(
    Path(sys.argv[1]),
    metadata={
        "artifact_class": "bpnet_production_preamble_regression",
        "model_id": "bpnet",
        "status": "pass",
        "production_outcomes_exposed": False,
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'BPNet production preamble regression: %s\n' "${TARGET}"
