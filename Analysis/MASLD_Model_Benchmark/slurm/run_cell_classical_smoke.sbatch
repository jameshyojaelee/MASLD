#!/usr/bin/env bash
#SBATCH --job-name=masld-cellsmoke
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export LC_ALL=C

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
module load python/3.11.5-GCCcore-13.2.0

CANDIDATE="${MASLD_CANDIDATE:?candidate path required}"
OUT="executions/classical-smoke/${SLURM_JOB_ID}"
mkdir -p "${OUT}"

# Five sequential runs in one job, matching the cpu_baseline_smoke profile
# (1 CPU / 16 GB each).  Sequential rather than an array: 5 short CPU runs do
# not justify an array, and arrays are disabled for this campaign.
python - "$CANDIDATE" <<'PY' > "${OUT}/run_ids.txt"
import json, sys
plan = json.load(open(f"{sys.argv[1]}/plan.json"))
for run in plan["runs"]:
    print(run["run_id"], run["model_id"])
PY

while read -r RUN_ID MODEL_ID; do
  echo "=== ${MODEL_ID} ${RUN_ID} ==="
  PYTHONPATH=src python -m masld_bench.cli run execute \
    --candidate "${CANDIDATE}" \
    --run-id "${RUN_ID}" \
    --output-root "${OUT}"
done < "${OUT}/run_ids.txt"

echo "ALL RUNS COMPLETE"
