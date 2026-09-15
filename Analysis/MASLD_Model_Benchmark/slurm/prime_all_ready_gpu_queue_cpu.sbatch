#!/usr/bin/env bash
#SBATCH --job-name=model-queue-prime
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
cd "$ROOT"
module load python/3.11.5-GCCcore-13.2.0

PYTHONPATH="src:." python -m unittest -v tests.unit.test_prime_all_ready_gpu_queue tests.unit.test_gpu_bundle_dispatcher
PYTHONPATH="src:." python scripts/prime_all_ready_gpu_queue.py \
  --root "$ROOT" \
  --queue "$ROOT/config/campaigns/gpu_bundle_queue" \
  --state "$ROOT/executions/gpu-bundle-dispatch-state" \
  --config "$ROOT/config/campaigns/gpu_queue_prime_all_ready_20260825.json" \
  --user "$USER" \
  --output "$ROOT/executions/gpu-queue-prime-all-ready-${SLURM_JOB_ID}"
