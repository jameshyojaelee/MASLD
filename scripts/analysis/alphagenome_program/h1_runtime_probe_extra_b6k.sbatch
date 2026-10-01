#!/bin/bash
#SBATCH --job-name=agp_agruntime_extrab6k
#SBATCH --partition=gpu
#SBATCH --gres=gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_agruntime_extrab6k_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_agruntime_extrab6k_%j.err

# Is 1,048,576 bp the ceiling because of the card or because of the model? Two lengths outside the
# five the SDK names: 786,432 (a multiple of 2,048 between two supported lengths) and 2,097,152
# (twice the largest supported length). Whatever happens is recorded, pass or fail.

set -uo pipefail
export PYTHONNOUSERSITE=1
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_CPP_MIN_LOG_LEVEL=3
export HF_HUB_OFFLINE=1

MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC=$PROJ/scripts/analysis/alphagenome_program
EXEC=$PROJ/Analysis/MASLD_Model_Benchmark/executions/alphagenome-runtime-probe-20260915T104500Z
mkdir -p "$EXEC/tables/b6k_extra"

echo "=== node $(hostname) job $SLURM_JOB_ID ==="
for L in 2097152; do
  echo "=== bench length $L (outside the SDK's named set) ==="
  $MM run -n alphagenome_local python "$SRC/h1_bench_length.py" \
      --length "$L" --repeats 2 \
      --out "$EXEC/tables/b6k_extra/bench_${L}.json"
  echo "  exit=$?"
done
echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
