#!/bin/bash
#SBATCH --job-name=agp_agruntime_l40s
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_agruntime_l40s_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_agruntime_l40s_%j.err

set -uo pipefail

export PYTHONNOUSERSITE=1
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_CPP_MIN_LOG_LEVEL=3
export HF_HUB_OFFLINE=1

MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
ENV=alphagenome_local
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC=$PROJ/scripts/analysis/alphagenome_program
EXEC=$PROJ/Analysis/MASLD_Model_Benchmark/executions/alphagenome-runtime-probe-20260915T104500Z
GPUTAG=l40s

mkdir -p "$EXEC/tables/$GPUTAG" "$EXEC/env" "$EXEC/logs"

echo "=== node $(hostname) job $SLURM_JOB_ID ==="
nvidia-smi | tee "$EXEC/env/nvidia_smi_${GPUTAG}.txt"
nvidia-smi --query-gpu=name,memory.total,driver_version,compute_cap --format=csv \
  | tee "$EXEC/env/gpu_${GPUTAG}.csv"

$MM run -n $ENV pip freeze > "$EXEC/env/pip_freeze_${GPUTAG}.txt" 2>/dev/null
$MM run -n $ENV python - <<'PY' > "$EXEC/env/versions_${GPUTAG}.txt" 2>&1
import sys, jax, jaxlib, numpy, platform
import alphagenome, alphagenome_research
print("python", sys.version.replace("\n", " "))
print("platform", platform.platform())
print("jax", jax.__version__)
print("jaxlib", jaxlib.__version__)
print("numpy", numpy.__version__)
print("alphagenome", getattr(alphagenome, "__version__", "n/a"))
print("alphagenome_research", getattr(alphagenome_research, "__version__", "n/a"))
print("jax devices", jax.devices())
print("backend", jax.default_backend())
PY
cat "$EXEC/env/versions_${GPUTAG}.txt"

# ---- stage 1: one process per input length, so each length gets its own peak-memory high-water mark
for L in 2048 16384 131072 524288 1048576; do
  echo "=== bench length $L ==="
  $MM run -n $ENV python "$SRC/h1_bench_length.py" \
      --length "$L" --repeats 3 \
      --out "$EXEC/tables/$GPUTAG/bench_${L}.json"
  echo "  exit=$?"
done

# ---- stage 2: variant-effect probes at 1 Mb
echo "=== probes ==="
$MM run -n $ENV python "$SRC/h1_variant_probes.py" \
    --length 1048576 --out "$EXEC/tables/$GPUTAG/variant_probes.json"
echo "  exit=$?"

# ---- stage 3: HSD17B13 anchor against the hosted value
echo "=== anchor ==="
$MM run -n $ENV python "$SRC/h1_anchor_hsd17b13.py" \
    --length 1048576 --out "$EXEC/tables/$GPUTAG/anchor_hsd17b13.json"
echo "  exit=$?"

echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
