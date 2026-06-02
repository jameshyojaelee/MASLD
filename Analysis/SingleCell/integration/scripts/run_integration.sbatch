#!/usr/bin/env bash
#SBATCH --job-name=scalesc_integration
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --output=Analysis/SingleCell/integration/logs/integration_%j.out
#SBATCH --error=Analysis/SingleCell/integration/logs/integration_%j.err

# ScaleSC single-cell integration pipeline
# Step 1: Convert Cell Ranger h5 → per-sample h5ad (CPU, ~10 min)
# Step 2: Run ScaleSC GPU integration (human, then mouse)

set -euo pipefail

echo "=== ScaleSC Integration Pipeline ==="
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Date: $(date)"
echo ""

# Activate environment first
eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

# Fix LD_LIBRARY_PATH for pip-installed NVIDIA libs (must be AFTER activation)
PIP_NVIDIA="${CONDA_PREFIX}/lib/python3.12/site-packages/nvidia"
if [ -d "${PIP_NVIDIA}" ]; then
    for pkg in cublas cusparse cusolver curand cufft nvjitlink cusparselt; do
        d="${PIP_NVIDIA}/${pkg}/lib"
        [ -d "$d" ] && export LD_LIBRARY_PATH="${d}:${LD_LIBRARY_PATH:-}"
    done
fi

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# --- Step 1: Prepare h5ad files (skip if already done) ---
H5AD_DIR="Analysis/SingleCell/integration/input_h5ad"
N_HUMAN=$(ls "${H5AD_DIR}/human/"*.h5ad 2>/dev/null | wc -l)
N_MOUSE=$(ls "${H5AD_DIR}/mouse/"*.h5ad 2>/dev/null | wc -l)
if [ "$N_HUMAN" -gt 0 ] && [ "$N_MOUSE" -gt 0 ]; then
    echo "=== Step 1: h5ad files already prepared (${N_HUMAN} human, ${N_MOUSE} mouse) ==="
else
    echo "=== Step 1: Converting Cell Ranger h5 → h5ad ==="
    python Analysis/SingleCell/integration/scripts/prepare_h5ad.py
fi

# --- Step 2: ScaleSC Integration (Human) ---
echo ""
echo "=== Step 2: ScaleSC Integration (Human) ==="
python Analysis/SingleCell/integration/scripts/run_scalesc_integration.py \
    --species human \
    --resolution 0.5 \
    --n-hvgs 4000 \
    --max-cell-batch 200000

# --- Step 3: ScaleSC Integration (Mouse) ---
echo ""
echo "=== Step 3: ScaleSC Integration (Mouse) ==="
python Analysis/SingleCell/integration/scripts/run_scalesc_integration.py \
    --species mouse \
    --resolution 0.3 \
    --n-hvgs 4000 \
    --max-cell-batch 200000

echo ""
echo "=== Integration Complete ==="
echo "End: $(date)"
