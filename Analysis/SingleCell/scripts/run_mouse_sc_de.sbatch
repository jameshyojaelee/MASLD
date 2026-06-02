#!/bin/bash
#SBATCH --job-name=mouse_sc_de
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/mouse_sc_de_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/mouse_sc_de_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SC_SCRIPTS="${BASE}/Analysis/SingleCell/scripts"

echo "=== Mouse scRNA-seq Pseudobulk DE + Species Effect Size Normalization ==="
echo "Job ID:   $SLURM_JOB_ID"
echo "Node:     $SLURMD_NODENAME"
echo "CPUs:     $SLURM_CPUS_PER_TASK"
echo "Start:    $(date)"

mkdir -p "${SC_SCRIPTS}/logs"

# ---------------------------------------------------------------------------
# Step 1: Generate pseudobulk CSVs from per-sample h5ad files
#   Uses rapids_singlecell env (has anndata 0.12 compatible with the h5ad format)
#   but only runs CPU code (no GPU ops needed)
# ---------------------------------------------------------------------------
echo ""
echo "--- Step 1: Mouse Pseudobulk Generation (Python) ---"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rapids_singlecell

# Set LD_LIBRARY_PATH for pip-installed NVIDIA libs
NVIDIA_LIB="/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/lib/python3.12/site-packages/nvidia"
if [ -d "$NVIDIA_LIB" ]; then
  for lib_dir in "$NVIDIA_LIB"/*/lib; do
    [ -d "$lib_dir" ] && export LD_LIBRARY_PATH="${lib_dir}:${LD_LIBRARY_PATH:-}"
  done
fi

python "${SC_SCRIPTS}/mouse_pseudobulk_de.py"
echo "Pseudobulk generation complete."

# ---------------------------------------------------------------------------
# Step 2: Run limma-voom DE (R) on the pseudobulk counts
# ---------------------------------------------------------------------------
echo ""
echo "--- Step 2: Mouse Pseudobulk DE (R/limma-voom) ---"

micromamba activate rnaseq

Rscript "${SC_SCRIPTS}/mouse_pseudobulk_de.R"
echo "Mouse pseudobulk DE complete."

# ---------------------------------------------------------------------------
# Step 3: Species effect size normalization + cross-species scatter
# ---------------------------------------------------------------------------
echo ""
echo "--- Step 3: Species Effect Size Normalization (R) ---"

Rscript "${SC_SCRIPTS}/species_effect_size_norm.R"
echo "Species normalization complete."

echo ""
echo "=== Done: $(date) ==="
echo "Outputs in: ${BASE}/Analysis/SingleCell/results_gpu_v2/mouse_sc/"
ls -lh "${BASE}/Analysis/SingleCell/results_gpu_v2/mouse_sc/" 2>/dev/null || true
