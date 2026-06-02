#!/usr/bin/env bash
#SBATCH --job-name=celltypist_annotation
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=12:00:00
#SBATCH --output=Analysis/SingleCell/integration/logs/annotation_%j.out
#SBATCH --error=Analysis/SingleCell/integration/logs/annotation_%j.err

set -euo pipefail

echo "=== CellTypist Annotation Pipeline ==="
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Date: $(date)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

# LD_LIBRARY_PATH fix (AFTER activation)
PIP_NVIDIA="${CONDA_PREFIX}/lib/python3.12/site-packages/nvidia"
if [ -d "${PIP_NVIDIA}" ]; then
    for pkg in cublas cusparse cusolver curand cufft nvjitlink cusparselt; do
        d="${PIP_NVIDIA}/${pkg}/lib"
        [ -d "$d" ] && export LD_LIBRARY_PATH="${d}:${LD_LIBRARY_PATH:-}"
    done
fi

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

python Analysis/SingleCell/integration/scripts/annotate_celltypist.py

echo ""
echo "=== Annotation Complete ==="
echo "End: $(date)"
