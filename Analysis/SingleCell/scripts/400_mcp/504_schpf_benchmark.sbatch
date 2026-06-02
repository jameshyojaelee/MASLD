#!/bin/bash
#SBATCH --job-name=alpha4_schpf
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=300G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/504_schpf_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/504_schpf_%j.err
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp
source "$HOME/.bashrc"
micromamba activate rapids_singlecell

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

# Install scHPF if missing (pure python, no compile)
if ! python -c "import schpf" 2>/dev/null; then
  echo "[504.sbatch] Installing scHPF from GitHub..."
  pip install --no-deps --quiet 'git+https://github.com/simslab/scHPF.git@master' || {
    echo "[504.sbatch] GitHub install failed; fallback to PyPI 'schpf'";
    pip install --quiet schpf || true
  }
fi

: "${SCHPF_K_LIST:=10,16}"
: "${SCHPF_N_TRIALS:=3}"
: "${SCHPF_N_CELLS:=0}"
export SCHPF_K_LIST SCHPF_N_TRIALS SCHPF_N_CELLS

python -u 504_schpf_benchmark.py
