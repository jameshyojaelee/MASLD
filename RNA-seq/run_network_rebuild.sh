#!/bin/bash
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --job-name=net_rebuild
#SBATCH --output=logs/net_rebuild_%j.out
#SBATCH --error=logs/net_rebuild_%j.err
# ===========================================================================
# Launcher for Top-5 RP10 deliverable: 293 + 297 sequential
# ===========================================================================

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Activate spatial env (has igraph + leidenalg)
export MAMBA_EXE='/gpfs/commons/home/jameslee/.local/bin/micromamba'
export MAMBA_ROOT_PREFIX='/gpfs/commons/home/jameslee/micromamba'
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate spatial

echo "=== $(date) START ==="
echo "Host: $(hostname)"
echo "JobID: ${SLURM_JOB_ID:-(no slurm)}"
echo "Python: $(which python)"
echo "PYTHONPATH: ${PYTHONPATH:-unset}"
echo

mkdir -p RNA-seq/logs

# ---- 293: continuous-stage Leiden ----
echo "=== $(date) 293 START ==="
python RNA-seq/293_continuous_stage_leiden.py
echo "=== $(date) 293 DONE ==="

# ---- 297: Maslov-Sneppen Q-null ----
# Use 50 rewires for first run (manuscript can re-run with N_REWIRES=200 later).
export N_REWIRES=${N_REWIRES:-50}
export REWIRE_PASSES=${REWIRE_PASSES:-5}
echo "=== $(date) 297 START (N_REWIRES=${N_REWIRES}, REWIRE_PASSES=${REWIRE_PASSES}) ==="
python RNA-seq/297_modularity_null.py
echo "=== $(date) 297 DONE ==="

echo "=== $(date) ALL DONE ==="
