#!/bin/bash
# 503b_mofa_pseudobulk.sbatch — pure-Python MOFA+ closure for showdown.
#
# Why bigmem + interactive: nslab QOS sat at 7TB; bigmem/interactive bypasses
# memory cap (CLAUDE.md HPC protocol). 200G is generous; smoke test fit easily
# under 16G.
#
# Track A reviewer-audit B1: replaces 4 failed R/MOFA2 retries (15495684,
# 15530928, 15539666 + 1 prior). See PROGRESS.md "MOFA+ closure" section.

#SBATCH --job-name=mofa_503b
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/503b_mofa_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/503b_mofa_%j.err

# Conda binutils unbound-var bug workaround — set +u around env activation.
set +u
ENV_PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
set -u

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT="$ROOT/Analysis/SingleCell/scripts/400_mcp/503b_mofa_pseudobulk_python.py"

mkdir -p "$ROOT/Analysis/SingleCell/scripts/400_mcp/logs"

echo "[$(date)] starting MOFA+ k=10 + k=16 full fit"
"$ENV_PY" "$SCRIPT" --k 10 --k 16 --seed 42

echo "[$(date)] DONE."
