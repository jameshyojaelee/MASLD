#!/bin/bash
#SBATCH --job-name=agp_ad1_smoke
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:40:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad1_smoke_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad1_smoke_%j.err
set -uo pipefail
export PYTHONNOUSERSITE=1
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_CPP_MIN_LOG_LEVEL=3
export HF_HUB_OFFLINE=1
MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC=$PROJ/scripts/analysis/alphagenome_program
RUN=$PROJ/GWAS/finemapping/results/alphagenome_program/ad1-gate-20260915T114615Z
echo "=== node $(hostname) job $SLURM_JOB_ID ==="
$MM run -n alphagenome_local python "$SRC/i1_01_score_local.py" \
    --variants "$RUN/inputs/sweep_variants.tsv" --length 16384 --limit 3 \
    --out "/scratch/claude-91825/smoke_16384.tsv" --flush-every 1
echo "exit=$?"
cat /scratch/claude-91825/smoke_16384.tsv 2>/dev/null | head -5
