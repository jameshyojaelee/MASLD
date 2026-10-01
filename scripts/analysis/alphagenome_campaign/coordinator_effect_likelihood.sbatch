#!/bin/bash
#SBATCH --job-name=masld_effect_likelihood
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_effect_likelihood_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_effect_likelihood_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false HF_HUB_OFFLINE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
campaign_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
campaign_python=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python
"$campaign_python" "$campaign_root/scripts/analysis/alphagenome_campaign/coordinator_effect_likelihood.py" \
 --out "$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model/effect_likelihood_$SLURM_JOB_ID" --max-hours 7.8
