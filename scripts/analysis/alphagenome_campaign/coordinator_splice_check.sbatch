#!/bin/bash
#SBATCH --job-name=masld_splice_reconstruct
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_splice_reconstruct_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_splice_reconstruct_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
campaign_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
 "$campaign_root/scripts/analysis/alphagenome_campaign/coordinator_splice_check.py" \
 --out "$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/coordinator-splice-$SLURM_JOB_ID"
