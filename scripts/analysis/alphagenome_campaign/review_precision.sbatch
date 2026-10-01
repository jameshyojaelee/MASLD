#!/bin/bash
#SBATCH --job-name=masld_ag_precision
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --dependency=afterok:21775091
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_precision_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_precision_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
precision_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
precision_base=$precision_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
precision_python=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$precision_python" "$precision_root/scripts/analysis/alphagenome_campaign/review_precision.py" \
  --final "$precision_base/model/final_21775091" \
  --out "$precision_base/review/precision-$SLURM_JOB_ID"
