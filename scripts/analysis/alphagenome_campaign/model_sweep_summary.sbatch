#!/bin/bash
#SBATCH --job-name=masld_ag_sweep_summary
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --dependency=afterok:21772726
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_sweep_summary_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_sweep_summary_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_summarize.py" \
 --sweep "$BASE/sweep_21772726" --native "$BASE/native" \
 --out "$BASE/sweep_summary_${SLURM_JOB_ID}"
