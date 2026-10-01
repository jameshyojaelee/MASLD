#!/bin/bash
#SBATCH --job-name=masld_f5background_reduce
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --dependency=afterok:21783065
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_f5background_reduce_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_f5background_reduce_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
cd "$PROJ/scripts/analysis/alphagenome_campaign"
"$PY" data_p4_extra_reduce.py --run "$PROJ/GWAS/finemapping/results/alphagenome_atlas/p4-f5-background-20260917T145358Z" \
  --out "$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/f5background-reduce-${SLURM_JOB_ID}" \
  --expect 18570 --workers 8
