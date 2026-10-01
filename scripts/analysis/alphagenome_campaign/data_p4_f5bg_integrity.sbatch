#!/bin/bash
#SBATCH --job-name=masld_f5bg_integrity
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_f5bg_integrity_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_f5bg_integrity_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/data_p4_archive_integrity.py" \
  --run "$PROJ/GWAS/finemapping/results/alphagenome_atlas/p4-f5-background-20260917T145358Z" \
  --out "$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/f5bg-integrity-${SLURM_JOB_ID}" \
  --workers 8
