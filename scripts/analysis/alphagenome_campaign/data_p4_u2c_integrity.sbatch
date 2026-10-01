#!/bin/bash
#SBATCH --job-name=masld_u2c_integrity
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_u2c_integrity_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_u2c_integrity_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/data_p4_archive_integrity.py" \
  --run "$PROJ/GWAS/finemapping/results/alphagenome_atlas/p4-u2c-genetic-20260917T151313Z" \
  --out "$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/u2c-integrity-${SLURM_JOB_ID}" \
  --workers 8
