#!/bin/bash
#SBATCH --job-name=masld_archive_integrity
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_archive_integrity_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_archive_integrity_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/data_p4_archive_integrity.py" \
  --run "$PROJ/GWAS/finemapping/results/alphagenome_atlas/p4-saturation-20260909T191917Z" \
  --out "$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/archive-integrity-${SLURM_JOB_ID}" \
  --workers 16
