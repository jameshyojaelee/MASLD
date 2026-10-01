#!/bin/bash
#SBATCH --job-name=masld_f5_background_universe
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_f5_background_universe_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_f5_background_universe_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/data_p4_f5_background_universe.py" \
  --out "$PROJ/GWAS/finemapping/results/alphagenome_atlas/p4-f5-background-20260917T145358Z"
