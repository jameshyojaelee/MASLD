#!/bin/bash
#SBATCH --job-name=masld_f5_f6_selftest
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_f5_f6_selftest_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_f5_f6_selftest_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
cd "$PROJ/scripts/analysis/alphagenome_campaign"
"$PY" data_p4_f5_f6_families.py --self-test \
  --out "$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/f5-f6-selftest-${SLURM_JOB_ID}"
