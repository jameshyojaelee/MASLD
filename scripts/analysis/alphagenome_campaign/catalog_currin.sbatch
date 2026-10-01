#!/bin/bash
#SBATCH --job-name=masld_catalog_currin
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_catalog_currin_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_catalog_currin_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
campaign_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
campaign_python=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$campaign_python" "$campaign_root/scripts/analysis/alphagenome_campaign/catalog_currin_checks.py"
"$campaign_python" "$campaign_root/scripts/analysis/alphagenome_campaign/catalog_currin.py" \
 --out "$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/catalog-currin-$SLURM_JOB_ID"
