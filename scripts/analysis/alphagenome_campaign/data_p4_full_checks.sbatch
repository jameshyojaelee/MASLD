#!/usr/bin/env bash
#SBATCH --job-name=masld_p4_full_checks
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/p4-checks-%j.log
set -euo pipefail
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_campaign
/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python data_p4_full_checks.py
