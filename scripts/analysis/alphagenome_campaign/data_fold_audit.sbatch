#!/usr/bin/env bash
#SBATCH --job-name=masld_reporter_folds
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/folds-%j.log
set -euo pipefail
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
campaign_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$campaign_root/scripts/analysis/alphagenome_campaign"
/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python data_fold_audit.py --endpoints "$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/endpoints-21771964" --out "$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/folds-$SLURM_JOB_ID"
/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python data_checks.py
/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python data_endpoints.py --out "$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/endpoints-revised-$SLURM_JOB_ID"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python review_axes_catalog.py --catalog "$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/catalog-21772153" --b2 "$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/review/b2-21771990" --out "$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/review/axes-catalog-$SLURM_JOB_ID"
