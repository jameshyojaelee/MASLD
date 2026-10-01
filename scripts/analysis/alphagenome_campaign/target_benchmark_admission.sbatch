#!/usr/bin/env bash
#SBATCH --job-name=masld_target_admission
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/target-admission-v1/admission-%j.log
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
campaign_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
"/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python" \
  "$campaign_root/scripts/analysis/alphagenome_campaign/target_benchmark_admission.py" \
  --manifest "$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/target-admission-v1/source_manifest.json" \
  --out "$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/target-admission-v1/admission-$SLURM_JOB_ID"
