#!/usr/bin/env bash
#SBATCH --job-name=masld_target_feasibility
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/target-admission-v1/followup-%j.log
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
campaign_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
target_data="$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "$campaign_root/scripts/analysis/alphagenome_campaign/target_benchmark_admission_followup.py" \
  --admission "$target_data/target-admission-v1/admission-21775202" \
  --f5-producer "$target_data/p4-background-21775203" \
  --activity-manifest "$target_data/target-admission-v1/activity_source_manifest.json" \
  --out "$target_data/target-admission-v1/followup-$SLURM_JOB_ID"
