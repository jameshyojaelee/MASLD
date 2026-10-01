#!/usr/bin/env bash
#SBATCH --job-name=masld_phase_census
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/read-phase-%j.log
set -euo pipefail
umask 077
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
campaign_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
crosswalk_out="$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/currin-donor-crosswalk-$SLURM_JOB_ID"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript \
  "$campaign_root/scripts/analysis/alphagenome_campaign/data_currin_donor_crosswalk.R" \
  "$campaign_root" "$crosswalk_out"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "$campaign_root/scripts/analysis/alphagenome_campaign/data_read_phase.py" \
  --out "$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/read-phase-$SLURM_JOB_ID"
