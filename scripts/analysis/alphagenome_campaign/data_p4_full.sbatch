#!/usr/bin/env bash
#SBATCH --job-name=masld_p4_full_reduce
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --dependency=afterok:21611740
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/p4-full-%j.log
set -euo pipefail
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
campaign_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
campaign_python=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
campaign_out="$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/p4-full-$SLURM_JOB_ID"
cd "$campaign_root/scripts/analysis/alphagenome_campaign"
"$campaign_python" data_checks.py
"$campaign_python" data_p4_full.py --out "$campaign_out" --workers 4
"$campaign_python" data_p4_families.py --features "$campaign_out" --out "$campaign_out/families"
