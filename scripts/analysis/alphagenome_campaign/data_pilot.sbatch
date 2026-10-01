#!/usr/bin/env bash
#SBATCH --job-name=masld_data_week1
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/data-%j.log
set -euo pipefail
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
campaign_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
campaign_python=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
campaign_out="$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data"
cd "$campaign_root/scripts/analysis/alphagenome_campaign"
"$campaign_python" data_checks.py
"$campaign_python" data_endpoints.py --out "$campaign_out/endpoints-$SLURM_JOB_ID"
"$campaign_python" data_p4.py --out "$campaign_out/p4-$SLURM_JOB_ID" --max-regions 100
"$campaign_python" -m pip freeze > "$campaign_out/environment-$SLURM_JOB_ID.txt"
