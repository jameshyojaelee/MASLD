#!/bin/bash
#SBATCH --job-name=masld_adapt_figure
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_adapt_figure_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_adapt_figure_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
campaign_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
campaign_base=$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
 "$campaign_root/scripts/analysis/alphagenome_campaign/coordinator_adaptation_figure.py" \
 --source "$campaign_base/review/adaptation-precision-21773558/ten_paired_contrasts.tsv" \
 --check "$campaign_base/coordinator-adaptation-21773575/checks.json" \
 --out "$campaign_root/figures/candidates/regulatory-adaptation-week1-20260915-v2"
