#!/bin/bash
#SBATCH --job-name=masld_review_models
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --dependency=afterok:21772726
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_models_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_models_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
campaign_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
campaign_base=$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
campaign_python=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$campaign_python" "$campaign_root/scripts/analysis/alphagenome_campaign/review_model_results.py" \
  --lx2 "$campaign_base/model/lx2_21772727" --sweep "$campaign_base/model/sweep_21772726" \
  --donor "$campaign_base/review/b2-21771990" \
  --figures "$campaign_root/figures/candidates/regulatory-week1-20260915-v3" \
  --require-sweep-complete --out "$campaign_base/review/model-results-$SLURM_JOB_ID"
python3 "$campaign_root/scripts/analysis/alphagenome_campaign/catalog_checks.py"
