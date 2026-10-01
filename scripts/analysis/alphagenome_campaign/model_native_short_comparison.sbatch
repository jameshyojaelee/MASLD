#!/bin/bash
#SBATCH --job-name=masld_native_short_compare
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --dependency=afterok:21775189
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_native_short_compare_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_native_short_compare_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
short_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
short_base=$short_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "$short_root/scripts/analysis/alphagenome_campaign/model_native_short_comparison.py" \
  --native "$short_base/model/native" --lengths "$short_base/model/native_lengths_21775189" \
  --adaptation "$short_base/model/sweep_summary_21773124/matched_predictions.tsv.gz" \
  --out "$short_base/model/native_short_compare_$SLURM_JOB_ID"
