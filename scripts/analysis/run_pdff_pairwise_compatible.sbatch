#!/bin/bash
#SBATCH --job-name=pdff_rg_pairs
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=90:00:00
#SBATCH --output=logs/pdff_rg_pairs_%j.out
#SBATCH --error=logs/pdff_rg_pairs_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
PDFF_RUN_ROOT="GWAS/ldsc/rg_pdff_compatible_$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}"
mkdir -p "$PDFF_RUN_ROOT"
echo "$PDFF_RUN_ROOT"
cp scripts/analysis/run_pdff_pairwise_compatible.sbatch "$PDFF_RUN_ROOT/"
cp scripts/analysis/pdff_pairwise_compatible.py "$PDFF_RUN_ROOT/"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python scripts/analysis/pdff_pairwise_compatible.py "$PDFF_RUN_ROOT"
