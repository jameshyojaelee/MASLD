#!/bin/bash
#SBATCH --job-name=ld_null_56b
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=../logs/56b_ld_null_%j.out
#SBATCH --error=../logs/56b_ld_null_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

# plink2 for MAF computation
module load PLINK/2.00a3.7-gfbf-2023a 2>/dev/null || module load PLINK/2.0a5.13 2>/dev/null || true

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs results/gwas_atac

echo "Running 56b_ld_aware_enrichment_null.R on $(hostname) at $(date)"
echo "SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK}"
which plink2 || true
which plink || true

Rscript src/56b_ld_aware_enrichment_null.R
echo "Done at $(date)"
