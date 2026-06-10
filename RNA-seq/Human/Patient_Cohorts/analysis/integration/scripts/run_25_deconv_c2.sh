#!/bin/bash
#SBATCH --job-name=limma
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/deconv_c2_%j.out
#SBATCH --error=logs/deconv_c2_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p logs
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
echo "[$(date)] Tier 4D — deconvolution attribution C2 recompute on $(hostname), CPUs=$SLURM_CPUS_PER_TASK"
Rscript 25_deconv_attribution_C2.R
echo "[$(date)] done"
