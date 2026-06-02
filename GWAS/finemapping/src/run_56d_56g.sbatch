#!/bin/bash
#SBATCH --job-name=motif_quality_regulon_audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/56d_56g_%j.out
#SBATCH --error=logs/56d_56g_%j.err

set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate motifbreakr  # needs MotifDb, TFBSTools, BSgenome, httr, jsonlite

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs results/gwas_atac

echo "================================================================"
echo "Job:     $SLURM_JOB_ID"
echo "Host:    $(hostname)"
echo "Started: $(date)"
echo "================================================================"

ts=$(date +%s)
echo ">>> Running 56d_motif_quality_filter.R"
Rscript src/56d_motif_quality_filter.R 2>&1 | tee logs/56d_${ts}.log

echo ""
echo ">>> Running 56g_regulon_threshold_audit.R"
Rscript src/56g_regulon_threshold_audit.R 2>&1 | tee logs/56g_${ts}.log

echo ""
echo "================================================================"
echo "Finished: $(date)"
echo "================================================================"
