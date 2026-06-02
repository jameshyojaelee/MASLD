#!/bin/bash
#SBATCH --job-name=fgsea
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/08b_unified_gsea_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/08b_unified_gsea_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00

# 08b_unified_pathway_gsea: fgsea across ALL contrasts x ALL collections
# Supersedes 08, 14c, 14f, 133, 115b
# Submit: sbatch run_08b_unified_gsea.sh
# Or with dependency: sbatch --dependency=afterok:<PHASE0_JOBID> run_08b_unified_gsea.sh

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

echo "=== 08b: Unified Pathway GSEA ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Started: $(date)"
echo ""

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

Rscript scripts/08b_unified_pathway_gsea.R 2>&1

echo ""
echo "=== Completed: $(date) ==="
