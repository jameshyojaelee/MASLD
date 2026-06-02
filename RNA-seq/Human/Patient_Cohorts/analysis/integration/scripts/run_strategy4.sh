#!/bin/bash
#SBATCH --job-name=liver_strategy4_sex
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/strategy4_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/strategy4_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlee@nygenome.org

set -euo pipefail

# Activate R environment
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration"

# Ensure log directory exists
mkdir -p "${SCRIPT_DIR}/logs"

echo "============================================================"
echo "  STRATEGY 4: SEX-STRATIFIED ANALYSIS"
echo "  Date: $(date)"
echo "  Job ID: ${SLURM_JOB_ID:-local}"
echo "  CPUs: $SLURM_CPUS_PER_TASK | RAM: 200G | Partition: cpu"
echo "============================================================"

# Use base dream model (Task 1 model selection pending)
export DREAM_MODEL=base

echo -e "\n=== 26: Sex-Stratified Analysis ==="
Rscript "${SCRIPT_DIR}/scripts/26_sex_stratified_analysis.R"

echo "============================================================"
echo "  STRATEGY 4 COMPLETE"
echo "  Finished: $(date)"
echo "============================================================"
