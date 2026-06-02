#!/bin/bash
#SBATCH --job-name=celltype_concordance
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs/celltype_concordance_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs/celltype_concordance_%j.err

set -euo pipefail

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/scripts"
LOG_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs"
mkdir -p "${LOG_DIR}"

echo "=== Cell-type concordance (Gap #17) ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript "${SCRIPT_DIR}/06_celltype_concordance.R"

echo ""
echo "Finished: $(date)"
