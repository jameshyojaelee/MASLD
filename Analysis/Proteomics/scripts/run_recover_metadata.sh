#!/bin/bash
#SBATCH --job-name=proteomics
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Proteomics/logs/recover_metadata_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Proteomics/logs/recover_metadata_%j.err

# Recover disease metadata for proteomics datasets (GSE276114 + PXD052937)
#
# Usage:
#   sbatch Analysis/Proteomics/scripts/run_recover_metadata.sh

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${BASE}/Analysis/Proteomics/scripts/recover_proteomics_metadata.R"

mkdir -p "${BASE}/Analysis/Proteomics/logs"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "============================================================"
echo "  Recover Proteomics Metadata"
echo "  Started: $(date)"
echo "  SLURM_JOB_ID: $SLURM_JOB_ID"
echo "  RAM: 32G | CPUs: $SLURM_CPUS_PER_TASK"
echo "============================================================"
echo ""

export MASLD_PROJECT_ROOT="$BASE"
Rscript "${SCRIPT}" 2>&1

echo ""
echo "============================================================"
echo "  COMPLETE: $(date)"
echo "  Results: ${BASE}/Analysis/Proteomics/results/"
echo "============================================================"
