#!/bin/bash
#SBATCH --job-name=limma
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/14b_lvqw_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/14b_lvqw_%j.err
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPTS="${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

# set strict mode AFTER activation (conda activate.d scripts reference unbound vars)
set -euo pipefail

export MASLD_PROJECT_ROOT="${BASE}"

echo "=== 14b/14d LVQW granular staging harmonization: $(date) ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo ""

Rscript "${SCRIPTS}/14b_nas_fibrosis_stage_lvqw.R" 2>&1

echo ""
echo "=== Done: $(date) ==="
