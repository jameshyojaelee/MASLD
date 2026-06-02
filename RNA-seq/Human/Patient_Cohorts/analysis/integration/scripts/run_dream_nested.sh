#!/bin/bash
#SBATCH --job-name=B2_dream_nested
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/B2_dream_nested_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/B2_dream_nested_%j.err

set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

WT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
MAIN=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_INPUT_ROOT="$MAIN"
export MASLD_OUTPUT_ROOT="$WT"

cd "$WT/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
mkdir -p ../logs

echo "=== B2 dream nested-RE ==="
echo "Start: $(date)"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-1}"

Rscript 05b_dream_nested_re.R

echo "Done: $(date)"
