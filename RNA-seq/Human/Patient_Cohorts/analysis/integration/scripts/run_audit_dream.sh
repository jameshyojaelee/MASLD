#!/bin/bash
#SBATCH --job-name=B2_audit_dream
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/B2_audit_dream_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/B2_audit_dream_%j.err

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

echo "=== B2 audit dream nested vs canonical ==="
echo "Start: $(date)"

Rscript audit_dream_nested_vs_canonical.R

echo "Done: $(date)"
