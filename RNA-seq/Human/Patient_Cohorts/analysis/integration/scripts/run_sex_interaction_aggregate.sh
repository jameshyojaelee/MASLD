#!/bin/bash
#SBATCH --job-name=sex_int_agg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_interaction_subsampling/logs/aggregate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_interaction_subsampling/logs/aggregate_%j.err

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/14.4d_aggregate_interaction_subsampling.R"

echo "=== Aggregating interaction subsampling results ==="
echo "Date: $(date)"

Rscript "$SCRIPT"

echo "=== Done: $(date) ==="
