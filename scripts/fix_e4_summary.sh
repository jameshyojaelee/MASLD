#!/bin/bash
#SBATCH --job-name=fix_e4_summary
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --qos=interactive
#SBATCH --output=scripts/logs/fix_e4_summary_%j.out
#SBATCH --error=scripts/logs/fix_e4_summary_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

micromamba run -n rnaseq Rscript RNA-seq/results/drug_repurposing/e4_tier_verification.R
