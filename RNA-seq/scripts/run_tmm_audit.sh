#!/bin/bash
#SBATCH --job-name=B4_tmm_audit
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=logs/B4/tmm_audit_%j.out
#SBATCH --error=logs/B4/tmm_audit_%j.err

set -eo pipefail
source /gpfs/commons/home/jameslee/.bashrc
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "Host: $(hostname)"
# Best run AFTER 14_5 produces jaccard_summary.csv (REPORT.md interpolates it),
# but the per-cohort audit half is independent — script handles missing input.
Rscript RNA-seq/scripts/tmm_directional_bias_audit.R
echo "Done."
