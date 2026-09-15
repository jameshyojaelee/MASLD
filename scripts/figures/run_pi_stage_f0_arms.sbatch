#!/bin/bash
#SBATCH --job-name=f0_arms
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/f0_arms_%j.out
#SBATCH --error=logs/f0_arms_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
source /etc/profile.d/modules.sh
module load R
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export STAGE_RELEASE_ROOT="$MASLD_PROJECT_ROOT/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-POOLED-REPRO"
export FIGURE_CANDIDATE_ROOT="$MASLD_PROJECT_ROOT/figures/candidates/f0-reference-arms-20260907T105813Z"
Rscript scripts/figures/build_pi_stage_f0_arms.R
echo "F0_ARMS_JOB_COMPLETE"
