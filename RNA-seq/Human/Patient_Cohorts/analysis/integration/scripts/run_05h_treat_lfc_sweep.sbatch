#!/bin/bash
#SBATCH --job-name=treatsweep
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/treatsweep_%j.out
#SBATCH --error=RNA-seq/logs/treatsweep_%j.err
#
# audit P2#18: TREAT effect-size-floor sweep {0.15,0.20,0.25,0.30}.
# Re-fits the C2 limma_voom_qw engine ONCE and evaluates treat() at each lfc.
# Writes ONLY results/audit_sensitivity/treat_lfc_sweep.csv; does NOT overwrite
# canonical_deg_results.csv.

set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

export TREAT_LFC_GRID="0.15,0.20,0.25,0.30"
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/05h_treat_lfc_sweep.R
