#!/bin/bash
#SBATCH --job-name=limma
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/per_study_de_qw_%j.log

# Re-run per-study limma-voom DE with quality weights enabled (voomWithQualityWeights),
# matching the canonical limma_voom_qw method. Triggered by config flip
# quality_weights: false -> true (2026-06-11).
set -eo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$MASLD_PROJECT_ROOT/RNA-seq/Human/Patient_Cohorts/analysis/integration"
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript "$MASLD_PROJECT_ROOT/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/02_per_study_de.R"
echo "[done] per-study QW DE complete"
