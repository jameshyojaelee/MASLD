#!/bin/bash
#SBATCH --job-name=liver_refactor_rerun
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/refactor_rerun_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/refactor_rerun_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/

echo "=== Running refactored scripts ==="
echo "Started at $(date)"

Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/07_consensus_degs.R
Rscript RNA-seq/Mouse/Unified_Integration/scripts/M04_mouse_consensus.R
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/08_pathway_analysis.R

# 16 needs 07 done first to build unified disease signatures
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/16_disease_signatures_consensus.R

# 12 is library intersection
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/12_library_intersection.R

# 19 and 20 depend on unified disease signatures from 16
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/19_mendelian_randomization.R
Rscript RNA-seq/20_drug_repurposing_v3.R

# 27a assembles multi-evidence atlas (replaces old weighted 27)
Rscript RNA-seq/27a_assemble_evidence_atlas.R

# 27b benchmarks preset configurations
Rscript RNA-seq/27b_benchmark_presets.R

# 33 acts upon multiple outputs limit to test assumptions
Rscript RNA-seq/33_audit_sensitivity_analyses.R

echo "=== Done at $(date) ==="
