#!/bin/bash
#SBATCH --job-name=liver_sex_analysis
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/sex_analysis_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/sex_analysis_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=3:00:00

# Targeted sex analysis re-run — 5 existing cohorts (no GSE213621, pending upstream).
#
# Runs:
#   01_sample_qc.R     — re-adds sex_source / sex_final columns to meta_matched.rds
#   02_per_study_de.R  — re-runs per-study DE; GSE135251 now includes inferred_sex covariate
#
# Does NOT re-run: 03 (counts unchanged), 04 (variance partition), 05 (dream,
# already uses inferred_sex), 07/08 (consensus/pathway — re-run after GSE213621).
#
# Re-run the full pipeline from script 00 once GSE213621 featureCounts finishes.

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "=== Job started: $(date) ==="
echo "SLURM_JOB_ID:        $SLURM_JOB_ID"
echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"
echo "Partition: cpu | RAM: 64G | CPUs: 4"
echo "Datasets:  5-cohort (excl. GSE213621 — upstream alignment pending)"
echo ""

# Script 01: Re-run QC to add sex_source / sex_final columns.
# GSE213621 counts absent → auto-skipped by the script (not in merged counts).
echo "=== 01: Sample QC + Sex Inference ==="
Rscript "$SCRIPTS/01_sample_qc.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 01_sample_qc.R"; exit 1; fi

# Script 02: Per-study DE.
# GSE135251 now uses '~ 0 + condition + inferred_sex' (was '~ 0 + condition').
# GSE213621 auto-skipped (no counts). All other datasets unchanged.
echo ""
echo "=== 02: Per-Study DE (all 5 cohorts) ==="
Rscript "$SCRIPTS/02_per_study_de.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 02_per_study_de.R"; exit 1; fi

echo ""
echo "=== Sex analysis complete: $(date) ==="
echo ""
echo "Next steps:"
echo "  1. Check logs for k-means concordance rate (script 01 output)"
echo "  2. Review sex_stratified_meta_results.csv for sex-differential DEGs"
echo "  3. Re-run full pipeline (scripts 00→08) after GSE213621 finishes"
