#!/bin/bash
#SBATCH --job-name=liver_10cohort
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/10cohort_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/10cohort_%j.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --time=12:00:00

# =============================================================================
# 10-cohort integration pipeline
#
# Runs the full integration pipeline (scripts 00-08, 25, 26) on all 10 human
# MASLD cohorts. This is the expanded re-run adding GSE162694, GSE174478,
# GSE193066, and GSE240729 to the original 6 cohorts.
#
# Mega-analysis cohorts (controls; canonical 5 after PRJNA512027 removed
#   2026-05-15): GSE126848, GSE130970, GSE135251, GSE213621, GSE162694
# Per-study DE only (no controls): GSE167523, GSE174478, GSE193066, GSE240729
#
# Usage:
#   sbatch run_10cohort_integration.sh
#   # Or with dependency on a prior job:
#   sbatch --dependency=afterok:JOBID run_10cohort_integration.sh
# =============================================================================

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

SCRIPTS="analysis/integration/scripts"
LOGDIR="analysis/integration/logs"
mkdir -p "${LOGDIR}"

echo "=== 10-Cohort Integration Pipeline ==="
echo "Job started: $(date)"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "Partition: bigmem | RAM: 400G | CPUs: 32"
echo ""

run_script() {
    local script="$1"
    local name=$(basename "$script" .R)
    echo "=== ${name} ==="
    echo "Started: $(date)"
    Rscript "${script}" 2>&1
    echo "Completed: $(date)"
    echo ""
}

# Phase 1: Core integration (sequential — each depends on prior)
run_script "${SCRIPTS}/00_harmonize_metadata.R"
run_script "${SCRIPTS}/01_sample_qc.R"
run_script "${SCRIPTS}/02_per_study_de.R"
run_script "${SCRIPTS}/03_integrate_counts.R"
run_script "${SCRIPTS}/04_variance_partition.R"
run_script "${SCRIPTS}/05_dream_mega_analysis.R"
run_script "${SCRIPTS}/07_consensus_degs.R"
run_script "${SCRIPTS}/08_pathway_analysis.R"

# Phase 2: Downstream analyses (depend on dream results from 05)
run_script "${SCRIPTS}/25_deconv_attribution.R"
run_script "${SCRIPTS}/26_sex_stratified_analysis.R"

# Phase 3: Visualization and annotation (depend on 05/07/08)
run_script "${SCRIPTS}/09_batch_correction_umap.R"
run_script "${SCRIPTS}/10_volcano_plots.R"
run_script "${SCRIPTS}/12_library_intersection.R"
run_script "${SCRIPTS}/17_annotate_gene_symbols.R"

# Phase 4: Disease subtype analyses (depend on 03/05)
run_script "${SCRIPTS}/13_nafl_vs_nash_de.R"
run_script "${SCRIPTS}/14_fibrosis_progression_de.R"
run_script "${SCRIPTS}/15_nas_component_de.R"
run_script "${SCRIPTS}/16_disease_signatures_consensus.R"

echo "=== All 10-cohort integration scripts completed ==="
echo "Finished: $(date)"
