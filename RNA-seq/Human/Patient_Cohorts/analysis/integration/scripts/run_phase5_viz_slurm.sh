#!/bin/bash
#SBATCH --job-name=liver_phase5_viz
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/phase5_viz_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/phase5_viz_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --dependency=afterok:13691948

# Phase 5: Visualization & Validation (runs after 03-08 pipeline completes)
# Dependencies: merged_dge.rds, dream_results.csv, consensus_degs.csv, per-study DE CSVs

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "=== Phase 5 started: $(date) ==="
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "Depends on: job 13691948 (03-08 pipeline)"
echo ""

# Script 09: Batch correction + UMAP
echo "=== 09: Batch Correction & UMAP ==="
Rscript analysis/integration/scripts/09_batch_correction_umap.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 09_batch_correction_umap.R"; exit 1; fi

# Script 10: Volcano plots
echo ""
echo "=== 10: Volcano Plots ==="
Rscript analysis/integration/scripts/10_volcano_plots.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 10_volcano_plots.R"; exit 1; fi

# Script 11: Ortholog mapping (requires internet for biomaRt)
echo ""
echo "=== 11: Ortholog Mapping ==="

# Script 12: Library intersection
echo ""
echo "=== 12: Library Intersection ==="
Rscript analysis/integration/scripts/12_library_intersection.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 12_library_intersection.R"; exit 1; fi

echo ""
echo "=== Phase 5 complete: $(date) ==="
