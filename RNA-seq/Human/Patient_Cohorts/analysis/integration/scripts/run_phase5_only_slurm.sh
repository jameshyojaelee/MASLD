#!/bin/bash
#SBATCH --job-name=liver_phase5
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/phase5_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/phase5_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

echo "=== 09: Batch Correction + UMAP ==="
Rscript scripts/09_batch_correction_umap.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 09"; exit 1; fi

echo ""
echo "=== 10: Volcano Plots ==="
Rscript scripts/10_volcano_plots.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 10"; exit 1; fi

echo ""
echo "=== 11: Ortholog Mapping — SKIPPED (11_ortholog_mapping.R removed; ortholog logic moved to Cross_Species_Concordance pipeline and script 17) ==="

echo ""
echo "=== 12: Library Intersection ==="
Rscript scripts/12_library_intersection.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 12"; exit 1; fi

echo ""
echo "=== Phase 5 complete: $(date) ==="
