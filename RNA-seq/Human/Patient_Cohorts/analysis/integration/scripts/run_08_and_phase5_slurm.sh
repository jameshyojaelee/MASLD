#!/bin/bash
#SBATCH --job-name=liver_08_and_phase5
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/08_and_phase5_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/08_and_phase5_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00

# Run script 08 (pathway analysis, fixed msigdbr API) then Phase 5 scripts (09-12)
set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

echo "=== 08: Pathway Analysis (fixed msigdbr) ==="
Rscript scripts/08_pathway_analysis.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 08"; exit 1; fi

echo ""
echo "=== 09: Batch Correction + UMAP ==="
Rscript scripts/09_batch_correction_umap.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 09"; exit 1; fi

echo ""
echo "=== 10: Volcano Plots ==="
Rscript scripts/10_volcano_plots.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 10"; exit 1; fi

echo ""
echo "=== 11: Ortholog Mapping ==="
if [ $? -ne 0 ]; then echo "FAILED: 11"; exit 1; fi

echo ""
echo "=== 12: Library Intersection ==="
Rscript scripts/12_library_intersection.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 12"; exit 1; fi

echo ""
echo "=== All scripts complete: $(date) ==="
