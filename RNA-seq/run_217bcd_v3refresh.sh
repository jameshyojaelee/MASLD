#!/bin/bash
#SBATCH --job-name=217bcd_v3
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/stratified_causal/logs/217bcd_v3_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/stratified_causal/logs/217bcd_v3_%j.err

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/results/stratified_causal/logs

echo "=== 217b/c/d hormone + motifbreakr v3 refresh ==="
echo "Started: $(date)"

echo ""
echo "--- 217b hormone_tf_coloc_enrichment_clean ---"
Rscript RNA-seq/217b_hormone_tf_coloc_enrichment_clean.R 2>&1 || echo "WARN: 217b failed"

echo ""
echo "--- 217c motifbreakr_sex_biased_clean ---"
Rscript RNA-seq/217c_motifbreakr_sex_biased_clean.R 2>&1 || echo "WARN: 217c failed"

echo ""
echo "--- 217d hormone_panel_annotation_clean ---"
Rscript RNA-seq/217d_hormone_panel_annotation_clean.R 2>&1 || echo "WARN: 217d failed"

echo ""
echo "Finished: $(date)"
