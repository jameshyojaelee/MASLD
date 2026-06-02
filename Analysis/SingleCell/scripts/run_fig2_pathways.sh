#!/bin/bash
#SBATCH --job-name=fig3_pathways
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=logs/fig3_pathways_%j.out
#SBATCH --error=logs/fig3_pathways_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell
mkdir -p logs

echo "=== Fig3 Cell-type Pathway Enrichment ==="
echo "Start: $(date)"

micromamba run -n rnaseq Rscript scripts/fig2_celltype_pathways.R

echo "=== Done: $(date) ==="
ls -la results_gpu_v2/fig2_data/celltype_pathway_enrichment.csv
