#!/bin/bash
#SBATCH --job-name=309_pb_retrofit
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/309_pb_retrofit_%j.out
#SBATCH --error=logs/309_pb_retrofit_%j.err

# S1 Pachter P0 retrofit: pseudobulk hepatocyte subtype markers
# Step 1 (Python, spatial env): re-run 309 to dump pseudobulk count matrices
#   + write VIZ_ONLY subtype_markers.csv from rank_genes_groups (legacy)
# Step 2 (R, rnaseq env): 309b_pseudobulk_hepatocyte_markers.R produces
#   subtype_{N}_de.csv via limma-voom on n=donors

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/Analysis/SingleCell
mkdir -p logs results_gpu_v2/hepatocyte_subtypes/pseudobulk_markers

# Step 1: Python — 309 (pseudobulk dump + VIZ_ONLY rank_genes_groups)
echo "[$(date)] Step 1: Python 309 (pseudobulk dump)"
micromamba activate spatial
python scripts/309_hepatocyte_subcluster_annotation.py
micromamba deactivate

# Step 2: R — 309b (limma-voom on pseudobulk)
echo "[$(date)] Step 2: R 309b (limma-voom)"
micromamba activate rnaseq
Rscript scripts/309b_pseudobulk_hepatocyte_markers.R
micromamba deactivate

echo "[$(date)] Done. Pseudobulk markers in results_gpu_v2/hepatocyte_subtypes/pseudobulk_markers/"
