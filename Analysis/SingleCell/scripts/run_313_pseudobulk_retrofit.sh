#!/bin/bash
#SBATCH --job-name=313_pb_retrofit
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/313_pb_retrofit_%j.out
#SBATCH --error=logs/313_pb_retrofit_%j.err

# Pachter P0 retrofit for Script 313: pseudobulk hepatocyte META-SUBTYPE markers
# Step 1 (Python, spatial env): re-run 313 to dump pseudobulk count matrices
#   (pseudobulk_meta_subtype/) + write VIZ_ONLY meta_subtype_markers.csv +
#   AUCell hepatocyte scores (unchanged).
# Step 2 (R, rnaseq env): 313b_pseudobulk_meta_subtype_markers.R produces
#   subtype_{META}_de.csv via limma-voom + duplicateCorrelation(block=donor),
#   n=donors design (Squair et al. 2021, Nat Commun).

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/Analysis/SingleCell
mkdir -p logs results_gpu_v2/hepatocyte_subtypes/crossmodal/bulk_sc_convergence/pseudobulk_meta_subtype

# Step 1: Python — 313 (pseudobulk dump + VIZ_ONLY rank_genes_groups + AUCell)
echo "[$(date)] Step 1: Python 313 (pseudobulk dump + AUCell)"
micromamba activate spatial
set -u
python scripts/313_bulk_sc_convergence.py
set +u
micromamba deactivate

# Step 2: R — 313b (limma-voom on pseudobulk by meta_subtype)
echo "[$(date)] Step 2: R 313b (limma-voom)"
micromamba activate rnaseq
set -u
Rscript scripts/313b_pseudobulk_meta_subtype_markers.R
set +u
micromamba deactivate

echo "[$(date)] Done. Pseudobulk meta-subtype markers in results_gpu_v2/hepatocyte_subtypes/crossmodal/bulk_sc_convergence/pseudobulk_meta_subtype/"
