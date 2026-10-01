#!/bin/bash
#SBATCH --job-name=masld_native_ATAC_panel
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --dependency=afterok:21773775
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_native_ATAC_panel_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_native_ATAC_panel_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
SRC=$BASE/data/currin-profile-21773775
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_native_retention_tests.py" --out "$BASE/model/native_retention_invariants_${SLURM_JOB_ID}.json"
"$R" "$PROJ/scripts/analysis/alphagenome_campaign/model_native_retention_prepare.R" \
 "$SRC/private_source/liver_ATAC_peakCounts_DESeqDatasets_WASP-filtered_consensusPeaks_autosomalOnly_refinedBoundaries.RData" \
 "$SRC/inspection.json" \
 "$PROJ/GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1/supplementalData1_liver_ATAC_peaks.bed.gz" \
 "$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_native_retention_panel.py" \
 "$BASE/model/native_retention_panel_${SLURM_JOB_ID}"
