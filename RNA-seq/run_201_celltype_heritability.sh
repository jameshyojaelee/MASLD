#!/bin/bash
#SBATCH --job-name=heritability
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/201_celltype_heritability_%j.out
#SBATCH --error=logs/201_celltype_heritability_%j.err
#
# Re-run 201_celltype_heritability.R with the :94 glob filter (mega-review A6.x
# FDR-family contamination fix) so the on-disk celltype_heritability_results.csv
# becomes the CLEAN 11-row per-cell-type family (was the OLD contaminated 63-row
# CSV that pulled in *_F{N}_vs_F{M}_bayesprism stage pseudo-celltypes). This also
# regenerates celltype_specificity_scores.csv (feeds 206 atlas sc_tau/sc_best_ct).

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$BASE/RNA-seq"

micromamba run -n rnaseq Rscript "$BASE/RNA-seq/201_celltype_heritability.R"
