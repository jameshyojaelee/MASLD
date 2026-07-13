#!/bin/bash
#SBATCH --job-name=edgeR
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/ATAC/Human_Multiome/logs/mash_vs_normal_%j.out
#SBATCH --error=Analysis/ATAC/Human_Multiome/logs/mash_vs_normal_%j.err

# MASH-vs-NORMAL documented endpoint edgeR (reuses the donor pseudobulk count
# matrices already on disk — no h5ad re-aggregation), then regenerate Fig S4g
# with both documented contrasts. Submitted with --dependency=afterok on the
# disease job so the MASLD-vs-control edgeR files exist before the figure runs.
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$(pwd)"
export PYTHONNOUSERSITE=1
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=== Step 1: MASH-vs-NORMAL edgeR (4 compartments) ==="
Rscript Analysis/ATAC/Human_Multiome/scripts/17d_edger_mash_vs_normal.R

echo ""
echo "=== Step 2: regenerate Fig S4g (MASLD-vs-Control + MASH-vs-Control) ==="
Rscript scripts/figures/figS05_scatac_disease_vs_control.R

echo ""
echo "=== Done: $(date) ==="
ls -la Analysis/ATAC/Human_Multiome/results/snapatac2/scatac_da_mashvsnormal_*_edger.csv || true
ls -la figures/main/fig4_validation/panels/figS4g.pdf || true
