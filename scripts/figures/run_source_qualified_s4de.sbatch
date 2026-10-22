#!/bin/bash
#SBATCH --job-name=qualified_S4DE
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=logs/qualified_S4DE_%j.out
#SBATCH --error=logs/qualified_S4DE_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export MASLD_PROJECT_ROOT="$PWD"
export S4DE_RENDER_ROOT="$PWD/figures/candidates/source-qualified-S4DE-20261014_${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/figures/render_source_qualified_s4de.R
cp scripts/figures/render_source_qualified_s4de.R scripts/figures/run_source_qualified_s4de.sbatch "$S4DE_RENDER_ROOT/"
for p in figs4d_communication_heatmap figs4e_program_linked_communication; do
 pdftotext "$S4DE_RENDER_ROOT/$p.pdf" "$S4DE_RENDER_ROOT/$p.txt"
 pdftoppm -r 150 -singlefile -png "$S4DE_RENDER_ROOT/$p.pdf" "$S4DE_RENDER_ROOT/$p"
done
