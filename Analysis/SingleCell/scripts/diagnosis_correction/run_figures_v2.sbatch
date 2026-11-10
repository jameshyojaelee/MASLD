#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/figures_v2_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/figures_v2_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
SCRIPT="$TASK/job_scripts/render_results_${SLURM_JOB_ID}.R"
cp Analysis/SingleCell/scripts/diagnosis_correction/render_results.R "$SCRIPT"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript "$SCRIPT"
mkdir "$TASK/figure_previews_v3"
for PDF in "$TASK"/figures_v3/*.pdf; do
  NAME="$(basename "$PDF" .pdf)"
  pdftotext -layout "$PDF" "$TASK/figure_previews_v3/$NAME.txt"
  pdftoppm -r 200 -singlefile -png "$PDF" "$TASK/figure_previews_v3/$NAME"
done
touch "$TASK/PANELS_V3_RENDERED_FOR_REVIEW"
