#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/review_panels_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/review_panels_%j.err
set -euo pipefail
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PY" Analysis/SingleCell/scripts/diagnosis_correction/qualify_ambient_states.py
mkdir "$TASK/figure_previews"
for PDF in "$TASK"/figures/*.pdf; do
  NAME="$(basename "$PDF" .pdf)"
  pdftotext -layout "$PDF" "$TASK/figure_previews/$NAME.txt"
  pdftoppm -r 200 -singlefile -png "$PDF" "$TASK/figure_previews/$NAME"
done
touch "$TASK/PANELS_RENDERED_FOR_REVIEW"
