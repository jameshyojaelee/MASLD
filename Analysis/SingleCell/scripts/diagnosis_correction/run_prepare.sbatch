#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/prepare_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/prepare_%j.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
TASK_ROOT="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
if ! test -e "$TASK_ROOT/metadata/library_metadata.tsv"; then
  "$PYTHON" Analysis/SingleCell/scripts/diagnosis_correction/build_metadata.py --output "$TASK_ROOT/metadata"
fi
"$PYTHON" Analysis/SingleCell/scripts/diagnosis_correction/build_atlas_view.py \
  --source Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad \
  --metadata "$TASK_ROOT/metadata/library_metadata.tsv" \
  --output "$TASK_ROOT/integrated_atlas_source_qualified_v2.h5ad"
touch "$TASK_ROOT/PREPARED"
