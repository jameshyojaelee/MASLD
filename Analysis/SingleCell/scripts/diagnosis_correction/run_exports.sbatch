#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/exports_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/exports_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PY" Analysis/SingleCell/scripts/diagnosis_correction/qualify_documented_F.py
"$PY" Analysis/SingleCell/scripts/diagnosis_correction/export_context.py
"$PY" scripts/portal/build_gene_detail_parquets.py --only pseudobulk_de --output-dir "$TASK/portal"
touch "$TASK/EXPORTS_COMPLETED"
