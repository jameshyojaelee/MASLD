#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/bf3_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/bf3_%j.err
set -euo pipefail
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
BF="$TASK/BF3"
S="$TASK/refit_scripts/BF3"
mkdir "$BF/primary"
"$PY" "$S/bf3_decompose.py" "$BF/build" "$BF/primary"
"$PY" "$S/bf3_m_slope.py" "$BF/primary" "$BF/primary"
for STEP in marker_split posthoc posthoc2; do
  mkdir "$BF/$STEP"
  "$PY" "$S/bf3_$STEP.py" "$BF/build" "$BF/primary" "$BF/$STEP"
done
mkdir "$BF/doublet_matched"
"$PY" "$S/bf3_doublet_matched.py" "$PWD/Analysis/Hypotheses/BF3-20261023T154841Z/doublet-20261023T232600Z" "$BF/build" "$BF/doublet_matched"
touch "$TASK/BF3_COMPLETED"
