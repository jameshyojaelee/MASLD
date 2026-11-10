#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/refits_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/refits_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
test -f "$TASK/REFITS_PREPARED"
export DIAGNOSIS_REFIT_ROOT="$TASK/OLS_ambient"
export QUALIFIED_DONOR_METADATA="$TASK/metadata/donor_metadata_stage_corrected_donor.tsv"
if ! test -f "$TASK/OLS_AMBIENT_COMPLETED"; then
  "$RS" "$TASK/refit_scripts/ols.R"
  "$PY" "$TASK/refit_scripts/ambient.py"
  touch "$TASK/OLS_AMBIENT_COMPLETED"
fi
export CCC_QUALIFIED_ROOT="$TASK/CCC"
if ! test -f "$TASK/CCC_COMPLETED"; then
  "$RS" "$TASK/refit_scripts/ccc.R"
  touch "$TASK/CCC_COMPLETED"
fi
export H17_QUALIFIED_ROSTER="$TASK/metadata/current_source_diagnosis_donor_roster.tsv"
mkdir -p "$TASK/H17"
for STEP in 01_de_and_concordance 02_resampling 03_sensitivity_and_call; do
  "$RS" "$TASK/refit_scripts/H17/$STEP.R" "$PWD/Analysis/Hypotheses/H17-20261009T015231Z" "$TASK/H17"
done
touch "$TASK/H17_COMPLETED"
"$PY" "$TASK/refit_scripts/h34.py" "$PWD" "$TASK/H34" "$PWD/Analysis/Hypotheses/H34-20261008T014810Z/intermediate"
touch "$TASK/H34_COMPLETED"
