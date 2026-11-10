#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/context_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/context_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PY" Analysis/SingleCell/scripts/diagnosis_correction/repair_library_census.py
"$PY" Analysis/SingleCell/integration/scripts/add_sample_metadata.py
"$RS" Analysis/SingleCell/scripts/diagnosis_correction/refit_composition.R
touch "$TASK/COMPOSITION_COMPLETED"
"$RS" Analysis/SingleCell/scripts/314_allcell_pseudobulk_de.R
touch "$TASK/ALLCELL_COMPLETED"
export SINGLECELL_ONLY=1 MULTI_EVIDENCE_ATLAS_OUT="$TASK/multi_evidence_atlas_sc_corrected.csv"
"$RS" RNA-seq/45a_integrate_all_sources.R
touch "$TASK/SC_ANNOTATIONS_COMPLETED"
