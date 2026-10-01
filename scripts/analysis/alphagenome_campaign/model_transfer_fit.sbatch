#!/bin/bash
#SBATCH --job-name=masld_ag_target_transfer
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --dependency=afterok:21772928
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_target_transfer_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_target_transfer_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_transfer_fit.py" \
 --features "$BASE/model/reporter_targets_21772928" \
 --reporters "$BASE/data/endpoints-revised-21772544/reporter_treatment_replicates.tsv.gz" \
 --strata-spec "$BASE/data/transfer_strata_spec.json" \
 --out "$BASE/model/transfer_${SLURM_JOB_ID}" --bootstrap 2000
