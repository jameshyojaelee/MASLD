#!/bin/bash
#SBATCH --job-name=masld_lx2_full_population
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --dependency=afterok:21772928
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_lx2_full_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_lx2_full_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_feature_agreement.py" \
 --archived "$PROJ/GWAS/finemapping/results/alphagenome_program/ad-arm2-mpra-20260915T114407Z/embeddings/alphagenome_frozen_trunk_2048bp.npz" \
 --current "$BASE/model/reporter_targets_21772928/reporters_2048.npz" \
 --out "$BASE/model/lx2_full_feature_agreement_${SLURM_JOB_ID}.json"
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_lx2.py" \
 --endpoints "$BASE/data/endpoints-revised-21772544/reporter_treatment_replicates.tsv.gz" \
 --embeddings "$BASE/model/reporter_targets_21772928/reporters_2048.npz" \
 --row-authority "$BASE/data/folds-21772544/reporter_fold_companion.tsv" \
 --population-rule "all4359_source_constructs_sequence_and_four_replicate_eligibility" \
 --out "$BASE/model/lx2_full_${SLURM_JOB_ID}" --bootstrap 2000
