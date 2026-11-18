#!/bin/bash
#SBATCH --job-name=eqtl_psd_D
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/logs/apply_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/logs/apply_%j.err
set -eo pipefail
export RUN_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118
export OUT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/evidence_passports/candidates/gene-catalog-conf1-eqtlpsd-20261118
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript ${RUN_DIR}/code/apply_refit_catalog.R
