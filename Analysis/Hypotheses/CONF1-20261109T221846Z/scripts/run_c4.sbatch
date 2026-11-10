#!/bin/bash
#SBATCH --job-name=CONF1c4
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z/logs/%j.err
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z
sha256sum "$RUN/PRESPEC_C3_C4.md" "$RUN/scripts/c4_catalog_candidate.R" > "$RUN/sha256_c4_${SLURM_JOB_ID}.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript "$RUN/scripts/c4_catalog_candidate.R" "$RUN/out_c4_catalog_candidate"
