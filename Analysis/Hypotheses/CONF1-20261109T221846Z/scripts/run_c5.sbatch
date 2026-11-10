#!/bin/bash
#SBATCH --job-name=CONF1c5
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=80G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z/logs/%j.err
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z
sha256sum "$RUN/PRESPEC_C5_C1b.md" "$RUN/scripts/c5_gse289173_extract_c1b.py" "$RUN/scripts/c5_trio_gse289173.R" > "$RUN/sha256_c5_${SLURM_JOB_ID}.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
$ENV/python -I "$RUN/scripts/c5_gse289173_extract_c1b.py" "$RUN/out_c5_extract_c1b"
$ENV/Rscript "$RUN/scripts/c5_trio_gse289173.R" "$RUN/out_c5_extract_c1b" "$RUN/out_c5_trio"
