#!/bin/bash
#SBATCH --job-name=DIR6_ph
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIR6-20261111T024808Z/logs/dir6_posthoc_%j.out
set -eo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIR6-20261111T024808Z
OUT=$RUN/out_posthoc_${SLURM_JOB_ID}
if [[ -e $OUT ]]; then echo "$OUT exists"; exit 1; fi
mkdir -p $OUT
sha256sum $RUN/scripts/03_posthoc_control_refit.R $RUN/out_23210381/instrument_effects.tsv > $OUT/sha256.txt
/gpfs/commons/home/jameslee/micromamba/envs/genomicsem/bin/Rscript --vanilla $RUN/scripts/03_posthoc_control_refit.R $RUN/out_23210381 $OUT
echo DONE
