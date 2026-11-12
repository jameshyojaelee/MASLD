#!/bin/bash
#SBATCH --job-name=ENZFAT1_S1
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/logs/S1_%j.out
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RUN=$ROOT/Analysis/Hypotheses/ENZFAT1-20261110T204401Z
OUT=$RUN/out_S1
if [[ -e $OUT ]]; then echo "$OUT exists; not overwriting"; exit 1; fi
mkdir -p $OUT
sha256sum $RUN/scripts/10_shares.R $RUN/scripts/run_S1.sbatch $RUN/out_S0/munged/*.sumstats.gz > $OUT/script_and_input_sha256.txt
/gpfs/commons/home/jameslee/micromamba/envs/genomicsem/bin/Rscript --vanilla $RUN/scripts/10_shares.R $RUN/out_S0 $OUT $ROOT/GWAS/ldsc_ref/eur_w_ld_chr/
