#!/bin/bash
#SBATCH --job-name=ENZFAT1_S1b
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/logs/S1b_%j.out
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RUN=$ROOT/Analysis/Hypotheses/ENZFAT1-20261110T204401Z
OUT=$RUN/out_S1b
if [[ -e $OUT ]]; then echo "$OUT exists; not overwriting"; exit 1; fi
mkdir -p $OUT
sha256sum $RUN/scripts/11_shares_subpd.R $RUN/scripts/run_S1b.sbatch $RUN/out_S1/ldsc_all.rds > $OUT/script_and_input_sha256.txt
/gpfs/commons/home/jameslee/micromamba/envs/genomicsem/bin/Rscript --vanilla $RUN/scripts/11_shares_subpd.R $RUN/out_S1 $OUT
