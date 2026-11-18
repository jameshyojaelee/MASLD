#!/bin/bash
#SBATCH --job-name=ENZFAT1_S1b_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:25:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/logs/S1b_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/logs/S1b_20261118_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RUN=$ROOT/Analysis/Hypotheses/ENZFAT1-20261110T204401Z
OUT=$RUN/out_S1b_20261118
if [[ -e $OUT ]]; then echo "$OUT exists; not overwriting"; exit 1; fi
mkdir -p $OUT
sha256sum $RUN/scripts/11_shares_subpd_20261118.R $RUN/scripts/run_S1b_20261118.sbatch $RUN/out_S1/ldsc_all.rds > $OUT/script_and_input_sha256.txt
/gpfs/commons/home/jameslee/micromamba/envs/genomicsem/bin/Rscript --vanilla $RUN/scripts/11_shares_subpd_20261118.R $RUN/out_S1 $OUT
