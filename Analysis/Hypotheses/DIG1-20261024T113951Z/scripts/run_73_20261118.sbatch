#!/bin/bash
#SBATCH --job-name=DIG1rel73_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/rel73_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/rel73_20261118_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RUN=$ROOT/Analysis/Hypotheses/DIG1-20261024T113951Z
OUT=$RUN/out_release_vs_disease_loci_20261118
echo "host $(hostname) | $(date -u)"
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I \
  $RUN/scripts/73_release_vs_disease_loci_20261118.py $OUT
sha256sum $RUN/scripts/73_release_vs_disease_loci_20261118.py \
  $ROOT/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/results/susie_coloc_all_gwas_refit_20261118.csv \
  $OUT/loci.tsv > $OUT/script_and_input_sha256.txt
echo "DONE $(date -u)"
