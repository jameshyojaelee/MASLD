#!/bin/bash
#SBATCH --job-name=ENZFAT1_S2_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/logs/S2_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/logs/S2_20261118_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RUN=$ROOT/Analysis/Hypotheses/ENZFAT1-20261110T204401Z
OUT=$RUN/out_S2_20261118
LOCI=$ROOT/Analysis/Hypotheses/DIG1-20261024T113951Z/out_release_vs_disease_loci_20261118/loci.tsv
if [[ -e $OUT ]]; then echo "$OUT exists; not overwriting"; exit 1; fi
mkdir -p $OUT
echo "host $(hostname) | $(date -u)"
sha256sum $RUN/scripts/20_extract_loci_20261118.py $RUN/scripts/21_locus_effects_20261118.R $RUN/scripts/run_S2_20261118.sbatch $LOCI > $OUT/script_and_input_sha256.txt
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I $RUN/scripts/20_extract_loci_20261118.py $ROOT $LOCI $OUT/locus_rows.tsv
/gpfs/commons/home/jameslee/micromamba/envs/genomicsem/bin/Rscript --vanilla $RUN/scripts/21_locus_effects_20261118.R $OUT/locus_rows.tsv $LOCI $OUT
echo "DONE $(date -u)"
