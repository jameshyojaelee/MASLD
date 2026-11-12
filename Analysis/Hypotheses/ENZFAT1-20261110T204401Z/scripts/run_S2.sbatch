#!/bin/bash
#SBATCH --job-name=ENZFAT1_S2
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/logs/S2_%j.out
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RUN=$ROOT/Analysis/Hypotheses/ENZFAT1-20261110T204401Z
OUT=$RUN/out_S2
LOCI=$ROOT/Analysis/Hypotheses/DIG1-20261024T113951Z/out_release_vs_disease_loci/loci.tsv
if [[ -e $OUT ]]; then echo "$OUT exists; not overwriting"; exit 1; fi
mkdir -p $OUT
sha256sum $RUN/scripts/20_extract_loci.py $RUN/scripts/21_locus_effects.R $RUN/scripts/run_S2.sbatch $LOCI > $OUT/script_and_input_sha256.txt
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I $RUN/scripts/20_extract_loci.py $ROOT $LOCI $OUT/locus_rows.tsv
/gpfs/commons/home/jameslee/micromamba/envs/genomicsem/bin/Rscript --vanilla $RUN/scripts/21_locus_effects.R $OUT/locus_rows.tsv $LOCI $OUT
