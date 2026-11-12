#!/bin/bash
#SBATCH --job-name=DIR6_s2
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIR6-20261111T024808Z/logs/dir6_step2_%j.out
# Rerun of step 2 only (summary-code fix); step-1 extraction copied read-only from out_23190132.
set -eo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RUN=$ROOT/Analysis/Hypotheses/DIR6-20261111T024808Z
S1=$RUN/out_23190132
OUT=$RUN/out_${SLURM_JOB_ID}
if [[ -e $OUT ]]; then echo "$OUT exists; not overwriting"; exit 1; fi
mkdir -p $OUT
for f in rows.tsv instruments.tsv controls_variants.tsv proxies.tsv orientation_asserts.tsv pip_freeze_spatial.txt; do cp $S1/$f $OUT/; done
sha256sum $RUN/PRESPEC_DIR6.md $RUN/scripts/01_extract.py $RUN/scripts/02_classify.R $RUN/scripts/run_dir6_step2.sbatch $S1/rows.tsv $S1/instruments.tsv > $OUT/script_and_prespec_sha256.txt
/gpfs/commons/home/jameslee/micromamba/envs/genomicsem/bin/Rscript --vanilla $RUN/scripts/02_classify.R $ROOT $OUT
echo "DIR6 DONE"
