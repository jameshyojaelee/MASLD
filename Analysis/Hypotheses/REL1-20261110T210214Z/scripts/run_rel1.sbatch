#!/bin/bash
#SBATCH --job-name=REL1
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REL1-20261110T210214Z/logs/rel1_%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REL1-20261110T210214Z
OUT=$RUN/out_${SLURM_JOB_ID}
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
cd $RUN
$R --vanilla $RUN/scripts/rel1.R $OUT
sha256sum $RUN/PRESPEC_REL1.md $RUN/scripts/rel1.R $RUN/scripts/run_rel1.sbatch > $OUT/script_sha256.tsv
echo "elapsed ${SECONDS}s"
