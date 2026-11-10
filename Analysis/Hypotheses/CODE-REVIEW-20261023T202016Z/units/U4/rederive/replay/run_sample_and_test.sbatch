#!/bin/bash
#SBATCH --job-name=u4_replay_test
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/logs/%j.out
set -eo pipefail
RP=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/replay
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla $RP/sample_pairs.R
# time the heaviest-looking task first: first EUR positive pair
K=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="ancestry") a=i} NR>1 && $2=="positive" && $a=="EUR"{print $1; exit}' $RP/pairs.tsv)
echo "test task = $K"
/usr/bin/time -v $RP/replay_task.sh $K
