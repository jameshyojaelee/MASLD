#!/bin/bash
#SBATCH --job-name=H36_int_test
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H36-20261009T180345Z/logs/slurm-codepath_test_integrate-%j.out

# Code-path test of h36_05_integrate.R on a NULL stand-in (observed := permutation 1, permuted := permutations 2..10
# of the first completed chunks). Not a result; output goes to codepath_test_null/.
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
H=$REPO/Analysis/Hypotheses/H36-20261009T180345Z
SC=$REPO/scripts/analysis/hypotheses_2026q4/H36_splicing
T=$H/codepath_test_null
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
test ! -e "$T"; mkdir -p "$T/root" "$T/out"
ln -s $H/build $T/root/build; ln -s $H/sqtl_coloc $T/root/sqtl_coloc
$RS $SC/tests/make_null_stage_standin.R $H/stage_test $T/root/stage_test 2 $H/build
$RS $SC/h36_05_integrate.R $REPO $T/root $T/out
