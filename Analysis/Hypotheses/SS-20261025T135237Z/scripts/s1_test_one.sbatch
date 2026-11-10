#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=interactive
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.err
# SS S1 end-to-end test of s1_one.sh on one library (usage: sbatch s1_test_one.sbatch <S>) (mate-order probe, STAR, strand gate, r8 counter, genotypes, Ig) on S0001, S0026 and
# S1169 (known per-run orientation from the S0a BAM; S0001/S0026 have superseded first alignments). Run before the array.
set -uo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
mkdir -p $SS/logs/test_unit
S=$1
/usr/bin/time -v -o $SS/logs/test_unit/${SLURM_JOB_ID}_$S.time bash $SS/scripts/s1_one.sh $S > $SS/logs/test_unit/${SLURM_JOB_ID}_$S.log 2>&1; echo "$S exit $?"
