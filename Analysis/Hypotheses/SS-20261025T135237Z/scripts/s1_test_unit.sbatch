#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.err
# SS S1 end-to-end test of s1_one.sh (mate-order probe, STAR, strand gate, r8 counter, genotypes, Ig) on S0001, S0026 and
# S1169 (known per-run orientation from the S0a BAM; S0001/S0026 have superseded first alignments). Run before the array.
set -uo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
mkdir -p $SS/logs/test_unit
printf 'S0001\nS0026\nS1169\n' | xargs -P 3 -I{} bash -c "/usr/bin/time -v -o $SS/logs/test_unit/${SLURM_JOB_ID}_{}.time bash $SS/scripts/s1_one.sh {} > $SS/logs/test_unit/${SLURM_JOB_ID}_{}.log 2>&1; echo {} exit \$?"
