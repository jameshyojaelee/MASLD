#!/bin/bash
#SBATCH --job-name=a4_zone
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A4_multiome_links_zonation/logs/%x_%j.log
# A4 job C: zonation validation, zone DA, LDSC peak sets and nulls, steatosis x zone with 1,000 permutations (04).
set -euo pipefail
A4=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A4_multiome_links_zonation
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$A4"
test -s work/zone/atac_peak_by_unit.f32
/usr/bin/time -v $ENV/bin/Rscript scripts/04_zonation_da.R "$A4" 1000 14 2> logs/time_04.txt
tail -n 40 logs/time_04.txt
