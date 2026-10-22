#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN6/gtex-area-20261021T124444Z/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN6/gtex-area-20261021T124444Z/%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 TMPDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN6/gtex-area-20261021T124444Z/tmp
mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN6/gtex-area-20261021T124444Z/tmp
/usr/bin/time -v /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN6/gtex-area-20261021T124444Z/scripts/h55_area_adjusted.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN6/gtex-area-20261021T124444Z/out
