#!/bin/bash
#SBATCH --job-name=a4_posthoc
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A4_multiome_links_zonation/logs/%x_%j.log
# A4 job E (post hoc, amendment 2026-10-25): per-replicate z of full-data hepatocyte links, same 100 bootstrap seeds.
set -euo pipefail
A4=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A4_multiome_links_zonation
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$A4"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/02_links.R "$A4" posthoc hepatocyte 100 7
