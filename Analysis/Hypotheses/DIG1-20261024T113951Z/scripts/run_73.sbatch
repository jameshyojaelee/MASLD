#!/bin/bash
#SBATCH --job-name=DIG1rel73
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
sha256sum $RUN/scripts/73_release_vs_disease_loci.py > $RUN/script_sha256_job73.txt
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I $RUN/scripts/73_release_vs_disease_loci.py $RUN/out_release_vs_disease_loci
