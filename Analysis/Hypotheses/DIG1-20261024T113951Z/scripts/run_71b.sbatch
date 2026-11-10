#!/bin/bash
#SBATCH --job-name=DIG1panx71b
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
sha256sum $RUN/scripts/71b_panx1_ase.R > $RUN/script_sha256_job71b.txt
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $RUN/scripts/71b_panx1_ase.R $RUN/out_panx1_ase_counts $RUN/out_panx1_ase
