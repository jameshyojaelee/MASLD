#!/bin/bash
#SBATCH --job-name=DIG1panx71a
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
sha256sum $RUN/scripts/71_panx1_pileup.sh $RUN/scripts/71_targets.tsv > $RUN/script_sha256_job71a.txt
bash $RUN/scripts/71_panx1_pileup.sh $RUN/out_panx1_ase_counts
