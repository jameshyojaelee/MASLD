#!/bin/bash
#SBATCH --job-name=DIG1null72b
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
sha256sum $RUN/scripts/72b_null_ase.R >> $RUN/script_sha256_job72.txt
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $RUN/scripts/72b_null_ase.R $RUN/out_ase_null_sites $RUN/out_panx1_ase/utr_het_counts.tsv
