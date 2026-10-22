#!/bin/bash
#SBATCH --job-name=6284
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H88-20261015T130436Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H88-20261015T130436Z/logs/%j.err
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H88-20261015T130436Z
export TMPDIR=$RUN/tmp
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $RUN/plot_h88.R $RUN/out $RUN/figures
