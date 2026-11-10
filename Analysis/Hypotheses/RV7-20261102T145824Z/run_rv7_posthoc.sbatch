#!/bin/bash
#SBATCH --job-name=rv7ph
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV7-20261102T145824Z/logs/%x_%j.out
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV7-20261102T145824Z
sha256sum scripts/rv7_posthoc_local_enzyme.R > out-20261102T150010Z/posthoc_script_sha256.txt
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/rv7_posthoc_local_enzyme.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV7-20261102T145824Z/out-20261102T150010Z
