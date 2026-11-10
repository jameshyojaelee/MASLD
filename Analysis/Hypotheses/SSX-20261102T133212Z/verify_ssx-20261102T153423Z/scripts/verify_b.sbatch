#!/bin/bash
#SBATCH --job-name=code-review
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/verify_ssx-20261102T153423Z/logs/verify_b_%j.log
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/verify_ssx-20261102T153423Z/scripts/verify_b.R
