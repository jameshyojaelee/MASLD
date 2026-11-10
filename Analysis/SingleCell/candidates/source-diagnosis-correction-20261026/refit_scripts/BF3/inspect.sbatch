#!/bin/bash
#SBATCH --job-name=bf3_inspect
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF3-20261023T154841Z/logs/%x_%j.out
export PYTHONNOUSERSITE=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF3-20261023T154841Z/scripts/inspect_atlas.py
