#!/bin/bash
#SBATCH --job-name=xmod1
#SBATCH --account=nslab --qos=nslab --partition=cpu --prefer=v6
#SBATCH --cpus-per-task=2 --mem=32G --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/XMOD1-20261112T202307Z/code/slurm-%j.log
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python Analysis/Hypotheses/XMOD1-20261112T202307Z/code/xmod.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/XMOD1-20261112T202307Z/out_$SLURM_JOB_ID
