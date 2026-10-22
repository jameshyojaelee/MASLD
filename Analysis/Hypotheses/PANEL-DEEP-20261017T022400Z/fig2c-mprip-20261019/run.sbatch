#!/bin/bash
#SBATCH --job-name=1019
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=mprip_%j.out
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript mprip.R
