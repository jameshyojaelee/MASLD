#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=16 --mem=96G --time=48:00:00
#SBATCH --job-name=masld_nmfsweep --output=RNA-seq/logs/masld_nmfsweep_%j.out --error=RNA-seq/logs/masld_nmfsweep_%j.err
set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript RNA-seq/94_nmf_nosex_ksweep.R
