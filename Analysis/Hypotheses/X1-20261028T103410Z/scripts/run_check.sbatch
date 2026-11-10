#!/bin/bash
#SBATCH --job-name=x1-labelcheck
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/X1-20261028T103410Z
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I scripts/check_labels.py \
  /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V7-20261021T235717Z/downloads_untrusted/scp3359/SCP3359 \
  out DNASE1L3,FCN2,CD36,PECAM1,VWF,PDGFRB,PDGFRA,LYVE1,HGF,ACTA2,MZB1,HPX
