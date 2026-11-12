#!/bin/bash
#SBATCH --job-name=mqtl2_ggt
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/MQTL2-20261111T000008Z/logs/mqtl2_ggt_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/MQTL2-20261111T000008Z/logs/mqtl2_ggt_%j.err
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/MQTL2-20261111T000008Z
module load PLINK/1.9b_6.21-x86_64
export PATH=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin:$PATH
sha256sum $RUN/scripts/04_ggt_gckr_control.R >> $RUN/out_23189934/script_and_input_sha256.txt
Rscript $RUN/scripts/04_ggt_gckr_control.R $RUN/out_23189934 $RUN/inputs/alt_loci_169_plus_ctrl.tsv
