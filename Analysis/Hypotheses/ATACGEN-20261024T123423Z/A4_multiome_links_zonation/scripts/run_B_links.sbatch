#!/bin/bash
#SBATCH --job-name=a4_links
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A4_multiome_links_zonation/logs/%x_%j.log
# A4 job B: links (02 full: K grid, Q1/Q2 calls, 100 donor bootstraps, split halves) for
# hepatocyte, fibroblast, macrophage; then Q2 agreement (03, 10,000 matched draws).
set -euo pipefail
A4=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A4_multiome_links_zonation
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$A4"
test -s work/variant_peak_cs.tsv
/usr/bin/time -v $ENV/bin/Rscript scripts/02_links.R "$A4" full hepatocyte,fibroblast,macrophage 100 14 2> logs/time_02_full.txt
tail -n 60 logs/time_02_full.txt
/usr/bin/time -v $ENV/bin/Rscript scripts/03_q2_agreement.R "$A4" hepatocyte,fibroblast,macrophage 10000 2> logs/time_03.txt
tail -n 30 logs/time_03.txt
