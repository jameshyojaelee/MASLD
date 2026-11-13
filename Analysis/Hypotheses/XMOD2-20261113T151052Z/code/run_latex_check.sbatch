#!/bin/bash
#SBATCH --job-name=XMOD2tex
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/XMOD2-20261113T151052Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/XMOD2-20261113T151052Z/logs/%j.err
set -euo pipefail
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/XMOD2-20261113T151052Z/latex_build_${SLURM_JOB_ID}
mkdir "$OUT"
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/docs/manuscript/latex
T=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/tectonic
$T -X compile main.tex --outdir "$OUT" > "$OUT/main-console.txt" 2>&1
$T -X compile supplement.tex --outdir "$OUT" > "$OUT/supplement-console.txt" 2>&1
