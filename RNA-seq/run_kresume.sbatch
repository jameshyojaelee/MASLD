#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=4 --mem=32G --time=12:00:00
#SBATCH --job-name=masld_kresume --output=RNA-seq/logs/masld_kresume_%j.out --error=RNA-seq/logs/masld_kresume_%j.err
set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== 217 ==="
Rscript RNA-seq/217_stratified_causal_atlas.R
echo "=== 93 sanity ==="
Rscript RNA-seq/93_kprogram_sanity.R
echo "=== figS08 ==="
Rscript scripts/figures/figS08_nmf_subtyping.R
echo "=== DONE ==="
