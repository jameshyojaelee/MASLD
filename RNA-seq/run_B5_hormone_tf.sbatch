#!/bin/bash
#SBATCH --job-name=B5_hormone_tf
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/B5_hormone_tf_%j.out
#SBATCH --error=logs/B5_hormone_tf_%j.err

# Agent: B5 (Sex-Hormone TF x COLOC) -- HEADLINE A7
# Chain: 217a -> 217b -> 217c -> 217d
# Outputs: RNA-seq/results/stratified_causal/{hormone_tf_sex_coloc,
#          sex_hormone_panel_expanded, motifbreakr_sex_biased_coloc_variants}.csv

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p logs

# Activate environment via micromamba (per CLAUDE.md)
# rnaseq env's binutils activate hook references unbound ADDR2LINE; relax `-u` around activate.
export MAMBA_NO_BANNER=1
set +u
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)" 2>/dev/null || \
  eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "[$(date '+%F %T')] === B5 launch ==="
echo "Host: $(hostname)  CPUs: ${SLURM_CPUS_PER_TASK:-?}  Job: ${SLURM_JOB_ID:-NA}"

echo "[$(date '+%F %T')] --- 217a: build TF target sets ---"
Rscript 217a_hormone_tf_targetsets.R

echo "[$(date '+%F %T')] --- 217b: Fisher + log-linear enrichment ---"
Rscript 217b_hormone_tf_coloc_enrichment.R

echo "[$(date '+%F %T')] --- 217c: motifbreakR cross-tab ---"
Rscript 217c_motifbreakr_sex_biased.R

echo "[$(date '+%F %T')] --- 217d: expanded panel annotation ---"
Rscript 217d_hormone_panel_annotation.R

echo "[$(date '+%F %T')] === B5 complete ==="
ls -lh /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/stratified_causal/hormone_tf*.csv \
       /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/stratified_causal/sex_hormone_panel_expanded.csv \
       /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/stratified_causal/motifbreakr_sex_biased*.csv 2>&1 | head -20
