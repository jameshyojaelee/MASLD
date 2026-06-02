#!/bin/bash
#SBATCH --job-name=B5b_hormone_tf_clean
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/B5b_hormone_tf_clean_%j.out
#SBATCH --error=logs/B5b_hormone_tf_clean_%j.err

# Agent: B5b (C5 critique fix -- non-circular hormone TF x COLOC)
# Removes disease_regulons.csv from target sets (was tautological).
# Adds permutation null, restricted universe, negative-control TFs.
# Chain: 217a_clean -> 217b_clean -> 217c_clean -> 217d_clean
#
# Outputs:
#   RNA-seq/results/stratified_causal/hormone_tf_sex_coloc_clean.csv
#   RNA-seq/results/stratified_causal/sex_hormone_panel_expanded_clean.csv
#   RNA-seq/results/stratified_causal/motifbreakr_sex_biased_coloc_variants_clean.csv
#   RNA-seq/results/stratified_causal/hormone_tf_negative_controls.csv
#   RNA-seq/results/stratified_causal/B5_vs_B5b_comparison.csv (headline diagnostic)

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p logs

# Activate environment (set +u around micromamba activate per ADDR2LINE bug)
export MAMBA_NO_BANNER=1
set +u
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)" 2>/dev/null || \
  eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export N_PERM=1000

echo "[$(date '+%F %T')] === B5b launch (non-circular sensitivity) ==="
echo "Host: $(hostname)  CPUs: ${SLURM_CPUS_PER_TASK:-?}  Job: ${SLURM_JOB_ID:-NA}"

echo "[$(date '+%F %T')] --- 217a_clean: build TF target sets (no disease_regulons) ---"
Rscript 217a_hormone_tf_targetsets_clean.R

echo "[$(date '+%F %T')] --- 217b_clean: Fisher + permutation null + B5 comparison ---"
Rscript 217b_hormone_tf_coloc_enrichment_clean.R

echo "[$(date '+%F %T')] --- 217c_clean: motifbreakR x sex-biased COLOC (+neg ctrls) ---"
Rscript 217c_motifbreakr_sex_biased_clean.R

echo "[$(date '+%F %T')] --- 217d_clean: expanded panel + dual_role_flag ---"
Rscript 217d_hormone_panel_annotation_clean.R

echo "[$(date '+%F %T')] === B5b complete ==="
ls -lh /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/stratified_causal/hormone_tf*_clean.csv \
       /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/stratified_causal/sex_hormone_panel_expanded_clean.csv \
       /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/stratified_causal/motifbreakr_sex_biased*_clean.csv \
       /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/stratified_causal/hormone_tf_negative_controls.csv \
       /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/stratified_causal/B5_vs_B5b_comparison.csv 2>&1 | head -20
