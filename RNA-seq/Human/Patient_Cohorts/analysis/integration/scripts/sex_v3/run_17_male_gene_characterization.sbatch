#!/bin/bash
#SBATCH --job-name=v6_8A
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=logs/M17_%j.out
#SBATCH --error=logs/M17_%j.err
# ---------------------------------------------------------------------------
# Pillar 8A wrapper -- characterize Male-biased genes from v5 classifier.
# Light-weight: read v5 csv + external xrefs, write characterization csv.
# ---------------------------------------------------------------------------
set -e

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "=== Pillar 8A start: $(date) ==="
echo "Host: $(hostname)  JID: ${SLURM_JOB_ID:-NA}"

Rscript 17_male_gene_characterization.R

echo "=== Pillar 8A end: $(date) ==="
