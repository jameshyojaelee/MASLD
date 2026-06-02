#!/bin/bash
#SBATCH --job-name=v6_consensus
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/M18_%j.out
#SBATCH --error=logs/M18_%j.err
# ---------------------------------------------------------------------------
# Consensus 18 + final manifest emit. Runs afterany on all v6 pillars.
# Graceful: writes whatever is available with tier_evidence_incomplete=TRUE
# for genes missing inputs.
# ---------------------------------------------------------------------------
set -e

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "=== Consensus 18 start: $(date) ==="
Rscript 18_consensus_v6.R
echo "--- Manifest emit ---"
Rscript v6_manifest_emit.R
echo "=== Consensus 18 end: $(date) ==="
