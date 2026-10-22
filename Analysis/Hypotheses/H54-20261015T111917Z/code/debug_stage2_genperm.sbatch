#!/bin/bash
# Debug run of stage 2 on GENE-PERMUTED data only (no real outcome). Needs H54_S1_DIR.
#SBATCH --job-name=5412
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z/logs/%x_%j.err
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z
: "${H54_S1_DIR:?H54_S1_DIR is required}"
D="${RUN}/debug/stage2_genperm_${SLURM_JOB_ID}"
mkdir "${D}"
H54_S2_OUT="${D}/stage2" H54_DEBUG_GENE_PERMUTED=1 H54_DEBUG_DIVISOR="${H54_DEBUG_DIVISOR:-10}" \
  H54_DEBUG_FORCE_OPEN_SECONDARY=1 bash "${RUN}/code/job3_stage2.sbatch"
