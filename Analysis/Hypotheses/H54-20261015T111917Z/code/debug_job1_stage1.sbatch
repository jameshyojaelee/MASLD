#!/bin/bash
# Debug run: composition into debug/, then stage 1 with draw counts divided by 20. Not the real run.
#SBATCH --job-name=5411
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z/logs/%x_%j.err
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z
D="${RUN}/debug/run_${SLURM_JOB_ID}"
mkdir "${D}"
if [ -z "${H54_COMP_DIR:-}" ]; then
  H54_COMP_OUT="${D}/composition" bash "${RUN}/code/job1_composition.sbatch"
  export H54_COMP_DIR="${D}/composition"
fi
SLURM_JOB_ID="${SLURM_JOB_ID}_s1" H54_S1_OUT="${D}/stage1" H54_DEBUG_DIVISOR=20 bash "${RUN}/code/job2_stage1.sbatch"
