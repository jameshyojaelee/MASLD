#!/bin/bash
#SBATCH --job-name=scdrs_nu25
#SBATCH --partition=cpu
#SBATCH --mem=160G
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --array=0-199%10
#SBATCH --output=scripts/figures/logs/scdrs_nu25_%A_%a.out
#SBATCH --error=scripts/figures/logs/scdrs_nu25_%A_%a.err

# Score 200 hep-expression-matched nulls on the 25% subsample atlas.
# cpu partition allows 10 concurrent under nslab QoS for fast wall time.

set -eo pipefail
mkdir -p scripts/figures/logs
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial

export SCDRS_GS=null
export SCDRS_TRAIT_INDEX=$SLURM_ARRAY_TASK_ID
python scripts/figures/run_scdrs_score_25pct.py
