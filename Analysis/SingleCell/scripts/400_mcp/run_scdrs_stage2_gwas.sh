#!/bin/bash
#SBATCH --job-name=scdrs_s2_gwas
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=700G
#SBATCH --cpus-per-task=4
#SBATCH --time=90:00:00
#SBATCH --array=0-16
#SBATCH --output=logs/scdrs_stage2_gwas_%A_%a.out
#SBATCH --error=logs/scdrs_stage2_gwas_%A_%a.err

# bigmem required: scDRS score_cell needs ~470GB peak (above cpu's 210GB cap).
# Per partition priority cpu→bigmem→io→dev: cpu doesn't fit, so bigmem is the
# correct first choice for this workload.

set -eo pipefail
mkdir -p logs

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

set -u
export SCDRS_TRAIT_INDEX=$SLURM_ARRAY_TASK_ID
python Analysis/SingleCell/scripts/400_mcp/run_scdrs_stage2_score_gwas.py
