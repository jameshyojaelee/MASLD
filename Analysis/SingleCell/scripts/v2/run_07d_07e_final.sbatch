#!/bin/bash
#SBATCH --job-name=07d_07e_final
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/07d_07e_final_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/07d_07e_final_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "[$(date)] 07d_chain_v3_final"
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/v2/07d_chain_v3_final.R
echo "[$(date)] 07e_figures_v3"
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/v2/07e_figures_v3.R
echo "[$(date)] DONE"
