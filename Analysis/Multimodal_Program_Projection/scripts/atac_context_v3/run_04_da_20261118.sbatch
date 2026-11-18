#!/usr/bin/env bash
#SBATCH --job-name=ws4a_t5_atacda
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:30:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/ws4a_t5_atacda_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/ws4a_t5_atacda_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

"$R_BIN" --vanilla Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/04_run_da_20261118.R
