#!/bin/bash
#SBATCH --job-name=jl_cf01
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/control_free_transfer/cft-v1-20260828T141625Z/slurm-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/control_free_transfer/cft-v1-20260828T141625Z/slurm-%j.err
set -euo pipefail
module load R/4.4.1-gfbf-2023b
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export CFT_OUT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/control_free_transfer/cft-v1-20260828T141625Z
Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/control_free_transfer/01_project_control_free_cohorts.R
python3 -m pip freeze > /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/control_free_transfer/cft-v1-20260828T141625Z/pip_freeze.txt 2>/dev/null || echo "no python env" > /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/control_free_transfer/cft-v1-20260828T141625Z/pip_freeze.txt
echo "EXIT_OK"
