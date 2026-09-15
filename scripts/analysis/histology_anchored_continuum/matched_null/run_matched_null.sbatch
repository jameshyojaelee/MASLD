#!/bin/bash
#SBATCH --job-name=jl_an01
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/matched_null/mnc-fib-v1-run-20260828T135030Z/slurm-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/matched_null/mnc-fib-v1-run-20260828T135030Z/slurm-%j.err
set -euo pipefail
module load R/4.4.1-gfbf-2023b
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export MNC_OUT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/matched_null/mnc-fib-v1-run-20260828T135030Z
Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/matched_null/02_run_matched_null.R
pip freeze > /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/matched_null/mnc-fib-v1-run-20260828T135030Z/pip_freeze.txt 2>/dev/null || python3 -m pip freeze > /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/matched_null/mnc-fib-v1-run-20260828T135030Z/pip_freeze.txt 2>/dev/null || echo "no python env" > /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/matched_null/mnc-fib-v1-run-20260828T135030Z/pip_freeze.txt
echo "EXIT_OK"
