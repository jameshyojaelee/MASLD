#!/bin/bash
#SBATCH --job-name=jl_cf02
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=6:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/control_free_transfer/cft-v1-20260828T140926Z/slurm-conv-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/control_free_transfer/cft-v1-20260828T140926Z/slurm-conv-%j.err
set -euo pipefail
module load R/4.4.1-gfbf-2023b
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export CFT_OUT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/control_free_transfer/cft-v1-20260828T140926Z
Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/control_free_transfer/02_counting_convention_impact.R
echo "EXIT_OK"
