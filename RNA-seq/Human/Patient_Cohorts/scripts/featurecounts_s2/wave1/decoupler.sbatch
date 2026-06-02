#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=decoupleR
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/decoupleR_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/decoupleR_%j.err
# Wave1 [C] TERMINAL: 51 decoupleR functional activity + 51b CollecTRI sensitivity. Reads -s2 dream.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== 51 decoupleR functional activity | $(date) ==="
Rscript RNA-seq/51_decoupler_functional_activity.R
if [ $? -ne 0 ]; then echo "FAILED: 51_decoupler_functional_activity"; exit 1; fi
echo "=== 51b decoupleR CollecTRI sensitivity | $(date) ==="
Rscript RNA-seq/51b_decoupler_collectri_sensitivity.R
if [ $? -ne 0 ]; then echo "FAILED: 51b_decoupler_collectri_sensitivity"; exit 1; fi
echo "=== done $(date) ==="
