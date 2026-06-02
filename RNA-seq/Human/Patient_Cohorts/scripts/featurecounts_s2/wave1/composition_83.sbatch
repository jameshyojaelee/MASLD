#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=TOAST
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/composition_83_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/composition_83_%j.err
# Wave1 [F]: 83 composition shifts. Reads MuSiC proportions + -s2 meta_matched.rds. Independent. NOT a feeder.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== 83 composition shifts | $(date) ==="
Rscript RNA-seq/83_composition_shifts.R
if [ $? -ne 0 ]; then echo "FAILED: 83_composition_shifts"; exit 1; fi
echo "=== done $(date) ==="
