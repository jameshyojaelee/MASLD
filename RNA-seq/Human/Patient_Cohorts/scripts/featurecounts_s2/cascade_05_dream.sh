#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=6:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/s2_05_dream_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/s2_05_dream_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "=== 05: Dream Mega-Analysis (post -s 2 recount, 32 CPUs, 400G) ==="
echo "Started: $(date)"

Rscript analysis/integration/scripts/05_dream_mega_analysis.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 05_dream_mega_analysis.R"; exit 1; fi

echo "Finished: $(date)"
