#!/bin/bash
#SBATCH --job-name=metafor
#SBATCH --partition=cpu,io,bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/s2_06_meta_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/s2_06_meta_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "=== 06: Meta-Analysis (post -s 2 recount) ==="
echo "Started: $(date)"

Rscript analysis/integration/scripts/06_meta_analysis.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 06_meta_analysis.R"; exit 1; fi

echo "Finished: $(date)"
