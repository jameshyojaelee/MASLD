#!/bin/bash
#SBATCH --job-name=qc
#SBATCH --partition=cpu,io,bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/s2_01_qc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/s2_01_qc_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "=== 01: Sample QC (post -s 2 recount) ==="
echo "Started: $(date)"

Rscript analysis/integration/scripts/01_sample_qc.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 01_sample_qc.R"; exit 1; fi

echo "Finished: $(date)"
