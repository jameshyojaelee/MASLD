#!/bin/bash
#SBATCH --job-name=consensus_07
#SBATCH --partition=cpu,io,bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/s2_07_consensus_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/s2_07_consensus_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "=== 07: Consensus DEGs (regenerate from STAR -s 2 dream_results_ashr.csv) ==="
echo "Started: $(date)"

Rscript analysis/integration/scripts/07_consensus_degs.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 07_consensus_degs.R"; exit 1; fi

echo "Finished: $(date)"
