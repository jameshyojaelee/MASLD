#!/bin/bash
#SBATCH --job-name=scripts_14_15_220
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=logs/scripts_14_15_220_%j.out
#SBATCH --error=logs/scripts_14_15_220_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=== Script 14: Fibrosis Progression DE (pass_technical fix) ==="
echo "Started: $(date)"
Rscript 14_fibrosis_progression_de.R 2>&1
echo "Completed 14: $(date)"

echo "=== Script 15: NAS Component DE (pass_technical fix) ==="
echo "Started: $(date)"
Rscript 15_nas_component_de.R 2>&1
echo "Completed 15: $(date)"

echo "=== Script 220: DIA-MS Processing (per-protein imputation fix) ==="
echo "Started: $(date)"
Rscript 220_process_pxd052937.R 2>&1
echo "Completed 220: $(date)"

echo "All done: $(date)"
