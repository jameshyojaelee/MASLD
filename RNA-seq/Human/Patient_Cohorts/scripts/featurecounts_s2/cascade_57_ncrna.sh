#!/bin/bash
#SBATCH --job-name=atlas_57
#SBATCH --partition=cpu,io,bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/s2_57_ncrna_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/s2_57_ncrna_%j.err
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "=== 57 ncRNA atlas integration (STAR -s2) $(date) ==="
Rscript 57_ncrna_atlas_integration.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 57_ncrna_atlas_integration"; exit 1; fi
echo "Finished: $(date)"
