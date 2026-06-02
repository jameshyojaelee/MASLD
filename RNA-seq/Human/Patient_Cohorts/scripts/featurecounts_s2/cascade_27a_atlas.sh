#!/bin/bash
#SBATCH --job-name=atlas_27a
#SBATCH --partition=cpu,io,bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_27a_atlas_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_27a_atlas_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

echo "=== 27a: Assemble Evidence Atlas (post -s 2 recount) ==="
echo "Started: $(date)"

Rscript 27a_assemble_evidence_atlas.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 27a_assemble_evidence_atlas.R"; exit 1; fi

echo "Finished: $(date)"
