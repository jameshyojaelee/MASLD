#!/bin/bash
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=dream
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/drug_cmap_v2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/drug_cmap_v2_%j.err
# Wave1 [B] TERMINAL: 20_drug_repurposing_cmap (LINCS/CGP reversal). Reads canonical -s2 dream.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export SLURM_CPUS=${SLURM_CPUS_PER_TASK:-8}
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p results/causal_inference/logs results/drug_repurposing
echo "=== drug cmap v2 | $(date) ==="
Rscript 20_drug_repurposing_cmap.R
if [ $? -ne 0 ]; then echo "FAILED: 20_drug_repurposing_cmap"; exit 1; fi
echo "=== done $(date) ==="
