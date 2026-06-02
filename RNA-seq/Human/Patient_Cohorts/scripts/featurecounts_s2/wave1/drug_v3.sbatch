#!/bin/bash
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=dream
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/drug_v3_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/drug_v3_%j.err
# Wave1 [B] ATLAS FEEDER: 20_drug_repurposing_v3 -> convergent_drug_targets.csv (read by 27a).
# Reads canonical -s2 dream_results.csv + consensus_degs.csv. Needs internet (CLUE/DGIdb/OpenTargets).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export SLURM_CPUS=${SLURM_CPUS_PER_TASK:-8}
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p results/drug_repurposing/logs
echo "=== drug v3 (convergent_drug_targets) | $(date) ==="
Rscript 20_drug_repurposing_v3.R
if [ $? -ne 0 ]; then echo "FAILED: 20_drug_repurposing_v3"; exit 1; fi
echo "=== done $(date) ==="
