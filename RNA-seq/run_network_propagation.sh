#!/bin/bash
#SBATCH --job-name=network_propagation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=RNA-seq/results/multi_evidence/logs/network_propagation_%j.out

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/results/multi_evidence/logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rapids_singlecell

python RNA-seq/46c_network_propagation.py
