#!/bin/bash
#SBATCH --job-name=pandas
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/269b_aggregate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/269b_aggregate_%j.err

set -eo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p "${BASE}/RNA-seq/logs"
cd "${BASE}/RNA-seq"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "=== 269b aggregate figure inputs start: $(date) on $(hostname) ==="
python 269b_aggregate_figure_inputs.py
echo "=== 269b end: $(date) ==="
