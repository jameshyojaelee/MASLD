#!/bin/bash
#SBATCH --job-name=figS_disease_sig
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=logs/figS_disease_sig_%j.out
#SBATCH --error=logs/figS_disease_sig_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs figures/figS_disease_signatures

echo "=== FigS Disease Signatures ==="
echo "Start: $(date)"

micromamba run -n rapids_singlecell python scripts/figures/figS_disease_signatures.py

echo "=== Done: $(date) ==="
ls -lh figures/figS_disease_signatures/
