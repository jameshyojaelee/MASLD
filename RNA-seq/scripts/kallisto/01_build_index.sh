#!/bin/bash
#SBATCH --job-name=B1_kallisto_index
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/index_%j.out
#SBATCH --error=logs/index_%j.err

# B1 - Build kallisto index from GENCODE v49 transcriptome
# Reused across all 5 cohorts. Writes worktree-local kidx to avoid touching shared ref dir.

set -euo pipefail
module load kallisto/0.51.1

TX_FA=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.transcripts.fa.gz
IDX_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/RNA-seq/results/kallisto/_index
mkdir -p "$IDX_DIR"
IDX="$IDX_DIR/gencode.v49.kidx"

if [[ -s "$IDX" ]]; then
  echo "[$(date)] Index already exists at $IDX -- skipping build"
  kallisto inspect "$IDX" || true
  exit 0
fi

echo "[$(date)] Building kallisto index: $IDX"
kallisto index -i "$IDX" -k 31 "$TX_FA"
echo "[$(date)] Index built. Inspecting:"
kallisto inspect "$IDX"
