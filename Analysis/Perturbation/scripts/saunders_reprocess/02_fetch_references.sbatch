#!/bin/bash
#SBATCH --job-name=prefetch
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Perturbation/scripts/saunders_reprocess/logs/02_refs_%j.log

# Fetch 10x Genomics Flex references for mouse:
# 1. refdata-gex-GRCm39-2024-A  (cellranger 8 transcriptome reference)
# 2. Chromium_Mouse_Transcriptome_Probe_Set_v1.1.0_GRCm39-2024-A.csv  (Flex probe set)
#
# The paper used cellranger 8.0.0 / "Assembly: mm10". For Flex on cellranger 8,
# the canonical mouse reference is GRCm39-2024-A (10x has dropped mm10 for Flex 8+).
# We'll attempt GRCm39 first; if that's mismatched, we fall back to mm10-2020-A
# (which has matching Flex probe set v1.0 for cellranger 7).

set -euo pipefail

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
REF_DIR=$PROJ/data/perturbation/datasets/saunders2025/refs
mkdir -p "$REF_DIR"
cd "$REF_DIR"

echo "[$(date +%T)] downloading refdata-gex-GRCm39-2024-A.tar.gz"
if [ ! -d "refdata-gex-GRCm39-2024-A" ]; then
    curl -sSL -o refdata-gex-GRCm39-2024-A.tar.gz \
        "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
    echo "[$(date +%T)] extracting..."
    tar -xzf refdata-gex-GRCm39-2024-A.tar.gz
    rm refdata-gex-GRCm39-2024-A.tar.gz
fi

echo "[$(date +%T)] downloading Flex Mouse Transcriptome Probe Set v1.1.0"
curl -sSL -o Chromium_Mouse_Transcriptome_Probe_Set_v1.1.0_GRCm39-2024-A.csv \
    "https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Mouse_Transcriptome_Probe_Set_v1.1.0_GRCm39-2024-A.csv"

# Also fetch the 10x Flex cell barcode whitelists (737K-fixed-rna-profiling)
# Cell Ranger 8 includes its own whitelists; bundled in CELLRANGER_DIR/lib/python/cellranger/barcodes/.
# We don't need to re-fetch, but let's note its path:
echo "[$(date +%T)] cellranger whitelist locations (from module-loaded cellranger):"
module load CellRanger/8.0.1
CR_DIR=$(dirname $(dirname $(which cellranger)))
echo "  $CR_DIR/lib/python/cellranger/barcodes/"
ls "$CR_DIR/lib/python/cellranger/barcodes/" 2>/dev/null | head -20 || true

echo "=== references ready in $REF_DIR ==="
ls -la "$REF_DIR"
