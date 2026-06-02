#!/bin/bash
#SBATCH --job-name=atlas_rebuild_v10
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/multi_evidence/logs/atlas_rebuild_v10_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/multi_evidence/logs/atlas_rebuild_v10_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
mkdir -p "$BASE/RNA-seq/results/multi_evidence/logs"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "============================================================"
echo "  ATLAS REBUILD v10 (GWAS Portfolio Expansion)"
echo "  Started: $(date)"
echo "  New sources: Ghouse Cirrhosis + HCC, Pan-UKBB AFR/CSA,"
echo "               GLGC TG, MAGIC HOMA-IR"
echo "============================================================"

cd "$BASE/RNA-seq"

# Step 1: Rebuild atlas with all new COLOC/TWAS layers
echo "--- Atlas rebuild (27a) ---"
Rscript 27a_assemble_evidence_atlas.R 2>&1

# Step 2: Cross-ancestry replication update
echo ""
echo "--- Cross-ancestry replication (48) ---"
Rscript 48_cross_ancestry_replication.R 2>&1

echo ""
echo "============================================================"
echo "  ATLAS v10 REBUILD COMPLETE"
echo "  Finished: $(date)"
echo "============================================================"
