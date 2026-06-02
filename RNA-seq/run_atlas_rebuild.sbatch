#!/bin/bash
#SBATCH --job-name=atlas_rebuild
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --time=48:00:00
#SBATCH --output=logs/atlas_rebuild_%j.out
#SBATCH --error=logs/atlas_rebuild_%j.err

set -euo pipefail

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$MASLD_PROJECT_ROOT/RNA-seq"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=== Step 1/4: Script 27a - Assemble Evidence Atlas ==="
echo "Started: $(date)"
Rscript 27a_assemble_evidence_atlas.R 2>&1
echo "Completed 27a: $(date)"
echo ""

echo "=== Step 2/4: Script 75 - Integrate Causal Overhaul ==="
echo "Started: $(date)"
Rscript 75_integrate_causal_overhaul.R 2>&1
echo "Completed 75: $(date)"
echo ""

echo "=== Step 3/4: Script 217 - Stratified Causal Atlas ==="
echo "Started: $(date)"
Rscript 217_stratified_causal_atlas.R 2>&1
echo "Completed 217: $(date)"
echo ""

echo "=== Step 4/4: Script 215 - Pharmacogenomic Targets ==="
echo "Started: $(date)"
Rscript 215_pharmacogenomic_targets.R 2>&1
echo "Completed 215: $(date)"
echo ""

echo "=== Validation ==="
Rscript -e '
library(data.table)
atlas <- fread("results/multi_evidence/multi_evidence_atlas.csv")
cat("Atlas dimensions:", nrow(atlas), "genes x", ncol(atlas), "columns\n")
cat("coloc_susie_best_pp4 max:", max(atlas$coloc_susie_best_pp4, na.rm=TRUE), "\n")
cat("coloc_susie PP4>0.5:", sum(atlas$coloc_susie_best_pp4 > 0.5, na.rm=TRUE), "\n")
cat("layers_active distribution:\n")
print(table(atlas$layers_active))
cat("l3_active TRUE:", sum(atlas$l3_active, na.rm=TRUE), "\n")
'
echo "All done: $(date)"
