#!/bin/bash
#SBATCH --job-name=mouse_rebuild
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=48:00:00
#SBATCH --output=logs/mouse_rebuild_%j.out
#SBATCH --error=logs/mouse_rebuild_%j.err

set -eo pipefail   # NOTE: do not add -u; micromamba activate hooks reference unbound vars

SCRIPTS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/scripts"
cd "$SCRIPTS"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=== M02: Per-Diet DE (with corrected HFD datasets) ==="
Rscript "$SCRIPTS/M02_mouse_per_diet_de.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: M02"; exit 1; fi

echo "=== M03: Meta-Analysis + Dream (re-enabled) ==="
Rscript "$SCRIPTS/M03_mouse_meta_analysis.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: M03"; exit 1; fi

echo "=== M04: Consensus DEGs ==="
Rscript "$SCRIPTS/M04_mouse_consensus.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: M04"; exit 1; fi

echo "=== M05: Volcano Plots ==="
Rscript "$SCRIPTS/M05_mouse_volcano_plots.R" 2>&1

echo "=== M06: Batch Correction UMAP ==="
Rscript "$SCRIPTS/M06_mouse_batch_correction_umap.R" 2>&1

echo "=== M07: Pathway Analysis ==="
Rscript "$SCRIPTS/M07_mouse_pathway_analysis.R" 2>&1

echo "Mouse pipeline complete: $(date)"
