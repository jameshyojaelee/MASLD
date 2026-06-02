#!/bin/bash
#SBATCH --job-name=liver_atlas_figs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/atlas_figs_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/atlas_figs_%j.err

set -euo pipefail

echo "=== Atlas Assembly + Figure Generation ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "--- Step 1: Re-running 27a (atlas assembly with UKBB ALT COLOC) ---"
Rscript RNA-seq/27a_assemble_evidence_atlas.R

echo ""
echo "--- Step 2: Re-running 27b (benchmark presets) ---"
Rscript RNA-seq/27b_benchmark_presets.R

echo ""
echo "--- Step 3: Generating Fig 5 (causal architecture) ---"
Rscript scripts/figures/fig5_causal_architecture.R

echo ""
echo "--- Step 4: Generating Fig 7 (multi-evidence) ---"
Rscript scripts/figures/fig7_multi_evidence.R

echo ""
echo "=== Complete ==="
echo "End: $(date)"
