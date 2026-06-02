#!/bin/bash
#SBATCH --job-name=integrate_overhaul
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/integrate_overhaul_%j.out
#SBATCH --error=RNA-seq/logs/integrate_overhaul_%j.err

# Phase 5A: Integrate all causal inference overhaul results into atlas
set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=== Phase 5A: Integrate Causal Overhaul ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"

# Step 1: Combine OTTERS results (if not already done)
for gwas in ghodsian ukbb_alt ukbb_ast ukbb_ggt; do
  combined="${BASE}/RNA-seq/results/causal_inference/otters_broadaway/${gwas}/otters_twas_combined.csv"
  if [ ! -f "$combined" ]; then
    echo "Combining OTTERS for ${gwas}..."
    python RNA-seq/63c_combine_otters_results.py --gwas "${gwas}"
  fi
done

# Step 2: Run OTTERS concordance analysis
echo ""
echo "Running OTTERS concordance..."
Rscript RNA-seq/63b_otters_concordance.R

# Step 3: Integrate into atlas
echo ""
echo "Integrating into atlas..."
Rscript RNA-seq/75_integrate_causal_overhaul.R

# Step 4: Rebuild atlas with new L4 columns
echo ""
echo "Rebuilding atlas assembly..."
Rscript RNA-seq/27a_assemble_evidence_atlas.R

echo ""
echo "=== Done: $(date) ==="
