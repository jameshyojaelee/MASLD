#!/bin/bash
#SBATCH --job-name=liver_9c_downstream
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/9cohort_downstream_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/9cohort_downstream_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00

# Downstream re-runs after 9-cohort dream: atlas rebuild + cross-species + figures
# Submit with: sbatch --dependency=afterok:$CORE_JOB scripts/run_9cohort_downstream.sh

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$MASLD_PROJECT_ROOT"
mkdir -p logs figures

echo "=========================================="
echo "  9-COHORT DOWNSTREAM REBUILD"
echo "  Job: $SLURM_JOB_ID | CPUs: $SLURM_CPUS_PER_TASK"
echo "  Started: $(date)"
echo "=========================================="

# Cross-species concordance (reads canonical_deg_results.csv)
echo ""
echo "=== Cross-Species Concordance (Phase 1: Gene-level) ==="
Rscript Analysis/Cross_Species_Concordance/scripts/01_corrected_gene_concordance.R 2>&1
if [ $? -ne 0 ]; then echo "WARNING: Cross-species concordance failed (non-fatal)"; fi

# Multi-evidence atlas rebuild
echo ""
echo "=== 27a: Assemble Evidence Atlas ==="
Rscript RNA-seq/27a_assemble_evidence_atlas.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 27a_assemble_evidence_atlas.R"; exit 1; fi

echo ""
echo "=== 27b: Benchmark Presets ==="
Rscript RNA-seq/27b_benchmark_presets.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 27b_benchmark_presets.R"; exit 1; fi

# Drug repurposing (reads dream t-statistics)
echo ""
echo "=== 20: Drug Repurposing ==="
Rscript RNA-seq/20_drug_repurposing_v3.R 2>&1
if [ $? -ne 0 ]; then echo "WARNING: Drug repurposing failed (non-fatal)"; fi

# GWAS overlay (reads canonical_deg_results.csv)
echo ""
echo "=== 30: GWAS Overlay ==="
Rscript RNA-seq/30_gwas_overlay_v2.R 2>&1
if [ $? -ne 0 ]; then echo "WARNING: GWAS overlay failed (non-fatal)"; fi

# Regenerate Figure 1
echo ""
echo "=== Figure 1: Atlas Overview ==="
Rscript scripts/figures/fig1_atlas_overview.R 2>&1
if [ $? -ne 0 ]; then echo "WARNING: Figure 1 failed (non-fatal)"; fi

echo ""
echo "=========================================="
echo "  Downstream rebuild complete: $(date)"
echo "=========================================="
