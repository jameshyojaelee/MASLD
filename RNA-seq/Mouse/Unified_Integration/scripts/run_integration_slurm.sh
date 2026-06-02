#!/bin/bash
#SBATCH --job-name=liver_mouse_integration
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/integration_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/integration_%j.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00

# Unified Mouse Integration Pipeline (M00 → M04)
# Runs after Phase A recount completes (dependency set at submission time)

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration

SCRIPTS="scripts"

echo "=== Unified Mouse Integration Pipeline ==="
echo "Start: $(date)"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo ""

# M00: Harmonize metadata
echo "=== M00: Harmonize Mouse Metadata ==="
Rscript "$SCRIPTS/M00_harmonize_mouse_metadata.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: M00"; exit 1; fi
echo ""

# M01: Sample QC + count merge
echo "=== M01: Sample QC ==="
Rscript "$SCRIPTS/M01_mouse_sample_qc.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: M01"; exit 1; fi
echo ""

# M02: Per-diet limma-voom DE
echo "=== M02: Per-Diet DE ==="
Rscript "$SCRIPTS/M02_mouse_per_diet_de.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: M02"; exit 1; fi
echo ""

# M02b: Pooled diet groups (Western + HFD three-dataset pools).
# MUST run after M02 — overwrites HFD_de_results.csv with a 3-cohort pooled result
# and requires M02's HFD_de_results.csv as a pre-condition assertion.
echo "=== M02b: Pooled Diet Groups ==="
Rscript "$SCRIPTS/M02b_pooled_diet_groups.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: M02b"; exit 1; fi
echo ""

# M02c: ashr shrinkage on all 4 per-diet results (MCD/CDAHFD/Western/HFD).
# MUST run after M02b — consumes the pooled Western/HFD files M02b wrote.
echo "=== M02c: ashr shrinkage ==="
Rscript "$SCRIPTS/M02c_ashr_shrinkage.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: M02c"; exit 1; fi
echo ""

# M03: Meta-analysis + Pooled dream
echo "=== M03: Meta-Analysis + Dream ==="
Rscript "$SCRIPTS/M03_mouse_meta_analysis.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: M03"; exit 1; fi
echo ""

# M04: Consensus DEGs
echo "=== M04: Consensus DEGs ==="
Rscript "$SCRIPTS/M04_mouse_consensus.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: M04"; exit 1; fi
echo ""

echo "=== Pipeline complete: $(date) ==="
