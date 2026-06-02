#!/bin/bash
#SBATCH --job-name=141_atac_progression
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression/141_atac_%j.out
#SBATCH --error=RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression/141_atac_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$(pwd)"
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2

SCRIPT="RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/141_atac_masl_vs_mash.py"
OUTDIR="RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression"

echo "=== 141: ATAC-seq MASL-vs-MASH DA ==="
echo "Started: $(date)"
echo "Node: $(hostname)"
echo ""

python "${SCRIPT}" \
    --input Analysis/ATAC/Human_Multiome/results/label_transfer/snapatac2_label_transferred.h5ad \
    --per-donor-dir Analysis/ATAC/Human_Multiome/results/snapatac2/per_donor \
    --output-dir "${OUTDIR}" \
    --gtf /gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz \
    --c2-degs RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/nafl_vs_nash_dream.csv

echo ""
echo "=== Done: $(date) ==="
echo "Output files:"
ls -lh "${OUTDIR}/atac_masl_vs_mash_da.csv" 2>/dev/null || echo "  (MASL-vs-MASH not generated)"
ls -lh "${OUTDIR}/atac_masl_vs_normal_da.csv" 2>/dev/null || echo "  (MASL-vs-NORMAL not generated)"
ls -lh "${OUTDIR}/atac_progression_overlap.csv" 2>/dev/null || echo "  (overlap not generated)"
