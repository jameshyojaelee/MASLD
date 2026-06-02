#!/bin/bash
# 62_otters_train_weights.sh
# ---------------------------------------------------------------------------
# Phase 2A Step 2: Train OTTERS TWAS weights from Broadaway summary eQTL
#
# Uses OTTERS (daiqile96/OTTERS) to train P+T and lassosum weights.
# SDPR and PRS-CS are optional (require separate compilation).
#
# Input:
#   - data/broadaway_eqtl/otters_format/chr{N}_broadaway.txt  (from Script 61)
#   - data/broadaway_eqtl/otters_format/gene_anno.txt
#   - data/1kg_eur/chr{N}_eur.{pgen,psam,pvar}  (LD reference)
#
# Output:
#   - data/broadaway_eqtl/otters_weights/chr{N}/
#
# Usage (SLURM array):
#   sbatch --array=1-22 RNA-seq/run_otters_train.sbatch
# ---------------------------------------------------------------------------

set -euo pipefail

# Configuration
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OTTERS_DIR="${BASE}/tools/OTTERS"
BROADAWAY_OTTERS="${BASE}/data/broadaway_eqtl/otters_format"
GENO_DIR="${BASE}/data/1kg_eur"
BED_DIR="${GENO_DIR}/bed"
OUT_DIR="${BASE}/data/broadaway_eqtl/otters_weights"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# Chromosome from SLURM array task ID
CHROM="${SLURM_ARRAY_TASK_ID:-1}"

echo "=== OTTERS Weight Training: chr${CHROM} ==="
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

# Verify inputs
SST_FILE="${BROADAWAY_OTTERS}/chr${CHROM}_broadaway.txt"
ANNO_FILE="${BROADAWAY_OTTERS}/gene_anno.txt"

if [[ ! -f "${SST_FILE}" ]]; then
    echo "ERROR: eQTL summary stats not found: ${SST_FILE}"
    echo "Run Script 61 (OTTERS format) first."
    exit 1
fi

if [[ ! -f "${ANNO_FILE}" ]]; then
    echo "ERROR: Gene annotation not found: ${ANNO_FILE}"
    exit 1
fi

# Count genes on this chromosome
N_GENES=$(awk -v c="${CHROM}" '$1 == c' "${ANNO_FILE}" | wc -l)
echo "Genes on chr${CHROM}: ${N_GENES}"
if [[ "${N_GENES}" -eq 0 ]]; then
    echo "No genes to process. Exiting."
    exit 0
fi

# Step 1: Convert pgen → bed if needed
BED_PREFIX="${BED_DIR}/chr${CHROM}_eur"
if [[ ! -f "${BED_PREFIX}.bed" ]]; then
    echo "Converting pgen → bed for chr${CHROM}..."
    mkdir -p "${BED_DIR}"
    PGEN_PREFIX="${GENO_DIR}/chr${CHROM}_eur"
    if [[ ! -f "${PGEN_PREFIX}.pgen" ]]; then
        echo "ERROR: No pgen file: ${PGEN_PREFIX}.pgen"
        exit 1
    fi
    # Use PLINK2 for conversion (full path — module names it 'plink')
    PLINK2="/nfs/sw/easybuild/software/PLINK/2.0a5.13/bin/plink"
    "${PLINK2}" --pfile "${PGEN_PREFIX}" \
           --make-bed \
           --out "${BED_PREFIX}" \
           --threads "${THREADS}" \
           --memory 8000 2>&1 | tail -5
    echo "BED conversion complete."
fi

# Verify BED files exist
for ext in bed bim fam; do
    if [[ ! -f "${BED_PREFIX}.${ext}" ]]; then
        echo "ERROR: Missing ${BED_PREFIX}.${ext}"
        exit 1
    fi
done
echo "BED prefix: ${BED_PREFIX}"
echo "Variants in BED: $(wc -l < "${BED_PREFIX}.bim")"

# Step 2: Create output directory
CHR_OUT="${OUT_DIR}/chr${CHROM}"
mkdir -p "${CHR_OUT}"

# Step 3: Run OTTERS training
# OTTERS expects plink to be available in PATH
PLINK19="/nfs/sw/easybuild/software/PLINK/1.9b_6.21-x86_64/plink"
if [[ ! -x "${PLINK19}" ]]; then
    echo "ERROR: PLINK 1.9 not found at ${PLINK19}"
    exit 1
fi

# Create temporary symlink so OTTERS finds 'plink'
TMPBIN=$(mktemp -d)
ln -s "${PLINK19}" "${TMPBIN}/plink"
export PATH="${TMPBIN}:${PATH}"

echo ""
echo "Running OTTERS training..."
echo "  OTTERS_DIR: ${OTTERS_DIR}"
echo "  SST: ${SST_FILE}"
echo "  GENO: ${BED_PREFIX}"
echo "  OUT: ${CHR_OUT}"
echo "  THREADS: ${THREADS}"
echo ""

python "${OTTERS_DIR}/training.py" \
    --OTTERS_dir="${OTTERS_DIR}" \
    --anno_dir="${ANNO_FILE}" \
    --geno_dir="${BED_PREFIX}" \
    --sst_file="${SST_FILE}" \
    --out_dir="${CHR_OUT}" \
    --chrom="${CHROM}" \
    --models=PT,lassosum \
    --lassosum_ld_blocks=EUR.hg19 \
    --r2=0.99 \
    --window=1000000 \
    --thread="${THREADS}" 2>&1

# Cleanup
rm -rf "${TMPBIN}"

# Report results
echo ""
echo "=== Training Results ==="
N_WEIGHT_FILES=$(find "${CHR_OUT}" -name "*.txt" -o -name "*.csv" 2>/dev/null | wc -l)
echo "Weight files generated: ${N_WEIGHT_FILES}"
ls -lh "${CHR_OUT}/"*.txt 2>/dev/null | head -20 || echo "  (no .txt weight files)"
echo ""
echo "=== Finished chr${CHROM}: $(date) ==="
