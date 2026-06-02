#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --partition=cpu
#SBATCH --array=1-22
#SBATCH --job-name=build_eas_ld
#SBATCH --output=logs/build_eas_ld_%A_%a.out
#SBATCH --error=logs/build_eas_ld_%A_%a.err

# build_eas_ld.sh
# Builds EAS LD matrices from 1000G Phase 3 VCFs in the same block structure
# as the UKBB EUR LD matrices (approx_LD_blocks.txt)
# One array task per chromosome (--array=1-22)

set -euo pipefail

module load PLINK/1.9

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"

CHR=${SLURM_ARRAY_TASK_ID}
echo "============================================================"
echo "Building EAS LD matrices for chr${CHR}"
echo "============================================================"

# Paths
VCF_DIR="/gpfs/commons/groups/sanjana_lab/jmorris/1000Genomes"
VCF="${VCF_DIR}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
BLOCKS_FILE="data/ld_ref/ukbb_eur/approx_LD_blocks.txt"
EAS_DIR="data/ld_ref/1kg_eas"
EAS_SAMPLES="${EAS_DIR}/eas_samples.txt"

# 1000G EAS population IDs (from integrated_call_samples_v3.20130502.ALL.panel)
# EAS super-population: CHB, JPT, CHS, CDX, KHV
PANEL_FILE="${VCF_DIR}/integrated_call_samples_v3.20130502.ALL.panel"

# Create EAS sample list if it doesn't exist
if [ ! -f "${EAS_SAMPLES}" ]; then
    echo "Creating EAS sample list..."
    mkdir -p "${EAS_DIR}"
    if [ -f "${PANEL_FILE}" ]; then
        awk '$3 == "EAS" {print $1, $1}' "${PANEL_FILE}" > "${EAS_SAMPLES}"
    else
        # Fallback: use known EAS population codes
        # Download panel file if needed
        echo "WARNING: Panel file not found at ${PANEL_FILE}"
        echo "Attempting to extract EAS samples from VCF header..."
        # This is a fallback — ideally the panel file should exist
        bcftools query -l "${VCF}" | head -504 > "${EAS_SAMPLES}.tmp"
        # For now, create placeholder — user should provide proper EAS sample list
        echo "ERROR: Cannot determine EAS samples without panel file"
        echo "Please provide EAS sample IDs in ${EAS_SAMPLES} (one per line, two columns: FID IID)"
        exit 1
    fi
    N_EAS=$(wc -l < "${EAS_SAMPLES}")
    echo "  Found ${N_EAS} EAS samples"
fi

if [ ! -f "${VCF}" ]; then
    echo "ERROR: VCF not found: ${VCF}"
    exit 1
fi

# First convert VCF to PLINK format for this chromosome (EAS only)
PLINK_PREFIX="${EAS_DIR}/chr${CHR}_eas"
if [ ! -f "${PLINK_PREFIX}.bed" ]; then
    echo "Converting VCF to PLINK format (EAS samples only)..."
    plink --vcf "${VCF}" \
        --keep "${EAS_SAMPLES}" \
        --make-bed \
        --snps-only just-acgt \
        --maf 0.01 \
        --geno 0.05 \
        --out "${PLINK_PREFIX}" \
        --allow-extra-chr \
        --chr "${CHR}" \
        2>/dev/null
    echo "  PLINK BED created: $(wc -l < ${PLINK_PREFIX}.bim) variants"
fi

# Process each LD block for this chromosome
mkdir -p "${EAS_DIR}/chr${CHR}"

N_BLOCKS=0
N_SUCCESS=0

while read -r block_chr block_start block_stop; do
    # Skip header and other chromosomes
    if [ "${block_chr}" != "${CHR}" ] || [ "${block_chr}" == "chr" ]; then
        continue
    fi

    N_BLOCKS=$((N_BLOCKS + 1))
    BLOCK_DIR="${EAS_DIR}/chr${CHR}/${block_start}.${block_stop}"
    BLOCK_PREFIX="${BLOCK_DIR}/${block_start}.${block_stop}"

    # Skip if already computed
    if [ -f "${BLOCK_PREFIX}.ld" ] && [ -f "${BLOCK_PREFIX}.bim" ]; then
        N_SUCCESS=$((N_SUCCESS + 1))
        continue
    fi

    mkdir -p "${BLOCK_DIR}"

    # Extract variants in this block and compute LD
    plink --bfile "${PLINK_PREFIX}" \
        --chr "${CHR}" \
        --from-bp "${block_start}" \
        --to-bp "${block_stop}" \
        --r square \
        --make-just-bim \
        --out "${BLOCK_PREFIX}" \
        --allow-extra-chr \
        2>/dev/null

    if [ -f "${BLOCK_PREFIX}.ld" ] && [ -f "${BLOCK_PREFIX}.bim" ]; then
        N_SUCCESS=$((N_SUCCESS + 1))
        N_VARS=$(wc -l < "${BLOCK_PREFIX}.bim")
        if [ $((N_BLOCKS % 10)) -eq 0 ]; then
            echo "  Block ${N_BLOCKS}: ${block_start}-${block_stop} (${N_VARS} variants)"
        fi
    else
        echo "  WARNING: Failed to compute LD for block ${block_start}-${block_stop}"
    fi
done < "${BLOCKS_FILE}"

echo ""
echo "chr${CHR}: Completed ${N_SUCCESS}/${N_BLOCKS} LD blocks"

# Copy the approx_LD_blocks.txt to EAS dir for consistency
if [ ! -f "${EAS_DIR}/approx_LD_blocks.txt" ]; then
    cp "${BLOCKS_FILE}" "${EAS_DIR}/approx_LD_blocks.txt"
fi

echo "============================================================"
echo "Done building EAS LD for chr${CHR}"
echo "============================================================"
