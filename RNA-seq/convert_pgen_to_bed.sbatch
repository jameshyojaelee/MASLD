#!/bin/bash
#SBATCH --job-name=pgen2bed
#SBATCH --partition=cpu
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/pgen2bed_%j.out
#SBATCH --error=RNA-seq/logs/pgen2bed_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
GENO_DIR="${BASE}/data/1kg_eur"
BED_DIR="${GENO_DIR}/bed"
PLINK2="/nfs/sw/easybuild/software/PLINK/2.0a5.13/bin/plink"

mkdir -p "${BED_DIR}"

echo "=== Converting pgen → bed for all chromosomes ==="
echo "Start: $(date)"

for CHR in $(seq 1 22); do
    BED_PREFIX="${BED_DIR}/chr${CHR}_eur"
    if [[ -f "${BED_PREFIX}.bed" ]]; then
        echo "chr${CHR}: already exists, skipping"
        continue
    fi

    PGEN_PREFIX="${GENO_DIR}/chr${CHR}_eur"
    if [[ ! -f "${PGEN_PREFIX}.pgen" ]]; then
        echo "chr${CHR}: no pgen file, skipping"
        continue
    fi

    echo "chr${CHR}: converting..."
    "${PLINK2}" --pfile "${PGEN_PREFIX}" \
        --make-bed \
        --out "${BED_PREFIX}" \
        --threads 4 \
        --memory 12000 2>&1 | tail -3
    echo "chr${CHR}: done ($(wc -l < "${BED_PREFIX}.bim") variants)"
done

echo ""
echo "=== All conversions complete: $(date) ==="
ls -lh "${BED_DIR}/"*.bed 2>/dev/null | wc -l
echo "BED files created"
