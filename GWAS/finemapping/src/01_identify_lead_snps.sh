#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --output=logs/lead_snps_%j.out
#SBATCH --error=logs/lead_snps_%j.err

# 01_identify_lead_snps.sh
# Identifies lead SNPs from reformatted GWAS using PLINK LD clumping
# Usage: sbatch 01_identify_lead_snps.sh
#
# For each reformatted GWAS in data/sumstats/*_reformatted_hg19.tsv:
#   1. Extract genome-wide significant variants (p < 5e-8)
#   2. LD-clump using 1000G reference panel
#   3. Output lead SNP file in CHR, BP, locus format

set -euo pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"

module load PLINK/1.9

# Reference panels for LD clumping
EUR_REF="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/1kg_eur/bed"
# EAS reference will be built separately

OUTDIR="data/lead_snps"
TMPDIR="data/lead_snps/tmp_clump"
mkdir -p "${OUTDIR}" "${TMPDIR}"

echo "============================================================"
echo "Lead SNP Identification via PLINK LD Clumping"
echo "============================================================"

for SUMSTATS in data/sumstats/*_reformatted_hg19.tsv; do
    BASENAME=$(basename "${SUMSTATS}" _reformatted_hg19.tsv)
    echo ""
    echo "--- Processing: ${BASENAME} ---"

    # Determine ancestry for LD reference
    if [[ "${BASENAME}" == BBJ_* ]]; then
        echo "  EAS ancestry — skipping PLINK clumping (no EAS reference yet)"
        echo "  Will use simple distance-based pruning instead"

        # Simple distance-based lead SNP identification for EAS
        # (no LD reference available yet — use 500kb window)
        awk -F'\t' 'NR>1 && $7 < 5e-8 {print $1, $2, $7}' "${SUMSTATS}" | \
            sort -k1,1n -k3,3g | \
            awk '{
                key = $1;
                if (!(key in best) || ($2 > last_pos[key] + 500000)) {
                    best[key] = $2;
                    last_pos[key] = $2;
                    print $1 "\t" $2 "\t" $1"."$2
                }
            }' > "${TMPDIR}/${BASENAME}_raw_leads.txt"

        # Add header and write output
        echo -e "CHR\tBP\tlocus" > "${OUTDIR}/${BASENAME}_leadSNPs.tsv"
        cat "${TMPDIR}/${BASENAME}_raw_leads.txt" >> "${OUTDIR}/${BASENAME}_leadSNPs.tsv"
        N_LEADS=$(wc -l < "${TMPDIR}/${BASENAME}_raw_leads.txt")
        echo "  Found ${N_LEADS} lead SNPs (distance-based)"
        continue
    fi

    # EUR ancestry — use PLINK LD clumping
    # Step 1: Create PLINK-compatible input (SNP P format)
    # Generate BOTH allele orderings so we match regardless of strand
    echo "  Preparing PLINK input (both allele orderings)..."
    awk -F'\t' 'NR>1 {
        id1 = $1":"$2":"$3":"$4
        id2 = $1":"$2":"$4":"$3
        print id1, $7
        print id2, $7
    }' "${SUMSTATS}" | sort -k1,1 -u | \
        sed '1i SNP P' > "${TMPDIR}/${BASENAME}_for_clump.txt"

    TOTAL_SIG=$(awk -F'\t' 'NR>1 && $7 < 5e-8' "${SUMSTATS}" | wc -l)
    echo "  Genome-wide significant variants: ${TOTAL_SIG}"

    if [ "${TOTAL_SIG}" -eq 0 ]; then
        echo "  No significant variants — writing empty lead SNP file"
        echo -e "CHR\tBP\tlocus" > "${OUTDIR}/${BASENAME}_leadSNPs.tsv"
        continue
    fi

    # Step 2: Find which chromosomes have BED files
    # The 1kg_eur/bed directory may have merged or per-chr files
    # Try per-chromosome approach first
    > "${TMPDIR}/${BASENAME}_all_clumped.txt"

    for CHR in $(seq 1 22); do
        # Check if this chromosome has significant hits
        CHR_SIG=$(awk -F'\t' -v chr="${CHR}" 'NR>1 && $1==chr && $7 < 5e-8' "${SUMSTATS}" | wc -l)
        if [ "${CHR_SIG}" -eq 0 ]; then
            continue
        fi

        BED_PREFIX="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/1kg_eur/chr${CHR}_eur"
        if [ ! -f "${BED_PREFIX}.bed" ]; then
            echo "  WARNING: No BED file for chr${CHR}, using distance-based for this chr"
            awk -F'\t' -v chr="${CHR}" 'NR>1 && $1==chr && $7 < 5e-8 {print chr":"$2":"$3":"$4, $7}' "${SUMSTATS}" | \
                sort -k2,2g | head -1 | awk -v chr="${CHR}" '{split($1,a,":"); print chr"\t"a[2]"\t"chr"."a[2]}' \
                >> "${TMPDIR}/${BASENAME}_all_clumped.txt"
            continue
        fi

        echo "  Clumping chr${CHR} (${CHR_SIG} sig variants)..."
        plink --bfile "${BED_PREFIX}" \
            --clump "${TMPDIR}/${BASENAME}_for_clump.txt" \
            --clump-p1 5e-8 \
            --clump-r2 0.1 \
            --clump-kb 500 \
            --chr "${CHR}" \
            --out "${TMPDIR}/${BASENAME}_chr${CHR}" \
            --allow-extra-chr \
            2>/dev/null || true

        # Extract clumped SNPs
        if [ -f "${TMPDIR}/${BASENAME}_chr${CHR}.clumped" ]; then
            awk 'NR>1 && $1!="" {split($3,a,":"); print a[1]"\t"a[2]"\t"a[1]"."a[2]}' \
                "${TMPDIR}/${BASENAME}_chr${CHR}.clumped" \
                >> "${TMPDIR}/${BASENAME}_all_clumped.txt"
        fi
    done

    # Step 3: Write final lead SNP file
    echo -e "CHR\tBP\tlocus" > "${OUTDIR}/${BASENAME}_leadSNPs.tsv"
    sort -k1,1n -k2,2n "${TMPDIR}/${BASENAME}_all_clumped.txt" \
        >> "${OUTDIR}/${BASENAME}_leadSNPs.tsv"

    N_LEADS=$(wc -l < "${TMPDIR}/${BASENAME}_all_clumped.txt")
    echo "  Final lead SNPs: ${N_LEADS}"
done

# Cleanup
rm -rf "${TMPDIR}"

echo ""
echo "============================================================"
echo "All lead SNP files written to ${OUTDIR}/"
echo "============================================================"
ls -la "${OUTDIR}"/*_leadSNPs.tsv
