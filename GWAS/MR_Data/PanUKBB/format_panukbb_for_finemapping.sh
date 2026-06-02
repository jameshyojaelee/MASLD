#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=fmt_panukbb_fm
#SBATCH --output=logs/format_panukbb_finemapping_%j.out
#SBATCH --error=logs/format_panukbb_finemapping_%j.err

# format_panukbb_for_finemapping.sh
# ---------------------------------------------------------------------------
# Reformat Pan-UKBB AFR & CSA liver enzyme GWAS for the finemapping pipeline,
# then identify lead SNPs via PLINK LD clumping with ancestry-matched panels.
#
# Usage: cd GWAS/MR_Data/PanUKBB && sbatch format_panukbb_for_finemapping.sh
# ---------------------------------------------------------------------------

set -eo pipefail

BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
FM_DIR="${BASE_DIR}/GWAS/finemapping"
PANUKBB_DIR="${BASE_DIR}/GWAS/MR_Data/PanUKBB"

cd "${PANUKBB_DIR}"
mkdir -p logs

echo "============================================================"
echo "Step 1: Reformat Pan-UKBB sumstats to hg19"
echo "============================================================"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

Rscript "${PANUKBB_DIR}/format_panukbb_for_finemapping.R"

echo ""
echo "============================================================"
echo "Step 2: Identify lead SNPs via PLINK LD clumping"
echo "============================================================"

module load PLINK/1.9

OUTDIR="${FM_DIR}/data/lead_snps"
TMPDIR="${FM_DIR}/data/lead_snps/tmp_clump_panukbb"
mkdir -p "${OUTDIR}" "${TMPDIR}"

# Ancestry → LD panel mapping
# AFR → 1kg_afr;  CSA (Pan-UKBB) → 1kg_sas (closest 1KG superpopulation)
declare -A LD_PANELS
LD_PANELS[AFR]="${FM_DIR}/data/ld_ref/1kg_afr"
LD_PANELS[CSA]="${FM_DIR}/data/ld_ref/1kg_sas"

declare -A LD_SUFFIX
LD_SUFFIX[AFR]="afr"
LD_SUFFIX[CSA]="sas"

for SUMSTATS in "${FM_DIR}/data/sumstats/PanUKBB_"*"_reformatted_hg19.tsv"; do
    BASENAME=$(basename "${SUMSTATS}" _reformatted_hg19.tsv)
    # Extract population from filename: PanUKBB_{POP}_{TRAIT}
    POP=$(echo "${BASENAME}" | sed 's/PanUKBB_\([A-Z]*\)_.*/\1/')

    LD_BASE="${LD_PANELS[${POP}]}"
    LD_SUF="${LD_SUFFIX[${POP}]}"

    echo ""
    echo "--- Processing: ${BASENAME} (LD panel: ${LD_SUF}) ---"

    # Create PLINK-compatible input: SNP (chr:pos:a1:a2) and P columns
    awk -F'\t' 'NR>1 {
        id1 = $1":"$2":"$3":"$4
        id2 = $1":"$2":"$4":"$3
        print id1, $7
        print id2, $7
    }' "${SUMSTATS}" | sort -k1,1 -u | sed '1i SNP P' > "${TMPDIR}/${BASENAME}_for_clump.txt"

    TOTAL_SIG=$(awk -F'\t' 'NR>1 && $7 < 5e-8' "${SUMSTATS}" | wc -l)
    echo "  Genome-wide significant variants: ${TOTAL_SIG}"

    if [ "${TOTAL_SIG}" -eq 0 ]; then
        echo "  No significant variants — writing empty lead SNP file"
        echo -e "CHR\tBP\tlocus" > "${OUTDIR}/${BASENAME}_leadSNPs.tsv"
        continue
    fi

    > "${TMPDIR}/${BASENAME}_all_clumped.txt"

    for CHR in $(seq 1 22); do
        CHR_SIG=$(awk -F'\t' -v chr="${CHR}" 'NR>1 && $1==chr && $7 < 5e-8' "${SUMSTATS}" | wc -l)
        if [ "${CHR_SIG}" -eq 0 ]; then continue; fi

        BED_PREFIX="${LD_BASE}/chr${CHR}_${LD_SUF}"
        if [ ! -f "${BED_PREFIX}.bed" ]; then
            echo "  WARNING: LD panel not found for chr${CHR}: ${BED_PREFIX}.bed"
            continue
        fi

        # Copy ref panel to tmpdir and rename BIM to chr:pos:a1:a2
        # (1KG BIM uses rsIDs; our clump file uses chr:pos:a1:a2)
        CHR_TMP="${TMPDIR}/ref_chr${CHR}_${LD_SUF}"
        cp "${BED_PREFIX}.bed" "${CHR_TMP}.bed"
        cp "${BED_PREFIX}.fam" "${CHR_TMP}.fam"
        awk -F'\t' 'BEGIN{OFS="\t"} {$2=$1":"$4":"$5":"$6; print}' \
            "${BED_PREFIX}.bim" > "${CHR_TMP}.bim"

        echo "  Clumping chr${CHR} (${CHR_SIG} sig)..."
        plink --bfile "${CHR_TMP}" \
            --clump "${TMPDIR}/${BASENAME}_for_clump.txt" \
            --clump-p1 5e-8 --clump-r2 0.1 --clump-kb 500 \
            --chr "${CHR}" \
            --out "${TMPDIR}/${BASENAME}_chr${CHR}" \
            --allow-extra-chr 2>/dev/null || true

        rm -f "${CHR_TMP}.bed" "${CHR_TMP}.bim" "${CHR_TMP}.fam"

        if [ -f "${TMPDIR}/${BASENAME}_chr${CHR}.clumped" ]; then
            awk 'NR>1 && NF>2 {split($3,a,":"); print a[1]"\t"a[2]"\t"a[1]"."a[2]}' \
                "${TMPDIR}/${BASENAME}_chr${CHR}.clumped" \
                >> "${TMPDIR}/${BASENAME}_all_clumped.txt"
        fi
    done

    echo -e "CHR\tBP\tlocus" > "${OUTDIR}/${BASENAME}_leadSNPs.tsv"
    sort -k1,1n -k2,2n "${TMPDIR}/${BASENAME}_all_clumped.txt" >> "${OUTDIR}/${BASENAME}_leadSNPs.tsv"
    N_LEADS=$(wc -l < "${TMPDIR}/${BASENAME}_all_clumped.txt")
    echo "  Final lead SNPs: ${N_LEADS}"
done

rm -rf "${TMPDIR}"

echo ""
echo "============================================================"
echo "All Pan-UKBB GWAS reformatted and lead SNPs identified"
echo "============================================================"
echo ""
echo "Next steps:"
echo "  1. Add registry entries to GWAS/finemapping/config/gwas_registry.tsv"
echo "  2. Run finemapping: sbatch src/02_run_susie.sh"
