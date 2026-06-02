#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=reformat_add

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=== Reformatting 4 additional GWAS ==="
Rscript src/00b_reformat_additional_gwas.R

echo ""
echo "=== Identifying lead SNPs ==="
module load PLINK/1.9

OUTDIR="data/lead_snps"
TMPDIR="data/lead_snps/tmp_clump"
mkdir -p "${OUTDIR}" "${TMPDIR}"

for SUMSTATS in data/sumstats/{Anstee,Pazoki,Ghouse}*_reformatted_hg19.tsv; do
    BASENAME=$(basename "${SUMSTATS}" _reformatted_hg19.tsv)
    echo "--- Processing: ${BASENAME} ---"

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

        BED_PREFIX="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/1kg_eur/chr${CHR}_eur"
        if [ ! -f "${BED_PREFIX}.bed" ]; then continue; fi

        echo "  Clumping chr${CHR} (${CHR_SIG} sig)..."
        plink --bfile "${BED_PREFIX}" \
            --clump "${TMPDIR}/${BASENAME}_for_clump.txt" \
            --clump-p1 5e-8 --clump-r2 0.1 --clump-kb 500 \
            --chr "${CHR}" \
            --out "${TMPDIR}/${BASENAME}_chr${CHR}" \
            --allow-extra-chr 2>/dev/null || true

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
echo "=== All done ==="
