#!/usr/bin/env bash
# Download the Hu et al. 2025 MASLD MPRA supplementary truth set.
#
# The processed DAV tables come from the Research Square/PMC supplement. The
# barcode map comes from GEO GSE281364 and defines the subset of oligos that
# passed barcode mapping. Large per-replicate DNA/RNA count files are optional.
set -euo pipefail

ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
OUT="${HU2025_MPRA_DIR:-${ROOT}/GWAS/finemapping/data/seqfunc_external/hu2025_mpra}"
RAW="${OUT}/raw"
mkdir -p "${RAW}"

download() {
    local url="$1"
    local dest="$2"
    if [[ -s "${dest}" ]]; then
        echo "[cached] ${dest}"
        return 0
    fi
    echo "[download] ${url}"
    curl -L --fail --retry 4 --retry-delay 5 --silent --show-error \
        -o "${dest}.partial" "${url}"
    mv "${dest}.partial" "${dest}"
}

RS_BASE="https://assets-eu.researchsquare.com/files/rs-6984670/v1"
download "${RS_BASE}/05c95998495329ab00c3344b.xlsx" "${RAW}/TableS1.xlsx"
download "${RS_BASE}/28ed16c000678dd7abbf266a.xlsx" "${RAW}/TableS2.xlsx"
download "${RS_BASE}/cdb18a66f49f0fa0e4980ef1.xlsx" "${RAW}/TableS3.xlsx"
download "${RS_BASE}/480c242b324269e5b55c5bbc.xlsx" "${RAW}/TableS4.xlsx"
download "${RS_BASE}/253bcc974497f7bd2049a765.xlsx" "${RAW}/TableS5.xlsx"
download "${RS_BASE}/bfc08cbd3cf616067ccf819a.xlsx" "${RAW}/TableS6.xlsx"

GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281364/suppl"
download "${GEO_BASE}/GSE281364_RAW.tar" "${RAW}/GSE281364_RAW.tar"

BARCODE_GZ="${RAW}/GSM8619256_bwa-CRS-bc_pair_filtered.txt.gz"
if [[ ! -s "${BARCODE_GZ}" ]]; then
    tar -xf "${RAW}/GSE281364_RAW.tar" -C "${RAW}" \
        GSM8619256_bwa-CRS-bc_pair_filtered.txt.gz
fi

# Reference and alternate oligos share the same interval name, with the alternate
# carrying a _Mut suffix. Collapse to the tested genomic interval universe.
TESTED="${OUT}/tested_oligo_intervals.tsv"
if [[ ! -s "${TESTED}" ]]; then
    gzip -cd "${BARCODE_GZ}" \
        | cut -f1 \
        | sed 's/_Mut$//' \
        | sort -u \
        > "${TESTED}.partial"
    mv "${TESTED}.partial" "${TESTED}"
fi

if [[ "${DOWNLOAD_COUNTS:-FALSE}" == "TRUE" ]]; then
    files=(
        GSE281364_bwa-CRS-bc_pair_filtered_G2_Ctrl_Rep1_DNA_RNA_counts.txt.gz
        GSE281364_bwa-CRS-bc_pair_filtered_G2_Ctrl_Rep2_DNA_RNA_counts.txt.gz
        GSE281364_bwa-CRS-bc_pair_filtered_G2_Ctrl_Rep3_DNA_RNA_counts.txt.gz
        GSE281364_bwa-CRS-bc_pair_filtered_G2_Ctrl_Rep4_DNA_RNA_counts.txt.gz
        GSE281364_bwa-CRS-bc_pair_filtered_G2_PAOA_Rep1_DNA_RNA_counts.txt.gz
        GSE281364_bwa-CRS-bc_pair_filtered_G2_PAOA_Rep2_DNA_RNA_counts.txt.gz
        GSE281364_bwa-CRS-bc_pair_filtered_G2_PAOA_Rep3_DNA_RNA_counts.txt.gz
        GSE281364_bwa-CRS-bc_pair_filtered_G2_PAOA_Rep4b_DNA_RNA_counts.txt.gz
        GSE281364_Lx2_Ctrl_pGL424_NAFLD-CRS-bc_pair_filtered_rep1_DNA_RNA_counts.txt.gz
        GSE281364_Lx2_Ctrl_pGL424_NAFLD-CRS-bc_pair_filtered_rep2_DNA_RNA_counts.txt.gz
        GSE281364_Lx2_Ctrl_pGL424_NAFLD-CRS-bc_pair_filtered_rep3_DNA_RNA_counts.txt.gz
        GSE281364_Lx2_Ctrl_pGL424_NAFLD-CRS-bc_pair_filtered_rep4_DNA_RNA_counts.txt.gz
        GSE281364_Lx2_TGFb_pGL424_NAFLD-CRS-bc_pair_filtered_rep1_DNA_RNA_counts.txt.gz
        GSE281364_Lx2_TGFb_pGL424_NAFLD-CRS-bc_pair_filtered_rep2_DNA_RNA_counts.txt.gz
        GSE281364_Lx2_TGFb_pGL424_NAFLD-CRS-bc_pair_filtered_rep3_DNA_RNA_counts.txt.gz
        GSE281364_Lx2_TGFb_pGL424_NAFLD-CRS-bc_pair_filtered_rep4_DNA_RNA_counts.txt.gz
    )
    for f in "${files[@]}"; do
        download "${GEO_BASE}/${f}" "${RAW}/${f}"
    done
fi

{
    printf 'file\tsha256\tbytes\n'
    for f in "${RAW}"/TableS{1,2,3,4,5,6}.xlsx "${RAW}/GSE281364_RAW.tar" "${BARCODE_GZ}" "${TESTED}"; do
        hash="$(sha256sum "${f}" | awk '{print $1}')"
        bytes="$(stat -c '%s' "${f}")"
        printf '%s\t%s\t%s\n' "${f#${ROOT}/}" "${hash}" "${bytes}"
    done
} > "${OUT}/source_manifest.tsv"

echo "[done] Hu2025 MPRA data: ${OUT}"
echo "[note] DOWNLOAD_COUNTS=TRUE downloads ~2.6 GB of per-replicate count files."

