#!/usr/bin/env bash
# Download the complete Hu et al. GSE281364 barcode-level count universe.
set -euo pipefail

ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
OUT=${ROOT}/GWAS/finemapping/data/seqfunc_external/hu2025_mpra/raw_counts_v2
BASE=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281364/suppl
PARALLEL=${HU_DOWNLOAD_PARALLEL:-6}
mkdir -p "${OUT}"

files=(
  GSE281364_Lx2_Ctrl_pGL424_NAFLD-CRS-bc_pair_filtered_rep1_DNA_RNA_counts.txt.gz
  GSE281364_Lx2_Ctrl_pGL424_NAFLD-CRS-bc_pair_filtered_rep2_DNA_RNA_counts.txt.gz
  GSE281364_Lx2_Ctrl_pGL424_NAFLD-CRS-bc_pair_filtered_rep3_DNA_RNA_counts.txt.gz
  GSE281364_Lx2_Ctrl_pGL424_NAFLD-CRS-bc_pair_filtered_rep4_DNA_RNA_counts.txt.gz
  GSE281364_Lx2_TGFb_pGL424_NAFLD-CRS-bc_pair_filtered_rep1_DNA_RNA_counts.txt.gz
  GSE281364_Lx2_TGFb_pGL424_NAFLD-CRS-bc_pair_filtered_rep2_DNA_RNA_counts.txt.gz
  GSE281364_Lx2_TGFb_pGL424_NAFLD-CRS-bc_pair_filtered_rep3_DNA_RNA_counts.txt.gz
  GSE281364_Lx2_TGFb_pGL424_NAFLD-CRS-bc_pair_filtered_rep4_DNA_RNA_counts.txt.gz
  GSE281364_bwa-CRS-bc_pair_filtered_G2_Ctrl_Rep1_DNA_RNA_counts.txt.gz
  GSE281364_bwa-CRS-bc_pair_filtered_G2_Ctrl_Rep2_DNA_RNA_counts.txt.gz
  GSE281364_bwa-CRS-bc_pair_filtered_G2_Ctrl_Rep3_DNA_RNA_counts.txt.gz
  GSE281364_bwa-CRS-bc_pair_filtered_G2_Ctrl_Rep4_DNA_RNA_counts.txt.gz
  GSE281364_bwa-CRS-bc_pair_filtered_G2_PAOA_Rep1_DNA_RNA_counts.txt.gz
  GSE281364_bwa-CRS-bc_pair_filtered_G2_PAOA_Rep2_DNA_RNA_counts.txt.gz
  GSE281364_bwa-CRS-bc_pair_filtered_G2_PAOA_Rep3_DNA_RNA_counts.txt.gz
  GSE281364_bwa-CRS-bc_pair_filtered_G2_PAOA_Rep4b_DNA_RNA_counts.txt.gz
)

download_one() {
  local name=$1
  wget --continue --tries=20 --waitretry=5 --retry-on-http-error=429,500,502,503,504 \
    --timeout=60 --directory-prefix "${OUT}" "${BASE}/${name}"
  gzip -t "${OUT}/${name}"
}
export -f download_one
export OUT BASE
printf '%s\n' "${files[@]}" | xargs -n 1 -P "${PARALLEL}" bash -c 'download_one "$0"'

if [[ $(find "${OUT}" -maxdepth 1 -type f -name '*DNA_RNA_counts.txt.gz' | wc -l) -ne 16 ]]; then
  echo "Expected 16 count files" >&2
  exit 1
fi
sha256sum "${OUT}"/*DNA_RNA_counts.txt.gz | sort -k2,2 > "${OUT}/sha256sum.txt"
printf 'source\tGSE281364\nfiles\t16\nstatus\tcomplete\n' > "${OUT}/download_contract.tsv"
