#!/usr/bin/env bash
set -euo pipefail

ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CANDIDATE="${ROOT}/Analysis/Multimodal_Program_Projection/candidates/multicellular-assembly-response-2026-08-09"
OUT="${CANDIDATE}/source_expression"
mkdir -p "${OUT}"

download_once() {
  local url="$1"
  local destination="$2"
  local partial="${destination}.part"
  if [[ -s "${destination}" ]]; then
    gzip -t "${destination}"
    return 0
  fi
  curl --fail --location --retry 5 --retry-delay 10 --output "${partial}" "${url}"
  gzip -t "${partial}"
  mv "${partial}" "${destination}"
}

download_once \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168285/suppl/GSE168285_raw_counts.txt.gz" \
  "${OUT}/GSE168285_raw_counts.txt.gz"
download_once \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168285/suppl/GSE168285_gene_annotation.txt.gz" \
  "${OUT}/GSE168285_gene_annotation.txt.gz"
download_once \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168285/suppl/GSE168285_sample_meta_data.txt.gz" \
  "${OUT}/GSE168285_sample_meta_data.txt.gz"
download_once \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175448/suppl/GSE175448_Centaur_rawReadCount.txt.gz" \
  "${OUT}/GSE175448_Centaur_rawReadCount.txt.gz"

manifest="${OUT}/expression_source_manifest.tsv"
temporary="${manifest}.tmp"
printf 'relative_path\tbytes\tsha256\tsource_url\tretrieval_date\n' > "${temporary}"
while IFS=$'\t' read -r filename url; do
  path="${OUT}/${filename}"
  bytes=$(stat -c '%s' "${path}")
  digest=$(sha256sum "${path}" | awk '{print $1}')
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "source_expression/${filename}" "${bytes}" "${digest}" "${url}" "2026-08-09" \
    >> "${temporary}"
done <<'EOF'
GSE168285_raw_counts.txt.gz	https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168285/suppl/GSE168285_raw_counts.txt.gz
GSE168285_gene_annotation.txt.gz	https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168285/suppl/GSE168285_gene_annotation.txt.gz
GSE168285_sample_meta_data.txt.gz	https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168285/suppl/GSE168285_sample_meta_data.txt.gz
GSE175448_Centaur_rawReadCount.txt.gz	https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175448/suppl/GSE175448_Centaur_rawReadCount.txt.gz
EOF
mv "${temporary}" "${manifest}"

printf 'PLAN44_EXPRESSION_ACQUISITION_COMPLETE\tfiles=4\tmanifest_sha256=%s\n' \
  "$(sha256sum "${manifest}" | awk '{print $1}')"
