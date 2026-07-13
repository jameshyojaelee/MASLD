#!/usr/bin/env bash
# Download the authoritative, compact Currin adult-liver caQTL truth files.
set -euo pipefail

ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
OUT=${ROOT}/GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1
BASE=https://zenodo.org/api/records/15025748/files
mkdir -p "${OUT}"

files=(
  background_variants_overlapping_non-caPeaks.bed.gz
  caQTL_variants_overlappingPeaks_LD-r2-0.8_withLead.bed.gz
  liver_significant_caQTL_leadVariants_1mb_analysis.bed.gz
  liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz
  supplementalData1_liver_ATAC_peaks.bed.gz
)
download_one() {
  local name=$1
  wget --continue --tries=8 --timeout=60 --output-document "${OUT}/${name}" \
    "${BASE}/${name}/content"
  gzip -t "${OUT}/${name}"
}
export -f download_one
export OUT BASE
printf '%s\n' "${files[@]}" | xargs -n 1 -P 5 bash -c 'download_one "$0"'
sha256sum "${OUT}"/*.gz | sort -k2,2 > "${OUT}/sha256sum.txt"
printf 'source\tCurrin_2025\nzenodo_record\t15025748\nrelease\tv1_2025-03-14\ngenome_build\tGRCh38\nstatus\tcomplete\n' > "${OUT}/source_contract.tsv"
