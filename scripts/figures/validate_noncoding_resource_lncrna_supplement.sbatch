#!/bin/bash
#SBATCH --job-name=nclnfigval
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nclnfigval_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nclnfigval_%j.err

set -euo pipefail
umask 027

SOURCE="figures/candidates/noncoding-resource-2026-08-11/lncrna_supplement_2026-08-11-r2"
REVIEW="figures/candidates/noncoding-resource-2026-08-11/lncrna_supplement_review_2026-08-11-r2"
test -d "${SOURCE}"
test ! -e "${REVIEW}"
mkdir "${REVIEW}"

awk -F '\t' 'NR > 1 {print $1 "\t" $2 "\t" $3}' "${SOURCE}/output_manifest.tsv" |
while IFS=$'\t' read -r relative_path size_bytes expected_sha; do
  test "$(stat -c %s "${SOURCE}/${relative_path}")" = "${size_bytes}"
  test "$(sha256sum "${SOURCE}/${relative_path}" | cut -d ' ' -f 1)" = "${expected_sha}"
done

for pdf in "${SOURCE}"/*.pdf; do
  test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
  pdftoppm -png -r 180 -singlefile "${pdf}" "${REVIEW}/$(basename "${pdf}" .pdf)"
done

{
  printf 'check\tstatus\n'
  printf 'manifest_size_and_sha256\tPASS\n'
  printf 'one_page_pdfs\tPASS\n'
  printf 'ghostscript_parse\tPASS\n'
  printf 'png_previews\tPASS\n'
} > "${REVIEW}/validation.tsv"
