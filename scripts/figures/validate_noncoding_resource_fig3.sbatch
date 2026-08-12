#!/bin/bash
#SBATCH --job-name=ncfigrev
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncfigrev_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncfigrev_%j.err

set -euo pipefail
umask 027

PDFINFO="/gpfs/commons/home/jameslee/.local/bin/pdfinfo"
PDFTOPPM="/gpfs/commons/home/jameslee/.local/bin/pdftoppm"

: "${NONCODING_FIGURE_ROOT:?Set NONCODING_FIGURE_ROOT to the immutable candidate}"
: "${NONCODING_FIGURE_REVIEW:?Set NONCODING_FIGURE_REVIEW to a new review directory}"

test -d "${NONCODING_FIGURE_ROOT}"
test ! -e "${NONCODING_FIGURE_REVIEW}"
mkdir "${NONCODING_FIGURE_REVIEW}"

for pdf in fig3_pooled_deg_biotype.pdf fig3_stage_deg_biotype.pdf; do
  test -f "${NONCODING_FIGURE_ROOT}/${pdf}"
  pages=$("${PDFINFO}" "${NONCODING_FIGURE_ROOT}/${pdf}" | awk -F: '$1 == "Pages" {gsub(/ /, "", $2); print $2}')
  test "${pages}" = "1"
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${NONCODING_FIGURE_ROOT}/${pdf}"
  "${PDFTOPPM}" -png -r 180 -singlefile \
    "${NONCODING_FIGURE_ROOT}/${pdf}" \
    "${NONCODING_FIGURE_REVIEW}/${pdf%.pdf}"
done

while IFS=$'\t' read -r relative_path size_bytes expected_sha; do
  test "${relative_path}" = "relative_path" && continue
  artifact="${NONCODING_FIGURE_ROOT}/${relative_path}"
  test -f "${artifact}"
  test "$(stat -c %s "${artifact}")" = "${size_bytes}"
  test "$(sha256sum "${artifact}" | cut -d ' ' -f1)" = "${expected_sha}"
done < "${NONCODING_FIGURE_ROOT}/output_manifest.tsv"

printf 'status\tPASS\nsource_manifest\tPASS\none_page_pdfs\tPASS\npdf_parse\tPASS\n' > \
  "${NONCODING_FIGURE_REVIEW}/validation.tsv"
