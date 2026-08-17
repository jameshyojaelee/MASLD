#!/bin/bash
#SBATCH --job-name=fig4_stage_validate
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4_stage_validate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4_stage_validate_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CANDIDATE="${PROJECT_ROOT}/figures/candidates/fig4-stage-terminology-corrected-2026-08-13-v5"
PANEL_ROOT="${CANDIDATE}/figure4/panels"
RNASEQ_RSCRIPT="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
export FIGURE_VALIDATE_ROOT="${CANDIDATE}"

test -s "${CANDIDATE}/sessionInfo.txt"
test ! -e "${CANDIDATE}/pdf_validation.tsv"
mapfile -t PDFS < <(find "${PANEL_ROOT}" -type f -name '*.pdf' | sort)
test "${#PDFS[@]}" -eq 7

printf 'relative_path\tpages\tghostscript\ttype3_marker\tstale_stage_label\n' > "${CANDIDATE}/pdf_validation.tsv"
mkdir -p "${CANDIDATE}/proofs/raster"
for pdf in "${PDFS[@]}"; do
  rel="${pdf#${CANDIDATE}/}"
  pages="$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')"
  test "${pages}" -eq 1
  gs -q -dBATCH -dNOPAUSE -sDEVICE=nullpage "${pdf}"
  type3='FALSE'
  if grep -aEq '/Subtype[[:space:]]*/Type3|/Type[[:space:]]*/Font[[:space:]]*/Subtype[[:space:]]*/Type3' "${pdf}"; then
    type3='TRUE'
  fi
  test "${type3}" = 'FALSE'
  text_file="${CANDIDATE}/proofs/raster/$(basename "${pdf}" .pdf).txt"
  pdftotext "${pdf}" "${text_file}"
  stale='FALSE'
  if grep -Eq '(^|[^A-Za-z])(MASL|MASH|NASH|Steat\.)([^A-Za-z]|$)' "${text_file}"; then
    stale='TRUE'
  fi
  test "${stale}" = 'FALSE'
  printf '%s\t%s\tPASS\t%s\t%s\n' "${rel}" "${pages}" "${type3}" "${stale}" >> "${CANDIDATE}/pdf_validation.tsv"
  pdftoppm -f 1 -singlefile -png -r 180 "${pdf}" \
    "${CANDIDATE}/proofs/raster/$(basename "${pdf}" .pdf)" >/dev/null 2>&1
done

"${RNASEQ_RSCRIPT}" --vanilla -e '
  root <- Sys.getenv("FIGURE_VALIDATE_ROOT")
  read <- function(name) utils::read.delim(file.path(root, "source_tables", name), check.names=FALSE)
  registry <- read("fig4c_all_117_programs.tsv")
  heatmap <- read("fig4c_hotspot_stage_heatmap.csv")
  profiles <- read("fig4c_selected_program_profiles.tsv")
  communication <- read("fig4d_fixed_13_pair_communication.tsv")
  projection <- read("fig4e_bulk_projection.tsv")
  tf_all <- read("fig4f_all_tf_activity.tsv")
  tf_display <- read("fig4f_display_tf_activity.tsv")
  stopifnot(nrow(registry) == 117L, length(unique(registry$program_uid)) == 117L)
  stopifnot(sum(registry$primary_selected) == 2L, sum(registry$hc3_supported) == 0L)
  stopifnot(nrow(heatmap) == 117L, nrow(profiles) == 6L)
  stopifnot(nrow(communication) == 39L, length(unique(communication$headline_label)) == 13L)
  stopifnot(nrow(projection) == 8L, length(unique(projection$program_uid)) == 2L)
  stopifnot(nrow(tf_all) == 1072L, all(table(tf_all$contrast) == 268L))
  stopifnot(nrow(tf_display) == 64L, length(unique(tf_display$tf)) == 16L)
  out <- data.frame(
    table=c("program_registry", "main_heatmap_source", "selected_profiles",
            "communication", "bulk_projection", "all_tf_activity", "display_tf_activity"),
    rows=c(nrow(registry), nrow(heatmap), nrow(profiles), nrow(communication),
           nrow(projection), nrow(tf_all), nrow(tf_display)),
    status="PASS"
  )
  utils::write.table(out, file.path(root, "source_cardinality_validation.tsv"),
                     sep="\t", quote=FALSE, row.names=FALSE)
' 

mapfile -t MAIN_PNGS < <(find "${CANDIDATE}/proofs/raster" -maxdepth 1 \
  -type f -name 'fig4[bcdef]*.png' ! -name 'fig4c_program_landscape.png' | sort)
test "${#MAIN_PNGS[@]}" -eq 5
montage "${MAIN_PNGS[@]}" -thumbnail '1200x900>' -background white \
  -geometry '+20+20' -tile 3x2 "${CANDIDATE}/proofs/figure4_candidate_proof.png"

find "${CANDIDATE}" -type f ! -path '*/proofs/raster/*' \
  ! -name 'output_manifest.sha256' -print0 | sort -z | \
  xargs -0 sha256sum > "${CANDIDATE}/output_manifest.sha256"
python3 scripts/manuscript/validate_resource_scope.py
