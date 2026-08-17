#!/bin/bash
#SBATCH --job-name=fig4ef_validate
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4ef_validate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4ef_validate_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
FIGURE_ROOT="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs"
MANIFEST="${FIGURE_ROOT}/manifests/current_candidate_manifest.tsv"
RNASEQ_R="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

while IFS=$'\t' read -r role relative_path expected_size expected_sha state; do
  [[ "${role}" == "role" ]] && continue
  target="${FIGURE_ROOT}/${relative_path}"
  [[ -f "${target}" ]]
  [[ "$(stat -c '%s' "${target}")" == "${expected_size}" ]]
  [[ "$(sha256sum "${target}" | cut -d' ' -f1)" == "${expected_sha}" ]]
done < "${MANIFEST}"

for pdf in \
  "${FIGURE_ROOT}/panels/fig4e_bulk_projection_pre_genetics.pdf" \
  "${FIGURE_ROOT}/panels/fig4f_tf_activity.pdf"; do
  [[ "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" == "1" ]]
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
  if gs -dNOSAFER -dNODISPLAY -dPDFDEBUG -f "${pdf}" -c quit 2>&1 | \
      grep -Eq '/Subtype[[:space:]]*/Type3'; then
    echo "Type 3 font found in ${pdf}" >&2
    exit 1
  fi
done

awk -F'\t' '
  NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
  { n++; program[$h["program_uid"]] = 1; contrast[$h["contrast"]] = 1 }
  END {
    if (n != 16 || length(program) != 4 || length(contrast) != 4) exit 1
  }
' "${FIGURE_ROOT}/source_tables/current_candidate/fig4e_bulk_projection.tsv"

awk -F'\t' '
  NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
  {
    n++; tf[$h["tf"]] = 1; contrast[$h["contrast"]] = 1
    if ($h["bh_family_size"] != 268) exit 1
  }
  END {
    if (n != 1072 || length(tf) != 268 || length(contrast) != 4) exit 1
  }
' "${FIGURE_ROOT}/source_tables/current_candidate/fig4f_all_tf_activity.tsv"

awk -F'\t' '
  NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
  {
    n++; tf[$h["tf"]] = 1; contrast[$h["contrast"]] = 1
    if ($h["bh_family_size"] != 268) exit 1
  }
  END {
    if (n != 64 || length(tf) != 16 || length(contrast) != 4) exit 1
  }
' "${FIGURE_ROOT}/source_tables/current_candidate/fig4f_display_tf_activity.tsv"

"${RNASEQ_R}" --vanilla -e \
  'invisible(parse(file="scripts/figures/render_pi_singlecell_panels.R")); invisible(parse(file="scripts/figures/tf_convergence_vsF0.R"))'
python3 scripts/manuscript/validate_resource_scope.py

echo "Figure 4E/4F validation passed"
