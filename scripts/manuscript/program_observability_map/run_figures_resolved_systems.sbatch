#!/bin/bash
#SBATCH --job-name=resolved-system-figs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/logs/resolved_system_figs_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/logs/resolved_system_figs_%j.err

set -euo pipefail

project_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
rscript=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
generator="$project_root/scripts/manuscript/program_observability_map/figures_resolved_systems.R"
candidate="$project_root/figures/candidates/program-systems-2026-08-12-v3"

cd "$project_root"

"$rscript" -e "parse(file='$generator')" >/dev/null
"$rscript" "$generator"

test -d "$candidate"
test "$(find "$candidate" -maxdepth 1 -type f -name '*.pdf' | wc -l)" -eq 4
test "$(awk 'END { print NR - 1 }' "$candidate/panel_t4_coordination_edges_source.tsv")" -eq 2426
test "$(awk 'END { print NR - 1 }' "$candidate/panel_t4_program_order_source.tsv")" -eq 113
test "$(awk 'END { print NR - 1 }' "$candidate/panel_t5_composition_coupling_source.tsv")" -eq 27
test "$(awk 'END { print NR - 1 }' "$candidate/panel_t6_sex_modification_source.tsv")" -eq 452

while IFS= read -r pdf; do
  test "$(pdfinfo "$pdf" | awk '/^Pages:/ { print $2 }')" -eq 1
  if grep -aEq '/Subtype[[:space:]]*/Type3' "$pdf"; then
    echo "Type 3 font found in $pdf" >&2
    exit 1
  fi
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "$pdf"
done < <(find "$candidate" -maxdepth 1 -type f -name '*.pdf' | sort)

python3 scripts/manuscript/validate_resource_scope.py

echo "RESOLVED_SYSTEM_FIGURES_VALIDATED $candidate"
