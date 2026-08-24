#!/bin/bash
#SBATCH --job-name=fig3d_gate
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3d_gate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3d_gate_%j.err

set -euo pipefail

candidate_root="figures/candidates/fig3d-endpoint-validation-2026-08-15-v4"
main_d="figures/main/fig3_bulk_transcriptomics/panels/fig3d_pca_fibrosis_gradient.pdf"
candidate_d="${candidate_root}/figure3/panels/fig3d_pca_fibrosis_gradient.pdf"
coordinates="figures/main/fig3_bulk_transcriptomics/source_tables/fig3d_endpoint_pca_coordinates.tsv"
census="figures/main/fig3_bulk_transcriptomics/source_tables/fig3d_endpoint_census.tsv"
figs3c_source="figures/supplementary/figS03_bulk_transcriptomics/source_tables/figs3c_cohort_robustness.tsv"
endpoint_metrics="figures/supplementary/figS03_deconvolution/source_tables/figs3_endpoint_vs_canonical_metrics.tsv"
main_manifest="figures/main/fig3_bulk_transcriptomics/manifests/fig3cd_endpoint_validation_manifest.tsv"
supp_all="figures/supplementary/figS03_deconvolution/panels/figs3d_pca_all_samples.pdf"
supp_concordance="figures/supplementary/figS03_deconvolution/panels/figs3_endpoint_deg_concordance.pdf"
supp_c="figures/supplementary/figS03_bulk_transcriptomics/panels/figs3c_cohort_robustness.pdf"
archive_all="figures/main/fig3_bulk_transcriptomics/panels/_superseded_all_samples_2026-08-14_v1/fig3d_pca_fibrosis_gradient.pdf"
archive_alluvial="figures/main/fig3_bulk_transcriptomics/panels/_superseded_replication_alluvial_2026-08-14_v1/fig3c_overlap_alluvial.pdf"
archive_strict="figures/main/fig3_bulk_transcriptomics/panels/_superseded_strict_extremes_2026-08-14_v1/fig3d_pca_fibrosis_gradient.pdf"
archive_magenta_blue="figures/main/fig3_bulk_transcriptomics/panels/_superseded_magenta_density_2026-08-14_v1/fig3d_pca_fibrosis_gradient.pdf"
archive_blue="figures/main/fig3_bulk_transcriptomics/panels/_superseded_blue_density_2026-08-14_v1/fig3d_pca_fibrosis_gradient.pdf"

for artifact in "${main_d}" "${candidate_d}" \
  "${coordinates}" "${census}" "${figs3c_source}" "${endpoint_metrics}" "${main_manifest}" \
  "${supp_all}" "${supp_concordance}" "${supp_c}" "${archive_all}" "${archive_alluvial}" \
  "${archive_strict}" "${archive_magenta_blue}" "${archive_blue}"; do
  test -f "${artifact}"
done
test "$(sha256sum "${main_d}" | cut -d' ' -f1)" = \
     "$(sha256sum "${candidate_d}" | cut -d' ' -f1)"
awk -F '\t' 'NR > 1 {n += $3; cohorts[$1] = 1; classes[$2] += $3}
  END {if (n != 418 || length(cohorts) != 4 || classes["Strict control"] != 115 || classes["Advanced disease"] != 303) exit 1}' "${census}"
test "$(awk 'END {print NR - 1}' "${coordinates}")" = "418"
awk -F '\t' 'NR > 1 {n++; if ($2 != 1347) exit 1}
  END {if (n != 5) exit 1}' "${figs3c_source}"
awk -F '\t' 'NR > 1 {v[$1]=$2}
  END {if (v["genes_tested"] != 23370 || v["canonical_degs"] != 1347 ||
           v["endpoint_degs"] != 2643 || v["overlap_degs"] != 1307) exit 1}' "${endpoint_metrics}"
for pdf in "${main_d}" "${supp_all}" "${supp_concordance}" "${supp_c}"; do
  test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
  if strings "${pdf}" | grep -q '/Subtype /Type3'; then
    echo "Type 3 font detected: ${pdf}" >&2
    exit 1
  fi
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
done
sha256sum -c "${candidate_root}/input_output_manifest.tsv"
awk -F '\t' 'NR > 1 {print $1 "  " $2}' "${main_manifest}" | sha256sum -c -
python3 scripts/manuscript/validate_resource_scope.py
echo "FIG3D_ENDPOINT_INDEX_GATE PASS"
