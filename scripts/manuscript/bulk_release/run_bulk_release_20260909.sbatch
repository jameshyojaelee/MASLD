#!/bin/bash
#SBATCH --job-name=bulk_release_r1
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/bulk_release_r1_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/bulk_release_r1_%j.err

# Regenerate the selected stage/NAS output and Figure 3A-E/S3 from the validated
# five-cohort substrate, then package only the bulk-side release.  No genetic
# membership, Gene Catalog, portal, or final Resource adoption is touched.
set -euo pipefail
umask 027

project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
rnaseq_env="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
# Defaults are the 2026-09-09 r1 run. The corrected-control chain
# (scripts/manuscript/resource_f_five/submit_corrected_controls_chain.sh) passes
# BULK_RELEASE_ID, POOLED_ROOT and F0_ROOT, and BULK_RELEASE_PACKAGE=0 until the
# genetics-side Figure 3G has been rebuilt on the same bulk fit.
release_id="${BULK_RELEASE_ID:-2026-09-09-bulk-r1}"
candidate_root="${project_root}/figures/candidates/bulk-release-${release_id}"
release_root="${project_root}/RNA-seq/results/manuscript_release/${release_id}"
pooled_root="${POOLED_ROOT:-${project_root}/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-POOLED-REPRO}"
f0_root="${F0_ROOT:-${project_root}/figures/candidates/f0-reference-arms-20260907T105813Z/analysis/f0_arms}"
package_release="${BULK_RELEASE_PACKAGE:-1}"
composition_path="${project_root}/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/precontract/composition_acceptance/accepted_composition.tsv.gz"
lncrna_root="${project_root}/RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/bulk_disease/run_canonical_gate_2026-08-12"
lncrna_biotype_root="${project_root}/RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/biotype_association/run_canonical_gate_2026-08-12"
fig3g_root="${project_root}/figures/candidates/fig3g-promoted-coloc-2026-08-17-v17-equal-points"
fig4_candidate="${project_root}/figures/candidates/pi-figure-redesign-2026-08-17-v13"

test ! -e "${candidate_root}"
test ! -e "${release_root}"
for path in "${pooled_root}/VALIDATED.json" "${f0_root}/f0_arm_results.tsv.gz" \
            "${composition_path}" "${lncrna_root}/output_manifest.tsv" \
            "${lncrna_biotype_root}/output_manifest.tsv" "${fig3g_root}/input_output_manifest.tsv"; do
  test -f "${path}"
done

export MASLD_PROJECT_ROOT="${project_root}"
export FIGURE_CANDIDATE_ROOT="${candidate_root}"
export STAGE_RELEASE_ROOT="${pooled_root}"
export STAGE_EXTENSION_ROOT="${candidate_root}/analysis/stage_extensions"
export COMPOSITION_RELEASE_PATH="${composition_path}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export PATH="${rnaseq_env}/bin:${PATH}"

"${rnaseq_env}/bin/Rscript" --vanilla "${project_root}/scripts/figures/build_pi_stage_extensions.R"

# These producers are scoped to Figure 3 and its supplement.  In particular,
# do not invoke the omnibus finalizer, which rewrites out-of-scope Figures 4-6.
"${rnaseq_env}/bin/Rscript" --vanilla "${project_root}/scripts/figures/nas_fib_grid.R"
"${rnaseq_env}/bin/Rscript" --vanilla "${project_root}/scripts/figures/fig3d_endpoint_validation.R"
"${rnaseq_env}/bin/Rscript" --vanilla "${project_root}/scripts/figures/progression_cascade.R"

for pdf in $(find "${candidate_root}/figure3" "${candidate_root}/supplementary/figureS3" -type f -name '*.pdf' | sort); do
  test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
  if strings "${pdf}" | grep -q '/Subtype /Type3'; then
    echo "Type 3 font detected: ${pdf}" >&2
    exit 1
  fi
  pdftotext -layout "${pdf}" - >/dev/null
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
done
mkdir -p "${candidate_root}/proofs"
gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r240 \
  -sOutputFile="${candidate_root}/proofs/fig3e_stage_remodeling.png" \
  "${candidate_root}/figure3/panels/fig3e_stage_remodeling.pdf"
gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r180 \
  -sOutputFile="${candidate_root}/proofs/figs3_stage_remodeling_full.png" \
  "${candidate_root}/supplementary/figureS3/panels/figs3_stage_remodeling_full.pdf"

if [[ "${package_release}" != "1" ]]; then
  echo "BULK_FIGURE_PRODUCERS_COMPLETE ${candidate_root} (release packaging skipped)"
  exit 0
fi
"${rnaseq_env}/bin/python3" "${project_root}/scripts/manuscript/bulk_release/build_bulk_release.py" \
  --release-id "${release_id}" \
  --output "${release_root}" \
  --pooled "${pooled_root}" \
  --stage "${candidate_root}/analysis/stage_extensions" \
  --f0 "${f0_root}" \
  --composition "${composition_path}" \
  --lncrna "${lncrna_root}" \
  --lncrna-biotype "${lncrna_biotype_root}" \
  --figure-candidate "${candidate_root}" \
  --fig3g "${fig3g_root}" \
  --fig4-candidate "${fig4_candidate}"
"${rnaseq_env}/bin/python3" "${project_root}/scripts/manuscript/bulk_release/validate_bulk_release.py" \
  --release-id "${release_id}" \
  --release "${release_root}"

echo "BULK_RELEASE_COMPLETE ${release_root}"
