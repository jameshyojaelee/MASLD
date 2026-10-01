#!/bin/bash
#SBATCH --job-name=fig3g_coloc_r2
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3g_coloc_r2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3g_coloc_r2_%j.err

# Figure 3G on the corrected bulk fit and the corrected COLOC release
# (review items 1-2, Phase 3). Interface recompute, then render, as
# run_recompute_genetics_state_interface.sbatch + run_fig3g_promoted.sbatch, with
# the adopted-release counts replaced by agreement between the two code paths.
#
#   sbatch --dependency=afterok:<coloc_r2 job>:<fcorr_repro job> --export=ALL,\
# GENETICS_RELEASE_ROOT=<FM>/results/susie_coloc_r2_<ts>,\
# BULK_ARM_ROOT=<RNA-seq>/results/manuscript_release/candidates/bulk-corrected-controls-<ts>-reproduction/primary_source_controls,\
# GENETICS_STATE_INTERFACE_OUT=<RNA-seq>/results/genetics_state_interface/<UTC ts>_coloc_r2_corrected_controls,\
# FIGURE_CANDIDATE_ROOT=<BASE>/figures/candidates/fig3g-coloc-r2-<UTC ts> \
#     scripts/figures/run_fig3g_coloc_r2.sbatch
#
# BULK_ARM_ROOT/model_input_manifest.tsv records the ORIGINAL F_five fit as
# accepted_reference_path, so the corrected deg_results.csv is passed explicitly.
set -euo pipefail
umask 027
: "${GENETICS_RELEASE_ROOT:?Set the corrected COLOC release (results/susie_coloc_r2_<ts>)}"
: "${BULK_ARM_ROOT:?Set the corrected primary_source_controls arm of the reproduction root}"
: "${GENETICS_STATE_INTERFACE_OUT:?Set a fresh interface output directory}"
: "${FIGURE_CANDIDATE_ROOT:?Set a fresh candidate path}"
test ! -e "${GENETICS_STATE_INTERFACE_OUT}"
test ! -e "${FIGURE_CANDIDATE_ROOT}"
project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
rnaseq_env="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
export PATH="/gpfs/commons/home/jameslee/.local/bin:${PATH}"
coloc="${GENETICS_RELEASE_ROOT}/gene_level_coloc_tier12.csv"
deg="${BULK_ARM_ROOT}/deg_results.csv"

test -s "${GENETICS_RELEASE_ROOT}/ASSEMBLY_COMPLETE.json"
sha256sum -c "${GENETICS_RELEASE_ROOT}/output_sha256.tsv"
test -s "${deg}"
test -s "${BULK_ARM_ROOT}/model_input_manifest.tsv"
# The reproduction job must have passed every check on this arm.
"${rnaseq_env}/bin/python3" - "${BULK_ARM_ROOT}/validation/validation_checks.tsv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
failed = [r["check"] for r in rows if r["pass"] != "TRUE"]
assert rows and not failed, f"bulk reproduction checks failed: {failed}"
PY

"${rnaseq_env}/bin/python3" \
  "${project_root}/scripts/manuscript/recompute_genetics_state_interface.py" \
  --out "${GENETICS_STATE_INTERFACE_OUT}" --coloc "${coloc}" --bulk-fragment "${deg}"
test -s "${GENETICS_STATE_INTERFACE_OUT}/matched_treat_by_arm.tsv"
test -s "${GENETICS_STATE_INTERFACE_OUT}/matched_pairs_by_arm.tsv"
test -s "${GENETICS_STATE_INTERFACE_OUT}/threshold_free_by_arm.tsv"

mkdir -p "${FIGURE_CANDIDATE_ROOT}"
export MASLD_PROJECT_ROOT="${project_root}"
export STAGE_RELEASE_ROOT="${BULK_ARM_ROOT}"
export GENETICS_STATE_INTERFACE_ROOT="${GENETICS_STATE_INTERFACE_OUT}"
export FIG3G_BULK_DEG="${deg}"
export FIG3G_EXPECT_ADOPTED_COUNTS=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
"${rnaseq_env}/bin/Rscript" --vanilla \
  "${project_root}/scripts/figures/fig3g_genetic_state_complementarity.R"
pdf="${FIGURE_CANDIDATE_ROOT}/figure3/panels/fig3g_genetic_state_complementarity.pdf"
supp_pdf="${FIGURE_CANDIDATE_ROOT}/figureS3/panels/figs3g_expression_matched_pairs.pdf"
for p in "${pdf}" "${supp_pdf}"; do
  test "$(pdfinfo "${p}" | awk '/^Pages:/ {print $2}')" = "1"
  if strings "${p}" | grep -q '/Subtype /Type3'; then exit 1; fi
  strings "${p}" | grep -q '/FontFamily (Helvetica)'
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${p}"
done

# The panel summary must agree with the interface tables it was drawn from, and
# the printed pair count must be the matched n.
"${rnaseq_env}/bin/python3" - "${FIGURE_CANDIDATE_ROOT}" "${GENETICS_STATE_INTERFACE_OUT}" "${supp_pdf}" <<'PY'
import csv, subprocess, sys
root, interface, supp = sys.argv[1:]
read = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
s = read(f"{root}/source_tables/fig3g_matched_treat_summary.tsv")[0]
ov = [r for r in read(f"{interface}/overlap_by_arm.tsv")
      if r["bulk_arm"] == "fragment_23370" and r["genetic_definition"] == "susie_pp4_0.5"][0]
m = [r for r in read(f"{interface}/matched_treat_by_arm.tsv")
     if r["bulk_arm"] == "fragment_23370" and r["genetic_definition"] == "susie_pp4_0.5"][0]
assert int(s["n_jointly_testable"]) == int(ov["n_joint_testable"])
assert int(s["n_genetically_anchored_joint"]) == int(ov["n_genetic_joint"])
for a, b in (("n_expression_matched_pairs", "n_pairs"), ("genetic_only_treat_pairs", "discordant_case_only"),
             ("matched_control_only_treat_pairs", "discordant_control_only"),
             ("neither_treat_pairs", "neither_treat"), ("both_treat_pairs", "both_treat")):
    assert int(s[a]) == int(m[b]), (a, s[a], m[b])
text = subprocess.check_output(["pdftotext", supp, "-"], text=True)
assert f"n = {int(m['n_pairs'])}" in text, "printed matched n differs from the interface"
print("\t".join(f"{k}={v}" for k, v in s.items()))
PY

mkdir -p "${FIGURE_CANDIDATE_ROOT}/proofs"
for p in "${pdf}" "${supp_pdf}"; do
  name=$(basename "${p}" .pdf)
  pdftotext "${p}" "${FIGURE_CANDIDATE_ROOT}/proofs/${name}.txt"
  gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r240 \
    -sOutputFile="${FIGURE_CANDIDATE_ROOT}/proofs/${name}.png" "${p}"
done
printf 'role\tpath\n' > "${FIGURE_CANDIDATE_ROOT}/inputs.tsv"
printf 'coloc_tier12\t%s\nbulk_deg\t%s\ninterface\t%s\n' "${coloc}" "${deg}" \
  "${GENETICS_STATE_INTERFACE_OUT}" >> "${FIGURE_CANDIDATE_ROOT}/inputs.tsv"
find "${FIGURE_CANDIDATE_ROOT}" -type f ! -name input_output_manifest.tsv -print0 | \
  sort -z | xargs -0 sha256sum > "${FIGURE_CANDIDATE_ROOT}/input_output_manifest.tsv"
sha256sum "${coloc}" "${deg}" >> "${FIGURE_CANDIDATE_ROOT}/input_output_manifest.tsv"
"${rnaseq_env}/bin/python3" \
  "${project_root}/scripts/manuscript/validate_resource_scope.py"
echo "FIG3G_COLOC_R2_COMPLETE ${FIGURE_CANDIDATE_ROOT}"
