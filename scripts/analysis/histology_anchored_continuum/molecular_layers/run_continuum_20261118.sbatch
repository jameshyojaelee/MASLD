#!/usr/bin/env bash
#SBATCH --job-name=hac_continuum_367_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/candidates/continuum_20261118_slurm_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/candidates/continuum_20261118_slurm_%j.err

set -euo pipefail
umask 027

ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
R_BIN="${RNASEQ_ENV}/bin/Rscript"
PY_BIN="${RNASEQ_ENV}/bin/python3"
export PATH="${RNASEQ_ENV}/bin:${PATH}"
export PYTHONNOUSERSITE=1

for name in $(compgen -e | grep -E '^(HAC_|KAMZOLAS_|MASLD_PROJECT_ROOT$)' || true); do
  unset "${name}"
done

ML_DIR="${ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers"
CAND_REL="${ROOT}/RNA-seq/results/manuscript_release/candidates"
REPRO="${CAND_REL}/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls"
STAGE_DIR="${ROOT}/figures/candidates/pi-stage-20261118/analysis/stage_extensions"
MANIFEST="${STAGE_DIR}/five_cohort_sample_manifest.tsv"
NMF_SCORES="${CAND_REL}/fixed-nmf-corrected-scores-20261009T141312Z"
SOURCE_CORE="${ROOT}/RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-corrected774-20261011T152106Z"
ATLAS="${ROOT}/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_atlas_20260824"

OUT_BASE="${ROOT}/scripts/analysis/histology_anchored_continuum/candidates/continuum_20261118"
MOL_OUT="${OUT_BASE}/molecular_layers"
SYSINF_OUT="${OUT_BASE}/system_inference"
STRESS_OUT="${OUT_BASE}/system_stress"
INP_DIR="${OUT_BASE}/inputs"
LOG_DIR="${OUT_BASE}/logs"
FIG_BASE="${ROOT}/figures/candidates/continuum-367-panels-20261118"
ATLAS_367="${FIG_BASE}/molecular_systems_atlas"
FIG_POOLED="${FIG_BASE}/molecular_systems_pooled"

for d in "${OUT_BASE}" "${FIG_BASE}"; do
  if [[ -e "${d}" ]]; then
    echo "Refusing to overwrite existing output dir: ${d}" >&2
    exit 1
  fi
done

mkdir -p "${MOL_OUT}" "${INP_DIR}" "${LOG_DIR}" "${FIG_BASE}"

status_file="${OUT_BASE}/stage_status.tsv"
printf 'stage\texit_code\tseconds\n' > "${status_file}"
stage() {
  local name="$1"; shift
  local t0=${SECONDS}
  echo "STAGE_START ${name} $(date -u +%FT%TZ)"
  "$@" > "${LOG_DIR}/${name}.log" 2>&1
  local rc=$?
  local dt=$(( SECONDS - t0 ))
  printf '%s\t%s\t%s\n' "${name}" "${rc}" "${dt}" >> "${status_file}"
  echo "STAGE_END ${name} rc=${rc} seconds=${dt}"
  return "${rc}"
}

# 1. Run default contract unit tests with zero override env vars set
stage tests_molecular "${R_BIN}" --vanilla "${ML_DIR}/tests/run_tests.R"

# 2. Build contract_override.json and system_context_roster.tsv from pi-stage-20261118
stage setup_inputs "${R_BIN}" --vanilla - "${REPRO}" "${MANIFEST}" "${STAGE_DIR}" "${INP_DIR}" <<'RS'
library(data.table)
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
repro_root <- normalizePath(args[[1L]], mustWork = TRUE)
manifest_path <- normalizePath(args[[2L]], mustWork = TRUE)
stage_dir <- normalizePath(args[[3L]], mustWork = TRUE)
out_dir <- normalizePath(args[[4L]], mustWork = TRUE)

expected <- fromJSON(file.path(repro_root, "expected_counts.json"), simplifyVector = TRUE)
model <- fread(file.path(repro_root, "model_input_manifest.tsv"))
validated <- fromJSON(file.path(repro_root, "VALIDATED.json"))
stopifnot(nrow(model) == 1L, isTRUE(as.logical(model$identity_qualified[[1L]])),
          identical(validated$status, "VALIDATED"),
          model$n_samples == 774L, model$n_genes == 24115L)

manifest <- fread(manifest_path, na.strings = c("", "NA"))
stopifnot(nrow(manifest) == 774L,
          manifest[dataset == "GSE162694" & !is.na(fibrosis_stage), .N] == 109L,
          manifest[dataset == "GSE213621" & !is.na(fibrosis_stage), .N] == 258L,
          manifest[dataset == "GSE213621" & group_binary == "Control" & !is.na(fibrosis_stage), .N] == 0L)

override <- list(
  dge_rds = model$dge_path[[1L]],
  canonical_deg_results = file.path(repro_root, "deg_results.csv"),
  sample_manifest = manifest_path,
  stage_results = file.path(stage_dir, "stage_all_gene_results.tsv"),
  stage_extension_results = file.path(stage_dir, "stage_extension_all_gene_results.tsv"),
  gene_family_size = as.integer(expected$n_genes),
  canonical_deg_size = as.integer(expected$canonical$n),
  participant_count = as.integer(expected$n_samples),
  signature_observed_size = 143L
)
stopifnot(all(file.exists(unlist(override[1:5]))))
write_json(override, file.path(out_dir, "contract_override.json"), auto_unbox = TRUE, pretty = TRUE)

roster <- manifest[, .(sample_id, dataset, analysis_unit_id, unit_evidence)]
stopifnot(!anyNA(roster), all(nzchar(roster$unit_evidence)),
          !anyDuplicated(roster$sample_id), !anyDuplicated(roster$analysis_unit_id))
fwrite(roster, file.path(out_dir, "system_context_roster.tsv"), sep = "\t")
cat("CONT_367_INPUTS_READY staged_eval=367\n")
RS

# 3. Export molecular layers env
export MASLD_PROJECT_ROOT="${ROOT}"
export HAC_ML_OUT_ROOT="${MOL_OUT}"
export HAC_SOURCE_ROOT="${SOURCE_CORE}"
export HAC_CONTRACT_OVERRIDE="${INP_DIR}/contract_override.json"
export HAC_FIXED_NMF_SCORES="${NMF_SCORES}"
export HAC_ML_WORKERS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

stage mol_01 "${R_BIN}" --vanilla "${ML_DIR}/01_freeze_inputs.R"

# Run bulk (10), pathway (20), hotspot (31) in parallel (respecting 8-CPU allocation)
stage mol_10 "${R_BIN}" --vanilla "${ML_DIR}/10_bulk_transcript_deg.R" & pid_10=$!
stage mol_20 "${R_BIN}" --vanilla "${ML_DIR}/20_pathway_continuum.R" & pid_20=$!
stage mol_31 "${R_BIN}" --vanilla "${ML_DIR}/31_hotspot_continuum.R" & pid_31=$!
for pid in "${pid_10}" "${pid_20}" "${pid_31}"; do wait "${pid}"; done

# Run TF (21), NMF with residualized diagnostics (30), paired (40) in parallel
stage mol_21 "${R_BIN}" --vanilla "${ML_DIR}/21_tf_continuum.R" & pid_21=$!
stage mol_30 "${R_BIN}" --vanilla "${ML_DIR}/30_nmf_continuum.R" & pid_30=$!
stage mol_40 "${R_BIN}" --vanilla "${ML_DIR}/40_paired_validation.R" & pid_40=$!
for pid in "${pid_21}" "${pid_30}" "${pid_40}"; do wait "${pid}"; done

stage mol_50 "${R_BIN}" --vanilla "${ML_DIR}/50_integrate_gates.R"

# 4. System inference (70) -> Stress test (73) + Comparison figure (62)
export HAC_SYSTEM_CONTEXT_ROSTER="${INP_DIR}/system_context_roster.tsv"
export HAC_CORRECTED_STAGE_ROOT="${REPRO}"
stage sys_70 "${R_BIN}" --vanilla "${ML_DIR}/70_infer_molecular_systems.R" "${MOL_OUT}" "${ATLAS}" "${SYSINF_OUT}"

stage stress_73 "${R_BIN}" --vanilla "${ML_DIR}/73_stress_test_molecular_systems.R" \
  "${MOL_OUT}" "${ATLAS}" "${SYSINF_OUT}" "${STRESS_OUT}" & pid_73=$!
stage fig_62 "${R_BIN}" --vanilla "${ML_DIR}/62_render_stage_continuum_comparison.R" \
  "${MOL_OUT}" "${FIG_BASE}/stage_continuum_comparison" & pid_62=$!
for pid in "${pid_73}" "${pid_62}"; do wait "${pid}"; done

# 5. Prepare frozen geometry copy under figures/candidates/ and run 67 + 68 (Fig S4U)
prepare_atlas() {
  mkdir -p "${ATLAS_367}/provenance"
  cp -a "${ATLAS}/geometry" "${ATLAS_367}/geometry"
  cp -a "${ATLAS}/source_tables" "${ATLAS_367}/source_tables"
  cp -a "${ATLAS}/provenance"/geometry_* "${ATLAS}/provenance/outcome_blind_audit.tsv" "${ATLAS_367}/provenance/"
  diff -rq "${ATLAS}/geometry" "${ATLAS_367}/geometry"
}
stage fig_prep_atlas prepare_atlas
stage fig_67 "${R_BIN}" --vanilla "${ML_DIR}/67_overlay_molecular_systems_atlas.R" \
  "${MOL_OUT}" "${FIG_BASE}/stage_continuum_comparison" "${ATLAS_367}"

write_checksums() {
  ( cd "$1" && printf 'path\tsha256\n' > promoted_source_checksums.tsv &&
    find . -type f ! -name promoted_source_checksums.tsv | sort | while read -r f; do
      printf '%s\t%s\n' "${1}/${f#./}" "$(sha256sum "${f}" | awk '{print $1}')"
    done >> promoted_source_checksums.tsv )
}
stage fig_67_checksums write_checksums "${ATLAS_367}"
stage fig_68_s4u "${R_BIN}" --vanilla "${ML_DIR}/68_render_pooled_molecular_systems.R" \
  "${ATLAS_367}" "${FIG_POOLED}"

# Mirror Fig S4U outputs cleanly inside continuum_20261118 and render PNG for visual verification
mkdir -p "${OUT_BASE}/fig_s4u_pooled"
cp -a "${FIG_POOLED}/." "${OUT_BASE}/fig_s4u_pooled/"
PANEL="${OUT_BASE}/fig_s4u_pooled/panels/s4_molecular_systems_pooled_continuum_windows.pdf"
pdfinfo "${PANEL}" > "${OUT_BASE}/fig_s4u_pooled/pdfinfo.txt"
pdftotext "${PANEL}" - > "${OUT_BASE}/fig_s4u_pooled/pdftotext.txt"
pdftoppm -png -r 150 "${PANEL}" "${OUT_BASE}/fig_s4u_pooled/s4u_preview"

echo "CONTINUUM_20261118_ALL_DONE ${OUT_BASE}"
