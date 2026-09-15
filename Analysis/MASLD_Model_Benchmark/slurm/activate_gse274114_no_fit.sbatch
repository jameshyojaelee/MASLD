#!/usr/bin/env bash
#SBATCH --job-name=model-data-080
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/bulk-expansion-source-audit-21065762
SOURCE_SHA256=883f870e2a76f78c39d257c24f8e795d2e74f1bd32aba785251b52aa7cb8229e
MATRIX=${ROOT}/executions/bulk-expansion-matrix-audit-21065853
MATRIX_SHA256=3f934ab9c6ee730c6bd0887b10d963f8210f274c36dc0878cb814c85026dc038
PRIOR=${ROOT}/executions/gse274114-activation-21066117
PRIOR_SHA256=f429021ebe8f742314a2073747267047447d07398a3e34eb799b3728c9d829d0
GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
GTF_SHA256=73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9
PYTHON=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python

TARGET=${ROOT}/executions/${SLURM_JOB_NAME}-${SLURM_JOB_ID}-gse274114
STAGE=${TARGET}.staging
SOURCE_STAGE=${TARGET}.source.staging
PROMOTION_PENDING=0
preserve_failed() {
  for path in "${STAGE}" "${SOURCE_STAGE}"; do
    if [[ -e "${path}" && ! -e "${path}.failed" ]]; then mv "${path}" "${path}.failed"; fi
  done
  if [[ "${PROMOTION_PENDING}" -eq 1 && -e "${TARGET}" && ! -e "${TARGET}.failed" ]]; then
    mv "${TARGET}" "${TARGET}.failed"
  fi
}
trap preserve_failed EXIT
for path in "${TARGET}" "${STAGE}" "${SOURCE_STAGE}"; do
  if [[ -e "${path}" || -e "${path}.failed" ]]; then
    printf 'Refusing to overwrite GSE274114 activation: %s\n' "${path}" >&2
    exit 1
  fi
done

cd "${ROOT}"
printf '%s  %s\n' \
  "${SOURCE_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  "${MATRIX_SHA256}" "${MATRIX}/ARTIFACTS.json" \
  "${PRIOR_SHA256}" "${PRIOR}/ARTIFACTS.json" \
  "${GTF_SHA256}" "${GTF}" \
  | sha256sum --check --strict
module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
python - "${SOURCE}" "${MATRIX}" "${PRIOR}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY

python -m py_compile scripts/activate_gse274114_no_fit.py
"${PYTHON}" -m unittest -v tests.unit.test_activate_gse274114_no_fit
python -m unittest -v tests.contract.test_gse274114_activation_contract

mkdir "${SOURCE_STAGE}"
curl -fsSL --retry 3 \
  'https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA1144975' \
  -o "${SOURCE_STAGE}/PRJNA1144975.runinfo.csv"
curl -fsSL --retry 3 \
  'https://pure.eur.nl/ws/files/180988305/Metabolic_dysfunction-associated_steatohepatitis_reduces_interferon_and_macrophage_liver_gene_signatures_in_patients_with_chronic_hepatitis_B.pdf' \
  -o "${SOURCE_STAGE}/article.pdf"
curl -fsSL --retry 3 \
  'https://ars.els-cdn.com/content/image/1-s2.0-S0168827824026539-mmc1.pdf' \
  -o "${SOURCE_STAGE}/supplement.pdf"
printf '%s  %s\n' \
  '090bd763671d3cecd1d96dfc3d9e8237944ac09e12744c69deff5d7bfc56ba7e' "${SOURCE_STAGE}/PRJNA1144975.runinfo.csv" \
  '809e671b6b23deaa52edc966fa9487d37b92e31e9c2f6e593d0fe0e637e25fb1' "${SOURCE_STAGE}/article.pdf" \
  'db1a82aca87974d6003f019014b80ae600882167b23ec50f74c56789cbd11953' "${SOURCE_STAGE}/supplement.pdf" \
  | sha256sum --check --strict
pdftotext -layout "${SOURCE_STAGE}/article.pdf" "${SOURCE_STAGE}/article.txt"

COMPARISON_ARGS=()
for path in config/datasets/*.toml; do
  if [[ "$(basename "${path}")" != 'gse274114_mash_hbv.toml' ]]; then
    COMPARISON_ARGS+=(--comparison-path "${path}")
  fi
done
for path in \
  executions/gse267145-authoritative-join-21064930 \
  executions/gse260666-activation-21066072 \
  executions/gse268273-activation-21066074 \
  executions/gse105127-authoritative-join-rights-21065552 \
  executions/model-data-072-21080651-participants \
  executions/gse296875-donor-lineage-fragments-21063889/fragment_manifest.tsv \
  executions/gse256398-qc-membership-21069669/retained_membership.tsv; do
  COMPARISON_ARGS+=(--comparison-path "${path}")
done

"${PYTHON}" scripts/activate_gse274114_no_fit.py \
  --source "${SOURCE}" --source-artifacts-sha256 "${SOURCE_SHA256}" \
  --matrix "${MATRIX}" --matrix-artifacts-sha256 "${MATRIX_SHA256}" \
  --prior-activation "${PRIOR}" --prior-activation-artifacts-sha256 "${PRIOR_SHA256}" \
  --gtf "${GTF}" --gtf-sha256 "${GTF_SHA256}" \
  --runinfo "${SOURCE_STAGE}/PRJNA1144975.runinfo.csv" \
  --runinfo-sha256 090bd763671d3cecd1d96dfc3d9e8237944ac09e12744c69deff5d7bfc56ba7e \
  --article-pdf "${SOURCE_STAGE}/article.pdf" \
  --article-pdf-sha256 809e671b6b23deaa52edc966fa9487d37b92e31e9c2f6e593d0fe0e637e25fb1 \
  --supplement-pdf "${SOURCE_STAGE}/supplement.pdf" \
  --supplement-pdf-sha256 db1a82aca87974d6003f019014b80ae600882167b23ec50f74c56789cbd11953 \
  --article-text "${SOURCE_STAGE}/article.txt" \
  "${COMPARISON_ARGS[@]}" \
  --output "${STAGE}"

mv "${SOURCE_STAGE}" "${STAGE}/source_evidence"
"${PYTHON}" -m pip freeze > "${STAGE}/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/activate_gse274114_no_fit.py" \
  "${ROOT}/tests/unit/test_activate_gse274114_no_fit.py" \
  "${ROOT}/tests/contract/test_gse274114_activation_contract.py" \
  "${ROOT}/slurm/activate_gse274114_no_fit.sbatch" \
  "${ROOT}/config/datasets/gse274114_mash_hbv.toml" \
  "${ROOT}/config/evaluation/gse274114_within_instrument_etiology_task.toml" \
  "${ROOT}/config/evaluation/gse274114_within_instrument_etiology_promotion_gate.json" \
  "${SOURCE}/ARTIFACTS.json" "${MATRIX}/ARTIFACTS.json" "${PRIOR}/ARTIFACTS.json" \
  "${GTF}" \
  > "${STAGE}/source.sha256"
(
  cd "${STAGE}"
  sha256sum \
    source_evidence/PRJNA1144975.runinfo.csv \
    source_evidence/article.pdf \
    source_evidence/article.txt \
    source_evidence/supplement.pdf \
    >> source.sha256
  sha256sum --check --strict source.sha256
)

python - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import gzip
from hashlib import sha256
import json
import math
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "activation_audit.json").read_text(encoding="utf-8"))
if (
    audit.get("participants") != 39
    or audit.get("sra_runs") != 71
    or audit.get("model_input_labels_separated") is not True
    or audit.get("four_class_performance_allowed") is not False
    or audit.get("MASH_vs_non_MASH_performance_allowed") is not False
    or audit.get("HBV_only_is_MASLD_negative_control") is not False
    or audit.get("normalization_run") is not False
    or audit.get("model_fit_or_scoring_run") is not False
):
    raise SystemExit("GSE274114 no-fit activation firewall differs")
with (root / "model_inputs" / "quant_manifest.tsv").open(newline="", encoding="utf-8") as handle:
    model_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(model_rows) != 39 or any("group" in key.lower() for key in model_rows[0]):
    raise SystemExit("GSE274114 model-input label separation differs")
if any(
    row["gene_axis"] != "GENCODE_v49_60324_unique_targets"
    or row["par_y_x_aggregation_applied"] != "True"
    or row["unmapped_source_rows_masked"] != "1230"
    or row["normalization_or_fit_applied"] != "False"
    or row["labels_present"] != "False"
    for row in model_rows
):
    raise SystemExit("GSE274114 materialized v49 model-input contract differs")
for row in model_rows:
    quant = (root / row["quant_path"]).resolve(strict=True)
    if root not in quant.parents:
        raise SystemExit("GSE274114 quant path escaped the frozen fixture")
    digest = sha256(quant.read_bytes()).hexdigest()
    if digest != row["quant_sha256"]:
        raise SystemExit("GSE274114 quant checksum differs")
    with gzip.open(quant, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        names = [record["Name"] for record in reader]
    if len(names) != 60324 or len(set(names)) != 60324 or any("_PAR_Y" in name for name in names):
        raise SystemExit("GSE274114 materialized v49 quant axis differs")
with (root / "authorities" / "quant_aggregation.tsv").open(newline="", encoding="utf-8") as handle:
    aggregation_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(aggregation_rows) != 39 or any(
    row["source_rows"] != "61598"
    or row["mapped_source_rows"] != "60368"
    or row["output_v49_genes"] != "60324"
    or not math.isclose(float(row["mapped_numreads"]), float(row["output_numreads"]), rel_tol=1e-12, abs_tol=1e-8)
    or not math.isclose(float(row["mapped_tpm"]), float(row["output_tpm"]), rel_tol=1e-12, abs_tol=1e-8)
    for row in aggregation_rows
):
    raise SystemExit("GSE274114 v49 aggregation conservation differs")
with (root / "authorities" / "participant_join.tsv").open(newline="", encoding="utf-8") as handle:
    participant_rows = list(csv.DictReader(handle, delimiter="\t"))
with (root / "authorities" / "run_join.tsv").open(newline="", encoding="utf-8") as handle:
    run_rows = list(csv.DictReader(handle, delimiter="\t"))
if (
    len(participant_rows) != 39
    or sum(int(row["technical_sra_run_count"]) for row in participant_rows) != 71
    or len(run_rows) != 71
    or any(row["technical_run_is_biological_replicate"] != "False" for row in run_rows)
):
    raise SystemExit("GSE274114 participant/technical-run topology differs")
with (root / "authorities" / "exact_identifier_overlap.tsv").open(newline="", encoding="utf-8") as handle:
    overlaps = list(csv.DictReader(handle, delimiter="\t"))
if len(overlaps) != 6 or any(int(row["exact_overlap_count"]) for row in overlaps):
    raise SystemExit("GSE274114 exact identifier overlap audit differs")
freeze_tree(root, {
    "artifact_class": "gse274114_no_fit_etiology_activation",
    "cohort_family_ids": ["gse274114_mash_hbv"],
    "participants": 39,
    "sra_technical_runs": 71,
    "supported_contrasts": ["hiseq_CTRL_vs_ENEG", "novaseq_NASH_vs_ENEG_NASH"],
    "four_class_performance_allowed": False,
    "MASH_vs_non_MASH_performance_allowed": False,
    "normalization_run": False,
    "model_fit_or_scoring_run": False,
    "status": "passed_label_separated_no_fit_fixture",
})
verify_frozen_tree(root)
PY
mv "${STAGE}" "${TARGET}"
PROMOTION_PENDING=1

(
  cd "${TARGET}"
  sha256sum --check --strict source.sha256
)
python - "${TARGET}" "${ROOT}/executions/model-data-075-21081936-gse274114" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree

current = Path(sys.argv[1]).resolve(strict=True)
previous = Path(sys.argv[2]).resolve(strict=True)
verify_frozen_tree(current)
current_manifest = json.loads((current / "ARTIFACTS.json").read_text(encoding="utf-8"))
previous_manifest = json.loads((previous / "ARTIFACTS.json").read_text(encoding="utf-8"))
if current_manifest["metadata"] != previous_manifest["metadata"]:
    raise SystemExit("GSE274114 scientific artifact metadata changed in provenance-only revision")
current_artifacts = {
    row["path"]: (row["sha256"], row["size_bytes"])
    for row in current_manifest["artifacts"]
    if row["path"] != "source.sha256"
}
previous_artifacts = {
    row["path"]: (row["sha256"], row["size_bytes"])
    for row in previous_manifest["artifacts"]
    if row["path"] != "source.sha256"
}
if current_artifacts != previous_artifacts:
    raise SystemExit("GSE274114 non-provenance artifacts changed in provenance-only revision")
print(json.dumps({
    "post_move_frozen_tree_verified": True,
    "post_move_source_lock_replayed": True,
    "non_provenance_artifacts_byte_identical_to_21081936": len(current_artifacts),
}, sort_keys=True))
PY
PROMOTION_PENDING=0

trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json"
