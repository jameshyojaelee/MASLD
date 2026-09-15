#!/usr/bin/env bash
#SBATCH --job-name=model-preflight
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
RUNTIME="${ROOT}/executions/environments/caduceus-native-torch2.2.0-21065466"
CHECKPOINTS="${ROOT}/executions/dna-lm-open-safetensors-21064556"
CODE="${ROOT}/executions/dna-lm-custom-code-21064626"
COMMON_PROBE="${ROOT}/executions/caduceus-native-probe-21065951"
LONG_PROBE="${ROOT}/executions/caduceus-131k-probe-21066018"
FIXTURE="${ROOT}/executions/gse281364-caduceus-131k-fixture-21071004"
FASTA="${ROOT}/executions/sequence-fasta-21063808"
TARGET="${ROOT}/executions/caduceus-131k-gse281364-preflight-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
RUNTIME_TMP=$(mktemp -d "${TMPDIR:-/tmp}/model-preflight-${SLURM_JOB_ID}.XXXXXX")
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
  if [[ -d "${RUNTIME_TMP}" ]]; then rm -r "${RUNTIME_TMP}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Caduceus preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for binding in \
  "3d59792d8a13ca10a2df3d6386bc266eb0e108def6279188793d1cfee4fc1167 ${RUNTIME}" \
  "6d21df999835624905e08573038dc72e893a5b25a8ca8f4ec1acea8f61d91c77 ${CHECKPOINTS}" \
  "70906faae158909a8972f674f680943d8a51467032ae3266a6a12f2152beb36c ${CODE}" \
  "6bc5288848ca6cd696cc90c36ceea4e789e97b7f11ceb70670526bd0ee696ceb ${COMMON_PROBE}" \
  "a9c531b8a044c8e7990e24cbafcb7b035d51d5861efccb7c2c326c66f41b5456 ${LONG_PROBE}" \
  "82b47ba2059193cae23e35db5401013659879024052c47558e602c03b4ee2c70 ${FIXTURE}" \
  "7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 ${FASTA}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" | sha256sum --check --strict
done
for binding in \
  "0c3e50fb5dc5a41c3686caf638416106b523a1e059a28c50fd5a75c2e37a0dbf ${ROOT}/config/gse281364_caduceus_131k_contract.json" \
  "e3d2cb3b7ccc209feb77d3cbd75e6d68e270be96d5f3cfcfef5de1dab1b618de ${ROOT}/scripts/caduceus_predict_gse281364_131k.py" \
  "6405c4dbc188e8d904b8e5d5afba68801e824773ac06becd462dd88b4327821f ${ROOT}/scripts/project_gse281364_caduceus_131k_embeddings.py" \
  "3aa9c27842c351af29d68f7e8642815e442eca26a1b470ccadffd53dd2d5876d ${ROOT}/tests/unit/test_caduceus_predict_gse281364_131k.py"; do
  read -r expected path <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RUNTIME}" "${CHECKPOINTS}" "${CODE}" \
  "${COMMON_PROBE}" "${LONG_PROBE}" "${FIXTURE}" "${FASTA}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY
python - "${ROOT}" "${STAGE}/validation/admission.json" <<'PY'
from hashlib import sha256
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
config = json.loads((root / "config/gse281364_caduceus_131k_contract.json").read_text())
checkpoint = json.loads((root / "config/artifacts/models/caduceus/checkpoints.json").read_text())
exposure = json.loads((root / "config/artifacts/models/caduceus/exposure_audit.json").read_text())
crosswalk = json.loads((root / "config/artifacts/models/caduceus/development_crosswalk.json").read_text())
if (
    checkpoint["terms_audit"]["result"] != "terms_pass_preflight"
    or checkpoint["checkpoint"]["sha256"] != config["checkpoint"]["sha256"]
    or checkpoint["checkpoint"]["size_bytes"] != 30937760
    or checkpoint["code"]["license"] != "Apache-2.0"
    or checkpoint["weight_license"] != "Apache-2.0"
    or exposure["checkpoint_finding"]["exposure_state"] != "target_label_unexposed"
    or exposure["checkpoint_finding"]["reference_sequence_exposure"] != "encoder_seen"
    or crosswalk["findings"]["long_range_eQTL_VEP_labels"]["exposure_state"] != "downstream_demo"
    or config["exposure"]["allowed_supervision"] != "locus-cross-fitted_assay-native_MPRA_activity_only"
    or config["exposure"]["prohibited_supervision"] != "cell-type_eQTL_or_ieQTL_effect_head"
    or any(config["firewall"].values())
    or not config["evo_family_disposition"].startswith("deferred_exact_Vortex_runtime_revision_unresolved")
):
    raise SystemExit("Caduceus admission, exposure, or Evo disposition differs")
result = {
    "status": "pass",
    "checkpoint_sha256": config["checkpoint"]["sha256"],
    "code_and_weight_license": "Apache-2.0",
    "exposure_state": "target_label_unexposed",
    "reference_sequence_exposure": "encoder_seen",
    "independent_eQTL_classifier_exposure": "downstream_demo",
    "allowed_head": config["exposure"]["allowed_supervision"],
    "prohibited_head": config["exposure"]["prohibited_supervision"],
    "evo_family": config["evo_family_disposition"],
    "project_outcomes_read": False,
    "sealed_outcomes_read": False,
  }
Path(sys.argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
PY
module purge
tar -xzf "${RUNTIME}/payload/caduceus-native-env.tar.gz" -C "${RUNTIME_TMP}"
cd "${ROOT}"
PYTHONPATH="${ROOT}" "${RUNTIME_TMP}/bin/python" -m unittest -v \
  tests.unit.test_caduceus_predict_gse281364_131k \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}" "${RUNTIME_TMP}/bin/python" scripts/caduceus_predict_gse281364_131k.py --help \
  > "${STAGE}/validation/predictor-help.stdout.txt" \
  2> "${STAGE}/validation/predictor-help.stderr.txt"
"${RUNTIME_TMP}/bin/python" -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  config/gse281364_caduceus_131k_contract.json \
  config/artifacts/models/caduceus/checkpoints.json \
  config/artifacts/models/caduceus/development_crosswalk.json \
  config/artifacts/models/caduceus/exposure_audit.json \
  scripts/caduceus_predict_gse281364_131k.py \
  scripts/project_gse281364_caduceus_131k_embeddings.py \
  tests/unit/test_caduceus_predict_gse281364_131k.py \
  slurm/preflight_gse281364_caduceus_131k_cpu.sbatch \
  > "${STAGE}/source.sha256"
rm -r "${RUNTIME_TMP}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
tests = (root / "validation/tests.stderr.txt").read_text()
admission = json.loads((root / "validation/admission.json").read_text())
if (
    "Ran 3 tests" not in tests
    or "OK" not in tests
    or admission["status"] != "pass"
    or admission["project_outcomes_read"]
    or admission["sealed_outcomes_read"]
    or (root / "validation/predictor-help.stderr.txt").stat().st_size != 0
):
    raise SystemExit("Caduceus campaign preflight differs")
freeze_tree(root, {
    "artifact_class": "gse281364_caduceus_131k_outcome_blind_preflight",
    "dataset_id": "gse281364",
    "model_id": "caduceus",
    "checkpoint_sha256": "a3e6976fe90460ff5d90d371b457d0263dc65e156ea5651b4a452e98228daef2",
    "license": "Apache-2.0",
    "exposure_state": "target_label_unexposed",
    "fixture_artifacts_sha256": "82b47ba2059193cae23e35db5401013659879024052c47558e602c03b4ee2c70",
    "synthetic_tests": 3,
    "project_outcomes_read": False,
    "sealed_outcomes_read": False,
    "evo_family": "deferred_runtime_unresolved",
    "status": "passed",
})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Caduceus GSE281364 preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
