#!/usr/bin/env bash
#SBATCH --job-name=model-check-073
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
target="${root}/executions/model-check-073-${SLURM_JOB_ID}"
stage="${target}.staging"
failed="${target}.failed"

preserve_failed_stage() {
  if [[ -d "${stage}" && ! -e "${failed}" ]]; then
    mv "${stage}" "${failed}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${target}" || -e "${stage}" || -e "${failed}" ]]; then
  printf 'Refusing to overwrite sequence terminal registry validation output.\n' >&2
  exit 1
fi
install -d -m 0750 "${stage}/validation"

sha256sum --check --strict <<EOF
8c6684a23ef40e6573c788fdfc9d697668f5f5fc978390118dc8f12f9c79ae1d  ${root}/config/models/regulatory_local.toml
cd7e67aa068c2a041fbf842564fc20d82be008020bb5226d70bc36c937e4a74d  ${root}/config/models/mandatory_baselines.toml
c779f06dedf69a15ec44280323a7701e0938b75bfa2660d11dd2a6f5de37868d  ${root}/config/evaluation/rna_conditioned_atac_capabilities.toml
54204d77794e2a1f3be78b848741bd3b2c6bf2aa136aca170dc3194dad6d054a  ${root}/config/evaluation/sequence_terminal_admission.json
d82c761e39e48c2488fbbd2fff45d27f0ac9782be2a80f2f7d0229f8d18428c2  ${root}/config/artifacts/models/bpnet/checkpoints.json
cc31df3f9172b261b9188edf3e88bd1326ec1e23a1202191a791fe0941ea3d51  ${root}/config/artifacts/models/bpnet/exposure_audit.json
75c435017ee2e7a3a7504cbd748e85f29a6430bf7820434832e10590b7e896d4  ${root}/config/artifacts/models/regulatory_native_controls/checkpoints.json
ae3e6866e8a4c209ed321ad1816fe035f419cd17d2e198e91663c176482cad4b  ${root}/config/artifacts/models/regulatory_native_controls/exposure_audit.json
f5a65ee9c0a60cdaed67615574a6c01fcf3bcf4bc70b570c0a1ad4dd926893ba  ${root}/config/artifacts/model-authorities/bpnet/cellxgene_uuids.json
4024e36e988b91aa6f87d06bbabe969aeb6224e49a35df48889b930317fa42f3  ${root}/config/artifacts/model-authorities/bpnet/code_license.json
6b74847e411737856fdb968d44b8bd6893d10a65650cb561fa251c2f6faa24e0  ${root}/config/artifacts/model-authorities/bpnet/declared_corpora.json
e032fb6285f224824ec736d60b0d88542f3f942c3f90a12fc740cc67c9e50997  ${root}/config/artifacts/model-authorities/bpnet/derivative_weights_license.json
c084eab18d910fc30093a61ba6857077c9394d3136b7cc0a8ab404019619447a  ${root}/config/artifacts/model-authorities/bpnet/exposure_audit.json
921ab84f00cc2944120661bc8adf4c847f61311a2873b41a2902f22cc63f0963  ${root}/config/artifacts/model-authorities/bpnet/preprocessing.json
ad40306a64ac78bbce1d39f17395a6640bb2962d45b609279085b32f6e09990a  ${root}/config/artifacts/model-authorities/bpnet/weights_license.json
34c72c5050db1ec4d8aad3ca4975afdf334a88767de169c6b2ffe4009c4a93d6  ${root}/tests/contract/test_registry_files.py
c648e18c69762287ea5863c2aef97a9496307def35793d38aa1b9513426db009  ${root}/tests/contract/test_regulatory_lane_preconditions.py
90f447bffa80248eb79908b5cfcebec6c36fd808de1dda023f21aaf774b68b7a  ${root}/tests/contract/test_sequence_terminal_admission.py
b379da6299a7fd9b42c77da7da0150403a6268683cc8ecd1a35dd3ad983a33a6  ${root}/executions/sequence_gpu_bundle_016_sequence_cnn_control-ready_36823fe74c-21076966/ARTIFACTS.json
14af86ac22e67bc29f8a7b85ee3a049b63f875ba5423210cc8e539d29b9fe9e8  ${root}/executions/sequence_gpu_bundle_004_bpnet-ready_71d22f6453-21077605/ARTIFACTS.json
EOF

module load python/3.11.5-GCCcore-13.2.0
python --version >"${stage}/validation/python.version.txt" 2>&1
python -m pip freeze >"${stage}/validation/pip-freeze.txt"
PYTHONPATH="${root}/src:${root}" python -m unittest -v \
  tests.contract.test_sequence_terminal_admission \
  tests.contract.test_registry_files.RegistryFileContractTest.test_model_roster_is_complete_unique_and_fail_closed \
  tests.contract.test_registry_files.RegistryFileContractTest.test_bpnet_preflight_is_matched_local_sequence_control \
  tests.contract.test_registry_files.RegistryFileContractTest.test_atac_profile_and_latent_capabilities_are_disjoint \
  tests.contract.test_regulatory_lane_preconditions \
  >"${stage}/validation/unittest.stdout.txt" \
  2>"${stage}/validation/unittest.stderr.txt"

PYTHONPATH="${root}/src" python - "${stage}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
summary = {
    "schema_version": "masld-bench-sequence-terminal-registry-validation-v1",
    "status": "pass",
    "tests": 10,
    "gpu_jobs_submitted": 0,
    "frozen_artifacts_mutated": False,
    "bpnet_development_profile_admitted": True,
    "cnn_registry_admitted": False,
    "external_evaluation_complete": False,
    "sealed_evaluation_complete": False,
    "champion_eligible": False,
}
(stage / "validation/summary.json").write_text(
    json.dumps(summary, sort_keys=True, indent=2) + "\n",
    encoding="utf-8",
)
freeze_tree(stage, {"artifact_class": "sequence_terminal_registry_validation", **summary})
PY

trap - EXIT
mv "${stage}" "${target}"
sha256sum "${target}/ARTIFACTS.json"
