#!/usr/bin/env bash
#SBATCH --job-name=masld-scooby-bounded-headers
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
REVIEW="${ROOT}/executions/scooby-cross-review-21064420"
REVIEW_ARTIFACTS_SHA256=27668bff0943a8baac6489f073e000c2ba181e729a6f3360ea1dc0a3cfd6aa7c
ARCHIVES="${ROOT}/executions/scooby-released-archive-inventory-21064432"
ARCHIVES_ARTIFACTS_SHA256=8cb983369ff43a93fbd2bc96bd6683613ce7d3cdf45efd6a4efe1c21715d2aa7
TARGET="${ROOT}/executions/scooby-bounded-headers-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite scooby bounded-header artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for binding in \
  "${REVIEW_ARTIFACTS_SHA256} ${REVIEW}/ARTIFACTS.json" \
  "${ARCHIVES_ARTIFACTS_SHA256} ${ARCHIVES}/ARTIFACTS.json"; do
  printf '%s\n' "${binding}" | sha256sum --check --strict
done

install -d -m 0750 "${STAGE}/metadata" "${STAGE}/headers" "${STAGE}/validation"
curl_flags=(--fail --location --silent --show-error --retry 5 --retry-all-errors)

ONEK_REPO=lauradmartens/onek1k-scooby
ONEK_REV=dd7a2b7be54c906b43de4bc0159b169caaf057c5
EPIC_REPO=lauradmartens/epicardioids-scooby
EPIC_REV=cefc3a626f5080baaf83878754ea924df6bd7e8d
NEURIPS_REPO=johahi/neurips-scooby
NEURIPS_REV=1591d15696c65f14053e4495e52c74812a86fd8d

for spec in \
  "onek1k ${ONEK_REPO} ${ONEK_REV}" \
  "epicardioids ${EPIC_REPO} ${EPIC_REV}" \
  "neurips ${NEURIPS_REPO} ${NEURIPS_REV}"; do
  read -r model repository revision <<< "${spec}"
  curl "${curl_flags[@]}" \
    "https://huggingface.co/api/models/${repository}/revision/${revision}?blobs=true" \
    --output "${STAGE}/metadata/${model}.api.json"
  curl "${curl_flags[@]}" \
    "https://huggingface.co/${repository}/raw/${revision}/model.safetensors" \
    --output "${STAGE}/metadata/${model}.model.safetensors.lfs_pointer"
  curl "${curl_flags[@]}" \
    "https://huggingface.co/${repository}/resolve/${revision}/config.json" \
    --output "${STAGE}/metadata/${model}.config.json"
done
for spec in \
  "onek1k ${ONEK_REPO} ${ONEK_REV}" \
  "epicardioids ${EPIC_REPO} ${EPIC_REV}"; do
  read -r model repository revision <<< "${spec}"
  curl "${curl_flags[@]}" \
    "https://huggingface.co/${repository}/resolve/${revision}/README.md" \
    --output "${STAGE}/metadata/${model}.README.md"
done
neurips_readme_status=$(curl --location --silent --show-error \
  --output "${STAGE}/validation/neurips.README.response" \
  --write-out '%{http_code}' \
  "https://huggingface.co/${NEURIPS_REPO}/resolve/${NEURIPS_REV}/README.md")
if [[ "${neurips_readme_status}" != 404 ]]; then
  printf 'NeurIPS immutable checkpoint card status changed; reassess terms\n' >&2
  exit 1
fi
rm -- "${STAGE}/validation/neurips.README.response"

module load python/3.11.5-GCCcore-13.2.0
python "${ROOT}/scripts/scooby_remote_safetensors_header.py" \
  --url "https://huggingface.co/${ONEK_REPO}/resolve/${ONEK_REV}/model.safetensors" \
  --expected-size 765092032 \
  --expected-parameters 191267309 \
  --reported-parameters 191267298 \
  --expected-sha256 2ae335c651182040462d965bb879118be40a6cb174550a3686fb10b04a74f7cb \
  --output "${STAGE}/headers/onek1k.header.json" \
  > "${STAGE}/validation/onek1k.header.stdout.json"
python "${ROOT}/scripts/scooby_remote_safetensors_header.py" \
  --url "https://huggingface.co/${EPIC_REPO}/resolve/${EPIC_REV}/model.safetensors" \
  --expected-size 772988628 \
  --expected-parameters 193241454 \
  --reported-parameters 193236835 \
  --expected-sha256 f945e1937a46ddd5cb31249ef71e721e7a7b1381aa207fb2de449f23ac75a9e8 \
  --output "${STAGE}/headers/epicardioids.header.json" \
  > "${STAGE}/validation/epicardioids.header.stdout.json"

python - "${STAGE}" "${REVIEW}" "${ARCHIVES}" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
review_root = Path(sys.argv[2]).resolve(strict=True)
archive_root = Path(sys.argv[3]).resolve(strict=True)
review = json.loads((review_root / "scooby_cross_review_receipt.json").read_text())
training = json.loads(
    (archive_root / "inventories/training_data.inventory.json").read_text()
)
submission = json.loads(
    (archive_root / "inventories/submission_data.inventory.json").read_text()
)
expected = {
    "onek1k": {
        "repository": "lauradmartens/onek1k-scooby",
        "revision": "dd7a2b7be54c906b43de4bc0159b169caaf057c5",
        "size": 765092032,
        "sha256": "2ae335c651182040462d965bb879118be40a6cb174550a3686fb10b04a74f7cb",
        "config_sha256": "e4802fb5c2d5c74a4f3fc852436c3cb5a5f41ded5c80a4e2b1dcc43a0d3f47d9",
        "readme_sha256": "d8d7a46d41a1a37fe4f0a5f637bf55c649310185329127d8a2204632e480be17",
        "license": "mit",
    },
    "epicardioids": {
        "repository": "lauradmartens/epicardioids-scooby",
        "revision": "cefc3a626f5080baaf83878754ea924df6bd7e8d",
        "size": 772988628,
        "sha256": "f945e1937a46ddd5cb31249ef71e721e7a7b1381aa207fb2de449f23ac75a9e8",
        "config_sha256": "f74372c396e1dd3d7d139eae6fed6fbfc6742130061e7e0cdaa11d93a00e700e",
        "readme_sha256": "d8d7a46d41a1a37fe4f0a5f637bf55c649310185329127d8a2204632e480be17",
        "license": "mit",
    },
    "neurips": {
        "repository": "johahi/neurips-scooby",
        "revision": "1591d15696c65f14053e4495e52c74812a86fd8d",
        "size": 772970196,
        "sha256": "8f2e2de3e86016378116c90accf94572e75963d5445beb6cabd585c8372170dd",
        "config_sha256": "26862280d15319c004a062e1823699d9e1150e753c1afcc6c40ffe9676e68eba",
        "license": None,
    },
}

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

metadata_receipts = {}
for model, contract in expected.items():
    api = json.loads((root / "metadata" / f"{model}.api.json").read_text())
    if api.get("sha") != contract["revision"]:
        raise SystemExit(f"{model} immutable API revision differs")
    siblings = {row["rfilename"]: row for row in api.get("siblings", [])}
    weight = siblings.get("model.safetensors", {})
    lfs = weight.get("lfs") or {}
    if int(weight.get("size", -1)) != contract["size"] or lfs.get("sha256") != contract["sha256"]:
        raise SystemExit(f"{model} API weight identity differs")
    card_license = (api.get("cardData") or {}).get("license")
    if card_license != contract["license"]:
        raise SystemExit(f"{model} immutable card license differs")
    config = root / "metadata" / f"{model}.config.json"
    if digest(config) != contract["config_sha256"]:
        raise SystemExit(f"{model} config bytes differ")
    pointer = (root / "metadata" / f"{model}.model.safetensors.lfs_pointer").read_text()
    expected_pointer = (
        "version https://git-lfs.github.com/spec/v1\n"
        f"oid sha256:{contract['sha256']}\n"
        f"size {contract['size']}\n"
    )
    if pointer != expected_pointer:
        raise SystemExit(f"{model} LFS pointer differs")
    if "readme_sha256" in contract and digest(root / "metadata" / f"{model}.README.md") != contract["readme_sha256"]:
        raise SystemExit(f"{model} README bytes differ")
    metadata_receipts[model] = {
        "revision": contract["revision"],
        "size_bytes": contract["size"],
        "declared_sha256": contract["sha256"],
        "card_license": card_license if card_license is not None else "UNDECLARED",
    }

headers = {
    model: json.loads((root / "headers" / f"{model}.header.json").read_text())
    for model in ("onek1k", "epicardioids")
}
if any(
    header["status"] != "pass"
    or header["checkpoint_payload_downloaded"]
    or header["checkpoint_deserialized"]
    or header["model_forward_executed"]
    or header["declared_file_sha256_verified_from_full_payload"]
    for header in headers.values()
):
    raise SystemExit("bounded safetensors header disposition differs")

training_members = {row["path"]: row for row in training["members"]}
submission_members = {row["path"]: row for row in submission["members"]}
required_context = {
    "onek1k_embedding": ("training_data/onek1k_training_data/embedding_no_val_genes.pq", 61200253),
    "epicardioids_embedding": ("training_data/epicardioids_training_data/embedding.pq", 7073807),
    "neurips_embedding": ("training_data/scooby_training_data/embedding_no_val_genes_new.pq", 4906504),
    "neurips_poisson_multivi": ("training_data/scooby_training_data/neurips_poissonmultivi_no_normo/model.pt", 249786484),
}
for label, (path, size) in required_context.items():
    if training_members.get(path, {}).get("uncompressed_size") != size:
        raise SystemExit(f"released context inventory differs: {label}")
for path in (
    "submission_data/model_weights/multiome_model/pytorch_model.bin",
    "submission_data/model_weights/rnaonly_model/pytorch_model.bin",
):
    if path not in submission_members:
        raise SystemExit(f"historical checkpoint inventory differs: {path}")

if (
    review["checkpoint_forward_allowed"]
    or review["rna_conditioned_atac_models"]
    or review["sealed_inference_allowed"]
):
    raise SystemExit("cross-review execution disposition differs")
receipt = {
    "schema_version": "masld-bench-scooby-bounded-header-admission-v1",
    "status": "pass",
    "metadata": metadata_receipts,
    "headers": {
        model: {
            "header_size_bytes": value["header_size_bytes"],
            "tensor_count": value["tensor_count"],
            "parameter_count": value["parameter_count"],
            "registry_reported_parameter_count": value["registry_reported_parameter_count"],
            "registry_to_header_element_difference": value["registry_to_header_element_difference"],
            "http_bytes_downloaded": value["http_bytes_downloaded"],
        }
        for model, value in headers.items()
    },
    "released_context_findings": {
        "scooby_onek1k": "exact_embedding_released_but_no_inductive_scPoli_encoder_or_liver_context",
        "scooby_epicardioids": "exact_embedding_released_but_no_inductive_scGLUE_encoder_and_context_used_observed_unpaired_ATAC",
        "scooby_neurips": "exact_embedding_and_PoissonMultiVI_artifact_released_but_registered_HF_checkpoint_terms_undeclared_and_context_used_observed_ATAC",
    },
    "archive_checkpoint_relationship": "historical_pickle_weights_are_noncanonical_until_byte_identity_source_runtime_and_safe_loading_are_separately_bound",
    "checkpoint_payload_downloaded": False,
    "checkpoint_deserialized": False,
    "model_forward_executed": False,
    "project_data_read": False,
    "outcomes_read": False,
    "rna_conditioned_atac_models": [],
    "checkpoint_forward_allowed": [],
    "terminal_disposition": "all_three_registered_scooby_checkpoints_have_terminal_admission_dispositions; no_forward_permitted",
}
(root / "scooby_bounded_header_admission_receipt.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
PY

sha256sum \
  "${ROOT}/scripts/scooby_remote_safetensors_header.py" \
  "${ROOT}/slurm/scooby_acquire_bounded_headers.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads(
    (root / "scooby_bounded_header_admission_receipt.json").read_text()
)
freeze_tree(
    root,
    {
        "artifact_class": "scooby_bounded_header_admission",
        "models": ["scooby_onek1k", "scooby_epicardioids", "scooby_neurips"],
        "checkpoint_payload_downloaded": False,
        "checkpoint_deserialized": False,
        "model_forward_executed": False,
        "project_data_read": False,
        "outcomes_read": False,
        "rna_conditioned_atac_models": [],
        "checkpoint_forward_allowed": [],
        "terminal_disposition": receipt["terminal_disposition"],
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scooby bounded-header admission: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
