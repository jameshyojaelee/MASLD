#!/usr/bin/env bash
#SBATCH --job-name=bench-dna-lm-native-preflight
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export TOKENIZERS_PARALLELISM=false
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
CONFIG="${ROOT}/config/gse281364_dna_lm_native_contract.json"
RUNTIME="${ROOT}/executions/environments/dna-lm-transformers-4.38.1-21064683"
TARGET="${ROOT}/executions/gse281364-dna-lm-native-contract-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite DNA-LM native preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
declare -a bindings=(
  "74a3a5c36bd1416daf3fe4b63fe317cfa7e2d008b96e2503d4ca43770b4357de executions/gse281364-dna-lm-common-fixture-21069076"
  "876c072889324661df0500e2f1eb8d09a8f1e3054deddb153d8296b069845a63 executions/environments/dna-lm-transformers-4.38.1-21064683"
  "5b64097862b0db5ad3bdb1d0de1429e932bc06f315e7106209696c702878e3fa executions/dnabert2-runtime-bundle-21064834"
  "25a1d8571cb6e7188f398ac820402ad9f4aea7571a99e5b794cfb92dea57b536 executions/environments/dnabert2-einops-0.6.1-21064904"
  "6c881e7f390d301a818c305bd67e6e305077de87daafd642340814dedbe5f02f executions/dna-lm-weights-only-dnabert2-21064812"
  "a0324f4b774c0af4d8bfbdc9bb4ce4514d62e704e4f0ae34bf8cd8d2c67f9545 executions/dnabert2-common-runtime-probe-21065134"
  "4e1703b03ff28924904053e842b89079d510a80d2a9ad082de5369adaeea75c8 executions/dna-lm-admission-sources-21064540"
  "9ec36ffea9bf2c94de1d1d6be6f7f6bc42822a6003f5cf9cd7ee342da2d6f395 executions/dna-lm-weights-only-nucleotide_transformer-21064811"
  "f0c65ebdb5880b9392dfd131e43b3f320df88cd8eccb7accf1f7cc9dd358e325 executions/nucleotide-transformer-common-runtime-probe-21065102"
  "70906faae158909a8972f674f680943d8a51467032ae3266a6a12f2152beb36c executions/dna-lm-custom-code-21064626"
  "6d21df999835624905e08573038dc72e893a5b25a8ca8f4ec1acea8f61d91c77 executions/dna-lm-open-safetensors-21064556"
  "9deb82d5faae1636cb6129ce6d30b70b8264c6cc30313caa8fb7506293aa12b6 executions/hyenadna-common-runtime-probe-21064769"
)
for binding in "${bindings[@]}"; do
  read -r expected relative <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}/ARTIFACTS.json" \
    | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${ROOT}" "${bindings[@]}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
for binding in sys.argv[2:]:
    _, relative = binding.split(" ", 1)
    verify_frozen_tree(root / relative)
print("all frozen runtime, checkpoint, license, probe, and fixture trees verified")
PY
module purge
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH="${ROOT}:${RUNTIME}/site-packages"
"${BASE}" -m unittest -v \
  tests.unit.test_gse281364_dna_lm_native_contract \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${BASE}" "${ROOT}/scripts/preflight_gse281364_dna_lm_native_contract.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/contract" \
  > "${STAGE}/validation/preflight.stdout.json" \
  2> "${STAGE}/validation/preflight.stderr.txt"
"${BASE}" -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/gse281364_dna_lm_native_contract.py" \
  "${ROOT}/scripts/preflight_gse281364_dna_lm_native_contract.py" \
  "${ROOT}/scripts/extract_gse281364_dna_lm_embeddings.py" \
  "${ROOT}/scripts/project_gse281364_dna_lm_embeddings.py" \
  "${ROOT}/tests/unit/test_gse281364_dna_lm_native_contract.py" \
  "${ROOT}/slurm/preflight_gse281364_dna_lm_native_contract.sbatch" \
  > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
contract = json.loads((root / "contract/contract.json").read_text())
with (root / "contract/tokenization_manifest.tsv").open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if (
    contract["status"] != "pass_outcome_blind_native_preflight"
    or contract["elements"] != 1033
    or contract["outer_locus_sequence_groups"] != 1033
    or contract["outer_folds"] != 5
    or contract["allele_sequences"] != 4132
    or contract["models"] != ["dnabert2", "nucleotide_transformer", "hyenadna"]
    or contract["nucleotide_transformer"]["tokens_per_phase_including_cls"] != 682
    or contract["hyenadna"]["tokens_including_terminal_sep"] != 4097
    or contract["checkpoint_bytes_loaded"]
    or contract["model_forward_executed"]
    or contract["reporter_counts_read"]
    or contract["reporter_outcomes_read"]
    or contract["sealed_labels_read"]
    or contract["downstream_head_fit"]
    or len(rows) != 1033
    or {row["outer_fold"] for row in rows} != {"0", "1", "2", "3", "4"}
):
    raise SystemExit("GSE281364 DNA-LM native preflight differs")
freeze_tree(root, {
    "artifact_class": "gse281364_dna_lm_native_token_coordinate_contract",
    "dataset_id": "gse281364",
    "models": ["dnabert2", "nucleotide_transformer", "hyenadna"],
    "elements": 1033,
    "outer_locus_sequence_groups": 1033,
    "outer_folds": 5,
    "checkpoint_bytes_loaded": False,
    "model_forward_executed": False,
    "reporter_outcomes_read": False,
    "downstream_head_fit": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 DNA-LM native contract: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
