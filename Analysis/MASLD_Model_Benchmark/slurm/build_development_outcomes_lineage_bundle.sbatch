#!/usr/bin/env bash
#SBATCH --job-name=seq-outcome-bundle
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${LINEAGE:?LINEAGE is required}"
case "${LINEAGE}" in cholangiocyte|fibroblast|hepatocyte|macrophage|t_cell) ;;
  *) printf 'Invalid primary lineage: %s\n' "${LINEAGE}" >&2; exit 1 ;;
esac

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
DONOR_BIGWIGS="${ROOT}/executions/gse296875-deduplicated-tn5-bigwigs-21063992"
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
BUNDLE_TARGET="${ROOT}/executions/sequence-development-outcomes-${LINEAGE}-bundle-${SLURM_JOB_ID}"
BUNDLE_STAGE="${BUNDLE_TARGET}.staging"
BUNDLE_FAILED="${BUNDLE_TARGET}.failed"
preserve_bundle() { if [[ -d "${BUNDLE_STAGE}" && ! -e "${BUNDLE_FAILED}" ]]; then mv "${BUNDLE_STAGE}" "${BUNDLE_FAILED}"; fi; }
trap preserve_bundle EXIT
if [[ -e "${BUNDLE_TARGET}" || -e "${BUNDLE_STAGE}" || -e "${BUNDLE_FAILED}" ]]; then
  printf 'Refusing to overwrite outcome bundle: %s\n' "${BUNDLE_TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${BUNDLE_STAGE}/validation" "${BUNDLE_STAGE}/receipts"
for pair in \
  "b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 ${BASE}/ARTIFACTS.json" \
  "96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa ${IMAGE}" \
  "61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 ${CORE}/ARTIFACTS.json" \
  "01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f ${DONOR_BIGWIGS}/ARTIFACTS.json" \
  "00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f ${SPLIT}/ARTIFACTS.json" \
  "4675b979b9b43aa9b75e6f2c6c0a2d2725292d4ab840e8c5e0665a2639f425a4 ${ROOT}/scripts/build_chrombpnet_development_outcomes.py" \
  "678cf4fadb8e8226c563ff91e594438417c7f5f4aad534a2be6b224bcb259f32 ${ROOT}/tests/unit/test_build_chrombpnet_development_outcomes.py"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SPLIT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]).resolve(strict=True))
PY
module purge
module load singularity/4.2.2
cd "${ROOT}"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m unittest -v tests.unit.test_build_chrombpnet_development_outcomes \
  > "${BUNDLE_STAGE}/validation/tests.stdout.txt" \
  2> "${BUNDLE_STAGE}/validation/tests.stderr.txt"

failures=0
for fold in 0 1 2 3 4; do
  split_id="donor${fold}_genomic${fold}"
  receipt="${BUNDLE_STAGE}/receipts/${split_id}.json"
  if [[ "${LINEAGE}" == hepatocyte && "${fold}" == 0 ]]; then
    existing="${ROOT}/executions/chrombpnet-hepatocyte-donor0-genomic0-outcomes-21064141"
    existing_sha=c890168b9bd8d11d233659baee4d6674c9cf470be8d8fea23ee653de865762c3
    printf '%s  %s\n' "${existing_sha}" "${existing}/ARTIFACTS.json" | sha256sum --check --strict
    module purge
    module load python/3.11.5-GCCcore-13.2.0
    PYTHONPATH="${ROOT}/src" python - "${existing}" "${receipt}" "${split_id}" "${existing_sha}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree
root, receipt, split_id, manifest_sha = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4]
verify_frozen_tree(root)
receipt.write_text(json.dumps({
    "schema_version": "masld-bench-development-outcome-bundle-task-v1",
    "status": "reused_frozen_artifact", "split_id": split_id,
    "artifact_root": str(root.resolve()), "artifacts_sha256": manifest_sha,
}, sort_keys=True, indent=2) + "\n")
PY
    module purge
    module load singularity/4.2.2
    continue
  fi
  target="${ROOT}/executions/chrombpnet-${LINEAGE}-${split_id}-outcomes-${SLURM_JOB_ID}"
  stage="${target}.staging"
  failed="${target}.failed"
  if [[ -e "${target}" || -e "${stage}" || -e "${failed}" ]]; then
    printf 'Outcome task target already exists: %s\n' "${target}" >&2
    failures=$((failures + 1))
    continue
  fi
  install -d -m 0750 "${stage}/validation"
  module purge
  module load singularity/4.2.2
  if singularity exec --cleanenv --containall --net --network none \
      --bind "${ROOT}:${ROOT}:rw" \
      --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
      --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
      "${IMAGE}" python "${ROOT}/scripts/build_chrombpnet_development_outcomes.py" \
        --donor-bigwigs "${DONOR_BIGWIGS}" \
        --donor-bigwigs-artifacts-sha256 01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f \
        --split-contract "${SPLIT}" \
        --split-artifacts-sha256 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
        --split-id "${split_id}" --lineage "${LINEAGE}" \
        --output "${stage}/outcomes" --workers 8 \
        > "${stage}/validation/build.stdout.txt" \
        2> "${stage}/validation/build.stderr.txt"; then
    singularity exec --cleanenv --containall --net --network none \
      --bind "${ROOT}:${ROOT}:rw" --env PYTHONNOUSERSITE=1 \
      "${IMAGE}" python - "${stage}" "${DONOR_BIGWIGS}" "${split_id}" "${LINEAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys
import h5py
import numpy as np

root, bigwigs, split_id, lineage = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4]
rows = list(csv.DictReader((bigwigs / "bigwig_manifest.tsv").open(), delimiter="\t"))
folds = {fold: 0 for fold in range(5)}
for row in rows:
    if row["lineage_id"] == lineage:
        folds[int(row["outer_fold"])] += 1
prefix = split_id.split("_")[0]
if not prefix.startswith("donor"):
    raise SystemExit("split prefix differs")
outer = int(prefix[len("donor"):])
expected = {"valid": folds[(outer + 1) % 5], "test": folds[outer]}
summary = json.loads((root / "outcomes/summary.json").read_text())
if summary["split_id"] != split_id or summary["lineage_id"] != lineage:
    raise SystemExit("outcome identity differs")
if summary["model_inference_input_eligible"] or summary["benchmark_metrics_calculated"]:
    raise SystemExit("outcome firewall differs")
with h5py.File(root / "outcomes/outcomes.h5", "r") as handle:
    for role, donors in expected.items():
        if handle[role]["counts"].shape != (donors, 16_000):
            raise SystemExit(f"{role} outcome shape differs")
        if np.any(np.asarray(handle[role]["counts"][:, :100]) < 0):
            raise SystemExit(f"{role} outcome is negative")
(root / "validation/probe.json").write_text(json.dumps({
    "schema_version": "masld-bench-development-outcome-matrix-probe-v1",
    "status": "pass", "split_id": split_id, "lineage_id": lineage,
    "role_donors": expected, "windows_per_role": 16000,
    "observed_atac_outcomes": True, "model_inference_input_eligible": False,
    "benchmark_metrics_calculated": False,
}, sort_keys=True, indent=2) + "\n")
PY
    sha256sum \
      "${ROOT}/scripts/build_chrombpnet_development_outcomes.py" \
      "${ROOT}/tests/unit/test_build_chrombpnet_development_outcomes.py" \
      "${ROOT}/slurm/build_development_outcomes_lineage_bundle.sbatch" \
      > "${stage}/source.sha256"
    module purge
    module load python/3.11.5-GCCcore-13.2.0
    PYTHONPATH="${ROOT}/src" python - "${stage}" "${split_id}" "${LINEAGE}" "${receipt}" <<'PY'
from pathlib import Path
from hashlib import sha256
import json
import sys
from masld_bench.artifacts import freeze_tree

root, split_id, lineage, receipt = Path(sys.argv[1]).resolve(strict=True), sys.argv[2], sys.argv[3], Path(sys.argv[4])
freeze_tree(root, {
    "artifact_class": "chrombpnet_development_outcomes",
    "dataset_id": "gse296875", "split_id": split_id, "lineage_id": lineage,
    "biological_unit": "donor", "observed_atac_outcomes": True,
    "model_inference_input_eligible": False, "benchmark_metrics_calculated": False,
    "selected_bigwig_members_verified": True, "status": "passed",
})
manifest_sha = sha256((root / "ARTIFACTS.json").read_bytes()).hexdigest()
receipt.write_text(json.dumps({
    "schema_version": "masld-bench-development-outcome-bundle-task-v1",
    "status": "passed", "split_id": split_id,
    "artifact_root": str(root).removesuffix(".staging"),
    "artifacts_sha256": manifest_sha,
}, sort_keys=True, indent=2) + "\n")
PY
    mv "${stage}" "${target}"
  else
    mv "${stage}" "${failed}"
    printf '{"schema_version":"masld-bench-development-outcome-bundle-task-v1","status":"failed","split_id":"%s"}\n' \
      "${split_id}" > "${receipt}"
    failures=$((failures + 1))
  fi
done

module purge
module load singularity/4.2.2
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${BUNDLE_STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/build_chrombpnet_development_outcomes.py" \
  "${ROOT}/tests/unit/test_build_chrombpnet_development_outcomes.py" \
  "${ROOT}/slurm/build_development_outcomes_lineage_bundle.sbatch" \
  > "${BUNDLE_STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BUNDLE_STAGE}" "${LINEAGE}" "${failures}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root, lineage, failures = Path(sys.argv[1]).resolve(strict=True), sys.argv[2], int(sys.argv[3])
receipts = [json.loads(path.read_text()) for path in sorted((root / "receipts").glob("*.json"))]
if len(receipts) != 5:
    raise SystemExit("bundle receipt count differs")
freeze_tree(root, {
    "artifact_class": "sequence_development_outcomes_lineage_bundle",
    "dataset_id": "gse296875", "lineage_id": lineage,
    "logical_tasks": 5, "failures": failures,
    "observed_atac_outcomes": True, "model_inference_input_eligible": False,
    "benchmark_metrics_calculated": False,
    "status": "passed" if failures == 0 else "partial_failure",
})
PY
trap - EXIT
mv "${BUNDLE_STAGE}" "${BUNDLE_TARGET}"
sha256sum "${BUNDLE_TARGET}/ARTIFACTS.json"
if [[ "${failures}" -ne 0 ]]; then exit 2; fi
