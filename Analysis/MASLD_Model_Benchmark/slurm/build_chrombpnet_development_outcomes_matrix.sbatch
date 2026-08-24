#!/usr/bin/env bash
#SBATCH --job-name=masld-dev-atac-outcomes
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${OUTER_FOLD:?OUTER_FOLD is required}"
: "${LINEAGE:?LINEAGE is required}"
: "${GENOMIC_FOLD:=${OUTER_FOLD}}"

case "${OUTER_FOLD}" in 0|1|2|3|4) ;; *) printf 'Invalid OUTER_FOLD\n' >&2; exit 1 ;; esac
case "${GENOMIC_FOLD}" in 0|1|2|3|4) ;; *) printf 'Invalid GENOMIC_FOLD\n' >&2; exit 1 ;; esac
case "${LINEAGE}" in b_cell|cholangiocyte|endothelial_cell|fibroblast|hepatocyte|macrophage|t_cell) ;;
  *) printf 'Invalid LINEAGE\n' >&2; exit 1 ;;
esac

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
DONOR_BIGWIGS="${ROOT}/executions/gse296875-deduplicated-tn5-bigwigs-21063992"
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
SPLIT_ID="donor${OUTER_FOLD}_genomic${GENOMIC_FOLD}"
TARGET="${ROOT}/executions/chrombpnet-${LINEAGE}-${SPLIT_ID}-outcomes-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite development outcome artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for pair in \
  "b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 ${BASE}/ARTIFACTS.json" \
  "96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa ${IMAGE}" \
  "61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 ${CORE}/ARTIFACTS.json" \
  "01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f ${DONOR_BIGWIGS}/ARTIFACTS.json" \
  "00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f ${SPLIT}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${DONOR_BIGWIGS}" "${SPLIT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY

module purge
module load singularity/4.2.2
cd "${ROOT}"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python "${ROOT}/tests/unit/test_build_chrombpnet_development_outcomes.py"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python "${ROOT}/scripts/build_chrombpnet_development_outcomes.py" \
    --donor-bigwigs "${DONOR_BIGWIGS}" \
    --donor-bigwigs-artifacts-sha256 01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f \
    --split-contract "${SPLIT}" \
    --split-artifacts-sha256 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
    --split-id "${SPLIT_ID}" --lineage "${LINEAGE}" \
    --output "${STAGE}/outcomes" --workers 8 \
    > "${STAGE}/validation/build.stdout.txt" \
    2> "${STAGE}/validation/build.stderr.txt"

singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python - "${STAGE}" "${DONOR_BIGWIGS}" "${SPLIT_ID}" "${LINEAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys
import h5py
import numpy as np

root, bigwigs, split_id, lineage = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4]
rows = list(csv.DictReader((bigwigs / "bigwig_manifest.tsv").open(), delimiter="\t"))
folds = {int(row["outer_fold"]): 0 for row in rows if row["lineage_id"] == lineage}
for row in rows:
    if row["lineage_id"] == lineage:
        folds[int(row["outer_fold"])] += 1
prefix = split_id.split("_")[0]
if not prefix.startswith("donor"):
    raise SystemExit("split prefix differs")
outer = int(prefix[len("donor"):])
expected = {"valid": folds[(outer + 1) % 5], "test": folds[outer]}
summary = json.loads((root / "outcomes" / "summary.json").read_text())
if summary["split_id"] != split_id or summary["lineage_id"] != lineage:
    raise SystemExit("outcome identity differs")
if summary["model_inference_input_eligible"] or summary["benchmark_metrics_calculated"]:
    raise SystemExit("outcome firewall differs")
with h5py.File(root / "outcomes" / "outcomes.h5", "r") as handle:
    for role, donors in expected.items():
        values = handle[role]["counts"]
        if values.shape != (donors, 16000) or np.any(np.asarray(values[:, :100]) < 0):
            raise SystemExit(f"{role} outcome matrix differs")
(root / "validation" / "probe.json").write_text(json.dumps({
    "schema_version": "masld-bench-development-outcome-matrix-probe-v1",
    "status": "pass", "split_id": split_id, "lineage_id": lineage,
    "role_donors": expected, "windows_per_role": 16000,
    "observed_atac_outcomes": True, "model_inference_input_eligible": False,
    "benchmark_metrics_calculated": False,
}, sort_keys=True, indent=2) + "\n")
PY
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum "${ROOT}/scripts/build_chrombpnet_development_outcomes.py" \
  "${ROOT}/tests/unit/test_build_chrombpnet_development_outcomes.py" \
  "${ROOT}/slurm/build_chrombpnet_development_outcomes_matrix.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${SPLIT_ID}" "${LINEAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
root, split_id, lineage = Path(sys.argv[1]).resolve(strict=True), sys.argv[2], sys.argv[3]
freeze_tree(root, {
    "artifact_class": "chrombpnet_development_outcomes",
    "dataset_id": "gse296875", "split_id": split_id, "lineage_id": lineage,
    "biological_unit": "donor", "observed_atac_outcomes": True,
    "model_inference_input_eligible": False, "benchmark_metrics_calculated": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Development outcome artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
