#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-preflight-101
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
BIGWIGS="${ROOT}/executions/gse296875-deduplicated-tn5-bigwigs-21063992"
TARGET="${ROOT}/executions/corgi-regional-evaluator-runtime-preflight-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi evaluator preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for pair in \
  "b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 ${BASE}/ARTIFACTS.json" \
  "96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa ${IMAGE}" \
  "61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 ${CORE}/ARTIFACTS.json" \
  "01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f ${BIGWIGS}/ARTIFACTS.json" \
  "8008d5a8b8803b61890db4823d4b27be4c472f44881ca3788ce728f01be8f731 ${ROOT}/scripts/evaluate_corgi_regional_smoke.py" \
  "f2725c66fae08cb4e2002008441a36e30b8591abc6cdc11731fff1fa3893fa16 ${ROOT}/tests/unit/test_evaluate_corgi_regional_smoke.py" \
  "5e9c60694f4af7687a2a88603a664f83d441a70bbb87607624c5be2ba13bc289 ${ROOT}/slurm/evaluate_corgi_regional_smoke.sbatch"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BASE}" "${CORE}" "${BIGWIGS}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY
module purge
module load singularity/4.2.2
container=(
  singularity exec --cleanenv --containall --net --network none
  --bind "${ROOT}:${ROOT}:rw"
  --pwd "${ROOT}"
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages"
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1
  "${IMAGE}"
)
"${container[@]}" python -m unittest -v tests.unit.test_evaluate_corgi_regional_smoke \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${container[@]}" python scripts/evaluate_corgi_regional_smoke.py --help \
  > "${STAGE}/validation/cli-help.txt"
"${container[@]}" python - "${BIGWIGS}/bigwig_manifest.tsv" <<'PY' \
  > "${STAGE}/validation/pybigwig-probe.json"
import csv
import json
from pathlib import Path
import sys
import pyBigWig

with Path(sys.argv[1]).open(newline="", encoding="utf-8") as handle:
    row = next(csv.DictReader(handle, delimiter="\t"))
path = Path(row["path"])
if not path.is_absolute():
    path = Path(sys.argv[1]).parent / path
with pyBigWig.open(str(path)) as track:
    contigs = track.chroms()
    contig = next(iter(contigs))
    value = track.stats(contig, 0, min(1000, contigs[contig]), type="sum")[0]
print(json.dumps({"pyBigWig": pyBigWig.__version__, "contig": contig, "finite_or_null": value is None or value == value}, sort_keys=True))
PY
"${container[@]}" python -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/evaluate_corgi_regional_smoke.py" \
  "${ROOT}/tests/unit/test_evaluate_corgi_regional_smoke.py" \
  "${ROOT}/slurm/evaluate_corgi_regional_smoke.sbatch" \
  "${ROOT}/slurm/preflight_corgi_regional_evaluator_runtime.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "validation/pybigwig-probe.json").read_text())
if not probe.get("finite_or_null") or not probe.get("contig"):
    raise SystemExit("Corgi evaluator runtime probe differs")
freeze_tree(root, {
    "artifact_class": "corgi_regional_evaluator_runtime_preflight",
    "outcomes_read": False,
    "held_validation_atac_read": False,
    "direct_cli_executed": True,
    "pybigwig_file_read": True,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi evaluator runtime preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
