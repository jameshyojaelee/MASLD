#!/usr/bin/env bash
#SBATCH --job-name=data-stage-901
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
LOCAL=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/PXD051911
FTP=ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2025/03/PXD051911
TARGET="${ROOT}/executions/pxd051911-activation-readiness-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"

cd "${ROOT}"
module purge
module load h5py/3.9.0-foss-2023a

mkdir -p "${STAGING}/source" "${STAGING}/audit" "${STAGING}/validation"

# The deposit's own checksum.txt is the provenance authority for every matrix.
wget -q -O "${STAGING}/source/checksum.txt" "${FTP}/checksum.txt"

# scWAT and oWAT are not held locally and must come from the deposit.
for f in scwat_protein_quant.txt owat_protein_quant.txt; do
  wget -q -O "${STAGING}/source/${f}" "${FTP}/${f}"
done

# Liver and metadata are already held locally; copy them in so the fixture is
# self-contained and every matrix is authenticated by the same code path.
for f in liver_protein_quant.txt meta_data.txt; do
  cp -n "${LOCAL}/${f}" "${STAGING}/source/${f}"
done

# Authenticate every staged matrix against the digests the deposit published.
python - "${STAGING}/source" <<'PY'
import hashlib, sys
from pathlib import Path

source = Path(sys.argv[1])
published = {}
for line in (source / "checksum.txt").read_text(encoding="utf-8").splitlines():
    if "\t" not in line or line.startswith("#"):
        continue
    path_text, digest = line.rsplit("\t", 1)
    published[path_text.replace("\\", "/").rsplit("/", 1)[-1]] = digest.strip()

for name in (
    "liver_protein_quant.txt",
    "meta_data.txt",
    "scwat_protein_quant.txt",
    "owat_protein_quant.txt",
):
    value = hashlib.sha1()
    with (source / name).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    actual = value.hexdigest()
    if actual != published[name]:
        raise SystemExit(f"FAIL {name}: {actual} != published {published[name]}")
    print(f"ok {name} {actual}")
PY

# Full contract suite, no test excluded. The suite carries pre-existing
# model-registry failures owned by other campaigns, so the check is that this
# campaign introduces no NEW failure and that its own tests ran and passed.
set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_build_pxd051911_activation_readiness \
  tests.unit.test_check_contract_suite_baseline \
  tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_build_pxd051911_activation_readiness \
  --require-module test_registry_files \
  | tee "${STAGING}/validation/contract_baseline.txt"

PYTHONPATH=${ROOT}:${ROOT}/src python scripts/build_pxd051911_activation_readiness.py \
  --source "${STAGING}/source" \
  --deposit-checksum "${STAGING}/source/checksum.txt" \
  --output "${STAGING}/audit/readiness" \
  > "${STAGING}/validation/build.stdout.txt" \
  2> "${STAGING}/validation/build.stderr.txt"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true

sha256sum \
  "${ROOT}/config/datasets/pxd051911.toml" \
  "${ROOT}/config/evaluation/family_native_tournament.toml" \
  "${ROOT}/scripts/build_pxd051911_activation_readiness.py" \
  "${ROOT}/tests/unit/test_build_pxd051911_activation_readiness.py" \
  "${ROOT}/scripts/check_contract_suite_baseline.py" \
  "${ROOT}/tests/unit/test_check_contract_suite_baseline.py" \
  "${ROOT}/slurm/stage_pxd051911_activation_readiness.sbatch" \
  > "${STAGING}/source.sha256"

PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree

digest = freeze_tree(
    Path("${STAGING}"),
    metadata={
        "artifact_class": "pxd051911_activation_readiness",
        "dataset_id": "pxd051911",
        "status": "activation_ready_development_only",
        "unit_of_replication": "participant",
        "sealed_outcomes_read": False,
        "models_fit": False,
        "models_scored": False,
        "rna_protein_pair_inferred": False,
        "champion_claim_allowed": False,
    },
)
print("ARTIFACTS.json sha256:", digest)
PY

mv "${STAGING}" "${TARGET}"
echo "FROZEN: ${TARGET}"
