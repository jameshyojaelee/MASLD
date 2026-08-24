#!/usr/bin/env bash
#SBATCH --job-name=masld-alpha-sei-fixture
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG="${ROOT}/config/alphagenome_sei_admission.json"
FASTA_ROOT="${ROOT}/executions/sequence-fasta-21063808"
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
TARGET="${ROOT}/executions/alphagenome-sei-sequence-fixture-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite AlphaGenome/Sei sequence fixture: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' 7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 \
  "${FASTA_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
  "${SPLIT}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${FASTA_ROOT}" "${SPLIT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("frozen reference and split verification passed")
PY
python "${ROOT}/scripts/alphagenome_sei_build_fixture.py" \
  --config "${CONFIG}" \
  --fasta "${FASTA_ROOT}/GRCh38.p14.sequence_model.fa" \
  --fai "${FASTA_ROOT}/GRCh38.p14.sequence_model.fa.fai" \
  --windows "${SPLIT}/ccre_evaluation_windows.tsv" \
  --output "${STAGE}/fixture" \
  > "${STAGE}/validation/fixture.stdout.txt" \
  2> "${STAGE}/validation/fixture.stderr.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/alphagenome_sei_build_fixture.py" \
  "${ROOT}/slurm/alphagenome_sei_build_sequence_fixture.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import gzip
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
fixture = root / "fixture"
summary = json.loads((fixture / "summary.json").read_text())
firewall = json.loads((fixture / "outcome_firewall.json").read_text())
with (fixture / "sequence_manifest.tsv").open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if (
    summary["status"] != "pass"
    or summary["anchors"] != 3
    or summary["sequence_rows"] != 6
    or summary["fasta_records_per_model"] != 12
    or summary["observed_outcomes_loaded"]
    or summary["checkpoint_bytes_loaded"]
    or summary["runtime_forward_executed"]
    or firewall["fixture_genomic_roles"] != ["train"]
    or firewall["observed_atac_loaded"]
    or firewall["rna_loaded"]
    or firewall["variant_or_regulatory_outcomes_loaded"]
    or firewall["sealed_labels_loaded"]
    or len(rows) != 6
):
    raise SystemExit("AlphaGenome/Sei sequence fixture differs")
for model_id, expected_length in (("alphagenome_all_folds", 1048576), ("sei", 4096)):
    model_rows = [row for row in rows if row["model_id"] == model_id]
    if (
        len(model_rows) != 3
        or {int(row["input_length"]) for row in model_rows} != {expected_length}
        or {int(row["variant_index0"]) for row in model_rows} != {expected_length // 2}
    ):
        raise SystemExit("AlphaGenome/Sei sequence geometry differs")
    with gzip.open(fixture / f"{model_id}.alleles.fa.gz", "rt", encoding="ascii") as handle:
        if sum(line.startswith(">") for line in handle) != 12:
            raise SystemExit("AlphaGenome/Sei FASTA census differs")
freeze_tree(
    root,
    {
        "artifact_class": "alphagenome_sei_sequence_fixture",
        "models": ["alphagenome_all_folds", "sei"],
        "reference_build": "GRCh38.p14",
        "eligible_genomic_folds": [2, 3, 4],
        "checkpoint_bytes_loaded": False,
        "observed_outcomes_loaded": False,
        "runtime_forward_executed": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'AlphaGenome/Sei sequence fixture: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
