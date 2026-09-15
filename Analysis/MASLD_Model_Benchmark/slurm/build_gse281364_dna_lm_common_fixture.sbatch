#!/usr/bin/env bash
#SBATCH --job-name=bench-mpra-dna-lm-fixture
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG="${ROOT}/config/gse281364_dna_lm_common_fixture.json"
SPLIT="${ROOT}/executions/gse281364-outcome-blind-splits-21068885"
SPLIT_SHA=26a20a18c7697c4e16eb645640428d3fd6290957efcef4df18219f9e16320de8
SEI="${ROOT}/executions/gse281364-sei-fixture-21068925"
SEI_SHA=620464c0ef7cf378eb2b44a1608cb6b9c273bd70b8770a43b2a95d21e725af4d
TARGET="${ROOT}/executions/gse281364-dna-lm-common-fixture-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite common DNA-LM fixture artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
printf '%s  %s\n' "${SPLIT_SHA}" "${SPLIT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${SEI_SHA}" "${SEI}/ARTIFACTS.json" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SPLIT}" "${SEI}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY
install -d -m 0750 "${STAGE}/validation"
PYTHONPATH="${ROOT}" python -m unittest -v \
  tests.unit.test_build_gse281364_dna_lm_common_fixture \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}" python "${ROOT}/scripts/build_gse281364_dna_lm_common_fixture.py" \
  --config "${CONFIG}" \
  --split-root "${SPLIT}" \
  --sei-root "${SEI}" \
  --output "${STAGE}/fixture" \
  > "${STAGE}/validation/build.stdout.json" \
  2> "${STAGE}/validation/build.stderr.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/build_gse281364_dna_lm_common_fixture.py" \
  "${ROOT}/tests/unit/test_build_gse281364_dna_lm_common_fixture.py" \
  "${ROOT}/slurm/build_gse281364_dna_lm_common_fixture.sbatch" \
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
with (fixture / "model_contracts.tsv").open(newline="", encoding="utf-8") as handle:
    models = list(csv.DictReader(handle, delimiter="\t"))
with gzip.open(fixture / "common_4096.alleles.fa.gz", "rt", encoding="ascii") as handle:
    fasta_records = sum(line.startswith(">") for line in handle)
if (
    summary["status"] != "pass_outcome_blind_common_fixture"
    or summary["elements"] != 1033
    or summary["outer_locus_sequence_groups"] != 1033
    or summary["outer_folds"] != 5
    or summary["sequences"] != 4132
    or summary["input_length_bp"] != 4096
    or summary["forward_variant_index0"] != 2048
    or summary["reverse_complement_variant_index0"] != 2047
    or summary["pool_length_bp"] != 1536
    or len(rows) != 1033
    or len(models) != 3
    or fasta_records != 4132
    or any(firewall[key] for key in (
        "reporter_counts_loaded",
        "reporter_outcomes_loaded",
        "sealed_labels_loaded",
        "checkpoint_bytes_loaded",
        "model_forward_executed",
        "downstream_head_fit",
    ))
):
    raise SystemExit("GSE281364 common DNA-LM fixture differs")
if {row["outer_locus_sequence_group_id"] for row in rows}.__len__() != 1033:
    raise SystemExit("locus-group assignment differs")
if {row["outer_fold"] for row in rows} != {"0", "1", "2", "3", "4"}:
    raise SystemExit("outer-fold coverage differs")
if {row["model_id"] for row in models} != {
    "dnabert2", "nucleotide_transformer", "hyenadna"
}:
    raise SystemExit("model contract differs")
freeze_tree(root, {
    "artifact_class": "gse281364_outcome_blind_dna_lm_common_fixture",
    "dataset_id": "gse281364",
    "models": ["dnabert2", "nucleotide_transformer", "hyenadna"],
    "elements": 1033,
    "outer_locus_sequence_groups": 1033,
    "outer_folds": 5,
    "sequences": 4132,
    "input_length_bp": 4096,
    "reporter_outcomes_read": False,
    "checkpoint_bytes_loaded": False,
    "model_forward_executed": False,
    "downstream_head_fit": False,
    "champion_eligible": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 common DNA-LM fixture: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
