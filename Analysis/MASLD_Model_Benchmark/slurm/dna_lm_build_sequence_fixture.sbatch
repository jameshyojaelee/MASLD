#!/usr/bin/env bash
#SBATCH --job-name=masld-dna-lm-fixture
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

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG="${ROOT}/config/dna_lm_fixture.json"
SOURCE="${ROOT}/executions/dna-lm-admission-sources-21064540"
UPSTREAM="${ROOT}/executions/alphagenome-sei-sequence-fixture-21064515"
TARGET="${ROOT}/executions/dna-lm-sequence-fixture-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite DNA-LM sequence fixture: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' 4e1703b03ff28924904053e842b89079d510a80d2a9ad082de5369adaeea75c8 \
  "${SOURCE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' e87468da0f64da1c022751834b624471221aec3ba03fcd17d6dbf756e3b8472e \
  "${UPSTREAM}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SOURCE}" "${UPSTREAM}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

source, upstream = (Path(value).resolve(strict=True) for value in sys.argv[1:])
verify_frozen_tree(source)
verify_frozen_tree(upstream)
admission = json.loads((source / "review" / "source_admission_receipt.json").read_text())
firewall = json.loads((upstream / "fixture" / "outcome_firewall.json").read_text())
if (
    admission["status"] != "pass"
    or admission["checkpoint_bytes_downloaded"]
    or admission["model_forward_executed"]
    or firewall["observed_atac_loaded"]
    or firewall["rna_loaded"]
    or firewall["variant_or_regulatory_outcomes_loaded"]
    or firewall["sealed_labels_loaded"]
    or firewall["fixture_genomic_roles"] != ["train"]
):
    raise SystemExit("DNA-LM source review or outcome firewall differs")
print("DNA-LM source review and upstream outcome firewall passed")
PY

python "${ROOT}/scripts/dna_lm_build_fixture.py" \
  --config "${CONFIG}" \
  --upstream "${UPSTREAM}" \
  --output "${STAGE}/fixture" \
  > "${STAGE}/validation/fixture.stdout.txt" \
  2> "${STAGE}/validation/fixture.stderr.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/dna_lm_build_fixture.py" \
  "${ROOT}/slurm/dna_lm_build_sequence_fixture.sbatch" \
  "${SOURCE}/ARTIFACTS.json" \
  "${UPSTREAM}/ARTIFACTS.json" \
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
expected_contexts = {
    "common_6000": (6000, 5),
    "hyenadna_native_1000000": (1000000, 1),
    "caduceus_native_131072": (131072, 1),
    "evo_native_1048576": (1048576, 1),
}
if (
    summary["status"] != "pass"
    or summary["anchors"] != 3
    or summary["sequence_manifest_rows"] != 24
    or summary["common_pool_length_bp"] != 1536
    or summary["embeddings_biologically_interpretable_without_head"]
    or any(firewall[key] for key in (
        "checkpoint_bytes_loaded",
        "observed_atac_loaded",
        "rna_loaded",
        "variant_or_regulatory_outcomes_loaded",
        "sealed_labels_loaded",
        "supervised_eqtl_or_ieqtl_head_fit",
        "model_forward_executed",
    ))
    or len(rows) != 24
):
    raise SystemExit("DNA-LM fixture summary or firewall differs")
for context_id, (length, model_count) in expected_contexts.items():
    context_rows = [row for row in rows if row["context_id"] == context_id]
    if (
        len(context_rows) != 3 * model_count
        or {int(row["input_length_bp"]) for row in context_rows} != {length}
        or {int(row["forward_variant_index0"]) for row in context_rows} != {length // 2}
        or {int(row["reverse_complement_variant_index0"]) for row in context_rows}
        != {length // 2 - 1}
    ):
        raise SystemExit("DNA-LM context geometry differs")
    with gzip.open(fixture / f"{context_id}.alleles.fa.gz", "rt", encoding="ascii") as handle:
        if sum(line.startswith(">") for line in handle) != 12:
            raise SystemExit("DNA-LM context FASTA census differs")
freeze_tree(
    root,
    {
        "artifact_class": "dna_lm_sequence_fixture",
        "models": [
            "dnabert2",
            "nucleotide_transformer",
            "hyenadna",
            "caduceus",
            "evo",
        ],
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
printf 'DNA-LM sequence fixture artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
