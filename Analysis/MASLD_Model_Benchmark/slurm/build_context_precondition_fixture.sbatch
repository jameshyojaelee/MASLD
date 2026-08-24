#!/usr/bin/env bash
#SBATCH --job-name=masld-context-fixture
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/context-conditioning-precondition-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
MODULE_SHA256=55c14bcf03cbd7c0edb976f6ad73c9d892c3327f7fe91622f11ba76edc00da70
CAMPAIGN_SHA256=3f4db25525fc8dcba655d214cf3841747c877e96321d98998246811138be785c

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite context precondition fixture: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${MODULE_SHA256}" \
  "${ROOT}/src/masld_bench/context_preconditions.py" | sha256sum --check --strict
printf '%s  %s\n' "${CAMPAIGN_SHA256}" \
  "${ROOT}/config/campaigns/v1_conditional_model.toml" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="src:." python - "${STAGE}/fixture.json" <<'PY'
from pathlib import Path
import json
import sys
import tomllib

from masld_bench.context_preconditions import (
    ContextMapperPrecondition,
    ContextRow,
    ContextVector,
    apply_film_lora_precondition,
    require_precondition_only_campaign,
)

with Path("config/campaigns/v1_conditional_model.toml").open("rb") as handle:
    campaign = tomllib.load(handle)
require_precondition_only_campaign(campaign)
observed = ContextVector(
    values=(3.0, 1.0, 2.0),
    observed_mask=(True, True, True),
    states=("observed", "observed", "observed"),
)
train = ContextRow(
    row_id="train-donor:train-block",
    dataset_id="GSE296875",
    donor_id="train-donor",
    study_id="GSE296875",
    split_role="train",
    genomic_block="train-block",
    pairing_level="same_nucleus",
    rna=observed,
)
query = ContextRow(
    row_id="held-donor:held-block",
    dataset_id="GSE296875",
    donor_id="held-donor",
    study_id="GSE296875",
    split_role="test",
    genomic_block="held-block",
    pairing_level="same_nucleus",
    rna=observed,
)
mapper = ContextMapperPrecondition.fit(
    [train],
    ordered_gene_ids=("g1", "g2", "g3"),
    reference_values=(10.0, 20.0, 30.0),
    mapper_arm="released_rank_quantile_from_counts",
    held_donor_ids=("held-donor",),
    held_genomic_blocks=("held-block",),
    receptive_field_bp=524_288,
    boundary_buffer_bp=524_288,
)
mapped = mapper.transform_query(query)
features = ((1.0, 2.0, 3.0),)
neutral = apply_film_lora_precondition(
    features,
    gamma_delta=(0.0, 0.0, 0.0),
    beta=(0.0, 0.0, 0.0),
    lora_a=((1.0, 0.0, 0.0),),
    lora_b=((0.0,), (0.0,), (0.0,)),
    alpha=1.0,
    context_observed=True,
    gate=1.0,
)
if mapped != (30.0, 10.0, 20.0) or neutral != features:
    raise SystemExit("context precondition fixture differs")
payload = {
    "schema_version": "masld-bench-context-conditioning-precondition-v1",
    "status": "pass",
    "mapper_fit_roles": ["train"],
    "training_donors": sorted(mapper.training_donor_ids),
    "query_donors": [query.donor_id],
    "donor_overlap": False,
    "training_genomic_blocks": sorted(mapper.training_genomic_blocks),
    "query_genomic_blocks": [query.genomic_block],
    "genomic_block_overlap": False,
    "boundary_buffer_bp": mapper.boundary_buffer_bp,
    "receptive_field_bp": mapper.receptive_field_bp,
    "observed_atac_input_exposed": False,
    "sealed_outcomes_exposed": False,
    "missing_context_numeric_zero_allowed": False,
    "neutral_film_identity": True,
    "zero_initialized_lora_identity": True,
    "architecture_built": False,
    "model_training_started": False,
    "conditional_campaign_status": campaign["status"],
    "conditional_submit_enabled": campaign["submit_enabled"],
    "conditional_model_authorized": False,
}
Path(sys.argv[1]).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(json.dumps(payload, sort_keys=True))
PY

python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ROOT}/src/masld_bench/context_preconditions.py" \
  "${ROOT}/tests/unit/test_context_preconditions.py" \
  "${ROOT}/config/campaigns/v1_conditional_model.toml" \
  "${ROOT}/slurm/build_context_precondition_fixture.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

freeze_tree(
    Path(sys.argv[1]),
    metadata={
        "artifact_class": "context_conditioning_precondition_fixture",
        "status": "pass",
        "conditional_model_authorized": False,
        "architecture_built": False,
        "model_training_started": False,
        "sealed_outcomes_exposed": False,
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Context conditioning precondition fixture: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
