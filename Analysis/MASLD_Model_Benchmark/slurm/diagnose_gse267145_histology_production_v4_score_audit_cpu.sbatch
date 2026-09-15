#!/usr/bin/env bash
#SBATCH --job-name=model-check-080
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src:${ROOT}"

"${PYTHON}" - <<'PY'
from pathlib import Path
import json

from masld_bench.artifacts import verify_frozen_tree
from masld_bench.contracts import MissingState, PredictionBundle
from masld_bench.hashing import canonical_sha256, sha256_file

import scripts.audit_gse267145_histology_production_v4_scores as auditor

root = Path("executions/model-scoring-078-21092073.failed/scores")
production = Path("executions/model-training-069-production-v4")
participants, outer_folds, _ = auditor._load_outcomes(
    Path("executions/model-data-061-21079623/activation"),
    Path("executions/model-data-064-21079902/fixture/folds"),
)
_, lineage = auditor._verify_five_seed_ensembles(
    prediction_root=production / "aggregate/predictions",
    participants=participants,
    outer_folds=outer_folds,
)
_, index_rows = auditor._read_tsv(
    root / "prediction_bundles/prediction_bundle_index.tsv"
)
record = {row["prediction_set_id"]: row for row in index_rows}[
    "training_stage_distribution--endpoint-stage3"
]
model_id = "training_stage_distribution"
endpoint_id = "stage3"
production_sha256 = "ebbe02b9cf0e11ea6d01ccf981652631801b382fbd12225710b9fe17013da59d"
bundle_root = root / "prediction_bundles/training_stage_distribution--endpoint-stage3"
manifest = verify_frozen_tree(bundle_root)
document = bundle_root / "prediction_bundle.json"
bundle = PredictionBundle.load_json(document)
bundle.validate_artifacts(bundle_root)
standardized_fields = auditor.STANDARDIZED_BASE_FIELDS + auditor.ENDPOINT_VALUE_FIELDS[
    endpoint_id
]
table_schema = canonical_sha256(
    {"format": "tsv", "fields": list(standardized_fields)}
)
join_key = canonical_sha256(
    {
        "task_id": auditor.TASK_ID,
        "dataset_ids": list(auditor.DATASET_IDS),
        "split_id": auditor.SPLIT_ID,
        "row_id_field": "prediction_row_id",
        "unit_id_field": "participant_id",
        "unit_id_namespace": auditor.UNIT_NAMESPACE,
        "biological_unit": "participant",
    }
)
identity = {
    "run_id": production_sha256,
    "task_id": auditor.TASK_ID,
    "model_id": model_id,
    "endpoint_id": endpoint_id,
    "prediction_aggregation": auditor.SEED_AGGREGATION,
    "dataset_ids": list(auditor.DATASET_IDS),
    "split_id": auditor.SPLIT_ID,
    "standardized_table_sha256": bundle.standardized_table.sha256,
    "row_ids_sha256": bundle.row_ids.sha256,
    "n_predictions": auditor.EXPECTED_PARTICIPANTS,
    "source_join_key_sha256": join_key,
}
checks = {
    "record_model": record.get("model_id") == model_id,
    "record_endpoint": record.get("endpoint_id") == endpoint_id,
    "record_aggregation": record.get("prediction_aggregation") == auditor.SEED_AGGREGATION,
    "record_source_filename": record.get("source_filename") == f"{model_id}.tsv",
    "record_source_sha": record.get("source_prediction_sha256")
    == lineage[model_id]["ensemble_prediction_sha256"],
    "record_bundle_id": record.get("bundle_id") == bundle.bundle_id,
    "record_document_sha": record.get("bundle_document_sha256") == sha256_file(document),
    "record_artifacts_sha": record.get("bundle_artifacts_sha256")
    == sha256_file(bundle_root / "ARTIFACTS.json"),
    "manifest_outcomes": manifest.get("metadata", {}).get("outcomes_read") is False,
    "manifest_metrics": manifest.get("metadata", {}).get("metrics_calculated") is False,
    "run_id": bundle.run_id == production_sha256,
    "bundle_identity": bundle.bundle_id == canonical_sha256(identity),
    "task_id": bundle.task_id == auditor.TASK_ID,
    "model_id": bundle.model_id == model_id,
    "dataset_ids": tuple(bundle.dataset_ids) == auditor.DATASET_IDS,
    "split_id": bundle.split_id == auditor.SPLIT_ID,
    "n_predictions": bundle.n_predictions == auditor.EXPECTED_PARTICIPANTS,
    "row_id_field": bundle.row_id_field == "prediction_row_id",
    "unit_id_field": bundle.unit_id_field == "participant_id",
    "unit_namespace": bundle.unit_id_namespace == auditor.UNIT_NAMESPACE,
    "biological_unit": bundle.biological_unit == "participant",
    "table_schema": bundle.table_schema_sha256 == table_schema,
    "join_key": bundle.source_join_key_sha256 == join_key,
    "missing_state": bundle.missing_state is MissingState.OBSERVED,
    "metadata_endpoint": bundle.metadata.get("endpoint_id") == endpoint_id,
    "metadata_aggregation": bundle.metadata.get("prediction_aggregation")
    == auditor.SEED_AGGREGATION,
    "metadata_seeds": bundle.metadata.get("model_seed_roster")
    == list(auditor.reference.SEEDS),
    "metadata_fold": bundle.metadata.get("outer_fold_aggregation")
    == "each_participant_once_from_held_fold",
    "metadata_seed_sha": bundle.metadata.get("source_seed_prediction_sha256")
    == lineage[model_id]["seed_prediction_sha256"],
    "metadata_membership": bundle.metadata.get("five_seed_membership_sha256")
    == lineage[model_id]["five_seed_membership_sha256"],
    "metadata_mutated": bundle.metadata.get("prediction_mutated") is False,
    "metadata_pre_outcome": bundle.metadata.get("frozen_before_outcomes_read") is True,
    "metadata_seed_replicates": bundle.metadata.get("seeds_are_biological_replicates") is False,
    "metadata_champion": bundle.metadata.get("champion_claim_allowed") is False,
    "metadata_external": bundle.metadata.get("external_claim_allowed") is False,
    "metadata_diagnostic": bundle.metadata.get("diagnostic_or_prognostic_claim_allowed") is False,
}
print(json.dumps({"failed": sorted(key for key, value in checks.items() if not value), "checks": checks}, indent=2, sort_keys=True))
PY
