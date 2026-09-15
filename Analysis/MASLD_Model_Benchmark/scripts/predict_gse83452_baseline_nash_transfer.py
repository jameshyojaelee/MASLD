#!/usr/bin/env python3
"""Predict GSE83452 baseline NASH status from frozen GSE135251 models, label-blind.

This step opens GSE83452 expression for the first time.  It never opens the
GSE83452 record table, so the deposited NASH status, age, sex and intervention
cannot reach it even by accident; row identity comes from the frozen label-blind
baseline mask.  Every transform applied to an array here was fit on GSE135251 or
is a function of that single array alone, and the shared representation code is
imported from the source-fit module so the two sides cannot drift apart.

Only the 152 baseline arrays are scored.  The 79 one-year follow-up arrays
belong to the separate frozen paired/intervention task; the mask that selects
them read ``timepoint`` and nothing else.

Per-seed scores are emitted alongside the seed-averaged prediction.  A linear
fit repeated under five schema seeds is one fitted prediction unless the seeds
actually move it, and the receipt records how many distinct seed vectors the
five seeds really produced rather than asserting five replicates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.contracts import ArtifactRef, MissingState, PredictionBundle
from masld_bench.hashing import canonical_sha256

from scripts.fit_gse135251_nash_transfer_source_models import (
    CLASSES,
    LEARNED_MODEL_IDS,
    PRIOR_MODEL_ID,
    apply_pipeline,
    build_representation,
)


TASK_ID = "gse83452_baseline_nash_transfer"
DATASET_IDS = ("antwerp_inserm_shared",)
SPLIT_ID = "train_study_held_gse83452_baseline_external_v1"
UNIT_NAMESPACE = "gse83452_outcome_blind_participant_id"
ENDPOINT_ID = "nash_status"
PREDICTION_FIELDS = (
    "prediction_row_id",
    "row_id",
    "endpoint_id",
    "prediction_state",
    "model_id",
    "probability_no_nash",
    "probability_nash",
    "predicted_nash_status",
)
SEED_SCORE_FIELDS = (
    "prediction_row_id",
    "row_id",
    "model_id",
    "model_seed",
    "probability_nash",
)
# A mask record may never carry these keys at all.  `nash_status`, `age`, `sex`
# and `intervention` are deliberately present as explicit nulls: the TaskSpec
# keeps them structurally withheld from preprocessing, and an explicit null is
# the record of that, not a leak.  A filled value is checked separately.
FORBIDDEN_MASK_KEYS = (
    "fibrosis_stage_group",
    "primary_transfer_evaluable",
    "paired_expression_stress_evaluable",
    "paired_nash_transition_evaluable",
)
STRUCTURALLY_WITHHELD_KEYS = ("nash_status", "age", "sex", "intervention")
# Any of these appearing as a value or a matrix column means the outcome leaked.
OUTCOME_TOKENS = ("nash_status", "nash", "no_nash", "undefined")


class ExternalPredictionError(RuntimeError):
    """Raised when label-blind external prediction would violate its selection record."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ExternalPredictionError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def load_outcome_blind_masks(path: Path) -> list[dict[str, Any]]:
    """Read the frozen baseline mask and refuse anything outcome-bearing."""

    payload = json.loads(path.read_text(encoding="utf-8"))
    if (
        payload.get("status") != "pass_label_blind_participant_masks"
        or payload.get("labels_read") is not False
        or payload.get("outcome_columns_read") is not False
        or payload.get("withheld_fields_imputed") is not False
        or payload.get("selection_is_label_blind") is not True
        or payload.get("arrays_per_participant_after_selection") != 1
    ):
        raise ExternalPredictionError("baseline participant mask is not label-blind")
    records = payload.get("mask_records")
    if not isinstance(records, list) or not records:
        raise ExternalPredictionError("baseline participant mask carries no records")
    for record in records:
        if set(FORBIDDEN_MASK_KEYS) & set(record):
            raise ExternalPredictionError("participant mask carries an outcome column")
        values = {value for value in record.values() if isinstance(value, str)}
        if set(OUTCOME_TOKENS) & values:
            raise ExternalPredictionError("participant mask leaked an outcome value")
        for key in STRUCTURALLY_WITHHELD_KEYS:
            if record.get(key, "sentinel") is not None or record.get(
                f"{key}_observed", "sentinel"
            ) is not False:
                raise ExternalPredictionError(
                    "a structurally withheld field was filled or declared observed"
                )
        if record.get("prediction_eligible") is not True:
            raise ExternalPredictionError("participant mask marks a row ineligible")
        if record.get("baseline_selected") is not True:
            raise ExternalPredictionError("a non-baseline record reached the predictor")
    return sorted(records, key=lambda record: record["participant_id"])


def read_external_matrix(
    path: Path, *, key_field: str, gene_ids: Sequence[str], accessions: Sequence[str]
) -> np.ndarray:
    """Read the external gene matrix onto the shared axis in shared-axis order.

    Only the declared baseline columns are decoded.  The follow-up columns sit in
    the same file and are never read, so no post-intervention array can influence
    a baseline score even through a per-array transform.
    """

    wanted = {value: index for index, value in enumerate(gene_ids)}
    matrix = np.full((len(accessions), len(gene_ids)), np.nan, dtype=np.float64)
    seen = 0
    with path.open(encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if not header or header[0] != key_field:
            raise ExternalPredictionError("external matrix row key differs")
        columns = header[1:]
        for token in OUTCOME_TOKENS:
            if token in columns:
                raise ExternalPredictionError("external matrix carries an outcome column")
        position = {value: index for index, value in enumerate(columns)}
        if len(position) != len(columns) or any(
            value not in position for value in accessions
        ):
            raise ExternalPredictionError("external matrix lacks a declared array")
        take = [position[value] for value in accessions]
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            row = wanted.get(fields[0])
            if row is None:
                continue
            if len(fields) != len(columns) + 1:
                raise ExternalPredictionError("external matrix row is ragged")
            matrix[:, row] = [float(fields[1 + index]) for index in take]
            seen += 1
    if seen != len(gene_ids) or not np.all(np.isfinite(matrix)):
        raise ExternalPredictionError("external matrix is incomplete on the shared axis")
    return matrix


def _prediction_row_id(run_id: str, model_id: str, row_id: str) -> str:
    return hashlib.sha256(
        f"{run_id}\0{model_id}\0{ENDPOINT_ID}\0{row_id}".encode("utf-8")
    ).hexdigest()


def _freeze_bundle(
    *,
    output: Path,
    model_id: str,
    row_ids: Sequence[str],
    probabilities: np.ndarray,
    seed_probabilities: np.ndarray,
    seeds: Sequence[int],
    run_id: str,
    arm_id: str,
    source_fit_artifacts_sha256: str,
    activation_artifacts_sha256: str,
    summarization_artifacts_sha256: str,
    distinct_seed_vectors: int,
    seed_dependent: bool,
) -> dict[str, Any]:
    target = output / model_id
    target.mkdir()
    rows: list[dict[str, Any]] = []
    row_id_rows: list[dict[str, str]] = []
    seed_rows: list[dict[str, Any]] = []
    for index, (row_id, probability) in enumerate(
        zip(row_ids, probabilities, strict=True)
    ):
        prediction_row_id = _prediction_row_id(run_id, model_id, row_id)
        rows.append(
            {
                "prediction_row_id": prediction_row_id,
                "row_id": row_id,
                "endpoint_id": ENDPOINT_ID,
                "prediction_state": "observed",
                "model_id": model_id,
                "probability_no_nash": format(float(1.0 - probability), ".17g"),
                "probability_nash": format(float(probability), ".17g"),
                "predicted_nash_status": CLASSES[1] if probability > 0.5 else CLASSES[0],
            }
        )
        row_id_rows.append({"prediction_row_id": prediction_row_id, "row_id": row_id})
        for seed_index, seed in enumerate(seeds):
            seed_rows.append(
                {
                    "prediction_row_id": prediction_row_id,
                    "row_id": row_id,
                    "model_id": model_id,
                    "model_seed": seed,
                    "probability_nash": format(
                        float(seed_probabilities[seed_index, index]), ".17g"
                    ),
                }
            )
    table = target / "predictions.tsv"
    row_inventory = target / "row_ids.tsv"
    seed_table = target / "seed_scores.tsv"
    write_tsv(table, PREDICTION_FIELDS, rows)
    write_tsv(row_inventory, ("prediction_row_id", "row_id"), row_id_rows)
    write_tsv(seed_table, SEED_SCORE_FIELDS, seed_rows)
    table_ref = ArtifactRef.from_path(
        table,
        relative_to=target,
        media_type="text/tab-separated-values",
        role=f"standardized_prediction_table:{TASK_ID}",
    )
    row_ref = ArtifactRef.from_path(
        row_inventory,
        relative_to=target,
        media_type="text/tab-separated-values",
        role=f"prediction_row_ids:{TASK_ID}",
    )
    seed_ref = ArtifactRef.from_path(
        seed_table,
        relative_to=target,
        media_type="text/tab-separated-values",
        role=f"seed_level_scores:{TASK_ID}",
    )
    join_key = canonical_sha256(
        {
            "task_id": TASK_ID,
            "dataset_ids": list(DATASET_IDS),
            "split_id": SPLIT_ID,
            "row_id_field": "prediction_row_id",
            "unit_id_field": "row_id",
            "unit_id_namespace": UNIT_NAMESPACE,
            "biological_unit": "participant",
        }
    )
    identity = {
        "run_id": run_id,
        "task_id": TASK_ID,
        "arm_id": arm_id,
        "model_id": model_id,
        "dataset_ids": list(DATASET_IDS),
        "split_id": SPLIT_ID,
        "table_sha256": table_ref.sha256,
        "row_ids_sha256": row_ref.sha256,
        "seed_scores_sha256": seed_ref.sha256,
        "n_predictions": len(row_ids),
        "source_join_key_sha256": join_key,
    }
    bundle = PredictionBundle(
        schema_version="masld-bench-prediction-bundle-v1",
        bundle_id=canonical_sha256(identity),
        run_id=run_id,
        task_id=TASK_ID,
        model_id=model_id,
        dataset_ids=DATASET_IDS,
        split_id=SPLIT_ID,
        artifacts=(table_ref, row_ref, seed_ref),
        standardized_table=table_ref,
        row_ids=row_ref,
        n_predictions=len(row_ids),
        row_id_field="prediction_row_id",
        unit_id_field="row_id",
        unit_id_namespace=UNIT_NAMESPACE,
        biological_unit="participant",
        table_schema_sha256=canonical_sha256(
            {"format": "tsv", "fields": list(PREDICTION_FIELDS)}
        ),
        source_join_key_sha256=join_key,
        format_version="gse83452-external-baseline-nash-status-prediction-tsv-v1",
        missing_state=MissingState.OBSERVED,
        metadata={
            "endpoint_id": ENDPOINT_ID,
            "arm_id": arm_id,
            "positive_class": CLASSES[1],
            "timepoint": "baseline",
            "prediction_role": "five_seed_source_fit_external_development_transfer",
            "prediction_aggregation": "mean_probability_across_five_source_model_seeds",
            "model_seed_roster": list(seeds),
            "seed_level_results_emitted": True,
            "seed_dependent": seed_dependent,
            "distinct_seed_prediction_vectors": distinct_seed_vectors,
            "source_fit_artifacts_sha256": source_fit_artifacts_sha256,
            "activation_artifacts_sha256": activation_artifacts_sha256,
            "summarization_artifacts_sha256": summarization_artifacts_sha256,
            "external_labels_read": False,
            "external_participant_table_read": False,
            "query_fit_or_calibration_performed": False,
            "query_participants_jointly_normalized": False,
            "prediction_frozen_before_evaluator_label_join": True,
            "project_sealed": False,
            "champion_claim_allowed": False,
            "diagnostic_or_prognostic_claim_allowed": False,
            "clinical_claim_allowed": False,
        },
    )
    document = target / "prediction_bundle.json"
    write_json_exclusive(document, bundle.to_dict(), mode=0o440)
    loaded = PredictionBundle.load_json(document)
    loaded.validate_artifacts(target)
    artifacts_sha256 = freeze_tree(
        target,
        {
            "artifact_class": "gse83452_baseline_nash_transfer_prediction_bundle",
            "bundle_id": bundle.bundle_id,
            "arm_id": arm_id,
            "model_id": model_id,
            "participants": len(row_ids),
            "biological_unit": "participant",
            "external_labels_read": False,
            "query_fit_or_calibration_performed": False,
            "metrics_calculated": False,
            "project_sealed": False,
            "status": "passed_unscored",
        },
    )
    return {
        "model_id": model_id,
        "bundle_path": f"{model_id}/prediction_bundle.json",
        "bundle_id": bundle.bundle_id,
        "bundle_document_sha256": sha256_file(document),
        "bundle_artifacts_sha256": artifacts_sha256,
        "prediction_table_sha256": table_ref.sha256,
        "seed_scores_sha256": seed_ref.sha256,
        "distinct_seed_prediction_vectors": distinct_seed_vectors,
    }


def _check_hash(path: Path, expected: str, label: str) -> None:
    if sha256_file(path) != expected:
        raise ExternalPredictionError(f"{label} SHA-256 differs")


def predict_external(
    *, benchmark_root: Path, contract_path: Path, contract_sha256: str, output: Path
) -> dict[str, Any]:
    from masld_bench.artifacts import verify_frozen_tree

    if output.exists():
        raise ExternalPredictionError(f"refusing to overwrite predictions: {output}")
    _check_hash(contract_path, contract_sha256, "external prediction contract")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if (
        contract.get("status") != "registered_prediction_only_pending"
        or contract.get("external_labels_read") is not False
        or contract.get("query_fit_or_calibration_allowed") is not False
        or contract.get("external_metrics_allowed") is not False
    ):
        raise ExternalPredictionError("external prediction contract differs")

    source_spec = contract["source_fit"]
    source_fit = benchmark_root / source_spec["path"]
    _check_hash(
        source_fit / "ARTIFACTS.json", source_spec["artifacts_sha256"], "source fit"
    )
    verify_frozen_tree(source_fit)
    fit_receipt = json.loads(
        (source_fit / "fit_receipt.json").read_text(encoding="utf-8")
    )
    if (
        fit_receipt.get("status") != source_spec["required_status"]
        or fit_receipt.get("external_labels_read") is not False
        or fit_receipt.get("external_expression_values_read") is not False
        or fit_receipt.get("arm_id") != contract["arm_id"]
        or fit_receipt.get("selected_source_model_id")
        != contract["selected_source_model_id"]
    ):
        raise ExternalPredictionError("source fit receipt differs from the contract")

    activation_spec = contract["activation"]
    activation = benchmark_root / activation_spec["path"]
    _check_hash(
        activation / "ARTIFACTS.json",
        activation_spec["artifacts_sha256"],
        "transfer activation",
    )
    verify_frozen_tree(activation)
    summary_spec = contract["external_expression"]
    summary = benchmark_root / summary_spec["path"]
    _check_hash(
        summary / "ARTIFACTS.json",
        summary_spec["artifacts_sha256"],
        "external summarization",
    )
    verify_frozen_tree(summary)

    records = load_outcome_blind_masks(activation / "baseline_participant_masks.json")
    if len(records) != summary_spec["baseline_participants"]:
        raise ExternalPredictionError("external baseline participant census differs")
    row_ids = [str(record["participant_id"]) for record in records]
    accessions = [str(record["sample_accession"]) for record in records]
    if len(set(row_ids)) != len(row_ids) or len(set(accessions)) != len(accessions):
        raise ExternalPredictionError("external row identity is duplicated")

    fitted_model_ids = tuple(contract["fitted_model_ids"])
    unfittable_model_ids = tuple(contract["unfittable_model_ids"])
    if (
        sorted(fitted_model_ids) != sorted(fit_receipt["fitted_model_ids"])
        or sorted(unfittable_model_ids) != sorted(fit_receipt["unfittable_model_ids"])
        or set(fitted_model_ids) & set(unfittable_model_ids)
    ):
        raise ExternalPredictionError(
            "the contract's fitted and unfittable rosters do not agree with the "
            "source fit receipt"
        )
    if not set(fitted_model_ids) <= set(LEARNED_MODEL_IDS):
        raise ExternalPredictionError("a fitted model is not a registered learned model")
    if not fitted_model_ids:
        raise ExternalPredictionError("no learned model was fitted on the source")
    model_ids = (PRIOR_MODEL_ID,) + tuple(sorted(fitted_model_ids))

    _, axis_rows = read_tsv(activation / "common_stable_gene_axis.tsv")
    shared_ids = [row["stable_gene_id"] for row in axis_rows]
    if len(shared_ids) != activation_spec["shared_genes"]:
        raise ExternalPredictionError("shared gene axis census differs")
    external = read_external_matrix(
        summary / summary_spec["gene_matrix"],
        key_field=summary_spec["key_field"],
        gene_ids=shared_ids,
        accessions=accessions,
    )
    representations = {
        name: build_representation(name=name, log2_shared=external)
        for name in sorted(
            {
                contract["models"][model_id]["representation"]
                for model_id in fitted_model_ids
            }
        )
    }
    spread = {
        name: {
            "external_mean_interquartile_range": float(
                np.mean(
                    np.quantile(values, 0.75, axis=1) - np.quantile(values, 0.25, axis=1)
                )
            ),
            "external_mean_standard_deviation": float(np.mean(np.std(values, axis=1))),
        }
        for name, values in representations.items()
    }

    seeds = [int(value) for value in source_spec["model_seeds"]]
    seed_probabilities: dict[str, np.ndarray] = {}
    for model_id in sorted(fitted_model_ids):
        config = contract["models"][model_id]
        stacked = []
        for seed in seeds:
            seed_root = source_fit / f"models/{model_id}/seed_{seed}"
            verify_frozen_tree(seed_root)
            with np.load(seed_root / "model_state.npz", allow_pickle=False) as loaded:
                state = {key: np.asarray(loaded[key]) for key in loaded.files}
            stacked.append(
                apply_pipeline(
                    config=config,
                    state=state,
                    representation=representations[config["representation"]],
                )
            )
        seed_probabilities[model_id] = np.stack(stacked, axis=0)
    prior = float(fit_receipt["training_prevalence_nash"])
    if not 0.0 < prior < 1.0:
        raise ExternalPredictionError("training prevalence is degenerate")
    seed_probabilities[PRIOR_MODEL_ID] = np.full(
        (len(seeds), len(row_ids)), prior, dtype=np.float64
    )
    if set(seed_probabilities) != set(model_ids):
        raise ExternalPredictionError("external prediction model roster differs")

    probabilities: dict[str, np.ndarray] = {}
    distinct_counts: dict[str, int] = {}
    for model_id, stack in seed_probabilities.items():
        if (
            stack.shape != (len(seeds), len(row_ids))
            or not np.all(np.isfinite(stack))
            or np.any(stack < 0.0)
            or np.any(stack > 1.0)
        ):
            raise ExternalPredictionError("external seed probabilities are invalid")
        probabilities[model_id] = stack.mean(axis=0)
        distinct_counts[model_id] = int(len(np.unique(np.round(stack, 12), axis=0)))

    output.mkdir(parents=True)
    index_rows = [
        _freeze_bundle(
            output=output,
            model_id=model_id,
            row_ids=row_ids,
            probabilities=probabilities[model_id],
            seed_probabilities=seed_probabilities[model_id],
            seeds=seeds,
            run_id=source_spec["artifacts_sha256"],
            arm_id=contract["arm_id"],
            source_fit_artifacts_sha256=source_spec["artifacts_sha256"],
            activation_artifacts_sha256=activation_spec["artifacts_sha256"],
            summarization_artifacts_sha256=summary_spec["artifacts_sha256"],
            distinct_seed_vectors=distinct_counts[model_id],
            seed_dependent=model_id != PRIOR_MODEL_ID,
        )
        for model_id in model_ids
    ]
    write_tsv(output / "prediction_bundle_index.tsv", tuple(index_rows[0]), index_rows)
    receipt = {
        "schema_version": "masld-bench-gse83452-external-prediction-only-v1",
        "status": "passed_frozen_prediction_only_no_label_access",
        "task_id": TASK_ID,
        "arm_id": contract["arm_id"],
        "gate_eligible_arm": bool(contract["gate_eligible_arm"]),
        "timepoint": "baseline",
        "participants": len(row_ids),
        "biological_unit": "participant",
        "cohort_family_id": "antwerp_inserm_shared",
        "followup_arrays_read": False,
        "shared_genes": len(shared_ids),
        "models": list(model_ids),
        "fitted_model_ids": sorted(fitted_model_ids),
        "unfittable_model_ids": sorted(unfittable_model_ids),
        "unfittable_models": fit_receipt["unfittable_models"],
        "prediction_bundles": len(index_rows),
        "model_seeds": seeds,
        "selected_source_model_id": contract["selected_source_model_id"],
        "training_prevalence_nash": prior,
        "prediction_aggregation": "mean_probability_across_five_source_model_seeds",
        "seed_level_results_emitted": True,
        "distinct_seed_prediction_vectors": dict(sorted(distinct_counts.items())),
        "representation_spread_diagnostics": spread,
        "external_labels_read": False,
        "external_participant_table_read": False,
        "query_fit_or_calibration_performed": False,
        "query_participants_jointly_normalized": False,
        "external_metrics_calculated": False,
        "prediction_frozen_before_evaluator_label_join": True,
        "project_sealed": False,
        "champion_claim_allowed": False,
        "diagnostic_or_prognostic_claim_allowed": False,
        "clinical_claim_allowed": False,
    }
    with (output / "prediction_receipt.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-root", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--contract-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = predict_external(
        benchmark_root=arguments.benchmark_root,
        contract_path=arguments.contract,
        contract_sha256=arguments.contract_sha256,
        output=arguments.output,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
