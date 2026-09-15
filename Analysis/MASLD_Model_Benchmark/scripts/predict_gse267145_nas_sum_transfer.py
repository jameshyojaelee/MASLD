#!/usr/bin/env python3
"""Predict the GSE267145 activity sum from frozen GSE135251 models, outcome-blind.

This step opens GSE267145 expression for the first time.  It never opens the
GSE267145 histology table, so the NASH-CRN component sum cannot reach it even by
accident; participant identity comes from the frozen label-blind roster.

Predictions are continuous and are never thresholded into classes.  Per-seed
scores are emitted alongside the seed mean so the evaluator can measure, rather
than assume, whether five schema seeds are five replicates or one fit.
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

from scripts.fit_gse267145_fibrosis_transfer_source_models import build_representation
from scripts.fit_gse135251_nas_sum_transfer_source_models import (
    LEARNED_MODEL_IDS,
    PRIOR_MODEL_ID,
    apply_pipeline,
    log2_cpm_complete_axis,
)

TASK_ID = "gse267145_nas_sum_transfer"
DATASET_IDS = ("gse267145_znf469_human_liver",)
SPLIT_ID = "train_gse135251_held_gse267145_nas_sum_external_v1"
UNIT_NAMESPACE = "gse267145_outcome_blind_participant_id"
ENDPOINT_ID = "nash_crn_component_sum"
PREDICTION_FIELDS = (
    "prediction_row_id", "row_id", "endpoint_id", "prediction_state",
    "model_id", "predicted_nas_sum",
)
SEED_FIELDS = ("prediction_row_id", "row_id", "model_id", "model_seed", "predicted_nas_sum")
FORBIDDEN_ROSTER_COLUMNS = (
    "nash_crn_component_sum", "steatosis", "ballooning", "lobular_inflammation",
    "lobular_necrosis", "fibrosis", "stage3", "stage5", "nas_score",
)


class NasPredictionError(RuntimeError):
    """Raised when outcome-blind prediction would violate its selection record."""


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
            raise NasPredictionError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def load_outcome_blind_roster(path: Path) -> list[dict[str, str]]:
    fields, rows = read_tsv(path)
    leaked = sorted(set(FORBIDDEN_ROSTER_COLUMNS) & set(fields))
    if leaked:
        raise NasPredictionError(f"target roster carries outcomes: {leaked}")
    if "participant_id" not in fields:
        raise NasPredictionError("target roster lacks participant_id")
    ids = [row["participant_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise NasPredictionError("target roster repeats a participant")
    return sorted(rows, key=lambda row: row["participant_id"])


def _prediction_row_id(run_id: str, model_id: str, row_id: str) -> str:
    return hashlib.sha256(
        f"{run_id}\0{model_id}\0{ENDPOINT_ID}\0{row_id}".encode("utf-8")
    ).hexdigest()


def _freeze_bundle(
    *, output: Path, model_id: str, row_ids: Sequence[str], predictions: np.ndarray,
    seed_predictions: np.ndarray, seeds: Sequence[int], run_id: str,
    source_fit_sha: str, activation_sha: str, target_sha: str,
    min_pairwise_correlation: float, seed_dependent: bool,
) -> dict[str, Any]:
    target = output / model_id
    target.mkdir()
    rows, id_rows, seed_rows = [], [], []
    for index, (row_id, value) in enumerate(zip(row_ids, predictions, strict=True)):
        prediction_row_id = _prediction_row_id(run_id, model_id, row_id)
        rows.append({
            "prediction_row_id": prediction_row_id, "row_id": row_id,
            "endpoint_id": ENDPOINT_ID, "prediction_state": "observed",
            "model_id": model_id, "predicted_nas_sum": format(float(value), ".17g"),
        })
        id_rows.append({"prediction_row_id": prediction_row_id, "row_id": row_id})
        for si, seed in enumerate(seeds):
            seed_rows.append({
                "prediction_row_id": prediction_row_id, "row_id": row_id,
                "model_id": model_id, "model_seed": seed,
                "predicted_nas_sum": format(float(seed_predictions[si, index]), ".17g"),
            })
    table, ids_table, seed_table = (
        target / "predictions.tsv", target / "row_ids.tsv", target / "seed_scores.tsv"
    )
    write_tsv(table, PREDICTION_FIELDS, rows)
    write_tsv(ids_table, ("prediction_row_id", "row_id"), id_rows)
    write_tsv(seed_table, SEED_FIELDS, seed_rows)
    refs = [
        ArtifactRef.from_path(p, relative_to=target,
                              media_type="text/tab-separated-values", role=role)
        for p, role in (
            (table, f"standardized_prediction_table:{TASK_ID}"),
            (ids_table, f"prediction_row_ids:{TASK_ID}"),
            (seed_table, f"seed_level_scores:{TASK_ID}"),
        )
    ]
    join_key = canonical_sha256({
        "task_id": TASK_ID, "dataset_ids": list(DATASET_IDS), "split_id": SPLIT_ID,
        "row_id_field": "prediction_row_id", "unit_id_field": "row_id",
        "unit_id_namespace": UNIT_NAMESPACE, "biological_unit": "participant",
    })
    identity = {
        "run_id": run_id, "task_id": TASK_ID, "model_id": model_id,
        "dataset_ids": list(DATASET_IDS), "split_id": SPLIT_ID,
        "table_sha256": refs[0].sha256, "row_ids_sha256": refs[1].sha256,
        "seed_scores_sha256": refs[2].sha256, "n_predictions": len(row_ids),
        "source_join_key_sha256": join_key,
    }
    bundle = PredictionBundle(
        schema_version="masld-bench-prediction-bundle-v1",
        bundle_id=canonical_sha256(identity), run_id=run_id, task_id=TASK_ID,
        model_id=model_id, dataset_ids=DATASET_IDS, split_id=SPLIT_ID,
        artifacts=tuple(refs), standardized_table=refs[0], row_ids=refs[1],
        n_predictions=len(row_ids), row_id_field="prediction_row_id",
        unit_id_field="row_id", unit_id_namespace=UNIT_NAMESPACE,
        biological_unit="participant",
        table_schema_sha256=canonical_sha256(
            {"format": "tsv", "fields": list(PREDICTION_FIELDS)}),
        source_join_key_sha256=join_key,
        format_version="gse267145-nas-sum-prediction-tsv-v1",
        missing_state=MissingState.OBSERVED,
        metadata={
            "endpoint_id": ENDPOINT_ID, "endpoint_is_continuous": True,
            "endpoint_binarised": False,
            "prediction_role": "five_seed_source_fit_external_development_transfer",
            "prediction_aggregation": "mean_prediction_across_five_source_model_seeds",
            "model_seed_roster": list(seeds), "seed_level_results_emitted": True,
            "seed_dependent": seed_dependent,
            "min_pairwise_correlation_across_seeds": min_pairwise_correlation,
            "source_fit_artifacts_sha256": source_fit_sha,
            "activation_artifacts_sha256": activation_sha,
            "target_molecular_artifacts_sha256": target_sha,
            "target_outcomes_read": False,
            "target_histology_table_read": False,
            "query_fit_or_calibration_performed": False,
            "recentred_or_rescaled_onto_the_target": False,
            "prediction_frozen_before_evaluator_label_join": True,
            "champion_claim_allowed": False,
        },
    )
    document = target / "prediction_bundle.json"
    write_json_exclusive(document, bundle.to_dict(), mode=0o440)
    PredictionBundle.load_json(document).validate_artifacts(target)
    artifacts_sha = freeze_tree(target, {
        "artifact_class": "gse267145_nas_sum_prediction_bundle",
        "bundle_id": bundle.bundle_id, "model_id": model_id,
        "participants": len(row_ids), "target_outcomes_read": False,
        "metrics_calculated": False, "status": "passed_unscored",
    })
    return {
        "model_id": model_id, "bundle_id": bundle.bundle_id,
        "bundle_artifacts_sha256": artifacts_sha,
        "prediction_table_sha256": refs[0].sha256,
        "seed_scores_sha256": refs[2].sha256,
        "min_pairwise_correlation_across_seeds": min_pairwise_correlation,
    }


def _min_pairwise_correlation(stack: np.ndarray) -> float:
    if stack.shape[0] < 2:
        return 1.0
    centred = stack - stack.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(centred, axis=1, keepdims=True)
    if np.any(norms <= 0):
        return 1.0
    unit = centred / norms
    corr = unit @ unit.T
    return float(np.min(corr[np.triu_indices(stack.shape[0], k=1)]))


def _check_hash(path: Path, expected: str, label: str) -> None:
    if sha256_file(path) != expected:
        raise NasPredictionError(f"{label} SHA-256 differs")


def predict(*, benchmark_root: Path, contract_path: Path, contract_sha256: str,
            output: Path) -> dict[str, Any]:
    from masld_bench.artifacts import verify_frozen_tree

    if output.exists():
        raise NasPredictionError(f"refusing to overwrite predictions: {output}")
    _check_hash(contract_path, contract_sha256, "prediction contract")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if (
        contract.get("status") != "registered_prediction_only_pending"
        or contract.get("target_outcomes_read") is not False
        or contract.get("query_fit_or_calibration_allowed") is not False
        or contract.get("target_metrics_allowed") is not False
    ):
        raise NasPredictionError("prediction contract differs")

    fit_spec = contract["source_fit"]
    source_fit = benchmark_root / fit_spec["path"]
    _check_hash(source_fit / "ARTIFACTS.json", fit_spec["artifacts_sha256"], "source fit")
    verify_frozen_tree(source_fit)
    fit_receipt = json.loads((source_fit / "fit_receipt.json").read_text(encoding="utf-8"))
    if (
        fit_receipt.get("status") != fit_spec["required_status"]
        or fit_receipt.get("target_outcomes_read") is not False
        or fit_receipt.get("selected_source_model_id") != contract["selected_source_model_id"]
        or sorted(fit_receipt["fitted_model_ids"]) != sorted(contract["fitted_model_ids"])
    ):
        raise NasPredictionError("source fit receipt differs from the contract")

    act_spec = contract["activation"]
    activation = benchmark_root / act_spec["path"]
    _check_hash(activation / "ARTIFACTS.json", act_spec["artifacts_sha256"], "activation")
    verify_frozen_tree(activation)
    tgt_spec = contract["target_expression"]
    target_root = benchmark_root / tgt_spec["path"]
    _check_hash(target_root / "ARTIFACTS.json", tgt_spec["artifacts_sha256"], "target molecular")
    verify_frozen_tree(target_root)

    roster = load_outcome_blind_roster(activation / "target_participant_roster.tsv")
    if len(roster) != tgt_spec["participants"]:
        raise NasPredictionError("target participant census differs")
    row_ids = [row["participant_id"] for row in roster]

    _, axis_rows = read_tsv(activation / "common_stable_gene_axis.tsv")
    shared_ids = [row["stable_gene_id"] for row in axis_rows]
    if len(shared_ids) != act_spec["shared_genes"]:
        raise NasPredictionError("shared gene axis census differs")
    _, target_axis = read_tsv(target_root / "rna_feature_axis.tsv")
    target_genes = [row["stable_gene_id"] for row in target_axis]
    lookup = {v: i for i, v in enumerate(target_genes)}
    if any(v not in lookup for v in shared_ids):
        raise NasPredictionError("a shared gene is absent from the target axis")
    _, target_participants = read_tsv(target_root / "participant_axis.tsv")
    order = {row["participant_id"]: index for index, row in enumerate(target_participants)}
    if set(order) != set(row_ids):
        raise NasPredictionError("target participant identities differ from the roster")

    values = np.load(target_root / "rna_values.npy", mmap_mode="r", allow_pickle=False)
    if values.shape != (len(order), len(target_genes)):
        raise NasPredictionError("target RNA matrix shape differs")
    rows_index = [order[value] for value in row_ids]
    log2_shared = log2_cpm_complete_axis(np.asarray(values))[
        np.ix_(rows_index, [lookup[v] for v in shared_ids])
    ]
    representations = {
        name: build_representation(name=name, log2_shared=log2_shared)
        for name in sorted({contract["models"][m]["representation"]
                            for m in contract["fitted_model_ids"]})
    }

    seeds = [int(v) for v in fit_spec["model_seeds"]]
    stacks: dict[str, np.ndarray] = {}
    for model_id in sorted(contract["fitted_model_ids"]):
        config = contract["models"][model_id]
        per_seed = []
        for seed in seeds:
            seed_root = source_fit / f"models/{model_id}/seed_{seed}"
            verify_frozen_tree(seed_root)
            with np.load(seed_root / "model_state.npz", allow_pickle=False) as loaded:
                state = {k: np.asarray(loaded[k]) for k in loaded.files}
            per_seed.append(apply_pipeline(
                state=state, representation=representations[config["representation"]]))
        stacks[model_id] = np.stack(per_seed, axis=0)
    prior = float(fit_receipt["training_mean_nas_sum"])
    stacks[PRIOR_MODEL_ID] = np.full((len(seeds), len(row_ids)), prior, dtype=np.float64)

    model_ids = (PRIOR_MODEL_ID,) + tuple(sorted(contract["fitted_model_ids"]))
    predictions, correlations = {}, {}
    for model_id, stack in stacks.items():
        if stack.shape != (len(seeds), len(row_ids)) or not np.all(np.isfinite(stack)):
            raise NasPredictionError("seed predictions are invalid")
        predictions[model_id] = stack.mean(axis=0)
        correlations[model_id] = _min_pairwise_correlation(stack)

    output.mkdir(parents=True)
    index_rows = [
        _freeze_bundle(
            output=output, model_id=model_id, row_ids=row_ids,
            predictions=predictions[model_id], seed_predictions=stacks[model_id],
            seeds=seeds, run_id=fit_spec["artifacts_sha256"],
            source_fit_sha=fit_spec["artifacts_sha256"],
            activation_sha=act_spec["artifacts_sha256"],
            target_sha=tgt_spec["artifacts_sha256"],
            min_pairwise_correlation=correlations[model_id],
            seed_dependent=model_id != PRIOR_MODEL_ID,
        )
        for model_id in model_ids
    ]
    write_tsv(output / "prediction_bundle_index.tsv", tuple(index_rows[0]), index_rows)
    receipt = {
        "schema_version": "masld-bench-gse267145-nas-sum-prediction-only-v1",
        "status": "passed_frozen_prediction_only_no_outcome_access",
        "task_id": TASK_ID, "participants": len(row_ids),
        "biological_unit": "participant", "shared_genes": len(shared_ids),
        "models": list(model_ids), "prediction_bundles": len(index_rows),
        "model_seeds": seeds,
        "selected_source_model_id": contract["selected_source_model_id"],
        "training_mean_nas_sum": prior,
        "endpoint_is_continuous": True, "endpoint_binarised": False,
        "prediction_aggregation": "mean_prediction_across_five_source_model_seeds",
        "seed_level_results_emitted": True,
        "min_pairwise_correlation_across_seeds": dict(sorted(correlations.items())),
        "predicted_range": {
            k: [float(np.min(v)), float(np.max(v))] for k, v in sorted(predictions.items())
        },
        "target_outcomes_read": False, "target_histology_table_read": False,
        "query_fit_or_calibration_performed": False,
        "recentred_or_rescaled_onto_the_target": False,
        "target_metrics_calculated": False,
        "prediction_frozen_before_evaluator_label_join": True,
        "champion_claim_allowed": False,
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
    a = parser.parse_args()
    print(json.dumps(predict(benchmark_root=a.benchmark_root, contract_path=a.contract,
                             contract_sha256=a.contract_sha256, output=a.output),
                     sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
