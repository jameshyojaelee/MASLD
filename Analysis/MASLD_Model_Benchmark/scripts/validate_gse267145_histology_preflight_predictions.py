#!/usr/bin/env python3
"""Validate a frozen one-fold GSE267145 timing preflight without scoring it."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


MODEL_IDS = (
    "training_stage_distribution",
    "rna_hvg_pca_elastic_net",
    "rna_hvg_pca_linear_svm",
    "rna_hvg_pca_nearest_centroid",
    "rna_hvg_pca_knn",
    "h3_variance_pca_elastic_net",
    "h3_variance_pca_linear_svm",
    "h3_variance_pca_nearest_centroid",
    "h3_variance_pca_knn",
    "block_pca_elastic_net",
    "calibrated_late_fusion",
)
PREDICTION_FIELDS = (
    "participant_id",
    "outer_fold",
    "probability_NOR",
    "probability_NAFL",
    "probability_NASH",
    "predicted_stage3",
    "predicted_nash_crn_component_sum",
    "predicted_fibrosis_cumulative_expected",
    "predicted_fibrosis_regression",
    "probability_fibrosis_F0",
    "probability_fibrosis_F1",
    "probability_fibrosis_F2_3",
    "predicted_fibrosis_group3",
)


class HistologyPreflightError(RuntimeError):
    """Raised when a timing preflight crosses the evaluator separation."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise HistologyPreflightError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def _validate_probability_rows(rows: list[dict[str, str]], fields: tuple[str, ...]) -> None:
    for row in rows:
        values = [float(row[field]) for field in fields]
        if any(value < 0 or value > 1 for value in values) or abs(sum(values) - 1) > 1e-6:
            raise HistologyPreflightError("preflight probabilities are invalid")


def validate_preflight(
    *,
    predictions: Path,
    outcomes: Path,
    folds: Path,
    output: Path,
    outer_fold: int,
    seed: int,
) -> dict[str, Any]:
    if output.exists():
        raise HistologyPreflightError(f"refusing to overwrite validation: {output}")
    expected_paths = {predictions / f"{model_id}.tsv" for model_id in MODEL_IDS}
    if set(predictions.glob("*.tsv")) != expected_paths:
        raise HistologyPreflightError("preflight prediction roster differs")
    prediction_ids: list[str] | None = None
    for path in sorted(expected_paths):
        fields, rows = read_tsv(path)
        if fields != PREDICTION_FIELDS or not rows:
            raise HistologyPreflightError("preflight prediction schema differs")
        row_ids = [row["participant_id"] for row in rows]
        if len(row_ids) != len(set(row_ids)) or any(
            int(row["outer_fold"]) != outer_fold for row in rows
        ):
            raise HistologyPreflightError("preflight participant/fold axis differs")
        if prediction_ids is None:
            prediction_ids = row_ids
        elif row_ids != prediction_ids:
            raise HistologyPreflightError("preflight model participant axes differ")
        _validate_probability_rows(
            rows, ("probability_NOR", "probability_NAFL", "probability_NASH")
        )
        _validate_probability_rows(
            rows,
            (
                "probability_fibrosis_F0",
                "probability_fibrosis_F1",
                "probability_fibrosis_F2_3",
            ),
        )
    # The evaluator-only join occurs only after every frozen prediction file passes.
    endpoint_fields, endpoint_rows = read_tsv(outcomes / "participant_endpoints.tsv")
    fold_fields, fold_rows = read_tsv(folds / "participant_outer_folds.tsv")
    if not {"participant_id", "outer_fold", "stage3", "fibrosis", "recorded_sex"} <= set(
        endpoint_fields
    ) or fold_fields != ("participant_id", "outer_fold"):
        raise HistologyPreflightError("evaluator input schema differs")
    expected_ids = [
        row["participant_id"] for row in fold_rows if int(row["outer_fold"]) == outer_fold
    ]
    if prediction_ids != expected_ids:
        raise HistologyPreflightError("preflight query IDs differ from the frozen fold")
    endpoint_by_id = {row["participant_id"]: row for row in endpoint_rows}
    held = [endpoint_by_id[participant] for participant in expected_ids]
    f3_total = sum(int(row["fibrosis"]) == 3 for row in endpoint_rows)
    male_rows = [row for row in endpoint_rows if row["recorded_sex"] == "M"]
    if f3_total != 4 or len(male_rows) != 14 or {row["stage3"] for row in male_rows} != {"NASH"}:
        raise HistologyPreflightError("machine-encoded evaluator limitations differ")
    output.mkdir(parents=True)
    disposition = {
        "schema_version": "masld-bench-gse267145-histology-preflight-disposition-v1",
        "status": "passed_firewall_validation_no_metrics",
        "outer_fold": outer_fold,
        "seed": seed,
        "query_participants": len(expected_ids),
        "held_f3_participants": sum(int(row["fibrosis"]) == 3 for row in held),
        "f3_limitation": {
            "code": "sparse_f3_inner_training_instability",
            "cohort_f3_participants": 4,
            "minimum_inner_training_f3": 1,
            "champion_gate_eligible": False,
        },
        "recorded_sex_limitation": {
            "code": "recorded_male_perfectly_confound_with_nash",
            "male_participants": 14,
            "male_stage3_classes": ["NASH"],
            "male_macro_f1_applicable": False,
            "sex_fairness_claim_allowed": False,
        },
        "metrics_calculated": False,
        "outcome_values_exported": False,
        "production_or_champion_claim_allowed": False,
    }
    (output / "disposition.json").write_text(
        json.dumps(disposition, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return disposition


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--predictions-artifacts-sha256", required=True)
    parser.add_argument("--outcomes", required=True, type=Path)
    parser.add_argument("--outcomes-artifacts-sha256", required=True)
    parser.add_argument("--folds", required=True, type=Path)
    parser.add_argument("--folds-artifacts-sha256", required=True)
    parser.add_argument("--outer-fold", required=True, type=int)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    for root, expected in (
        (arguments.predictions, arguments.predictions_artifacts_sha256),
        (arguments.outcomes, arguments.outcomes_artifacts_sha256),
        (arguments.folds, arguments.folds_artifacts_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise HistologyPreflightError(f"input ARTIFACTS SHA differs: {root}")
    result = validate_preflight(
        predictions=arguments.predictions,
        outcomes=arguments.outcomes,
        folds=arguments.folds,
        output=arguments.output,
        outer_fold=arguments.outer_fold,
        seed=arguments.seed,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
