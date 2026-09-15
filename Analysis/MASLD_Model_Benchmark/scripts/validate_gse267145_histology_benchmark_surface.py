#!/usr/bin/env python3
"""Validate and freeze the no-fit GSE267145 histology benchmark design."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping


class HistologySurfaceError(RuntimeError):
    """Raised when the benchmark surface is unsafe or inconsistent."""


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
            raise HistologySurfaceError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def _seeded_order(seed: int, outer_fold: int, namespace: str, participant: str) -> bytes:
    return sha256(
        f"{seed}\0{outer_fold}\0{namespace}\0{participant}".encode("utf-8")
    ).digest()


def _fibrosis_group(value: int) -> str:
    return "F0" if value == 0 else "F1" if value == 1 else "F2_3"


def _partition_census(rows: list[dict[str, str]]) -> dict[str, Any]:
    fibrosis = [int(row["fibrosis"]) for row in rows]
    nas = [int(row["nash_crn_component_sum"]) for row in rows]
    return {
        "participants": len(rows),
        "stage3": dict(sorted(Counter(row["stage3"] for row in rows).items())),
        "fibrosis_exact": {
            str(key): value for key, value in sorted(Counter(fibrosis).items())
        },
        "fibrosis_group3": dict(
            sorted(Counter(_fibrosis_group(value) for value in fibrosis).items())
        ),
        "fibrosis_thresholds": {
            f"ge_{threshold}": {
                "negative": sum(value < threshold for value in fibrosis),
                "positive": sum(value >= threshold for value in fibrosis),
            }
            for threshold in (1, 2, 3)
        },
        "nash_crn_component_sum": {
            "distinct_values": len(set(nas)),
            "value_counts": {
                str(key): value for key, value in sorted(Counter(nas).items())
            },
        },
        "source_stage5_descriptive_only": dict(
            sorted(Counter(row["stage5"] for row in rows).items())
        ),
    }


def build_endpoint_safe_inner_contract(
    endpoints: list[dict[str, str]], *, inner_folds: int = 4, seed: int = 20260824
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if inner_folds != 4:
        raise HistologySurfaceError("this surface requires four inner folds")
    participant_order = {row["participant_id"]: index for index, row in enumerate(endpoints)}
    roster: list[dict[str, Any]] = []
    outer_census: dict[str, Any] = {}
    minimum_f3_training = 10**9
    minimum_threshold_class = 10**9
    for outer_fold in range(5):
        training = [row for row in endpoints if int(row["outer_fold"]) != outer_fold]
        if not training:
            raise HistologySurfaceError("outer training partition is empty")
        assignment: dict[str, int] = {}
        validation_sizes = Counter()
        stage_sizes: dict[str, Counter[int]] = {
            stage: Counter() for stage in ("NOR", "NAFL", "NASH")
        }
        fibrosis3 = sorted(
            (row for row in training if int(row["fibrosis"]) == 3),
            key=lambda row: _seeded_order(
                seed, outer_fold, "fibrosis3", row["participant_id"]
            ),
        )
        if len(fibrosis3) < 2:
            raise HistologySurfaceError("outer training has fewer than two fibrosis3 participants")
        start = int.from_bytes(
            _seeded_order(seed, outer_fold, "fibrosis3_start", "roster")[:4], "big"
        ) % inner_folds
        for index, row in enumerate(fibrosis3):
            inner_fold = (start + index) % inner_folds
            participant = row["participant_id"]
            assignment[participant] = inner_fold
            validation_sizes[inner_fold] += 1
            stage_sizes[row["stage3"]][inner_fold] += 1
        for stage in ("NOR", "NAFL", "NASH"):
            remaining = sorted(
                (
                    row
                    for row in training
                    if row["stage3"] == stage and row["participant_id"] not in assignment
                ),
                key=lambda row: (
                    -int(row["fibrosis"]),
                    -int(row["nash_crn_component_sum"]),
                    _seeded_order(seed, outer_fold, stage, row["participant_id"]),
                ),
            )
            for row in remaining:
                inner_fold = min(
                    range(inner_folds),
                    key=lambda candidate: (
                        stage_sizes[stage][candidate],
                        validation_sizes[candidate],
                        _seeded_order(
                            seed,
                            outer_fold,
                            f"{stage}_fold_tie",
                            str(candidate),
                        ),
                    ),
                )
                participant = row["participant_id"]
                assignment[participant] = inner_fold
                validation_sizes[inner_fold] += 1
                stage_sizes[stage][inner_fold] += 1
        if len(assignment) != len(training):
            raise HistologySurfaceError("inner-fold assignment is incomplete")
        inner_census: dict[str, Any] = {}
        for inner_fold in range(inner_folds):
            inner_validation = [
                row for row in training if assignment[row["participant_id"]] == inner_fold
            ]
            inner_training = [
                row for row in training if assignment[row["participant_id"]] != inner_fold
            ]
            training_census = _partition_census(inner_training)
            validation_census = _partition_census(inner_validation)
            if set(training_census["stage3"]) != {"NOR", "NAFL", "NASH"}:
                raise HistologySurfaceError("inner training lacks a stage3 class")
            if set(training_census["fibrosis_group3"]) != {"F0", "F1", "F2_3"}:
                raise HistologySurfaceError("inner training lacks a fibrosis group")
            for counts in training_census["fibrosis_thresholds"].values():
                if counts["negative"] == 0 or counts["positive"] == 0:
                    raise HistologySurfaceError("inner training lacks a threshold class")
                minimum_threshold_class = min(
                    minimum_threshold_class, counts["negative"], counts["positive"]
                )
            if training_census["nash_crn_component_sum"]["distinct_values"] < 2:
                raise HistologySurfaceError("inner training NAS endpoint is constant")
            f3_training = training_census["fibrosis_exact"].get("3", 0)
            if f3_training < 1:
                raise HistologySurfaceError("inner training lacks fibrosis3")
            minimum_f3_training = min(minimum_f3_training, f3_training)
            inner_census[str(inner_fold)] = {
                "training": training_census,
                "validation": validation_census,
            }
        for row in sorted(training, key=lambda item: participant_order[item["participant_id"]]):
            roster.append(
                {
                    "outer_fold": outer_fold,
                    "participant_id": row["participant_id"],
                    "inner_validation_fold": assignment[row["participant_id"]],
                }
            )
        outer_census[str(outer_fold)] = {
            "outer_training": _partition_census(training),
            "inner_folds": inner_census,
        }
    census = {
        "schema_version": "masld-bench-gse267145-endpoint-safe-inner-folds-v1",
        "status": "passed",
        "algorithm_uses": ["stage3", "fibrosis", "nash_crn_component_sum"],
        "algorithm_does_not_use": ["source_stage5", "recorded_sex"],
        "source_stage5_role": "descriptive_only_not_selection",
        "inner_folds": inner_folds,
        "minimum_inner_training_fibrosis3": minimum_f3_training,
        "minimum_inner_training_threshold_class_count": minimum_threshold_class,
        "every_inner_training_endpoint_safe": True,
        "outer_folds": outer_census,
    }
    return roster, census


def validate_surface(
    surface: Mapping[str, Any], root: Path, *, verify_hashes: bool = True
) -> dict[str, Any]:
    if (
        surface.get("schema_version")
        != "masld-bench-gse267145-histology-surface-v1"
        or surface.get("status") != "design_review_required_no_production_fit"
        or surface.get("task_id") != "paired_bulk_histology_state"
        or surface.get("biological_unit") != "participant"
        or surface.get("production_fit_authorized") is not False
    ):
        raise HistologySurfaceError("surface identity or no-fit gate differs")
    inputs = surface.get("inputs")
    if not isinstance(inputs, Mapping) or set(inputs) != {
        "molecular_fixture",
        "participant_folds",
        "evaluator_outcomes",
    }:
        raise HistologySurfaceError("input roles differ")
    resolved: dict[str, Path] = {}
    for role, record in inputs.items():
        if not isinstance(record, Mapping):
            raise HistologySurfaceError(f"{role} input is not an object")
        relative = Path(str(record.get("path", "")))
        if relative.is_absolute() or ".." in relative.parts:
            raise HistologySurfaceError(f"{role} path is unsafe")
        path = root / relative
        if verify_hashes and sha256_file(path / "ARTIFACTS.json") != record.get(
            "artifacts_sha256"
        ):
            raise HistologySurfaceError(f"{role} ARTIFACTS hash differs")
        resolved[role] = path
    participant_fields, participants = read_tsv(
        resolved["molecular_fixture"] / "participant_axis.tsv"
    )
    fold_fields, folds = read_tsv(
        resolved["participant_folds"] / "participant_outer_folds.tsv"
    )
    endpoint_fields, endpoints = read_tsv(
        resolved["evaluator_outcomes"] / "participant_endpoints.tsv"
    )
    if participant_fields != (
        "participant_index",
        "participant_id",
        "rna_source_sample_accession",
        "h3k27ac_source_sample_accession",
        "pairing",
        "rna_observation_state",
        "h3k27ac_observation_state",
    ) or fold_fields != ("participant_id", "outer_fold"):
        raise HistologySurfaceError("outcome-free participant or fold fields differ")
    endpoint_required = {
        "participant_id",
        "outer_fold",
        "stage3",
        "stage5",
        "nash_crn_component_sum",
        "fibrosis",
        "recorded_sex",
    }
    if not endpoint_required <= set(endpoint_fields):
        raise HistologySurfaceError("evaluator endpoint fields differ")
    participant_ids = [row["participant_id"] for row in participants]
    if (
        len(participant_ids) != 99
        or len(set(participant_ids)) != 99
        or participant_ids != [row["participant_id"] for row in folds]
        or participant_ids != [row["participant_id"] for row in endpoints]
        or [row["outer_fold"] for row in folds]
        != [row["outer_fold"] for row in endpoints]
    ):
        raise HistologySurfaceError("participant/fold/outcome axes differ")
    firewall = surface.get("firewall")
    if not isinstance(firewall, Mapping) or any(
        firewall.get(field) is not expected
        for field, expected in (
            ("single_cell_methods_allowed", False),
            ("published_outcome_selected_features_allowed", False),
            ("h3_coordinate_inference_allowed", False),
            ("missing_encoded_as_zero", False),
        )
    ):
        raise HistologySurfaceError("feature or modality firewall differs")
    if (
        firewall.get("participant_age_state") != "structurally_missing"
        or "recorded_sex" not in firewall.get("fit_worker_never_receives", [])
        or "outer_test_endpoint_values"
        not in firewall.get("fit_worker_never_receives", [])
    ):
        raise HistologySurfaceError("metadata or outer-test firewall differs")
    resampling = surface.get("resampling")
    selection = resampling.get("inner_hyperparameter_selection") if isinstance(resampling, Mapping) else None
    if (
        not isinstance(resampling, Mapping)
        or resampling.get("inner_split") != "four_fold_endpoint_safe_within_outer_training_only"
        or resampling.get("inner_split_fail_closed") is not True
        or not isinstance(selection, Mapping)
        or selection.get("rule") != "one_standard_error"
        or selection.get("exact_tie_rule_only") is not False
        or selection.get("applies_to_every_primary_secondary_and_calibration_hyperparameter_choice") is not True
        or selection.get("source_stage5_allowed_to_drive_selection") is not False
        or selection.get("candidate_eligibility_requires") != [
            "convergence_on_all_four_inner_training_fits",
            "convergence_on_the_exact_outer_training_refit_without_held_outcomes",
        ]
        or selection.get("invalid_candidate_action")
        != "audit_and_exclude_before_one_standard_error_selection"
        or selection.get("all_candidates_invalid_action") != "fail_closed"
        or selection.get("selected_prediction_fit")
        != "reuse_exact_candidate_inner_and_outer_training_fits_without_stochastic_rerun"
        or selection.get("candidate_audit_persisted_before_prediction_bundle") is not True
    ):
        raise HistologySurfaceError("inner-fold or one-standard-error contract differs")
    heads = surface.get("shared_task_native_heads")
    regression = heads.get("fibrosis_exact_regression") if isinstance(heads, Mapping) else None
    if (
        not isinstance(regression, Mapping)
        or regression.get("heads") != ["ridge_linear_regression", "elastic_net_linear_regression"]
        or regression.get("prediction_clipping") != [0.0, 3.0]
        or regression.get("reported_separately_from_cumulative_threshold_head") is not True
        or heads.get("source_stage5", {}).get("selection_metric")
        != "not_applicable_descriptive_only"
    ):
        raise HistologySurfaceError("fibrosis regression or stage5 descriptive contract differs")
    lanes = surface.get("lanes")
    if not isinstance(lanes, Mapping) or set(lanes) != {
        "rna_only",
        "h3k27ac_only",
        "observed_pair_multimodal",
    }:
        raise HistologySurfaceError("lane roster differs")
    expected_modalities = {
        "rna_only": ["bulk_rna"],
        "h3k27ac_only": ["h3k27ac_cutrun"],
        "observed_pair_multimodal": ["bulk_rna", "h3k27ac_cutrun"],
    }
    for lane, modalities in expected_modalities.items():
        record = lanes[lane]
        if not isinstance(record, Mapping) or record.get("query_modalities") != modalities:
            raise HistologySurfaceError(f"{lane} modality boundary differs")
        baselines = record.get("baselines")
        if not isinstance(baselines, Mapping) or len(baselines) < 2:
            raise HistologySurfaceError(f"{lane} lacks two task-native baselines")
    for lane_name, prefix in (("rna_only", "rna_hvg_pca"), ("h3k27ac_only", "h3_variance_pca")):
        baseline_ids = set(lanes[lane_name]["baselines"])
        required_ids = {
            f"{prefix}_elastic_net",
            f"{prefix}_linear_svm",
            f"{prefix}_nearest_centroid",
            f"{prefix}_knn",
        }
        if baseline_ids != required_ids:
            raise HistologySurfaceError(f"{lane_name} baseline roster differs")
        for model_id in (f"{prefix}_nearest_centroid", f"{prefix}_knn"):
            spec = lanes[lane_name]["baselines"][model_id]
            if spec.get("probability_temperature") != [0.5, 1.0, 2.0]:
                raise HistologySurfaceError(f"{model_id} probability calibration differs")
        if lanes[lane_name]["baselines"][f"{prefix}_knn"].get("neighbor_counts") != [3, 5, 7, 11]:
            raise HistologySurfaceError(f"{lane_name} kNN roster differs")
    forbidden_text = (
        "scvi",
        "scanvi",
        "scbasset",
        "published_differential_genes_as_features",
        "published_differential_regions_as_features",
    )
    serialized = json.dumps(surface, sort_keys=True).lower()
    if any(token in serialized for token in forbidden_text):
        raise HistologySurfaceError("surface admits a blocked method or feature source")
    evaluator = surface.get("independent_evaluator")
    if (
        not isinstance(evaluator, Mapping)
        or evaluator.get("evaluator_id") != "gse267145_participant_histology_v1"
        or evaluator.get("model_source_import_allowed") is not False
        or evaluator.get("prediction_bundle_mutation_allowed") is not False
        or evaluator.get("primary") != ["stage3_macro_f1"]
        or evaluator.get("age_error_audit")
        != "not_applicable_participant_age_structurally_missing"
    ):
        raise HistologySurfaceError("independent evaluator boundary differs")
    sex_audit = evaluator.get("sex_error_audit")
    if (
        not isinstance(sex_audit, Mapping)
        or sex_audit.get("male_stage3_macro_f1")
        != "not_applicable_missing_NOR_and_NAFL"
        or sex_audit.get("between_sex_global_performance_gap_claim_allowed") is not False
    ):
        raise HistologySurfaceError("sex audit does not acknowledge class confounding")
    observed_sex_stage: dict[str, dict[str, int]] = {}
    for row in endpoints:
        observed_sex_stage.setdefault(row["recorded_sex"], {}).setdefault(row["stage3"], 0)
        observed_sex_stage[row["recorded_sex"]][row["stage3"]] += 1
    expected_census = {"F": {"NAFL": 28, "NASH": 33, "NOR": 24}, "M": {"NASH": 14}}
    if observed_sex_stage != expected_census:
        raise HistologySurfaceError("recorded-sex/stage census differs")
    _, inner_census = build_endpoint_safe_inner_contract(
        endpoints,
        inner_folds=4,
        seed=int(resampling["inner_split_seed"]),
    )
    prediction = surface.get("prediction_contract")
    if (
        not isinstance(prediction, Mapping)
        or prediction.get("rows") != 99
        or prediction.get("artifact_frozen_before_evaluator_join") is not True
        or prediction.get("outer_test_predictions_only") is not True
    ):
        raise HistologySurfaceError("prediction contract differs")
    return {
        "schema_version": "masld-bench-gse267145-histology-surface-validation-v1",
        "status": "passed_design_only_no_fit",
        "surface_sha256": sha256_file(
            root / "config/evaluation/gse267145_histology_benchmark_surface.json"
        ),
        "participants": 99,
        "outer_folds": 5,
        "lanes": list(lanes),
        "baseline_count": 1
        + sum(len(record["baselines"]) for record in lanes.values()),
        "primary_endpoint": "stage3_macro_f1",
        "secondary_endpoint_family": [
            "fibrosis_group3",
            "fibrosis_ordinal",
            "nash_crn_component_sum",
        ],
        "recorded_sex_model_input": False,
        "male_stage3_macro_f1_applicable": False,
        "participant_age_state": "structurally_missing",
        "one_standard_error_selection": True,
        "endpoint_safe_inner_splits": True,
        "minimum_inner_training_fibrosis3": inner_census[
            "minimum_inner_training_fibrosis3"
        ],
        "minimum_inner_training_threshold_class_count": inner_census[
            "minimum_inner_training_threshold_class_count"
        ],
        "production_fit_authorized": False,
        "model_fitted": False,
        "metrics_calculated": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--surface", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    root = arguments.root.resolve(strict=True)
    surface_path = arguments.surface.resolve(strict=True)
    if root not in surface_path.parents:
        raise HistologySurfaceError("surface is outside the benchmark root")
    surface = json.loads(surface_path.read_text(encoding="utf-8"))
    result = validate_surface(surface, root)
    _, endpoints = read_tsv(
        root / surface["inputs"]["evaluator_outcomes"]["path"] / "participant_endpoints.tsv"
    )
    inner_roster, inner_census = build_endpoint_safe_inner_contract(
        endpoints,
        inner_folds=4,
        seed=int(surface["resampling"]["inner_split_seed"]),
    )
    arguments.output.mkdir(parents=True, exist_ok=False)
    (arguments.output / "validation_receipt.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with (arguments.output / "inner_fold_roster.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("outer_fold", "participant_id", "inner_validation_fold"),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(inner_roster)
    (arguments.output / "inner_fold_endpoint_census.json").write_text(
        json.dumps(inner_census, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
