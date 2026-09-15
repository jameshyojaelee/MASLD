#!/usr/bin/env python3
"""Freeze a participant-level GSE267145 histology development task."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from masld_bench.registry import load_task_spec


STAGE5 = ("NOR", "NAFL", "NASH_F0", "NASH_F1", "NASH_F23")
STAGE3 = ("NOR", "NAFL", "NASH")
EXPECTED_STAGE5 = {
    "NAFL": 28,
    "NASH_F0": 19,
    "NASH_F1": 15,
    "NASH_F23": 13,
    "NOR": 24,
}
EXPECTED_SEX = {"F": 85, "M": 14}


class HistologyActivationError(RuntimeError):
    """Raised when one task, participant, or endpoint requirement differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_artifact(root: Path, expected_sha256: str) -> None:
    observed = sha256_file(root / "ARTIFACTS.json")
    if observed != expected_sha256:
        raise HistologyActivationError(
            f"input ARTIFACTS SHA-256 differs for {root}: {observed}"
        )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise HistologyActivationError(f"TSV lacks a header: {path}")
        return list(reader)


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise HistologyActivationError("refusing to write an empty endpoint table")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=tuple(rows[0]),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def integer_field(row: Mapping[str, str], field: str, minimum: int, maximum: int) -> int:
    try:
        value = int(row[field])
    except (KeyError, ValueError) as error:
        raise HistologyActivationError(f"{field} is not an integer") from error
    if value < minimum or value > maximum:
        raise HistologyActivationError(f"{field} is outside its source-native scale")
    return value


def stage3(stage5: str) -> str:
    if stage5 == "NOR":
        return "NOR"
    if stage5 == "NAFL":
        return "NAFL"
    if stage5 in {"NASH_F0", "NASH_F1", "NASH_F23"}:
        return "NASH"
    raise HistologyActivationError("source stage is outside the frozen ontology")


def build_endpoint_records(
    rows: Sequence[Mapping[str, str]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if len(rows) != 99:
        raise HistologyActivationError("paired participant count differs")
    identifiers = [row.get("participant_id", "") for row in rows]
    if len(set(identifiers)) != 99 or any(not value for value in identifiers):
        raise HistologyActivationError("participant identity is missing or duplicated")
    for field in ("rna_gsm", "h3k27ac_gsm"):
        values = [row.get(field, "") for row in rows]
        if len(set(values)) != 99 or any(not value for value in values):
            raise HistologyActivationError(f"{field} is missing or duplicated")

    records: list[dict[str, Any]] = []
    fold_stage: dict[int, Counter[str]] = defaultdict(Counter)
    for row in rows:
        source_stage = row.get("stage", "")
        if source_stage not in STAGE5:
            raise HistologyActivationError("source stage differs")
        steatosis = integer_field(row, "steatosis", 0, 3)
        ballooning = integer_field(row, "ballooning", 0, 2)
        inflammation = integer_field(row, "lobular_inflammation", 0, 3)
        necrosis = integer_field(row, "lobular_necrosis", 0, 2)
        fibrosis = integer_field(row, "fibrosis", 0, 3)
        try:
            fold = int(row["outer_fold"])
        except (KeyError, ValueError) as error:
            raise HistologyActivationError("outer fold is not an integer") from error
        if fold not in range(5):
            raise HistologyActivationError("outer fold is outside the frozen roster")
        if (
            row.get("pairing") != "same_sample_different_aliquot"
            or row.get("rna_measurement_state")
            != "observed_unit_semantics_unresolved"
            or row.get("h3k27ac_measurement_state") != "observed_integer_counts"
            or row.get("sex") not in {"F", "M"}
        ):
            raise HistologyActivationError("pairing, assay state, or recorded sex differs")
        expected_fibrosis = {
            "NOR": {0},
            "NAFL": {0},
            "NASH_F0": {0},
            "NASH_F1": {1},
            "NASH_F23": {2, 3},
        }[source_stage]
        if fibrosis not in expected_fibrosis:
            raise HistologyActivationError("source stage and fibrosis are inconsistent")
        if source_stage == "NOR" and any(
            value != 0 for value in (steatosis, ballooning, inflammation, fibrosis)
        ):
            raise HistologyActivationError("NOR histology components differ")
        collapsed = stage3(source_stage)
        fold_stage[fold][collapsed] += 1
        records.append(
            {
                "participant_id": row["participant_id"],
                "outer_fold": fold,
                "stage3": collapsed,
                "stage5": source_stage,
                "nash_crn_component_sum": steatosis + ballooning + inflammation,
                "steatosis": steatosis,
                "ballooning": ballooning,
                "lobular_inflammation": inflammation,
                "lobular_necrosis": necrosis,
                "fibrosis": fibrosis,
                "recorded_sex": row["sex"],
                "histology_state": "observed",
                "sex_state": "observed",
            }
        )
    stage5_counts = Counter(record["stage5"] for record in records)
    sex_counts = Counter(record["recorded_sex"] for record in records)
    if dict(stage5_counts) != EXPECTED_STAGE5 or dict(sex_counts) != EXPECTED_SEX:
        raise HistologyActivationError("stage or recorded-sex census differs")
    if set(fold_stage) != set(range(5)) or any(
        set(fold_stage[fold]) != set(STAGE3) for fold in range(5)
    ):
        raise HistologyActivationError("a held participant fold lacks a primary class")
    nas_counts = Counter(record["nash_crn_component_sum"] for record in records)
    fibrosis_counts = Counter(record["fibrosis"] for record in records)
    summary = {
        "schema_version": "masld-bench-gse267145-histology-endpoints-v1",
        "status": "passed",
        "participants": 99,
        "biological_unit": "participant",
        "source_stage5_counts": dict(sorted(stage5_counts.items())),
        "stage3_counts": dict(
            sorted(Counter(record["stage3"] for record in records).items())
        ),
        "nash_crn_component_sum_counts": {
            str(key): value for key, value in sorted(nas_counts.items())
        },
        "fibrosis_counts": {
            str(key): value for key, value in sorted(fibrosis_counts.items())
        },
        "recorded_sex_counts": dict(sorted(sex_counts.items())),
        "fold_stage3_counts": {
            str(fold): dict(sorted(fold_stage[fold].items())) for fold in range(5)
        },
        "all_endpoint_states": "observed",
        "nas_derivation": "steatosis_plus_ballooning_plus_lobular_inflammation",
        "lobular_necrosis_in_nas": False,
        "cohort_summary_age_or_bmi_assigned_to_participants": False,
    }
    return records, summary


def validate_task_contract(task_path: Path, gate_path: Path) -> dict[str, Any]:
    task = load_task_spec(task_path)
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    if (
        task.task_id != "paired_bulk_histology_state"
        or task.unit_of_inference != "participant"
        or task.primary_metric != "participant_macro_f1_stage3"
        or task.datasets_train != ("gse267145_znf469_human_liver",)
        or task.datasets_development
        or task.datasets_sealed
        or tuple(value.value for value in task.required_pairing_levels)
        != ("same_sample_different_aliquot",)
    ):
        raise HistologyActivationError("histology TaskSpec identity differs")
    if task.promotion_gate_config_sha256 != sha256_file(gate_path):
        raise HistologyActivationError("promotion-gate hash differs")
    if (
        gate.get("champion_eligible") is not False
        or gate.get("claim_mode") != "development_only"
        or gate.get("promotion_gate_id") != task.promotion_gate_id
    ):
        raise HistologyActivationError("single-cohort promotion gate is unsafe")
    parameters = dict(task.evaluator_parameters)
    if (
        tuple(parameters.get("primary_class_roster", ())) != STAGE3
        or tuple(parameters.get("source_stage_roster", ())) != STAGE5
        or parameters.get("participants") != 99
        or parameters.get("rna_only_query_h3") is not False
        or parameters.get("observed_pair_query_h3") is not True
    ):
        raise HistologyActivationError("histology evaluator parameters differ")
    text = "\n".join((task.endpoint, *task.admission_gates, task.claim_gate)).lower()
    required = (
        "differential regions",
        "outer-training",
        "fractional expression",
        "opaque source feature",
        "recorded sex",
        "development",
        "derivative weights",
    )
    if any(token not in text for token in required):
        raise HistologyActivationError("TaskSpec omits a leakage, modality, or claim guard")
    return {
        "task_id": task.task_id,
        "task_spec_sha256": sha256_file(task_path),
        "promotion_gate_id": task.promotion_gate_id,
        "promotion_gate_sha256": sha256_file(gate_path),
        "primary_metric": task.primary_metric,
        "champion_eligible": False,
        "external_claim_eligible": False,
    }


def activation_contract(task: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": "masld-bench-gse267145-histology-task-activation-v1",
        "status": "participant_histology_development_active",
        "cohort_family_id": "gse267145_znf469_human_liver",
        "task_id": task["task_id"],
        "participants": 99,
        "biological_unit": "participant",
        "pairing": "same_sample_different_aliquot",
        "task_specific_model_fitting_allowed": True,
        "dataset_wide_activation": False,
        "allowed_lanes": {
            "rna_only": {
                "query_modalities": ["bulk_rna"],
                "held_h3_access": "evaluator_only_not_model_input",
            },
            "h3_only": {
                "query_modalities": ["h3k27ac_cutrun"],
                "held_rna_access": "evaluator_only_not_model_input",
            },
            "observed_pair_multimodal": {
                "query_modalities": ["bulk_rna", "h3k27ac_cutrun"],
                "requires_both_observed_assays": True,
                "cannot_support_rna_only_claim": True,
            },
        },
        "endpoint_contract": {
            "primary": "stage3_NOR_NAFL_NASH",
            "secondary": [
                "nash_crn_component_sum",
                "fibrosis_0_to_3",
                "source_stage5",
            ],
            "stage5_is_uniformly_spaced_ordinal": False,
            "histology_allowed_as_supervised_target": True,
            "histology_allowed_in_unsupervised_preprocessing": False,
        },
        "metadata_contract": {
            "recorded_sex": "adjustment_and_error_audit_only_v1",
            "participant_age": "structurally_missing",
            "participant_bmi": "structurally_missing",
            "missing_encoded_as_zero": False,
        },
        "mandatory_baselines": [
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
        ],
        "eligible_development_families": [
            "assay_native_regularized_linear",
            "partial_least_squares",
            "mofa_plus",
            "generalized_low_rank_multiview",
            "regularized_paired_bulk_autoencoder_screening_only",
            "calibrated_late_fusion",
        ],
        "feature_contract": {
            "complete_deposited_assay_axes_only": True,
            "feature_selection_fit_inside_outer_training": True,
            "published_outcome_selected_features_allowed": False,
            "h3_region_identity": "opaque_source_feature_key",
            "sequence_or_coordinate_features_allowed": False,
        },
        "claim_contract": {
            "cross_sectional_histology_associated_remodeling_only": True,
            "champion_eligible": False,
            "external_claim_eligible": False,
            "clinical_claim_eligible": False,
        },
        "rights": {
            "internal_nonclinical_research": "allowed_without_new_DUA",
            "source_data_redistribution": "prohibited_by_project_policy",
            "derivative_weight_release": "requires_model_specific_terms_and_legal_review",
        },
        "blocked_actions": [
            "published_differential_region_feature_selection",
            "published_differential_gene_or_signature_feature_selection",
            "sequence_extraction_from_unresolved_h3_coordinates",
            "normalization_or_feature_selection_fit_on_held_participants",
            "histology_or_sex_as_unsupervised_preprocessing_input",
            "task_champion_or_external_transfer_claim",
            "diagnostic_prognostic_progression_or_treatment_response_claim",
            "raw_or_processed_source_data_redistribution",
            "derivative_weight_release_without_review",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--join", required=True, type=Path)
    parser.add_argument("--join-artifacts-sha256", required=True)
    parser.add_argument("--measurement", required=True, type=Path)
    parser.add_argument("--measurement-artifacts-sha256", required=True)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--reference-artifacts-sha256", required=True)
    parser.add_argument("--qc-rights", required=True, type=Path)
    parser.add_argument("--qc-rights-artifacts-sha256", required=True)
    parser.add_argument("--coordinate-activation", required=True, type=Path)
    parser.add_argument("--coordinate-activation-artifacts-sha256", required=True)
    parser.add_argument("--task-spec", required=True, type=Path)
    parser.add_argument("--task-spec-sha256", required=True)
    parser.add_argument("--promotion-gate", required=True, type=Path)
    parser.add_argument("--promotion-gate-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    for root, expected in (
        (arguments.join, arguments.join_artifacts_sha256),
        (arguments.measurement, arguments.measurement_artifacts_sha256),
        (arguments.reference, arguments.reference_artifacts_sha256),
        (arguments.qc_rights, arguments.qc_rights_artifacts_sha256),
        (
            arguments.coordinate_activation,
            arguments.coordinate_activation_artifacts_sha256,
        ),
    ):
        verify_artifact(root, expected)
    if (
        sha256_file(arguments.task_spec) != arguments.task_spec_sha256
        or sha256_file(arguments.promotion_gate)
        != arguments.promotion_gate_sha256
    ):
        raise HistologyActivationError("TaskSpec or promotion-gate file hash differs")
    coordinate = json.loads(
        (arguments.coordinate_activation / "activation_contract.json").read_text(
            encoding="utf-8"
        )
    )
    if (
        coordinate.get("coordinate_semantics_resolved") is not False
        or coordinate.get("sequence_extraction_allowed") is not False
        or coordinate.get("participants") != 99
    ):
        raise HistologyActivationError("upstream coordinate/task boundary differs")
    records, endpoint_summary = build_endpoint_records(
        read_tsv(arguments.join / "participant_join.tsv")
    )
    task_receipt = validate_task_contract(
        arguments.task_spec, arguments.promotion_gate
    )
    arguments.output.mkdir(parents=True, exist_ok=False)
    write_tsv(arguments.output / "participant_endpoints.tsv", records)
    for name, value in (
        ("endpoint_summary.json", endpoint_summary),
        ("task_spec_receipt.json", task_receipt),
        ("activation_contract.json", activation_contract(task_receipt)),
    ):
        (arguments.output / name).write_text(
            json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(json.dumps(endpoint_summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
