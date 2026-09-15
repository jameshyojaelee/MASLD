#!/usr/bin/env python3
"""Validate the content-free scBasset queue-search incident record."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "masld-bench-protocol-exposure-incident-v1"
INCIDENT_ID = "scbasset_broad_rg_prediction_line_exposure_20260825"
STATUS = "acknowledged_non_scientific_input_not_used"
SEARCH_PATTERNS = ("model-training-605", '"priority": 64', "064-")
SEARCH_ROOTS = ("config/campaigns/gpu_bundle_queue", "executions")


class ScBassetSearchIncidentError(ValueError):
    """Raised when the content-free incident requirement differs."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def validate_incident(path: Path) -> dict[str, Any]:
    """Validate structure and exclusions without opening referenced predictions."""

    path = path.resolve(strict=True)
    try:
        incident = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ScBassetSearchIncidentError("incident JSON differs") from error
    if not isinstance(incident, dict):
        raise ScBassetSearchIncidentError("incident JSON is not an object")

    scope = incident.get("search_scope", {})
    exposure = incident.get("exposure_scope", {})
    exclusions = incident.get("explicit_exclusions", {})
    impact = incident.get("impact_assessment", {})
    containment = incident.get("containment", {})
    expected_false = (
        "scbasset_prediction_values_emitted",
        "scbasset_metrics_emitted",
        "raw_outcomes_emitted",
        "sealed_data_emitted",
        "sealed_labels_or_outcomes_emitted",
        "incident_values_used_for_rectangle_design",
        "incident_values_used_for_queue_design",
        "incident_values_used_for_model_design",
        "incident_values_used_for_training",
        "incident_values_used_for_model_selection",
        "incident_values_used_for_calibration_or_thresholding",
    )
    if (
        incident.get("schema_version") != SCHEMA_VERSION
        or incident.get("incident_id") != INCIDENT_ID
        or incident.get("timestamp_date") != "2026-08-25"
        or incident.get("incident_type")
        != "overbroad_queue_identity_search_under_executions"
        or incident.get("status") != STATUS
        or tuple(scope.get("search_roots", ())) != SEARCH_ROOTS
        or tuple(scope.get("search_patterns", ())) != SEARCH_PATTERNS
        or tuple(scope.get("file_classes_requested", ())) != ("json", "tsv")
        or scope.get("exact_shell_command_text_retained") is not False
        or scope.get("scope_reconstructed_from_incident_time_observation") is not True
        or exposure.get("unrelated_prediction_tsv_lines_emitted") is not True
        or exposure.get("emitted_line_count_recounted") is not False
        or exposure.get("emitted_line_count") is not None
        or exposure.get("prediction_values_copied_into_incident") is not False
        or exposure.get("prediction_values_reinspected_after_incident_recognition") is not False
        or any(exclusions.get(field) is not False for field in expected_false)
        or impact.get("scbasset_admission_affected") is not False
        or impact.get("scbasset_five_seed_rectangle_affected") is not False
        or impact.get("scbasset_prediction_view_affected") is not False
        or impact.get("sealed_test_compromise") is not False
        or impact.get("scientific_input_from_incident") is not False
        or impact.get("disposition") != STATUS
        or not all(
            containment.get(field) is True
            for field in (
                "broad_search_stopped",
                "future_queue_searches_restricted_to_queue_manifests_and_exact_claim_paths",
                "no_further_prediction_value_inspection_authorized",
                "incident_disclosed_to_parent_agent",
                "incident_artifact_contains_no_prediction_or_outcome_value",
            )
        )
    ):
        raise ScBassetSearchIncidentError("incident structure or containment differs")

    return {
        "schema_version": "masld-bench-protocol-exposure-incident-validation-v1",
        "status": "pass_structure_only_incident_validation",
        "incident_id": INCIDENT_ID,
        "incident_sha256": _digest(path),
        "unrelated_prediction_tsv_lines_emitted": True,
        "scbasset_prediction_values_emitted": False,
        "scbasset_metrics_emitted": False,
        "raw_outcomes_emitted": False,
        "sealed_data_emitted": False,
        "incident_values_used": False,
        "referenced_prediction_files_reopened_by_validator": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--incident", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    receipt = validate_incident(arguments.incident)
    rendered = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
