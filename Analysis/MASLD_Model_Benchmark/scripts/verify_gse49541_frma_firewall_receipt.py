#!/usr/bin/env python3
"""Adjudicate the GSE49541 single-array fRMA separation receipt, fail-closed.

The fixture compares each target array across two families of context, and this
module treats them differently on purpose.

Pipeline contexts are the separation.  Every one of them reads exactly one CEL per
``frma()`` call, which is what ``microarray_transfer_preprocessing.toml`` means
by ``held_array_application = "one_array_at_a_time_with_frozen_training_object"``.
They differ only in what else the same run processed and in what order.  All of
them must be bit-identical, with zero differing elements and a maximum absolute
difference of exactly zero.

Estimator contexts are multi-array ``frma()`` calls, retained as negative
controls.  They are *expected* to diverge, because fRMA's default
``robust_weighted_average`` re-estimates its robustness weights over every array
in the call.  This module requires that divergence to still be measured and
recorded, so that batching cannot be quietly reintroduced without the fixture
noticing, and it refuses any receipt claiming a multi-array call was used as a
preprocessing path.

Three further conditions must hold together, and none is sufficient alone:

* the perturbation control changed the digest, so bit-identity is a real
  measurement rather than a degenerate comparison;
* distinct targets carry distinct digests, so per-array independence is not
  every array collapsing onto one vector;
* no label, series matrix, all-sample RMA, or across-array normalization appears
  anywhere in the path.

A receipt that merely says ``pass`` is never trusted; the underlying per-target
counts are re-derived here.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

REQUIRED_PIPELINE_CONTEXTS = (
    "solo_only",
    "solo_after_P",
    "solo_before_P",
    "solo_interleaved_Q",
)
REQUIRED_ESTIMATOR_CONTEXTS = ("batch_P_first", "batch_P_last", "batch_Q_first")
REQUIRED_FALSE_FLAGS = (
    "all_sample_RMA_run",
    "across_array_quantile_normalization_run",
    "GEO_series_matrix_read",
    "expression_values_exported",
    "labels_read",
    "model_training_activated",
    "sealed_outcomes_read",
    "selection_used_labels",
    "selection_used_intensity",
    "multi_array_frma_call_used_as_preprocessing_path",
)
EXPECTED_PACKAGE_VERSIONS = {
    "affy": "1.88.0",
    "affyio": "1.80.0",
    "frma": "1.62.0",
    "hgu133plus2cdf": "2.18.0",
    "hgu133plus2frmavecs": "1.5.0",
}


class FirewallError(RuntimeError):
    """Raised when the single-array isolation evidence does not hold."""


def _scalar(value: Any) -> Any:
    """Unwrap the length-one lists that R's JSON writer emits."""
    if isinstance(value, list) and len(value) == 1:
        return value[0]
    return value


def _mapping(record: Mapping[str, Any], key: str, contexts: tuple[str, ...], label: str) -> dict[str, Any]:
    values = {name: _scalar(value) for name, value in (record.get(key) or {}).items()}
    if tuple(sorted(values)) != tuple(sorted(contexts)):
        raise FirewallError(f"{label} lacks a value for every required {key}")
    return values


def verify_receipt(receipt: Mapping[str, Any]) -> dict[str, Any]:
    if receipt.get("schema_version") != "masld-bench-gse49541-frma-single-array-firewall-v2":
        raise FirewallError("firewall receipt schema differs")
    if receipt.get("series") != "GSE49541" or receipt.get("platform_id") != "GPL570":
        raise FirewallError("firewall receipt names a different source")
    if receipt.get("method") != "frma_with_exact_hgu133plus2frmavecs":
        raise FirewallError("firewall receipt names a method the contract does not admit")
    if receipt.get("held_array_application") != "one_array_at_a_time_with_frozen_training_object":
        raise FirewallError("firewall receipt does not declare one-array-at-a-time application")
    if _scalar(receipt.get("platform_features")) != 54_675:
        raise FirewallError("firewall receipt platform axis differs")
    if _scalar(receipt.get("cel_records_available")) != 72:
        raise FirewallError("firewall receipt CEL record count differs")

    versions = {key: _scalar(value) for key, value in (receipt.get("package_versions") or {}).items()}
    if versions != EXPECTED_PACKAGE_VERSIONS:
        raise FirewallError(f"summarization package versions differ: {versions!r}")

    for flag in REQUIRED_FALSE_FLAGS:
        if _scalar(receipt.get(flag)) is not False:
            raise FirewallError(f"firewall receipt does not clear {flag}")

    per_target = receipt.get("per_target") or []
    if len(per_target) < 2:
        raise FirewallError("firewall evidence needs at least two target arrays")

    alone_digests: list[str] = []
    summaries: list[dict[str, Any]] = []
    estimator_divergence_seen = False
    for record in per_target:
        accession = _scalar(record.get("sample_accession"))

        if _scalar(record.get("pipeline_arrays_per_frma_call")) != 1:
            raise FirewallError(f"{accession} pipeline context read more than one array per call")

        digests = _mapping(record, "pipeline_digests", REQUIRED_PIPELINE_CONTEXTS, accession)
        if len(set(digests.values())) != 1:
            raise FirewallError(f"{accession} summary depends on what else the run processed")

        differing = _mapping(
            record, "pipeline_differing_elements", REQUIRED_PIPELINE_CONTEXTS, accession
        )
        maximum = {
            key: float(value)
            for key, value in _mapping(
                record,
                "pipeline_maximum_absolute_difference",
                REQUIRED_PIPELINE_CONTEXTS,
                accession,
            ).items()
        }
        if any(int(value) != 0 for value in differing.values()):
            raise FirewallError(f"{accession} carries a differing feature across pipeline contexts")
        if any(value != 0.0 for value in maximum.values()):
            raise FirewallError(f"{accession} carries a nonzero pipeline-context difference")
        if _scalar(record.get("pipeline_all_bitwise_identical")) is not True:
            raise FirewallError(f"{accession} is not bitwise identical across pipeline contexts")

        # The multi-array hazard must still be measured every run.
        estimator_differing = _mapping(
            record, "estimator_differing_elements", REQUIRED_ESTIMATOR_CONTEXTS, accession
        )
        estimator_maximum = {
            key: float(value)
            for key, value in _mapping(
                record,
                "estimator_maximum_absolute_difference",
                REQUIRED_ESTIMATOR_CONTEXTS,
                accession,
            ).items()
        }
        if any(int(value) != 0 for value in estimator_differing.values()):
            estimator_divergence_seen = True

        if _scalar(record.get("perturbation_control_digest_changed")) is not True:
            raise FirewallError(f"{accession} perturbation control did not change the digest")
        if _scalar(record.get("features")) != 54_675:
            raise FirewallError(f"{accession} summary axis differs")
        if _scalar(record.get("finite_values")) != 54_675:
            raise FirewallError(f"{accession} summary carries a nonfinite value")

        alone_digests.append(digests["solo_only"])
        summaries.append(
            {
                "sample_accession": accession,
                "pipeline_contexts": len(digests),
                "digest": digests["solo_only"],
                "pipeline_maximum_absolute_difference": max(maximum.values()),
                "estimator_maximum_absolute_difference": max(estimator_maximum.values()),
                "estimator_maximum_differing_elements": max(
                    int(value) for value in estimator_differing.values()
                ),
            }
        )

    if len(set(alone_digests)) != len(alone_digests):
        raise FirewallError("two target arrays produced the same summary digest")
    if not estimator_divergence_seen:
        raise FirewallError(
            "multi-array negative control no longer diverges; the hazard measurement is broken"
        )

    return {
        "schema_version": "masld-bench-gse49541-frma-firewall-verdict-v2",
        "status": "pass_inductive_firewall_on_real_CEL",
        "series": "GSE49541",
        "platform_id": "GPL570",
        "method": "frma_with_exact_hgu133plus2frmavecs",
        "held_array_application": "one_array_at_a_time_with_frozen_training_object",
        "targets_verified": len(summaries),
        "pipeline_contexts_per_target": len(REQUIRED_PIPELINE_CONTEXTS),
        "pipeline_contexts": list(REQUIRED_PIPELINE_CONTEXTS),
        "estimator_contexts": list(REQUIRED_ESTIMATOR_CONTEXTS),
        "per_target": summaries,
        "held_array_summary_invariant_to_other_arrays_in_the_run": True,
        "held_array_summary_invariant_to_processing_order": True,
        "arrays_per_frma_call": 1,
        "perturbation_control_detected_for_every_target": True,
        "distinct_targets_have_distinct_digests": True,
        "perturbation_control_magnitude": _scalar(
            per_target[0].get("perturbation_control_magnitude")
        ),
        "multi_array_frma_call_is_batch_dependent": True,
        "multi_array_frma_call_used_as_preprocessing_path": False,
        "batching_would_break_the_firewall": True,
        "random_effect_probe": receipt.get("random_effect_probe"),
        "frma_vector_object_sha256": _scalar(receipt.get("frma_vector_object_sha256")),
        "package_versions": versions,
        "all_sample_RMA_run": False,
        "across_array_quantile_normalization_run": False,
        "GEO_series_matrix_read": False,
        "labels_read": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    receipt = json.loads(arguments.receipt.read_text(encoding="utf-8"))
    verdict = verify_receipt(receipt)
    arguments.output.write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(verdict, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
