#!/usr/bin/env python3
"""Adjudicate the GSE83452 GPL16686 single-array separation receipt, fail-closed.

Deliberately a separate module from the GPL570 adjudicator.  The GPL570 pass is
evidence about ``affy``/``frma``/``rwaFit2`` and carries no weight for ``oligo``
or ``SCAN.UPC``, which are different packages and different code paths.  Sharing
an implementation would risk smuggling a GPL570 assumption into a GPL16686
verdict, so the standards are restated here rather than imported.

The structure mirrors GPL570 exactly, and for the same reasons:

* pipeline contexts read one CEL per summarization call and must be
  bit-identical -- that is the separation;
* multi-array estimator contexts are negative controls that must keep diverging,
  so batching cannot be reintroduced silently;
* the 1e-12 perturbation control must fire, so bit-identity is a measurement
  rather than a degenerate comparison;
* distinct targets must carry distinct digests, so per-array independence is not
  per-array collapse.

Unlike GPL570 the feature axis is not asserted against a frozen platform row
count.  The GEO GPL16686 table carries 53981 rows while the ``pd.hugene.2.0.st``
core target returns a different number of transcript clusters, and the
authoritative vendor mapping is still unavailable.  The axis is therefore
required only to be self-consistent and non-degenerate, and its size is
recorded, not judged.
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
    "multi_array_call_used_as_preprocessing_path",
    "gpl570_evidence_carried_over",
)
EXPECTED_PACKAGE_VERSIONS = {
    "oligo": "1.74.0",
    "affyio": "1.80.0",
    "affxparser": "1.82.0",
    "pd.hugene.2.0.st": "3.14.1",
    "SCAN.UPC": "2.52.0",
}
MINIMUM_CORE_FEATURES = 1_000


class FirewallError(RuntimeError):
    """Raised when the single-array isolation evidence does not hold."""


def _scalar(value: Any) -> Any:
    if isinstance(value, list) and len(value) == 1:
        return value[0]
    return value


def _mapping(record: Mapping[str, Any], key: str, contexts: tuple[str, ...], label: str) -> dict[str, Any]:
    values = {name: _scalar(value) for name, value in (record.get(key) or {}).items()}
    if tuple(sorted(values)) != tuple(sorted(contexts)):
        raise FirewallError(f"{label} lacks a value for every required {key}")
    return values


def verify_receipt(receipt: Mapping[str, Any]) -> dict[str, Any]:
    if receipt.get("schema_version") != "masld-bench-gse83452-oligo-single-array-firewall-v1":
        raise FirewallError("firewall receipt schema differs")
    if receipt.get("series") != "GSE83452" or receipt.get("platform_id") != "GPL16686":
        raise FirewallError("firewall receipt names a different source")
    if receipt.get("method") != "single_array_core_summary_normalize_false":
        raise FirewallError("firewall receipt names a method the contract does not admit")
    if receipt.get("held_array_application") != "one_array_at_a_time_with_frozen_training_object":
        raise FirewallError("firewall receipt does not declare one-array-at-a-time application")
    if _scalar(receipt.get("cel_records_available")) != 231:
        raise FirewallError("firewall receipt CEL record count differs")

    core_features = _scalar(receipt.get("core_features"))
    if not isinstance(core_features, int) or core_features < MINIMUM_CORE_FEATURES:
        raise FirewallError("core feature axis is absent or degenerate")

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

        if _scalar(record.get("pipeline_arrays_per_summarization_call")) != 1:
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
        if _scalar(record.get("features")) != core_features:
            raise FirewallError(f"{accession} summary axis differs from the core axis")
        if _scalar(record.get("finite_values")) != core_features:
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
        "schema_version": "masld-bench-gse83452-oligo-firewall-verdict-v1",
        "status": "pass_inductive_firewall_on_real_CEL",
        "series": "GSE83452",
        "platform_id": "GPL16686",
        "cohort_family_id": "antwerp_inserm_shared",
        "method": "single_array_core_summary_normalize_false",
        "held_array_application": "one_array_at_a_time_with_frozen_training_object",
        "core_features": core_features,
        "targets_verified": len(summaries),
        "pipeline_contexts_per_target": len(REQUIRED_PIPELINE_CONTEXTS),
        "pipeline_contexts": list(REQUIRED_PIPELINE_CONTEXTS),
        "estimator_contexts": list(REQUIRED_ESTIMATOR_CONTEXTS),
        "per_target": summaries,
        "held_array_summary_invariant_to_other_arrays_in_the_run": True,
        "held_array_summary_invariant_to_processing_order": True,
        "arrays_per_summarization_call": 1,
        "perturbation_control_detected_for_every_target": True,
        "distinct_targets_have_distinct_digests": True,
        "perturbation_control_magnitude": _scalar(
            per_target[0].get("perturbation_control_magnitude")
        ),
        "multi_array_call_is_batch_dependent": True,
        "multi_array_call_used_as_preprocessing_path": False,
        "batching_would_break_the_firewall": True,
        "gpl570_evidence_carried_over": False,
        "scan_probe": receipt.get("scan_probe"),
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
    verdict = verify_receipt(json.loads(arguments.receipt.read_text(encoding="utf-8")))
    arguments.output.write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(verdict, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
