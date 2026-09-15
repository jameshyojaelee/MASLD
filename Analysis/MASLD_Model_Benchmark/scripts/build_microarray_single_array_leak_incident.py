#!/usr/bin/env python3
"""Freeze the incident record for the fRMA multi-array summarization leak.

model-check-085 found that fRMA's default ``summarize = "robust_weighted_average"``
is not invariant to which arrays share its call.  The frozen requirements had
classified ``frma_with_exact_hgu133plus2frmavecs`` as an
``exact_frozen_single_array_method``; that classification is true only when the
call carries exactly one array.

This output file exists because the same hazard class is not specific to fRMA or to
GPL570.  Any summarizer that fits a robust or shared model across the columns it
is handed will borrow across held arrays, however frozen its reference vectors
are.  ``oligo`` RMA and SCAN on GPL16686/GSE83452 sit in that class and must be
measured the same way before that lane is trusted.

The record binds the failing evidence, the stage that isolated the mechanism,
and the measured cost of the alternatives, so the conclusion can be re-derived
rather than taken on trust.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


class IncidentError(RuntimeError):
    """Raised when the incident evidence is absent or does not say what is claimed."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _scalar(value: Any) -> Any:
    if isinstance(value, list) and len(value) == 1:
        return value[0]
    return value


def summarize_failure(receipt: dict[str, Any]) -> dict[str, Any]:
    """Re-derive the leak magnitudes from the failing v1 receipt."""
    if _scalar(receipt.get("every_target_bitwise_identical_across_contexts")) is not False:
        raise IncidentError("the recorded incident receipt does not show a failure")
    per_target = receipt.get("per_target") or []
    if not per_target:
        raise IncidentError("the recorded incident receipt has no per-target evidence")
    rows = []
    for record in per_target:
        differing = {
            key: int(_scalar(value))
            for key, value in (record.get("differing_elements") or {}).items()
        }
        maximum = {
            key: float(_scalar(value))
            for key, value in (record.get("maximum_absolute_difference") or {}).items()
        }
        digests = {key: _scalar(value) for key, value in (record.get("digests") or {}).items()}
        if _scalar(record.get("perturbation_control_digest_changed")) is not True:
            raise IncidentError("the failing run's perturbation control did not fire")
        rows.append(
            {
                "sample_accession": _scalar(record.get("sample_accession")),
                "features": _scalar(record.get("features")),
                "differing_elements": differing,
                "maximum_absolute_difference": maximum,
                "batch_position_changed_the_result": digests.get("batch_P_first")
                != digests.get("batch_P_last"),
                "batch_membership_changed_the_result": digests.get("batch_P_first")
                != digests.get("batch_Q_first"),
            }
        )
    return {
        "targets": rows,
        "maximum_differing_elements": max(
            value for row in rows for value in row["differing_elements"].values()
        ),
        "maximum_absolute_difference_log2": max(
            value for row in rows for value in row["maximum_absolute_difference"].values()
        ),
        "batch_position_matters": any(row["batch_position_changed_the_result"] for row in rows),
        "batch_membership_matters": any(
            row["batch_membership_changed_the_result"] for row in rows
        ),
    }


def summarize_mechanism(diagnostic: dict[str, Any]) -> dict[str, Any]:
    """Re-derive which fRMA argument carries the dependence."""
    if _scalar(diagnostic.get("raw_intensities_identical_across_contexts")) is not True:
        raise IncidentError("raw intensities differed; the mechanism is not isolated to fRMA")
    invariant, dependent = [], []
    for configuration in diagnostic.get("configurations") or []:
        if _scalar(configuration.get("status")) != "pass":
            continue
        identifier = _scalar(configuration.get("id"))
        entry = {
            "id": identifier,
            "background": _scalar(configuration.get("background")),
            "normalize": _scalar(configuration.get("normalize")),
            "summarize": _scalar(configuration.get("summarize")),
            "batch_P_differing": int(_scalar(configuration.get("batch_P_differing"))),
            "batch_Q_differing": int(_scalar(configuration.get("batch_Q_differing"))),
            "batch_P_max_absolute": float(_scalar(configuration.get("batch_P_max_absolute"))),
            "batch_Q_max_absolute": float(_scalar(configuration.get("batch_Q_max_absolute"))),
        }
        if entry["batch_P_differing"] == 0 and entry["batch_Q_differing"] == 0:
            invariant.append(entry)
        else:
            dependent.append(entry)
    if not invariant or not dependent:
        raise IncidentError("the diagnostic did not separate invariant from dependent settings")
    dependent_rules = sorted({entry["summarize"] for entry in dependent})
    invariant_rules = sorted({entry["summarize"] for entry in invariant})
    if set(dependent_rules) & set(invariant_rules):
        raise IncidentError("a summarize rule appears both invariant and dependent")
    return {
        "raw_intensities_identical_across_contexts": True,
        "batch_dependent_summarize_rules": dependent_rules,
        "batch_invariant_summarize_rules": invariant_rules,
        "background_setting_changes_the_dependence": False,
        "normalize_setting_changes_the_dependence": False,
        "configurations_batch_dependent": dependent,
        "configurations_batch_invariant": invariant,
        "implementation": (
            "frmaAffyBatch calls rwaFit2(pms[s, , drop = FALSE], ...) on the whole "
            "probe-by-array submatrix, so the robustness weights are estimated over "
            "every array present in the call"
        ),
    }


def summarize_alternative_cost(agreement: dict[str, Any]) -> dict[str, Any]:
    """Re-derive what switching to a batch-invariant rule would cost."""
    rows = []
    for record in agreement.get("per_array") or []:
        for comparison in record.get("comparisons") or []:
            rows.append(
                {
                    "sample_accession": _scalar(record.get("sample_accession")),
                    "reference_rule": _scalar(record.get("reference_rule")),
                    "rule": _scalar(comparison.get("rule")),
                    "reference_iqr": float(_scalar(record.get("reference_iqr"))),
                    "median_absolute_difference": float(
                        _scalar(comparison.get("median_absolute_difference"))
                    ),
                    "max_absolute_difference": float(
                        _scalar(comparison.get("max_absolute_difference"))
                    ),
                    "features_differing_by_more_than_0_10": int(
                        _scalar(comparison.get("features_differing_by_more_than_0_10"))
                    ),
                    "spearman_rho": float(_scalar(comparison.get("spearman_rho"))),
                }
            )
    if not rows:
        raise IncidentError("the agreement diagnostic has no comparisons")
    weighted = [row for row in rows if row["rule"] == "weighted_average"]
    if not weighted:
        raise IncidentError("the agreement diagnostic did not compare weighted_average")
    median_shift = max(row["median_absolute_difference"] for row in weighted)
    iqr = min(row["reference_iqr"] for row in weighted)
    return {
        "comparisons": rows,
        "weighted_average_median_absolute_shift_log2": median_shift,
        "weighted_average_median_shift_as_fraction_of_iqr": median_shift / iqr,
        "weighted_average_features_moving_more_than_0_10": max(
            row["features_differing_by_more_than_0_10"] for row in weighted
        ),
    }


def build_incident(
    *,
    failing_receipt: Path,
    mechanism_diagnostic: Path,
    agreement_diagnostic: Path,
    passing_verdict: Path,
    output: Path,
) -> dict[str, Any]:
    failure = summarize_failure(json.loads(failing_receipt.read_text(encoding="utf-8")))
    mechanism = summarize_mechanism(
        json.loads(mechanism_diagnostic.read_text(encoding="utf-8"))
    )
    cost = summarize_alternative_cost(
        json.loads(agreement_diagnostic.read_text(encoding="utf-8"))
    )
    resolution = json.loads(passing_verdict.read_text(encoding="utf-8"))
    if resolution.get("status") != "pass_inductive_firewall_on_real_CEL":
        raise IncidentError("the resolution verdict is not a pass")
    if resolution.get("arrays_per_frma_call") != 1:
        raise IncidentError("the resolution verdict does not hold the batch at one array")
    if resolution.get("multi_array_frma_call_is_batch_dependent") is not True:
        raise IncidentError("the resolution verdict stopped recording the hazard")

    # The leak is small; the alternative rule is not. Both numbers must survive
    # into the record, because the second is why the first was not "fixed" by
    # swapping the estimator.
    ratio = (
        cost["weighted_average_median_absolute_shift_log2"]
        / failure["maximum_absolute_difference_log2"]
    )

    incident = {
        "schema_version": "masld-bench-microarray-single-array-leak-incident-v1",
        "status": "resolved_by_holding_the_frma_batch_at_one_array",
        "incident_id": "frma_multi_array_summarization_leak_20260825",
        "series": "GSE49541",
        "platform_id": "GPL570",
        "cohort_family_id": "gse31803_gse49541_fibrosis_array",
        "method": "frma_with_exact_hgu133plus2frmavecs",
        "finding": (
            "fRMA's default summarize = robust_weighted_average is not invariant to "
            "which arrays share its call. The frozen contract classified this method "
            "as an exact_frozen_single_array_method; that holds only at one array per "
            "call."
        ),
        "detected_by": "gse49541_gpl570_single_array_inductive_firewall",
        "detection_would_have_been_missed_by": (
            "a capability probe that reads one array per platform, which is what the "
            "runtime audit ran and why the hazard survived until a batch was formed"
        ),
        "failure_evidence": failure,
        "mechanism": mechanism,
        "alternative_rule_cost": cost,
        "alternative_rule_cost_relative_to_leak": ratio,
        "resolution": {
            "action": "hold the frma() batch at one array by construction",
            "rationale": (
                "With one column rwaFit2 has nothing to borrow from, so the leak is "
                "structurally impossible rather than merely small. Switching to a "
                "batch-invariant summarize rule would instead change the measurement "
                "itself, by far more than the leak it removes."
            ),
            "summarize_rule_retained": "robust_weighted_average",
            "guarantee_rests_on": "pipeline discipline, machine-checked every run",
            "enforced_by": [
                "the firewall fixture's pipeline contexts, which must be bit-identical",
                "the firewall fixture's multi-array negative controls, which must keep diverging",
                "the summarization receipt's per-array arrays_read_in_this_call, which must be 1",
            ],
            "reproducibility_caveat": (
                "a standard batched fRMA rerun of this cohort will not reproduce these "
                "values exactly; disclose in the model card"
            ),
            "resolution_verdict": resolution,
        },
        "generalization": {
            "hazard_class": (
                "any summarizer that fits a robust or otherwise shared model across the "
                "columns it is handed borrows across held arrays, however frozen its "
                "reference vectors are"
            ),
            "frozen_reference_vectors_are_not_sufficient": True,
            "single_array_capability_probe_is_not_sufficient": True,
            "must_be_tested_before_use": [
                "GPL16686/GSE83452 oligo core summary with normalize=FALSE",
                "GPL16686/GSE83452 SCAN",
                "GPL570 SCAN",
            ],
            "required_test": (
                "summarize the same array alone and alongside others, compare the raw "
                "IEEE-754 doubles, and carry a perturbation control so a pass cannot be "
                "a broken comparator"
            ),
        },
        "contract_language_now_falsified": [
            "allowed_primary_paths.exact_frozen_single_array_method_after_platform_fixture",
            "platform.GPL570.candidate_single_array_methods.frma_with_exact_hgu133plus2frmavecs",
        ],
        "bound_evidence": {
            "failing_receipt": {
                "path": str(failing_receipt),
                "sha256": sha256_file(failing_receipt),
            },
            "mechanism_diagnostic": {
                "path": str(mechanism_diagnostic),
                "sha256": sha256_file(mechanism_diagnostic),
            },
            "agreement_diagnostic": {
                "path": str(agreement_diagnostic),
                "sha256": sha256_file(agreement_diagnostic),
            },
            "passing_verdict": {
                "path": str(passing_verdict),
                "sha256": sha256_file(passing_verdict),
            },
        },
        "labels_read": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
    }

    output.mkdir(parents=True, exist_ok=False)
    (output / "single_array_leak_incident.json").write_text(
        json.dumps(incident, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return incident


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--failing-receipt", type=Path, required=True)
    parser.add_argument("--mechanism-diagnostic", type=Path, required=True)
    parser.add_argument("--agreement-diagnostic", type=Path, required=True)
    parser.add_argument("--passing-verdict", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    incident = build_incident(
        failing_receipt=arguments.failing_receipt,
        mechanism_diagnostic=arguments.mechanism_diagnostic,
        agreement_diagnostic=arguments.agreement_diagnostic,
        passing_verdict=arguments.passing_verdict,
        output=arguments.output,
    )
    print(json.dumps({key: incident[key] for key in ("status", "incident_id", "finding")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
