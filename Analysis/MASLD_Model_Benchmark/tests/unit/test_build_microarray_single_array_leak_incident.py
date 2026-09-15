from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.build_microarray_single_array_leak_incident import (
    IncidentError,
    build_incident,
    summarize_alternative_cost,
    summarize_failure,
    summarize_mechanism,
)

BATCH_CONTEXTS = ("alone", "batch_P_first", "batch_P_last", "batch_Q_first")


def failing_receipt(bitwise: bool = False) -> dict[str, object]:
    return {
        "every_target_bitwise_identical_across_contexts": bitwise,
        "per_target": [
            {
                "sample_accession": "GSM789110",
                "features": 54_675,
                "digests": {
                    "alone": "a" * 64,
                    "batch_P_first": "b" * 64,
                    "batch_P_last": "b" * 64,
                    "batch_Q_first": "c" * 64,
                },
                "differing_elements": {
                    "alone": 0,
                    "batch_P_first": 38_463,
                    "batch_P_last": 38_463,
                    "batch_Q_first": 39_273,
                },
                "maximum_absolute_difference": {
                    "alone": 0.0,
                    "batch_P_first": 8.99e-4,
                    "batch_P_last": 8.99e-4,
                    "batch_Q_first": 8.99e-4,
                },
                "perturbation_control_digest_changed": True,
            }
        ],
    }


def mechanism_diagnostic() -> dict[str, object]:
    def configuration(identifier, summarize, differing, maximum):
        return {
            "id": identifier,
            "status": "pass",
            "background": "rma",
            "normalize": "quantile",
            "summarize": summarize,
            "batch_P_differing": differing,
            "batch_Q_differing": differing,
            "batch_P_max_absolute": maximum,
            "batch_Q_max_absolute": maximum,
        }

    return {
        "raw_intensities_identical_across_contexts": True,
        "configurations": [
            configuration("default", "robust_weighted_average", 38_463, 8.99e-4),
            configuration("avg", "average", 0, 0.0),
            configuration("med", "median", 0, 0.0),
            configuration("wavg", "weighted_average", 0, 0.0),
            {"id": "broken", "status": "failed", "error": "nope"},
        ],
    }


def agreement_diagnostic() -> dict[str, object]:
    return {
        "per_array": [
            {
                "sample_accession": "GSM789110",
                "reference_rule": "robust_weighted_average",
                "reference_iqr": 2.7068,
                "comparisons": [
                    {
                        "rule": "weighted_average",
                        "median_absolute_difference": 0.29186,
                        "max_absolute_difference": 3.35141,
                        "features_differing_by_more_than_0_10": 44_446,
                        "spearman_rho": 0.977437,
                    }
                ],
            }
        ]
    }


def passing_verdict() -> dict[str, object]:
    return {
        "status": "pass_inductive_firewall_on_real_CEL",
        "arrays_per_frma_call": 1,
        "multi_array_frma_call_is_batch_dependent": True,
    }


class SummarizeFailureTests(unittest.TestCase):
    def test_leak_magnitudes_are_rederived(self) -> None:
        summary = summarize_failure(failing_receipt())
        self.assertEqual(summary["maximum_differing_elements"], 39_273)
        self.assertAlmostEqual(summary["maximum_absolute_difference_log2"], 8.99e-4)

    def test_batch_position_and_membership_are_separated(self) -> None:
        summary = summarize_failure(failing_receipt())
        self.assertIs(summary["batch_position_matters"], False)
        self.assertIs(summary["batch_membership_matters"], True)

    def test_a_passing_receipt_is_not_an_incident(self) -> None:
        with self.assertRaisesRegex(IncidentError, "does not show a failure"):
            summarize_failure(failing_receipt(bitwise=True))

    def test_unfired_perturbation_control_invalidates_the_evidence(self) -> None:
        payload = failing_receipt()
        payload["per_target"][0]["perturbation_control_digest_changed"] = False  # type: ignore[index]
        with self.assertRaisesRegex(IncidentError, "perturbation control did not fire"):
            summarize_failure(payload)


class SummarizeMechanismTests(unittest.TestCase):
    def test_dependent_and_invariant_rules_are_separated(self) -> None:
        mechanism = summarize_mechanism(mechanism_diagnostic())
        self.assertEqual(
            mechanism["batch_dependent_summarize_rules"], ["robust_weighted_average"]
        )
        self.assertEqual(
            mechanism["batch_invariant_summarize_rules"],
            ["average", "median", "weighted_average"],
        )

    def test_differing_raw_intensities_break_the_attribution(self) -> None:
        payload = mechanism_diagnostic()
        payload["raw_intensities_identical_across_contexts"] = False
        with self.assertRaisesRegex(IncidentError, "not isolated to fRMA"):
            summarize_mechanism(payload)

    def test_a_rule_cannot_be_both_invariant_and_dependent(self) -> None:
        payload = mechanism_diagnostic()
        payload["configurations"].append(  # type: ignore[union-attr]
            {
                "id": "contradiction",
                "status": "pass",
                "background": "none",
                "normalize": "none",
                "summarize": "average",
                "batch_P_differing": 10,
                "batch_Q_differing": 10,
                "batch_P_max_absolute": 1.0,
                "batch_Q_max_absolute": 1.0,
            }
        )
        with self.assertRaisesRegex(IncidentError, "both invariant and dependent"):
            summarize_mechanism(payload)

    def test_a_diagnostic_with_no_contrast_is_rejected(self) -> None:
        payload = mechanism_diagnostic()
        payload["configurations"] = [payload["configurations"][0]]  # type: ignore[index]
        with self.assertRaisesRegex(IncidentError, "did not separate"):
            summarize_mechanism(payload)


class AlternativeCostTests(unittest.TestCase):
    def test_switching_cost_is_rederived_against_the_iqr(self) -> None:
        cost = summarize_alternative_cost(agreement_diagnostic())
        self.assertAlmostEqual(cost["weighted_average_median_absolute_shift_log2"], 0.29186)
        self.assertAlmostEqual(
            cost["weighted_average_median_shift_as_fraction_of_iqr"], 0.29186 / 2.7068
        )
        self.assertEqual(cost["weighted_average_features_moving_more_than_0_10"], 44_446)

    def test_missing_weighted_average_comparison_is_rejected(self) -> None:
        payload = agreement_diagnostic()
        payload["per_array"][0]["comparisons"][0]["rule"] = "median"  # type: ignore[index]
        with self.assertRaisesRegex(IncidentError, "did not compare weighted_average"):
            summarize_alternative_cost(payload)


class BuildIncidentTests(unittest.TestCase):
    def _write(self, base: Path) -> dict[str, Path]:
        paths = {}
        for name, payload in (
            ("failing_receipt", failing_receipt()),
            ("mechanism_diagnostic", mechanism_diagnostic()),
            ("agreement_diagnostic", agreement_diagnostic()),
            ("passing_verdict", passing_verdict()),
        ):
            path = base / f"{name}.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            paths[name] = path
        return paths

    def test_incident_binds_evidence_and_records_the_cost_ratio(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            paths = self._write(base)
            incident = build_incident(output=base / "incident", **paths)
            self.assertEqual(
                incident["status"], "resolved_by_holding_the_frma_batch_at_one_array"
            )
            self.assertEqual(
                incident["resolution"]["summarize_rule_retained"],  # type: ignore[index]
                "robust_weighted_average",
            )
            # The alternative rule costs far more than the leak it would remove.
            self.assertGreater(incident["alternative_rule_cost_relative_to_leak"], 100)
            self.assertEqual(len(incident["bound_evidence"]), 4)  # type: ignore[arg-type]
            for record in incident["bound_evidence"].values():  # type: ignore[union-attr]
                self.assertEqual(len(record["sha256"]), 64)
            self.assertIn(
                "GPL16686/GSE83452 oligo core summary with normalize=FALSE",
                incident["generalization"]["must_be_tested_before_use"],  # type: ignore[index]
            )

    def test_a_resolution_that_still_batches_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            paths = self._write(base)
            paths["passing_verdict"].write_text(
                json.dumps(passing_verdict() | {"arrays_per_frma_call": 6}), encoding="utf-8"
            )
            with self.assertRaisesRegex(IncidentError, "hold the batch at one array"):
                build_incident(output=base / "incident", **paths)

    def test_a_resolution_that_stopped_measuring_the_hazard_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            paths = self._write(base)
            paths["passing_verdict"].write_text(
                json.dumps(
                    passing_verdict() | {"multi_array_frma_call_is_batch_dependent": False}
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(IncidentError, "stopped recording the hazard"):
                build_incident(output=base / "incident", **paths)

    def test_a_failing_resolution_verdict_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            paths = self._write(base)
            paths["passing_verdict"].write_text(
                json.dumps(passing_verdict() | {"status": "fail"}), encoding="utf-8"
            )
            with self.assertRaisesRegex(IncidentError, "not a pass"):
                build_incident(output=base / "incident", **paths)


if __name__ == "__main__":
    unittest.main()
