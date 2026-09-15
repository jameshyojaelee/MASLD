from __future__ import annotations

import copy
import unittest

from scripts.verify_gse83452_oligo_firewall_receipt import (
    EXPECTED_PACKAGE_VERSIONS,
    REQUIRED_ESTIMATOR_CONTEXTS,
    REQUIRED_FALSE_FLAGS,
    REQUIRED_PIPELINE_CONTEXTS,
    FirewallError,
    verify_receipt,
)

CORE_FEATURES = 53_617


def target(accession: str, digest: str) -> dict[str, object]:
    return {
        "sample_accession": accession,
        "pipeline_arrays_per_summarization_call": 1,
        "pipeline_digests": {context: digest for context in REQUIRED_PIPELINE_CONTEXTS},
        "pipeline_differing_elements": {context: 0 for context in REQUIRED_PIPELINE_CONTEXTS},
        "pipeline_maximum_absolute_difference": {
            context: 0.0 for context in REQUIRED_PIPELINE_CONTEXTS
        },
        "pipeline_all_bitwise_identical": True,
        "estimator_digests": {
            context: f"{index}" * 64 for index, context in enumerate(REQUIRED_ESTIMATOR_CONTEXTS)
        },
        "estimator_differing_elements": {
            context: 51_000 for context in REQUIRED_ESTIMATOR_CONTEXTS
        },
        "estimator_maximum_absolute_difference": {
            context: 0.42 for context in REQUIRED_ESTIMATOR_CONTEXTS
        },
        "estimator_all_bitwise_identical": False,
        "perturbation_control_digest_changed": True,
        "perturbation_control_magnitude": 1e-12,
        "features": CORE_FEATURES,
        "finite_values": CORE_FEATURES,
    }


def receipt() -> dict[str, object]:
    payload: dict[str, object] = {
        "schema_version": "masld-bench-gse83452-oligo-single-array-firewall-v1",
        "series": "GSE83452",
        "platform_id": "GPL16686",
        "method": "single_array_core_summary_normalize_false",
        "held_array_application": "one_array_at_a_time_with_frozen_training_object",
        "core_features": CORE_FEATURES,
        "cel_records_available": 231,
        "package_versions": dict(EXPECTED_PACKAGE_VERSIONS),
        "scan_probe": {"method": "SCAN_with_pd.hugene.2.0.st"},
        "per_target": [target("GSM1", "a" * 64), target("GSM2", "b" * 64)],
    }
    for flag in REQUIRED_FALSE_FLAGS:
        payload[flag] = False
    return payload


class PipelineFirewallTests(unittest.TestCase):
    def test_clean_receipt_passes(self) -> None:
        verdict = verify_receipt(receipt())
        self.assertEqual(verdict["status"], "pass_inductive_firewall_on_real_CEL")
        self.assertEqual(verdict["platform_id"], "GPL16686")
        self.assertEqual(verdict["arrays_per_summarization_call"], 1)
        self.assertEqual(verdict["core_features"], CORE_FEATURES)
        self.assertIs(verdict["gpl570_evidence_carried_over"], False)

    def test_pipeline_divergence_is_rejected(self) -> None:
        payload = receipt()
        payload["per_target"][0]["pipeline_digests"]["solo_after_P"] = "c" * 64  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "what else the run processed"):
            verify_receipt(payload)

    def test_nonzero_pipeline_difference_is_rejected(self) -> None:
        payload = receipt()
        payload["per_target"][1]["pipeline_maximum_absolute_difference"]["solo_only"] = 1e-14  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "nonzero pipeline-context difference"):
            verify_receipt(payload)

    def test_multi_array_pipeline_call_is_rejected(self) -> None:
        payload = receipt()
        payload["per_target"][0]["pipeline_arrays_per_summarization_call"] = 4  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "more than one array per call"):
            verify_receipt(payload)


class HazardMeasurementTests(unittest.TestCase):
    def test_silent_negative_control_is_rejected(self) -> None:
        payload = receipt()
        for record in payload["per_target"]:  # type: ignore[union-attr]
            record["estimator_differing_elements"] = {
                context: 0 for context in REQUIRED_ESTIMATOR_CONTEXTS
            }
        with self.assertRaisesRegex(FirewallError, "hazard measurement is broken"):
            verify_receipt(payload)

    def test_multi_array_call_declared_as_preprocessing_path_is_rejected(self) -> None:
        payload = receipt()
        payload["multi_array_call_used_as_preprocessing_path"] = True
        with self.assertRaisesRegex(
            FirewallError, "multi_array_call_used_as_preprocessing_path"
        ):
            verify_receipt(payload)


class CrossPlatformIsolationTests(unittest.TestCase):
    """A GPL570 pass must never be usable as GPL16686 evidence."""

    def test_carried_over_gpl570_evidence_is_rejected(self) -> None:
        payload = receipt()
        payload["gpl570_evidence_carried_over"] = True
        with self.assertRaisesRegex(FirewallError, "gpl570_evidence_carried_over"):
            verify_receipt(payload)

    def test_a_gpl570_receipt_is_not_accepted_here(self) -> None:
        payload = receipt()
        payload["schema_version"] = "masld-bench-gse49541-frma-single-array-firewall-v2"
        with self.assertRaisesRegex(FirewallError, "schema differs"):
            verify_receipt(payload)

    def test_gpl570_packages_are_rejected(self) -> None:
        payload = receipt()
        payload["package_versions"] = {"affy": "1.88.0", "frma": "1.62.0"}
        with self.assertRaisesRegex(FirewallError, "package versions differ"):
            verify_receipt(payload)

    def test_frma_method_is_rejected_on_this_platform(self) -> None:
        payload = receipt()
        payload["method"] = "frma_with_exact_hgu133plus2frmavecs"
        with self.assertRaisesRegex(FirewallError, "contract does not admit"):
            verify_receipt(payload)


class AxisTests(unittest.TestCase):
    """The GPL16686 axis is recorded, not judged against the GEO row count."""

    def test_a_core_axis_unlike_the_geo_row_count_is_accepted(self) -> None:
        payload = receipt()
        self.assertNotEqual(CORE_FEATURES, 53_981)
        self.assertEqual(verify_receipt(payload)["core_features"], CORE_FEATURES)

    def test_a_degenerate_axis_is_rejected(self) -> None:
        payload = receipt()
        payload["core_features"] = 3
        with self.assertRaisesRegex(FirewallError, "absent or degenerate"):
            verify_receipt(payload)

    def test_per_target_axis_must_match_the_core_axis(self) -> None:
        payload = receipt()
        payload["per_target"][0]["features"] = CORE_FEATURES - 1  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "differs from the core axis"):
            verify_receipt(payload)

    def test_nonfinite_summary_is_rejected(self) -> None:
        payload = receipt()
        payload["per_target"][0]["finite_values"] = CORE_FEATURES - 5  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "nonfinite value"):
            verify_receipt(payload)


class ControlTests(unittest.TestCase):
    def test_undetected_perturbation_control_is_rejected(self) -> None:
        payload = receipt()
        payload["per_target"][1]["perturbation_control_digest_changed"] = False  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "perturbation control"):
            verify_receipt(payload)

    def test_identical_digests_across_targets_are_rejected(self) -> None:
        payload = receipt()
        payload["per_target"] = [target("GSM1", "a" * 64), target("GSM2", "a" * 64)]
        with self.assertRaisesRegex(FirewallError, "same summary digest"):
            verify_receipt(payload)

    def test_every_forbidden_flag_must_be_false(self) -> None:
        for flag in REQUIRED_FALSE_FLAGS:
            payload = copy.deepcopy(receipt())
            payload[flag] = True
            with self.assertRaisesRegex(FirewallError, flag):
                verify_receipt(payload)

    def test_wrong_cel_record_count_is_rejected(self) -> None:
        payload = receipt()
        payload["cel_records_available"] = 72
        with self.assertRaisesRegex(FirewallError, "CEL record count differs"):
            verify_receipt(payload)


if __name__ == "__main__":
    unittest.main()
