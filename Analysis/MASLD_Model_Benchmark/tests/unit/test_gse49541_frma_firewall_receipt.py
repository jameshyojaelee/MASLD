from __future__ import annotations

import copy
import unittest

from scripts.verify_gse49541_frma_firewall_receipt import (
    EXPECTED_PACKAGE_VERSIONS,
    REQUIRED_ESTIMATOR_CONTEXTS,
    REQUIRED_FALSE_FLAGS,
    REQUIRED_PIPELINE_CONTEXTS,
    FirewallError,
    verify_receipt,
)


def target(accession: str, digest: str) -> dict[str, object]:
    """A target whose pipeline contexts agree and whose batch contexts do not."""
    return {
        "sample_accession": accession,
        "pipeline_arrays_per_frma_call": 1,
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
            context: 38_463 for context in REQUIRED_ESTIMATOR_CONTEXTS
        },
        "estimator_maximum_absolute_difference": {
            context: 8.99e-4 for context in REQUIRED_ESTIMATOR_CONTEXTS
        },
        "estimator_all_bitwise_identical": False,
        "perturbation_control_digest_changed": True,
        "perturbation_control_magnitude": 1e-12,
        "features": 54_675,
        "finite_values": 54_675,
    }


def receipt() -> dict[str, object]:
    payload: dict[str, object] = {
        "schema_version": "masld-bench-gse49541-frma-single-array-firewall-v2",
        "series": "GSE49541",
        "platform_id": "GPL570",
        "method": "frma_with_exact_hgu133plus2frmavecs",
        "held_array_application": "one_array_at_a_time_with_frozen_training_object",
        "platform_features": 54_675,
        "cel_records_available": 72,
        "package_versions": dict(EXPECTED_PACKAGE_VERSIONS),
        "frma_vector_object_sha256": "f" * 64,
        "random_effect_probe": {"rule": "random_effect"},
        "per_target": [target("GSM1", "a" * 64), target("GSM2", "b" * 64)],
    }
    for flag in REQUIRED_FALSE_FLAGS:
        payload[flag] = False
    return payload


class PipelineFirewallTests(unittest.TestCase):
    def test_clean_receipt_passes(self) -> None:
        verdict = verify_receipt(receipt())
        self.assertEqual(verdict["status"], "pass_inductive_firewall_on_real_CEL")
        self.assertEqual(verdict["targets_verified"], 2)
        self.assertEqual(verdict["pipeline_contexts_per_target"], 4)
        self.assertEqual(verdict["arrays_per_frma_call"], 1)
        self.assertIs(verdict["held_array_summary_invariant_to_other_arrays_in_the_run"], True)
        self.assertIs(verdict["held_array_summary_invariant_to_processing_order"], True)

    def test_pipeline_context_divergence_is_rejected(self) -> None:
        payload = receipt()
        payload["per_target"][0]["pipeline_digests"]["solo_after_P"] = "c" * 64  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "what else the run processed"):
            verify_receipt(payload)

    def test_nonzero_pipeline_difference_is_rejected_even_when_digests_agree(self) -> None:
        payload = receipt()
        payload["per_target"][0]["pipeline_maximum_absolute_difference"]["solo_before_P"] = 1e-15  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "nonzero pipeline-context difference"):
            verify_receipt(payload)

    def test_differing_pipeline_element_count_is_rejected(self) -> None:
        payload = receipt()
        payload["per_target"][1]["pipeline_differing_elements"]["solo_interleaved_Q"] = 3  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "differing feature across pipeline contexts"):
            verify_receipt(payload)

    def test_multi_array_pipeline_call_is_rejected(self) -> None:
        payload = receipt()
        payload["per_target"][0]["pipeline_arrays_per_frma_call"] = 6  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "more than one array per call"):
            verify_receipt(payload)

    def test_missing_pipeline_context_is_rejected(self) -> None:
        payload = receipt()
        del payload["per_target"][0]["pipeline_digests"]["solo_before_P"]  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "every required pipeline_digests"):
            verify_receipt(payload)


class HazardMeasurementTests(unittest.TestCase):
    """The multi-array negative control must keep diverging, or it is broken."""

    def test_silent_negative_control_is_rejected(self) -> None:
        payload = receipt()
        for record in payload["per_target"]:  # type: ignore[union-attr]
            record["estimator_differing_elements"] = {
                context: 0 for context in REQUIRED_ESTIMATOR_CONTEXTS
            }
        with self.assertRaisesRegex(FirewallError, "hazard measurement is broken"):
            verify_receipt(payload)

    def test_divergence_on_one_target_is_enough_to_keep_the_control_alive(self) -> None:
        payload = receipt()
        payload["per_target"][0]["estimator_differing_elements"] = {  # type: ignore[index]
            context: 0 for context in REQUIRED_ESTIMATOR_CONTEXTS
        }
        verdict = verify_receipt(payload)
        self.assertIs(verdict["multi_array_frma_call_is_batch_dependent"], True)
        self.assertIs(verdict["batching_would_break_the_firewall"], True)

    def test_missing_estimator_context_is_rejected(self) -> None:
        payload = receipt()
        del payload["per_target"][0]["estimator_differing_elements"]["batch_Q_first"]  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "every required estimator_differing_elements"):
            verify_receipt(payload)

    def test_multi_array_call_declared_as_preprocessing_path_is_rejected(self) -> None:
        payload = receipt()
        payload["multi_array_frma_call_used_as_preprocessing_path"] = True
        with self.assertRaisesRegex(
            FirewallError, "multi_array_frma_call_used_as_preprocessing_path"
        ):
            verify_receipt(payload)


class ControlTests(unittest.TestCase):
    def test_undetected_perturbation_control_is_rejected(self) -> None:
        payload = receipt()
        payload["per_target"][0]["perturbation_control_digest_changed"] = False  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "perturbation control"):
            verify_receipt(payload)

    def test_identical_digests_across_targets_are_rejected(self) -> None:
        payload = receipt()
        payload["per_target"] = [target("GSM1", "a" * 64), target("GSM2", "a" * 64)]
        with self.assertRaisesRegex(FirewallError, "same summary digest"):
            verify_receipt(payload)

    def test_single_target_is_insufficient(self) -> None:
        payload = receipt()
        payload["per_target"] = [target("GSM1", "a" * 64)]
        with self.assertRaisesRegex(FirewallError, "at least two target arrays"):
            verify_receipt(payload)

    def test_nonfinite_summary_is_rejected(self) -> None:
        payload = receipt()
        payload["per_target"][0]["finite_values"] = 54_674  # type: ignore[index]
        with self.assertRaisesRegex(FirewallError, "nonfinite value"):
            verify_receipt(payload)


class ContractTests(unittest.TestCase):
    def test_every_forbidden_flag_must_be_false(self) -> None:
        for flag in REQUIRED_FALSE_FLAGS:
            payload = copy.deepcopy(receipt())
            payload[flag] = True
            with self.assertRaisesRegex(FirewallError, flag):
                verify_receipt(payload)

    def test_package_version_drift_is_rejected(self) -> None:
        payload = receipt()
        payload["package_versions"] = dict(EXPECTED_PACKAGE_VERSIONS) | {"frma": "1.60.0"}
        with self.assertRaisesRegex(FirewallError, "package versions differ"):
            verify_receipt(payload)

    def test_foreign_source_is_rejected(self) -> None:
        payload = receipt()
        payload["series"] = "GSE83452"
        with self.assertRaisesRegex(FirewallError, "different source"):
            verify_receipt(payload)

    def test_blocked_method_is_rejected(self) -> None:
        payload = receipt()
        payload["method"] = "all_sample_RMA"
        with self.assertRaisesRegex(FirewallError, "contract does not admit"):
            verify_receipt(payload)

    def test_missing_one_array_at_a_time_declaration_is_rejected(self) -> None:
        payload = receipt()
        payload["held_array_application"] = "whole_cohort"
        with self.assertRaisesRegex(FirewallError, "one-array-at-a-time"):
            verify_receipt(payload)

    def test_v1_schema_is_rejected(self) -> None:
        payload = receipt()
        payload["schema_version"] = "masld-bench-gse49541-frma-single-array-firewall-v1"
        with self.assertRaisesRegex(FirewallError, "schema differs"):
            verify_receipt(payload)

    def test_r_length_one_lists_are_unwrapped(self) -> None:
        payload = receipt()
        payload["platform_features"] = [54_675]
        payload["cel_records_available"] = [72]
        payload["labels_read"] = [False]
        payload["package_versions"] = {
            name: [value] for name, value in EXPECTED_PACKAGE_VERSIONS.items()
        }
        first = payload["per_target"][0]  # type: ignore[index]
        first["pipeline_all_bitwise_identical"] = [True]
        first["pipeline_arrays_per_frma_call"] = [1]
        first["features"] = [54_675]
        self.assertEqual(verify_receipt(payload)["targets_verified"], 2)


if __name__ == "__main__":
    unittest.main()
