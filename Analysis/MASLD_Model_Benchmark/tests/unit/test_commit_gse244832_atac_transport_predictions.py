from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.commit_gse244832_atac_transport_predictions import (
    ATACPredictionCommitError,
    resolve_model_manifest,
    validate_array_contract,
    validate_evidence_contract,
)


class CommitGSE244832ATACTransportPredictionsTests(unittest.TestCase):
    def _manifest(self) -> dict[str, object]:
        return {
            "input_regime": "observed_atac",
            "rotation_id": "valid_context_test_target",
            "prediction_layout": "donor_lineage",
            "evidence_contract": {
                "condition_labels_consumed": False,
                "phenotype_values_consumed": False,
                "rna_assay_consumed": False,
                "cross_assay_join_consumed": False,
                "sequence_features_consumed": False,
                "source_native_reference_bundle_assumed_resolved": False,
                "development_outcomes_consumed": False,
                "sealed_data_consumed": False,
                "target_role_atac_consumed": False,
                "target_role_contigs_entered_query_transform": False,
                "missing_as_zero": False,
                "scored_target_role": "test",
                "observed_atac_context_role": "valid",
                "query_atac_consumed": True,
                "receptive_field_bp": "not_applicable_global_observed_atac_encoder",
                "query_transform_scope": "per_donor_context_only",
                "per_donor_context_operation": (
                    "deterministic_normalization_without_fitted_target_state"
                ),
            },
        }

    def test_observed_atac_evidence_contract_passes(self) -> None:
        validate_evidence_contract(self._manifest())

    def test_target_transductive_query_fit_is_prohibited(self) -> None:
        manifest = self._manifest()
        manifest["evidence_contract"][
            "query_transform_scope"
        ] = "target_context_transductive_separate_lane"
        with self.assertRaises(ATACPredictionCommitError):
            validate_evidence_contract(manifest)

    def test_source_training_scope_has_no_target_operation(self) -> None:
        manifest = self._manifest()
        manifest["evidence_contract"]["query_transform_scope"] = "source_training_only"
        manifest["evidence_contract"]["per_donor_context_operation"] = "not_applicable"
        validate_evidence_contract(manifest)
        manifest["evidence_contract"]["per_donor_context_operation"] = "learned_target_fit"
        with self.assertRaises(ATACPredictionCommitError):
            validate_evidence_contract(manifest)

    def test_model_manifest_absolute_and_traversal_paths_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            for escaped in ("/tmp/model.json", "../model.json", "nested/../../model.json"):
                with self.assertRaises(ATACPredictionCommitError):
                    resolve_model_manifest(root, {"path": escaped, "sha256": "0" * 64})

    def test_sequence_or_rna_consumption_fails_closed(self) -> None:
        for field in (
            "sequence_features_consumed",
            "rna_assay_consumed",
            "cross_assay_join_consumed",
            "source_native_reference_bundle_assumed_resolved",
        ):
            manifest = self._manifest()
            manifest["evidence_contract"][field] = True
            with self.assertRaises(ATACPredictionCommitError):
                validate_evidence_contract(manifest)

    def test_sequence_input_regime_is_prohibited(self) -> None:
        manifest = self._manifest()
        manifest["input_regime"] = "sequence_only"
        with self.assertRaises(ATACPredictionCommitError):
            validate_evidence_contract(manifest)

    def test_count_geometry_and_missing_nan_are_enforced(self) -> None:
        values = np.ones((18, 4, 16000), dtype=np.float32)
        states = np.zeros((18, 4, 16000), dtype=np.uint8)
        unit_states = np.zeros((18, 4), dtype=np.uint8)
        unit_states[12:] = 3
        states[12:] = 3
        values[12:] = np.nan
        states[:, :, 0] = 5
        values[:, :, 0] = np.nan
        states[12:] = 3
        receipt = validate_array_contract(
            values,
            states,
            unit_states,
            output_family="masked_accessibility_count",
        )
        self.assertGreaterEqual(receipt["minimum_coverage"], 0.95)
        values[:, :, 0] = 0
        with self.assertRaises(ATACPredictionCommitError):
            validate_array_contract(
                values,
                states,
                unit_states,
                output_family="masked_accessibility_count",
            )

    def test_below_qc_unit_must_be_explicitly_masked(self) -> None:
        values = np.ones((18, 4, 16000), dtype=np.float32)
        states = np.zeros((18, 4, 16000), dtype=np.uint8)
        unit_states = np.zeros((18, 4), dtype=np.uint8)
        unit_states[12:] = 3
        states[12:] = 3
        values[12:] = np.nan
        states[12, 0, 0] = 0
        values[12, 0, 0] = 1
        with self.assertRaises(ATACPredictionCommitError):
            validate_array_contract(
                values,
                states,
                unit_states,
                output_family="masked_accessibility_count",
            )

    def test_profile_geometry_is_enforced(self) -> None:
        values = np.ones((18, 4, 16000, 20), dtype=np.float32)
        states = np.zeros((18, 4, 16000), dtype=np.uint8)
        unit_states = np.zeros((18, 4), dtype=np.uint8)
        unit_states[12:] = 1
        states[12:] = 1
        values[12:] = np.nan
        validate_array_contract(
            values,
            states,
            unit_states,
            output_family="functional_track_profile",
        )
        with self.assertRaises(ATACPredictionCommitError):
            validate_array_contract(
                values[:12],
                states[:12],
                unit_states,
                output_family="functional_track_profile",
            )


if __name__ == "__main__":
    unittest.main()
