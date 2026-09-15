from __future__ import annotations

import unittest

import numpy as np

from scripts.commit_gse281367_atac_transport_predictions import (
    ATACPredictionCommitError,
    validate_array_contract,
    validate_evidence_contract,
)


class CommitGSE281367ATACTransportPredictionsTests(unittest.TestCase):
    def _manifest(self, regime: str = "sequence_only") -> dict[str, object]:
        observed = regime == "observed_atac"
        return {
            "input_regime": regime,
            "rotation_id": "valid_context_test_target",
            "prediction_layout": "donor_lineage" if observed else "lineage_broadcast",
            "evidence_contract": {
                "condition_labels_consumed": False,
                "phenotype_values_consumed": False,
                "development_outcomes_consumed": False,
                "sealed_data_consumed": False,
                "target_role_atac_consumed": False,
                "target_role_contigs_entered_query_transform": False,
                "missing_as_zero": False,
                "scored_target_role": "test",
                "observed_atac_context_role": "valid" if observed else "none",
                "query_atac_consumed": observed,
                "receptive_field_bp": 1000,
                "query_transform_scope": "per_donor_context_only" if observed else "source_training_only",
            },
        }

    def test_sequence_and_observed_atac_evidence_contracts_pass_separately(self) -> None:
        validate_evidence_contract(self._manifest("sequence_only"))
        validate_evidence_contract(self._manifest("observed_atac"))

    def test_sequence_only_cannot_consume_query_atac(self) -> None:
        manifest = self._manifest("sequence_only")
        manifest["evidence_contract"]["query_atac_consumed"] = True
        with self.assertRaises(ATACPredictionCommitError):
            validate_evidence_contract(manifest)

    def test_target_atac_or_condition_consumption_fails(self) -> None:
        for field in ("target_role_atac_consumed", "condition_labels_consumed"):
            manifest = self._manifest("observed_atac")
            manifest["evidence_contract"][field] = True
            with self.assertRaises(ATACPredictionCommitError):
                validate_evidence_contract(manifest)

    def test_missing_prediction_is_nan_not_zero(self) -> None:
        values = np.ones((4, 16000), dtype=np.float32)
        states = np.zeros((4, 16000), dtype=np.uint8)
        states[:, 0] = 5
        values[:, 0] = np.nan
        receipt = validate_array_contract(
            values,
            states,
            output_family="functional_track_count",
            layout="lineage_broadcast",
        )
        self.assertGreaterEqual(receipt["minimum_coverage"], 0.95)
        values[:, 0] = 0
        with self.assertRaises(ATACPredictionCommitError):
            validate_array_contract(
                values,
                states,
                output_family="functional_track_count",
                layout="lineage_broadcast",
            )

    def test_profile_geometry_and_coverage_are_enforced(self) -> None:
        values = np.ones((4, 16000, 20), dtype=np.float32)
        states = np.zeros((4, 16000), dtype=np.uint8)
        validate_array_contract(
            values,
            states,
            output_family="functional_track_profile",
            layout="lineage_broadcast",
        )
        states[:, :1000] = 5
        values[:, :1000] = np.nan
        with self.assertRaises(ATACPredictionCommitError):
            validate_array_contract(
                values,
                states,
                output_family="functional_track_profile",
                layout="lineage_broadcast",
            )


if __name__ == "__main__":
    unittest.main()
