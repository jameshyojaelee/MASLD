from __future__ import annotations

import copy
import unittest

from scripts.aggregate_lsgkm_gse281364_dinucleotide_scores import (
    CONTEXTS,
    LSGKMAggregationError,
    READOUTS,
    SEEDS,
    SPLIT_SEEDS,
    bind_manifest_record_only,
    expected_remaining_fits,
    validate_config,
)


def valid_config() -> dict:
    return {
        "schema_version": "masld-bench-lsgkm-gse281364-dinucleotide-aggregation-input-v1",
        "status": "aggregate_25_frozen_fits_without_control_or_outcome_values",
        "design_id": "lsgkm_hepatocyte_atac_peak_exact_dinucleotide_null_v1",
        "fit_grid": {
            "fixed_seeds": list(SEEDS),
            "total_fit_count": 25,
            "reused_fit": {"split_id": "donor0_genomic0", "model_seed": 1103},
            "remaining_fit_count": 24,
            "expected_elements_per_seed": 1033,
            "expected_rows_per_readout": 10330,
        },
        "bundle_fits": [
            {
                "split_id": split_id,
                "model_seeds": list(seeds),
                "artifacts_path": f"executions/{split_id}/ARTIFACTS.json",
                "artifacts_sha256": "0" * 64,
            }
            for split_id, seeds in SPLIT_SEEDS.items()
        ],
        "action_firewall": {
            "lsgkm_prediction_value_read_authorized": True,
            "row_universe_metadata_read_authorized": True,
            "raw_score_aggregation_authorized": True,
            "control_manifest_identity_bind_authorized": True,
            "control_prediction_member_open_authorized": False,
            "outcome_access_authorized": False,
            "reporter_count_access_authorized": False,
            "benchmark_metric_calculation_authorized": False,
            "calibration_or_head_fit_authorized": False,
            "candidate_selection_authorized": False,
            "sealed_asset_access_authorized": False,
        },
        "output_contract": {
            "readout_ids": list(READOUTS),
            "contexts": list(CONTEXTS),
            "rows_per_readout": 10330,
            "total_rows": 20660,
            "static_score_is_context_specific": False,
            "static_native_score_directly_stackable": False,
            "calibration_status": "uncalibrated_static_rank_score",
        },
    }


class LSGKMAggregationTests(unittest.TestCase):
    def test_valid_24_fit_rectangle_plus_reused_probe(self) -> None:
        config = valid_config()
        validate_config(config)
        self.assertEqual(len(expected_remaining_fits()), 24)
        self.assertNotIn(("donor0_genomic0", 1103), expected_remaining_fits())

    def test_control_member_binding_uses_manifest_record_only(self) -> None:
        record = bind_manifest_record_only(
            {"oof_predictions.tsv.gz": {"path": "oof_predictions.tsv.gz", "sha256": "a" * 64}},
            "oof_predictions.tsv.gz",
            "a" * 64,
        )
        self.assertEqual(record["sha256"], "a" * 64)

    def test_control_open_permission_is_rejected(self) -> None:
        changed = copy.deepcopy(valid_config())
        changed["action_firewall"]["control_prediction_member_open_authorized"] = True
        with self.assertRaises(LSGKMAggregationError):
            validate_config(changed)

    def test_metric_or_outcome_permission_is_rejected(self) -> None:
        for field in ("outcome_access_authorized", "benchmark_metric_calculation_authorized"):
            with self.subTest(field=field):
                changed = copy.deepcopy(valid_config())
                changed["action_firewall"][field] = True
                with self.assertRaises(LSGKMAggregationError):
                    validate_config(changed)

    def test_missing_bundle_seed_is_rejected(self) -> None:
        changed = copy.deepcopy(valid_config())
        changed["bundle_fits"][1]["model_seeds"].pop()
        with self.assertRaises(LSGKMAggregationError):
            validate_config(changed)


if __name__ == "__main__":
    unittest.main()
