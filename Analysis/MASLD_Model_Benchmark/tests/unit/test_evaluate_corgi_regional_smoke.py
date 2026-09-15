from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.evaluate_corgi_regional_smoke import (
    CONTEXT_ARMS,
    CorgiEvaluationError,
    multinomial_deviance_per_insertion,
    spearman_or_zero,
    validate_prediction_tables,
    verify_selected_member,
)


class CorgiRegionalSmokeEvaluationTests(unittest.TestCase):
    def test_perfect_profile_beats_uniform_profile(self) -> None:
        observed = [1.0, 4.0, 15.0, 80.0]
        perfect = multinomial_deviance_per_insertion(observed, observed)
        uniform = multinomial_deviance_per_insertion(observed, [1.0, 1.0, 1.0, 1.0])
        self.assertLess(perfect, uniform)
        self.assertAlmostEqual(spearman_or_zero(observed, observed), 1.0)

    def test_tied_spearman_is_finite(self) -> None:
        value = spearman_or_zero([0.0, 0.0, 1.0, 2.0], [0.0, 0.0, 2.0, 1.0])
        self.assertGreater(value, 0.0)
        self.assertLess(value, 1.0)

    def test_prediction_contract_rejects_test_role(self) -> None:
        records = []
        for index, arm in enumerate(CONTEXT_ARMS):
            records.append(
                {
                    "prediction_index": str(index),
                    "unit_index": "7",
                    "donor_id": "D1",
                    "lineage_id": "hepatocyte",
                    "outer_fold": "1",
                    "donor_fold": "1",
                    "context_arm": arm,
                    "context_source_unit": "7",
                    "outcome_role": "test" if index == 0 else "valid",
                }
            )
        windows = [
            {
                "contig": "chr1",
                "output_start": str(index * 1_000),
                "output_end": str(index * 1_000 + 1_000),
                "window_id": f"w{index}",
                "genomic_fold": "1",
                "selection_hash": f"h{index}",
            }
            for index in range(3)
        ]
        receipt = {
            "schema_version": "masld-bench-corgi-outcome-aligned-prediction-v3",
            "status": "pass_outcome_free_prediction",
            "outcome_role": "valid",
            "lineage_id": "hepatocyte",
            "mapper_outer_fold": 0,
            "valid_donor_fold": 1,
            "valid_genomic_fold": 1,
            "held_ATAC_or_other_outcomes_used": False,
            "test_or_sealed_features_or_labels_read": False,
            "model_fitted_or_adapted": False,
            "context_arms": list(CONTEXT_ARMS),
            "regional_score_names": ["strand_tta_orientation_mean_softplus_regional_sum"],
            "valid_donors": 1,
            "scoreable_windows": 3,
        }
        with self.assertRaisesRegex(CorgiEvaluationError, "firewall"):
            validate_prediction_tables(records, windows, receipt)

    def test_deviance_rejects_missing_or_negative_signal(self) -> None:
        with self.assertRaises(CorgiEvaluationError):
            multinomial_deviance_per_insertion([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
        with self.assertRaises(CorgiEvaluationError):
            multinomial_deviance_per_insertion([1.0, 2.0, 3.0], [1.0, -1.0, 1.0])

    def test_selected_prediction_member_rejects_post_freeze_tampering(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            member = root / "fold0" / "receipt.json"
            member.parent.mkdir()
            member.write_text('{"status":"frozen"}\n', encoding="utf-8")
            checksum = sha256(member.read_bytes()).hexdigest()
            (root / "ARTIFACTS.json").write_text(
                json.dumps(
                    {
                        "schema_version": "masld-bench-artifacts-v1",
                        "artifacts": [
                            {
                                "path": "fold0/receipt.json",
                                "sha256": checksum,
                                "size_bytes": member.stat().st_size,
                            }
                        ],
                        "metadata": {},
                    }
                ),
                encoding="utf-8",
            )
            self.assertEqual(
                verify_selected_member(root, "fold0/receipt.json", checksum), member
            )
            member.write_text('{"status":"tampered"}\n', encoding="utf-8")
            with self.assertRaisesRegex(CorgiEvaluationError, "differs"):
                verify_selected_member(root, "fold0/receipt.json", checksum)


if __name__ == "__main__":
    unittest.main()
