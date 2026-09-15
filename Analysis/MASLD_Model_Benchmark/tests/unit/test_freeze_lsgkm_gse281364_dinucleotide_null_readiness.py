from __future__ import annotations

import copy
from pathlib import Path
import unittest

from scripts.freeze_lsgkm_gse281364_dinucleotide_null_readiness import (
    DinucleotideNullReadinessError,
    SEEDS,
    canonical_sequence,
    dinucleotide_counts,
    exact_dinucleotide_shuffle,
    five_seed_shuffles,
    load_json,
    monomer_counts,
    validate_config_boundary,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/lsgkm_gse281364_dinucleotide_null_readiness.json"


class LSGKMDinucleotideNullReadinessTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = load_json(CONFIG)
        self.sequence = ("ACGTTGCAAGTCGATCGGATCCGATGCTAGCTAGGCTA" * 8)[:300]

    def test_exact_shuffle_is_deterministic_and_preserves_native_statistics(self) -> None:
        kwargs = {
            "design_id": self.config["design_id"],
            "split_id": "donor0_genomic0",
            "positive_id": "peak-1",
            "model_seed": 1103,
            "attempt": 0,
        }
        first = exact_dinucleotide_shuffle(self.sequence, **kwargs)
        second = exact_dinucleotide_shuffle(self.sequence, **kwargs)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 300)
        self.assertEqual(first[0], self.sequence[0])
        self.assertEqual(first[-1], self.sequence[-1])
        self.assertEqual(monomer_counts(first), monomer_counts(self.sequence))
        self.assertEqual(dinucleotide_counts(first), dinucleotide_counts(self.sequence))

    def test_five_seed_shuffle_roster_is_canonically_distinct(self) -> None:
        result = five_seed_shuffles(
            self.sequence,
            design_id=self.config["design_id"],
            split_id="donor0_genomic0",
            positive_id="peak-1",
            max_attempts=256,
        )
        self.assertEqual(list(result), SEEDS)
        canonical = {canonical_sequence(value) for value in result.values()}
        self.assertEqual(len(canonical), 5)
        self.assertNotIn(canonical_sequence(self.sequence), canonical)

    def test_low_complexity_sequence_fails_closed(self) -> None:
        with self.assertRaises(DinucleotideNullReadinessError):
            five_seed_shuffles(
                "A" * 300,
                design_id=self.config["design_id"],
                split_id="donor0_genomic0",
                positive_id="degenerate",
                max_attempts=8,
            )

    def test_ambiguous_sequence_fails_closed(self) -> None:
        with self.assertRaises(DinucleotideNullReadinessError):
            exact_dinucleotide_shuffle(
                "A" * 299 + "N",
                design_id=self.config["design_id"],
                split_id="donor0_genomic0",
                positive_id="ambiguous",
                model_seed=1103,
                attempt=0,
            )

    def test_fixed_10000_pair_floor_cannot_be_relaxed(self) -> None:
        changed = copy.deepcopy(self.config)
        changed["negative_design"]["positive_count_per_fit"] = 9999
        changed["negative_design"]["minimum_positive_count_per_fit"] = 9999
        with self.assertRaises(DinucleotideNullReadinessError):
            validate_config_boundary(changed)

    def test_old_genomic_matching_tolerance_cannot_be_reintroduced(self) -> None:
        changed = copy.deepcopy(self.config)
        changed["negative_design"]["gc_tolerance"] = 0.05
        with self.assertRaises(DinucleotideNullReadinessError):
            validate_config_boundary(changed)

    def test_failed_campaign_is_history_only_and_not_evidence(self) -> None:
        validate_config_boundary(self.config)
        relation = self.config["campaign_relation"]
        self.assertTrue(relation["independent_redesign"])
        self.assertFalse(relation["failed_campaign_evidence_reused"])
        self.assertFalse(relation["failed_campaign_thresholds_reinterpreted"])
        self.assertIn("not-scientific-feasibility-evidence", relation["historical_boundary"]["role"])

    def test_outcomes_predictions_and_training_are_firewalled(self) -> None:
        validate_config_boundary(self.config)
        self.assertTrue(all(value is False for value in self.config["action_firewall"].values()))
        self.assertFalse(self.config["comparison"]["prediction_values_read_during_design"])
        self.assertFalse(self.config["comparison"]["reporter_outcomes_read_during_design"])
        self.assertFalse(self.config["production_gate"]["production_training_currently_authorized"])

    def test_claim_is_peak_sequence_order_not_accessibility_classification(self) -> None:
        validate_config_boundary(self.config)
        negative = self.config["negative_design"]
        claims = self.config["claim_limits"]
        self.assertIn("sequence_order", negative["estimand"])
        self.assertIn("open_versus_closed_chromatin", negative["not_an_estimand"])
        self.assertFalse(claims["accessibility_classifier_claim_supported"])
        self.assertFalse(claims["condition_specific_effect_supported"])

    def test_exact_control_identity_and_row_universe_are_fixed(self) -> None:
        validate_config_boundary(self.config)
        comparison = self.config["comparison"]
        scoring = self.config["scoring_universe"]
        self.assertEqual(comparison["control_model_id"], "available_simple_controls")
        self.assertEqual(comparison["control_head_id"], "allele_identity_ridge")
        self.assertEqual(comparison["canonical_control_rows_sha256"], "7e76a9d8030ee8ed763223c01a660911aacbd47285d4c522a83297a1300ae42a")
        self.assertEqual(scoring["element_count"], 1033)
        self.assertEqual(scoring["long_range_block_count"], 239)
        self.assertEqual(scoring["rows_per_readout"], 10330)


if __name__ == "__main__":
    unittest.main()
