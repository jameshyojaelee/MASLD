from __future__ import annotations

import unittest

import json
from pathlib import Path

from scripts.run_gse296875_observed_multiome_factorized_seed_bundle import SEEDS, SeedBundleError, planned_surfaces, validate


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_factorized_full_rectangle_contract_20260825.json"
CONTRACT = ROOT / "executions/model-check-321-21099900"


class ObservedMultiomeFactorizedSeedBundleTests(unittest.TestCase):
    def test_real_contract_and_seed_bundle_parents_pass_without_evaluator_values(self) -> None:
        config = json.loads(CONFIG.read_text(encoding="utf-8"))
        inputs = validate(ROOT, config, CONTRACT, SEEDS[0])
        self.assertEqual(len(inputs["labels"]), 25)
        self.assertEqual(len(inputs["evaluator"]), 5)

    def test_each_seed_has_exact_complete_surface_rectangle(self) -> None:
        all_rows = []
        for seed in SEEDS:
            rows = planned_surfaces(seed)
            self.assertEqual(len(rows), 25)
            self.assertEqual({(row["held_genomic_fold"], row["held_donor_fold"]) for row in rows}, {(genomic, donor) for genomic in range(5) for donor in range(5)})
            all_rows.extend((row["held_genomic_fold"], row["held_donor_fold"], row["seed"]) for row in rows)
        self.assertEqual(len(set(all_rows)), 125)

    def test_unregistered_seed_is_rejected(self) -> None:
        with self.assertRaises(SeedBundleError):
            planned_surfaces(7)


if __name__ == "__main__":
    unittest.main()
