from __future__ import annotations

import json
from pathlib import Path
import tomllib
import unittest

from scripts.build_gse244832_reference_guarded_atac_exchange import (
    LabelFreeATACExchangeError,
    target_family_eligibility,
    validate_guard,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/evaluation/gse244832_label_free_atac_exchange.toml"
GUARD = ROOT / "config/evaluation/gse244832_atac_reference_guard_20260825.json"
BUILDER = ROOT / "scripts/build_gse244832_reference_guarded_atac_exchange.py"
TEST = ROOT / "tests/unit/test_build_gse244832_reference_guarded_atac_exchange.py"


class GSE244832ReferenceGuardedATACExchangeTests(unittest.TestCase):
    def test_real_guard_and_implementation_are_bound(self) -> None:
        contract, guard, resolved = validate_guard(ROOT, CONTRACT, GUARD, BUILDER, TEST)
        self.assertEqual(contract["target_dataset_id"], "gse244832")
        self.assertFalse(guard["sequence_or_coordinate_family_execution_authorized"])
        self.assertTrue(guard["observed_atac_query_materialization_authorized"])
        self.assertEqual(len(resolved["family"]), 25)

    def test_only_observed_atac_regime_remains_potentially_executable(self) -> None:
        with CONTRACT.open("rb") as handle:
            template = tomllib.load(handle)
        authority = ROOT / template["family_eligibility_authority"]["path"]
        with authority.open("rb") as handle:
            family = tomllib.load(handle)["family_eligibility"]
        guard = json.loads(GUARD.read_text(encoding="utf-8"))
        guarded = target_family_eligibility(family, guard)
        for row in guarded:
            if row["input_regime"] != "observed_atac":
                self.assertTrue(row["eligibility"].startswith("ineligible_"))
                self.assertFalse(row["query_atac_allowed"])

    def test_unclassified_regime_fails_closed(self) -> None:
        guard = json.loads(GUARD.read_text(encoding="utf-8"))
        with self.assertRaises(LabelFreeATACExchangeError):
            target_family_eligibility(
                [{"model_id": "x", "input_regime": "unknown", "eligibility": "eligible"}],
                guard,
            )


if __name__ == "__main__":
    unittest.main()
