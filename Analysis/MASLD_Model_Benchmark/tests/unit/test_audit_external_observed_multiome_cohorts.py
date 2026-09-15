from __future__ import annotations

import copy
from pathlib import Path
import tomllib
import unittest

from scripts.audit_external_observed_multiome_cohorts import (
    CohortAuditError,
    build_summary,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/evaluation/external_observed_multiome_cohort_expansion.toml"


class ExternalObservedMultiomeCohortAuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = tomllib.loads(CONTRACT.read_text(encoding="utf-8"))

    def test_roster_separates_topology_and_biological_role(self) -> None:
        candidates = validate_contract(self.contract)
        by_id = {row["candidate_id"]: row for row in candidates}
        self.assertEqual(by_id["fnih_86_liver"]["topology"], "same_nucleus")
        self.assertEqual(by_id["fnih_86_liver"]["exact_donor_n"], 86)
        self.assertEqual(
            by_id["gse281574_alcohol_multiome"]["contamination_role"],
            "etiology_ood_development_stress_test_not_masld_negative_control",
        )
        self.assertEqual(
            by_id["gse244832"]["topology"], "same_sample_different_aliquot"
        )
        self.assertIn("reference_only", by_id["encode_gse286690_liver_multiome"]["contamination_role"])
        self.assertIn(
            "EGAD00001010049", by_id["ega_hepatoblastoma_multiome"]["accessions"]
        )
        self.assertNotIn(
            "EGAD00001015262", by_id["ega_hepatoblastoma_multiome"]["accessions"]
        )
        self.assertIn("access_right=restricted", by_id["fnih_86_liver"]["permissions"])

    def test_no_candidate_can_claim_external_champion(self) -> None:
        candidates = validate_contract(self.contract)
        self.assertTrue(
            all(not row["eligible_for_masld_external_champion"] for row in candidates)
        )
        summary = build_summary(
            candidates, [{"source_id": "synthetic", "available": True}]
        )
        self.assertFalse(summary["paired_masld_external_seal_available_now"])
        self.assertFalse(summary["observed_multiome_external_champion_available_now"])

    def test_sealed_outcome_and_matrix_firewalls_are_binding(self) -> None:
        for field in ("sealed_outcomes_read", "biological_matrices_downloaded"):
            changed = copy.deepcopy(self.contract)
            changed[field] = True
            with self.assertRaises(CohortAuditError):
                validate_contract(changed)

    def test_missing_category_is_rejected(self) -> None:
        changed = copy.deepcopy(self.contract)
        changed["candidate"].pop()
        with self.assertRaises(CohortAuditError):
            validate_contract(changed)


if __name__ == "__main__":
    unittest.main()
