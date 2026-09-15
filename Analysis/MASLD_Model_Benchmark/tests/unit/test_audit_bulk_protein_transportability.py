#!/usr/bin/env python3
"""Unit checks for bulk/protein transportability reconciliation."""

from __future__ import annotations

import copy
import importlib.util
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts/audit_bulk_protein_transportability.py"
SPEC = importlib.util.spec_from_file_location("bulk_protein_transportability", MODULE)
assert SPEC and SPEC.loader
transport_audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(transport_audit)


class BulkProteinTransportabilityAuditTest(unittest.TestCase):
    @staticmethod
    def _census() -> dict:
        return json.loads(
            (ROOT / transport_audit.CENSUS_RELATIVE).read_text(encoding="utf-8")
        )

    def test_validate_census_accepts_fail_closed_authority(self) -> None:
        cohorts = transport_audit.validate_census(self._census())
        self.assertEqual(len(cohorts), 10)
        self.assertEqual(
            cohorts["pxd051911"]["cross_cohort_rna_pairing"],
            "same_study_unpaired",
        )

    def test_validate_census_rejects_missing_as_zero_roster(self) -> None:
        census = copy.deepcopy(self._census())
        census["global_contract"]["allowed_missingness_states"].append("zero")
        with self.assertRaisesRegex(
            transport_audit.BulkProteinAuditError, "missingness state roster"
        ):
            transport_audit.validate_census(census)

    def test_validate_census_rejects_sealed_metadata_escape(self) -> None:
        census = copy.deepcopy(self._census())
        row = next(
            value
            for value in census["cohorts"]
            if value["cohort_family_id"] == "gse267031"
        )
        row["metadata"]["fibrosis_or_stage"]["state"] = "observed"
        with self.assertRaisesRegex(
            transport_audit.BulkProteinAuditError, "escaped the seal"
        ):
            transport_audit.validate_census(census)

    def test_validate_census_rejects_false_rna_protein_pair(self) -> None:
        census = copy.deepcopy(self._census())
        row = next(
            value
            for value in census["cohorts"]
            if value["cohort_family_id"] == "pxd051911"
        )
        row["cross_cohort_rna_pairing"] = "same_donor_different_tissue"
        with self.assertRaisesRegex(
            transport_audit.BulkProteinAuditError, "PXD051911 topology"
        ):
            transport_audit.validate_census(census)

    def test_validate_census_rejects_silent_submission_authority(self) -> None:
        census = copy.deepcopy(self._census())
        census["next_outcome_blind_campaign"][
            "new_submission_authorized_by_this_reconciliation"
        ] = True
        with self.assertRaisesRegex(
            transport_audit.BulkProteinAuditError,
            "next outcome-blind campaign boundary",
        ):
            transport_audit.validate_census(census)

    def test_index_rows_rejects_duplicate_ids(self) -> None:
        value = {"cohorts": [{"id": "same"}, {"id": "same"}]}
        with self.assertRaisesRegex(
            transport_audit.BulkProteinAuditError, "duplicate or malformed"
        ):
            transport_audit.index_rows(value, "cohorts", "id")


if __name__ == "__main__":
    unittest.main()
