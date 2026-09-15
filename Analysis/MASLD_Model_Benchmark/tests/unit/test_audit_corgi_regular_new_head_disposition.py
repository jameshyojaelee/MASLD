from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.audit_corgi_regular_new_head_disposition import (
    CorgiNewHeadDispositionError,
    audit_disposition,
)


ROOT = Path(__file__).resolve().parents[2]
AUTHORITY = (
    ROOT
    / "config/artifacts/models/corgi/new_assay_head_development_disposition_20260825.json"
)


class CorgiRegularNewHeadDispositionTests(unittest.TestCase):
    def test_negative_development_result_is_frozen(self) -> None:
        receipt = audit_disposition(ROOT, AUTHORITY)
        self.assertEqual(receipt["status"], "pass_frozen_negative_development_result")
        self.assertEqual(receipt["selected_alpha_by_fold"], [0.0] * 5)
        self.assertEqual(receipt["donors"], 39)
        self.assertEqual(receipt["donor_by_chromosome_units_per_context_arm"], 191)
        self.assertIs(receipt["numeric_residual_is_biological_gain"], False)
        self.assertIs(receipt["incremental_corgi_profile_contribution"], False)
        self.assertIs(receipt["film_plus_head_execution_authorized"], False)
        self.assertIs(receipt["sealed_or_test_data_opened_by_audit"], False)

    def test_rejects_recasting_numeric_residual_as_gain(self) -> None:
        authority = json.loads(AUTHORITY.read_text(encoding="utf-8"))
        authority["independent_evaluation"]["reported_numeric_residual_interpretation"] = (
            "biological_gain"
        )
        with tempfile.TemporaryDirectory() as temporary:
            mutated = Path(temporary) / "authority.json"
            mutated.write_text(json.dumps(authority), encoding="utf-8")
            with self.assertRaises(CorgiNewHeadDispositionError):
                audit_disposition(ROOT, mutated)

    def test_rejects_open_execution_gate(self) -> None:
        authority = json.loads(AUTHORITY.read_text(encoding="utf-8"))
        authority["next_rung"]["training_or_prediction_authorized_by_this_disposition"] = True
        with tempfile.TemporaryDirectory() as temporary:
            mutated = Path(temporary) / "authority.json"
            mutated.write_text(json.dumps(authority), encoding="utf-8")
            with self.assertRaises(CorgiNewHeadDispositionError):
                audit_disposition(ROOT, mutated)


if __name__ == "__main__":
    unittest.main()
