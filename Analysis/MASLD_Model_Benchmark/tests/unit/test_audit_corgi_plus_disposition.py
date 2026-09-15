from __future__ import annotations

from pathlib import Path
import unittest

from scripts.audit_corgi_plus_disposition import audit_disposition


ROOT = Path(__file__).resolve().parents[2]
AUTHORITY = (
    ROOT / "config/artifacts/models/corgi_plus/release_disposition_20260824.json"
)


class CorgiPlusDispositionAuditTests(unittest.TestCase):
    def test_current_release_is_terminal_blocked_and_regular_corgi_is_not_relabelled(self) -> None:
        receipt = audit_disposition(ROOT, AUTHORITY)
        self.assertEqual(receipt["status"], "pass_frozen_terminal_blocker")
        self.assertIs(receipt["corgi_plus_executable_for_registered_task"], False)
        self.assertIs(receipt["corgi_plus_production_sbatch_authorized"], False)
        self.assertEqual(receipt["nearest_legitimate_comparator"], "corgi_regular")
        self.assertIs(receipt["outcomes_opened"], False)
        self.assertIs(receipt["sealed_features_opened"], False)


if __name__ == "__main__":
    unittest.main()
