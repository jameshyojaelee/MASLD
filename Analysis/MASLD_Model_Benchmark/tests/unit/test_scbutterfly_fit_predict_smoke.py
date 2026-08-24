from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "scbutterfly_fit_predict_smoke.py"
SPEC = importlib.util.spec_from_file_location("scbutterfly_safe_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
safe = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(safe)


class SCButterflySafeTests(unittest.TestCase):
    def test_cli_exposes_independent_fold_selection(self) -> None:
        source = MODULE.read_text(encoding="utf-8")
        self.assertIn('parser.add_argument("--fold", type=int, action="append")', source)
        self.assertIn('for fold in sorted(folds):', source)

    def test_chromosome_contract_accepts_contiguous_blocks(self) -> None:
        peaks = [
            {"chromosome": "chr1"},
            {"chromosome": "chr1"},
            {"chromosome": "chr2"},
        ]
        self.assertEqual(safe.chromosome_contract(peaks), [2, 1])

    def test_chromosome_contract_rejects_reentry(self) -> None:
        peaks = [
            {"chromosome": "chr1"},
            {"chromosome": "chr2"},
            {"chromosome": "chr1"},
        ]
        with self.assertRaisesRegex(safe.SCButterflyError, "not contiguous"):
            safe.chromosome_contract(peaks)


if __name__ == "__main__":
    unittest.main()
