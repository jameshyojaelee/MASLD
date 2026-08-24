from __future__ import annotations

import unittest

from scripts.build_corgi_gse296875_context_counts import CorgiContextError, resolve_roster


class CorgiContextCountTests(unittest.TestCase):
    def test_resolution_prefers_exact_id_then_unique_symbol(self) -> None:
        roster_ids = [f"ENSG{i:011d}" for i in range(2_891)]
        roster_symbols = [f"G{i}" for i in range(2_891)]
        source_ids = roster_ids[:-2] + ["ENSG99999999998"]
        source_symbols = roster_symbols[:-2] + [roster_symbols[-2]]
        rows = resolve_roster(roster_ids, roster_symbols, source_ids, source_symbols)
        self.assertEqual(rows[0]["mapping_state"], "exact_stable_id")
        self.assertEqual(rows[-2]["mapping_state"], "unique_symbol_rescue")
        self.assertEqual(rows[-1]["mapping_state"], "absent_source_annotation_structurally_missing")
        self.assertFalse(rows[-1]["observed_feature"])

    def test_resolution_rejects_duplicate_roster(self) -> None:
        roster_ids = ["ENSG00000000001"] * 2_891
        symbols = [f"G{i}" for i in range(2_891)]
        with self.assertRaisesRegex(CorgiContextError, "duplicated"):
            resolve_roster(roster_ids, symbols, ["ENSG00000000001"], ["G0"])


if __name__ == "__main__":
    unittest.main()
