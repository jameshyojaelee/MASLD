from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "summarize_stabmap_seed_stability.py"
SPEC = importlib.util.spec_from_file_location("stabmap_stability_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
stability = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stability)


class StabMapStabilityTests(unittest.TestCase):
    def test_parser_requires_three_unique_seeds(self) -> None:
        values = stability.parse_evaluations(
            ["4817=/tmp/a", "5839=/tmp/b", "6857=/tmp/c"]
        )
        self.assertEqual(set(values), {4817, 5839, 6857})
        with self.assertRaises(stability.StabMapStabilityError):
            stability.parse_evaluations(["4817=/tmp/a"])

    def test_seeds_are_not_declared_replicates(self) -> None:
        source = MODULE.read_text(encoding="utf-8")
        self.assertIn("seeds are not biological replicates", source)
        self.assertIn('"cells_used_as_biological_replicates": False', source)


if __name__ == "__main__":
    unittest.main()
