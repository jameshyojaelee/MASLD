import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("summarize_candidate", Path(__file__).resolve().parents[1] / "scripts/summarize_candidate.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class FamilyTests(unittest.TestCase):
    def test_missing_methods_stay_in_family(self):
        family = list(module.planned_comparisons())
        self.assertEqual(len(family), 48)
        self.assertEqual(len(set(family)), len(family))
        self.assertTrue(all(row[2] != "rna" for row in family if row[0] != "GSE296875"))

    def test_holm_known_values(self):
        self.assertEqual(module.holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])
        self.assertEqual(module.holm([0.001] + [1.0] * 47)[0], 0.048)

    def test_pilot_status_comes_from_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "completion.json").write_text(json.dumps({"interpretation": "technical_pilot_no_research_calls"}))
            (root / "run_settings.json").write_text(json.dumps({"pilot": False}))
            self.assertTrue(module.source_is_technical(root))
            (root / "completion.json").write_text(json.dumps({"interpretation": "internal_prediction"}))
            self.assertFalse(module.source_is_technical(root))
            (root / "run_settings.json").write_text(json.dumps({"pilot": False, "analysis_scope": "bounded_gene_development"}))
            self.assertTrue(module.source_is_technical(root))
            (root / "run_settings.json").write_text(json.dumps({"pilot": True}))
            self.assertTrue(module.source_is_technical(root))

    def test_invalidated_results_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "run_summary.json").write_text("{}")
            (root / "completion.json").write_text('{"status":"completed"}')
            (root / "INVALIDATED_EDGE_IDENTITY.json").write_text('{"reason":"wrong peak identity"}')
            with patch("sys.argv", ["summary", "--da", str(root), "--network", str(root), "--output", str(root / "output")]):
                with self.assertRaisesRegex(ValueError, "invalidated"):
                    module.main()


if __name__ == "__main__":
    unittest.main()
