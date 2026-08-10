from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "14_prepare_guide_design_worklist.py"
sys.path.insert(0, str(SCRIPT.parent))
SPEC = importlib.util.spec_from_file_location("plan45_guide_worklist", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class GuideDesignTests(unittest.TestCase):
    def test_forward_allele_resolution(self):
        self.assertEqual(
            MODULE.forward_alleles("A", "A", "G", "G", False),
            ("G", "A", "direct_forward_match"),
        )
        self.assertEqual(
            MODULE.forward_alleles("T", "A", "G", "G", False),
            ("C", "T", "complemented_to_forward_match"),
        )
        self.assertEqual(
            MODULE.forward_alleles("A", "A", "T", "T", True),
            ("", "", "unresolved_palindromic_forward_strand"),
        )

    def test_indexed_fasta_across_line_boundaries(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "test.fa"
            fai = root / "test.fa.fai"
            fasta.write_text(">chr1\nACGTA\nCGTAC\nGT\n", encoding="ascii")
            fai.write_text("chr1\t12\t6\t5\t6\n", encoding="ascii")
            indexed = MODULE.IndexedFasta(fasta, fai)
            self.assertEqual(indexed.fetch("chr1", 4, 10), "TACGTAC")

    def test_base_editor_route(self):
        self.assertEqual(
            MODULE.base_editor_route("A", "G"),
            "ABE_candidate_reference_to_alternate",
        )
        self.assertEqual(
            MODULE.base_editor_route("G", "A"),
            "CBE_candidate_reference_to_alternate",
        )
        self.assertEqual(
            MODULE.base_editor_route("A", "C"),
            "prime_edit_required_reference_to_alternate",
        )

    def test_stage_a_mode_follows_oriented_target_expression(self):
        self.assertEqual(
            MODULE.stage_a_mode("risk_increases_expression"), "CRISPRa"
        )
        self.assertEqual(
            MODULE.stage_a_mode("risk_decreases_expression"), "CRISPRi"
        )
        with self.assertRaises(RuntimeError):
            MODULE.stage_a_mode("unresolved")


if __name__ == "__main__":
    unittest.main()
