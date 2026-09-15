"""The family must be derived on one quantification, and the ragged read asserted.

The property worth pinning is that pooling across quantifications is refused
rather than caveated: the rule evaluated on the two kallisto arms and the rule
evaluated on all three give different families, and the code must derive the
first. A test that only checked "a family was produced" would pass either way.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "build_composition_substrate.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cs = _load(MODULE, "composition_substrate_tested")


def _write(text: str) -> Path:
    handle = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
    handle.write(text)
    handle.close()
    return Path(handle.name)


class RaggedReadTests(unittest.TestCase):
    def test_the_leading_unnamed_key_is_taken_as_the_sample(self) -> None:
        path = _write("A\tB\nS1\t0.6\t0.4\nS2\t0.7\t0.3\n")
        header, keys, values, audit = cs.read_ragged(path)
        self.assertEqual(header, ["A", "B"])
        self.assertEqual(keys, ["S1", "S2"])
        self.assertEqual(audit["measured_offset"], 1)
        self.assertEqual(values.shape, (2, 2))

    def test_an_aligned_table_is_refused_not_silently_reinterpreted(self) -> None:
        path = _write("A\tB\n0.6\t0.4\n")
        with self.assertRaises(cs.CompositionError):
            cs.read_ragged(path)

    def test_rows_that_do_not_sum_to_one_are_refused(self) -> None:
        path = _write("A\tB\nS1\t0.6\t0.1\n")
        with self.assertRaises(cs.CompositionError):
            cs.read_ragged(path)

    def test_a_row_sum_off_by_float_noise_is_accepted(self) -> None:
        """Tolerance, not equality - most real rows are not bitwise 1.0."""

        path = _write("A\tB\nS1\t0.6000000000000001\t0.3999999999999999\n")
        _, _, _, audit = cs.read_ragged(path)
        self.assertLess(audit["max_abs_deviation_from_1"], cs.ROW_SUM_TOLERANCE)


class FamilyDerivationTests(unittest.TestCase):
    @staticmethod
    def _tables():
        return {
            "GSE135251": {"Hep": 1.00, "Mac": 1.00, "Mono": 1.00, "T": 0.99, "Fib": 0.23},
            "GSE130970": {"Hep": 1.00, "Mac": 1.00, "Mono": 0.97, "T": 0.94, "Fib": 1.00},
            "GSE193066": {"Hep": 1.00, "Mac": 1.00, "Mono": 0.16, "T": 0.45, "Fib": 0.99},
        }

    def test_the_kallisto_family_differs_from_the_pooled_family(self) -> None:
        tables = self._tables()
        kallisto = cs.family_from(tables, cs.KALLISTO_ARMS)
        pooled = cs.family_from(tables, cs.ALL_ARMS)
        self.assertEqual(kallisto, ["Hep", "Mac", "Mono", "T"])
        self.assertEqual(pooled, ["Hep", "Mac"])
        self.assertNotEqual(kallisto, pooled)

    def test_the_pooled_rule_loses_lineages_to_the_star_arm(self) -> None:
        tables = self._tables()
        lost = set(cs.family_from(tables, cs.KALLISTO_ARMS)) - set(
            cs.family_from(tables, cs.ALL_ARMS))
        self.assertEqual(lost, {"Mono", "T"})

    def test_a_lineage_failing_either_kallisto_arm_is_excluded(self) -> None:
        tables = self._tables()
        self.assertNotIn("Fib", cs.family_from(tables, cs.KALLISTO_ARMS))

    def test_the_threshold_is_the_declared_ninety_percent(self) -> None:
        self.assertEqual(cs.SAMPLE_FRACTION, 0.90)
        self.assertEqual(cs.DETECTION_FLOOR, 1e-4)

    def test_a_lineage_exactly_at_ninety_percent_is_included(self) -> None:
        tables = {a: {"X": 0.90} for a in cs.ALL_ARMS}
        self.assertEqual(cs.family_from(tables, cs.KALLISTO_ARMS), ["X"])


class QuantificationTests(unittest.TestCase):
    def test_the_arms_are_declared_non_uniform(self) -> None:
        self.assertEqual(cs.QUANTIFICATION["GSE193066"], "STAR")
        self.assertEqual(cs.QUANTIFICATION["GSE135251"], "kallisto")
        self.assertEqual(cs.QUANTIFICATION["GSE130970"], "kallisto")
        self.assertNotIn("GSE193066", cs.KALLISTO_ARMS)

    def test_the_superseded_five_are_recorded_for_comparison(self) -> None:
        self.assertEqual(len(cs.SUPERSEDED_FIVE), 5)
        self.assertIn("Plasma cells", cs.SUPERSEDED_FIVE)
        self.assertIn("cDC1s", cs.SUPERSEDED_FIVE)


if __name__ == "__main__":
    unittest.main()
