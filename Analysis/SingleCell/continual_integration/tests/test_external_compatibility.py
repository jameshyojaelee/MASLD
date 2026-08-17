from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from masld_cl.external_compatibility import (
    ExternalCompatibilityError,
    gene_mapping_diagnostics,
    read_label_map,
    resolve_gene_indices,
)


class TestExternalCompatibility(unittest.TestCase):
    def setUp(self):
        self.root = Path(__file__).resolve().parents[1]

    def test_production_label_map_uses_only_frozen_vocabulary(self):
        mapping = read_label_map(self.root / "external_label_map_v1.tsv")
        self.assertEqual(len(mapping), 48)
        self.assertEqual(mapping["Cycling"]["include_reference"], "False")
        self.assertEqual(mapping["DCmac"]["include_reference"], "False")

    def test_gene_resolution_prefers_id_then_unique_name(self):
        indices, summary = resolve_gene_indices(
            ["ENSG1", "B"], ["ENSG1", "ENSG2"], ["A", "B"]
        )
        self.assertEqual(indices, [0, 1])
        self.assertEqual(
            summary["mapping_methods"],
            {"external_feature_name": 1, "external_var_id": 1},
        )

    def test_missing_and_ambiguous_genes_fail(self):
        with self.assertRaisesRegex(ExternalCompatibilityError, "missing"):
            resolve_gene_indices(["C"], ["ENSG1"], ["A"])
        with self.assertRaisesRegex(ExternalCompatibilityError, "ambiguous"):
            resolve_gene_indices(["A"], ["ENSG1", "ENSG2"], ["A", "A"])
        diagnostics = gene_mapping_diagnostics(
            ["A", "MISSING"], ["ENSG1", "ENSG2"], ["A", "A"]
        )
        self.assertFalse(diagnostics["eligible"])
        self.assertEqual(diagnostics["ambiguous_genes"], ["A"])
        self.assertEqual(diagnostics["missing_genes"], ["MISSING"])

    def test_required_field_set_operations_accept_pandas_index(self):
        required = {"STUDY", "donor_id"}
        self.assertEqual(required - set(pd.Index(["STUDY"])), {"donor_id"})

    def test_invalid_label_map_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "map.tsv"
            path.write_text(
                "author_cell_type\tfrozen_broad_label\tinclude_reference\trationale\n"
                "X\tInvented\tTrue\tbad\n"
            )
            with self.assertRaisesRegex(ExternalCompatibilityError, "non-frozen"):
                read_label_map(path)
