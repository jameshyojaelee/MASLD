"""Structural tests for the donor-by-lineage pseudobulk fixture builder.

These run on the scanpy environment (Python 3.10) alongside the builder, and
deliberately do not import the benchmark package: the fitting and fixture path
must stay unable to reach the scoring module.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "build_pseudobulk", ROOT / "scripts/build_gse296875_donor_lineage_pseudobulk.py"
)
assert _spec is not None and _spec.loader is not None
builder = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(builder)


class FrozenConstantTests(unittest.TestCase):
    def test_partition_is_the_frozen_membership_partition(self) -> None:
        self.assertEqual(
            builder.PRIMARY_LINEAGES,
            ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell"),
        )
        self.assertEqual(
            builder.SECONDARY_LINEAGES, ("endothelial_cell", "b_cell")
        )
        self.assertEqual(len(builder.LINEAGES), 7)

    def test_threshold_and_cohort_shape(self) -> None:
        self.assertEqual(builder.MINIMUM_CELLS_PER_UNIT, 20)
        self.assertEqual(builder.EXPECTED_DONORS, 39)
        self.assertEqual(builder.EXPECTED_NUCLEI, 68_398)
        self.assertEqual(builder.EXPECTED_RNA_FEATURES, 36_601)

    def test_builder_imports_nothing_from_the_benchmark_package(self) -> None:
        """The fixture path must not be able to reach the scoring module.

        Checked on the import graph rather than on the file text, so a mention
        in prose does not pass or fail the assertion.
        """

        import ast

        tree = ast.parse(
            (ROOT / "scripts/build_gse296875_donor_lineage_pseudobulk.py").read_text()
        )
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])
        self.assertNotIn("masld_bench", imported)
        self.assertEqual(
            imported & {"sklearn", "masld_bench"}, set()
        )

    def test_builder_never_opens_an_outcome_table(self) -> None:
        source = (
            ROOT / "scripts/build_gse296875_donor_lineage_pseudobulk.py"
        ).read_text()
        for outcome_artifact in (
            "donor_endpoints",
            "endpoint_masks",
            "steatosis_numeric",
            "fibrosis_any",
        ):
            self.assertNotIn(outcome_artifact, source)


class BoundInputTests(unittest.TestCase):
    def test_bound_fixtures_still_verify(self) -> None:
        corgi = ROOT / "executions/corgi-gse296875-context-counts-21066278"
        lock = Path(
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis"
            "/SingleCell/continual_integration/reference/gse296875_source_lock_v44.json"
        )
        builder.verify_bound(ROOT, corgi, lock)

    def test_a_changed_fixture_is_rejected(self) -> None:
        corgi = ROOT / "executions/corgi-gse296875-context-counts-21066278"
        lock = Path(
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis"
            "/SingleCell/continual_integration/reference/gse296875_source_lock_v44.json"
        )
        original = builder.BOUND["donor_folds"]
        builder.BOUND["donor_folds"] = "0" * 64
        try:
            with self.assertRaises(builder.PseudobulkError):
                builder.verify_bound(ROOT, corgi, lock)
        finally:
            builder.BOUND["donor_folds"] = original


class AggregationTests(unittest.TestCase):
    """Aggregate a tiny synthetic well and check the sums land per unit."""

    def _write_h5(self, path: Path) -> None:
        # Three features, two of them Gene Expression; four barcodes.
        counts = np.array(
            [
                [1, 2, 3, 4],  # gene A
                [5, 6, 7, 8],  # gene B
                [9, 9, 9, 9],  # a Peaks feature, must be excluded
            ],
            dtype=np.int32,
        )
        from scipy import sparse

        matrix = sparse.csc_matrix(counts)
        with h5py.File(path, "w") as handle:
            group = handle.create_group("matrix")
            group.create_dataset("data", data=matrix.data)
            group.create_dataset("indices", data=matrix.indices)
            group.create_dataset("indptr", data=matrix.indptr)
            group.create_dataset("shape", data=np.array([3, 4], dtype=np.int64))
            group.create_dataset(
                "barcodes",
                data=np.array([b"BC1-1", b"BC2-1", b"BC3-1", b"BC4-1"]),
            )
            features = group.create_group("features")
            features.create_dataset(
                "feature_type",
                data=np.array([b"Gene Expression", b"Gene Expression", b"Peaks"]),
            )
            features.create_dataset(
                "id", data=np.array([b"ENSG1.1", b"ENSG2.1", b"PEAK1"])
            )
            features.create_dataset("name", data=np.array([b"A", b"B", b"P"]))

    def test_units_sum_their_own_nuclei_only(self) -> None:
        builder.EXPECTED_RNA_FEATURES = 2
        try:
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "well.h5"
                self._write_h5(path)
                rows = [
                    {"raw_barcode": "BC1-1", "cell_id": "c1", "donor_id": "D1", "lineage_id": "hepatocyte"},
                    {"raw_barcode": "BC2-1", "cell_id": "c2", "donor_id": "D1", "lineage_id": "hepatocyte"},
                    {"raw_barcode": "BC3-1", "cell_id": "c3", "donor_id": "D2", "lineage_id": "t_cell"},
                ]
                unit_index = {("D1", "hepatocyte"): 0, ("D2", "t_cell"): 1}
                counts, _ = builder.aggregate_well(path, rows, unit_index, 2)
                # Gene A: barcodes 1+2 = 1+2 = 3 for D1; barcode 3 = 3 for D2.
                # Gene B: 5+6 = 11 for D1; 7 for D2.  The Peaks row is excluded.
                self.assertEqual(counts.shape, (2, 2))
                np.testing.assert_array_equal(counts[0], np.array([3, 11]))
                np.testing.assert_array_equal(counts[1], np.array([3, 7]))
                # Barcode 4 belongs to no unit and contributes nowhere.
                self.assertEqual(int(counts.sum()), 3 + 11 + 3 + 7)
        finally:
            builder.EXPECTED_RNA_FEATURES = 36_601

    def test_an_unmapped_barcode_is_an_error_not_a_default_donor(self) -> None:
        builder.EXPECTED_RNA_FEATURES = 2
        try:
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "well.h5"
                self._write_h5(path)
                rows = [
                    {"raw_barcode": "MISSING-1", "cell_id": "c9", "donor_id": "D1", "lineage_id": "hepatocyte"}
                ]
                with self.assertRaises(builder.PseudobulkError):
                    builder.aggregate_well(path, rows, {("D1", "hepatocyte"): 0}, 1)
        finally:
            builder.EXPECTED_RNA_FEATURES = 36_601


class MaskSemanticTests(unittest.TestCase):
    def test_sub_threshold_units_are_missing_not_zero(self) -> None:
        source = (
            ROOT / "scripts/build_gse296875_donor_lineage_pseudobulk.py"
        ).read_text()
        self.assertIn('"insufficient_cells"', source)
        self.assertIn("masked_units_retain_observed_cell_count", source)
        self.assertIn('"missing_encoded_as_biological_zero": False', source)

    def test_threshold_is_a_greater_or_equal_comparison(self) -> None:
        for observed, expected in ((19, False), (20, True), (21, True), (0, False)):
            self.assertEqual(
                observed >= builder.MINIMUM_CELLS_PER_UNIT, expected
            )


if __name__ == "__main__":
    unittest.main()
