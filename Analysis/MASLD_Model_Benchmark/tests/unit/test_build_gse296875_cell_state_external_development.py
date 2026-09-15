from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse

from scripts.build_gse296875_cell_state_external_development import (
    CLASSES,
    ExternalCellStateError,
    SOURCE_TO_BROAD,
    identifier,
    select_records,
    write_feature_h5,
)


class ExternalCellStateViewTests(unittest.TestCase):
    def joined(self) -> list[dict[str, str]]:
        rows: list[dict[str, str]] = []
        for donor in ("d1", "d2", "d3"):
            for source_label in SOURCE_TO_BROAD:
                for replicate in range(3):
                    rows.append(
                        {
                            "cell_id": f"{donor}:{source_label}:{replicate}",
                            "donor_id": donor,
                            "well_id": "well1",
                            "raw_barcode": f"ACGT{replicate}-1",
                            "source_label": source_label,
                        }
                    )
        return rows

    def test_selection_is_deterministic_balanced_and_donor_complete(self) -> None:
        first = select_records(self.joined(), cells_per_class=6, selection_seed=17)
        second = select_records(self.joined(), cells_per_class=6, selection_seed=17)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 30)
        self.assertEqual(len({row["cell_id"] for row in first}), 30)
        for broad in CLASSES:
            rows = [row for row in first if row["broad_label"] == broad]
            self.assertEqual(len(rows), 6)
            self.assertEqual({row["donor_id"] for row in rows}, {"d1", "d2", "d3"})
        immune_sources = {
            row["source_label"] for row in first if row["broad_label"] == "immune"
        }
        self.assertTrue(immune_sources.issubset({"B cells", "Kupffer", "NK-T"}))

    def test_selection_fails_when_one_class_is_absent_from_one_donor(self) -> None:
        joined = [
            row
            for row in self.joined()
            if not (row["donor_id"] == "d3" and row["source_label"] == "LSEC")
        ]
        with self.assertRaisesRegex(ExternalCellStateError, "not observed in every donor"):
            select_records(joined, cells_per_class=3)

    def test_feature_artifact_excludes_labels_and_donor_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "features.h5"
            matrix = sparse.csr_matrix(np.asarray([[1, 0, 2], [0, 3, 0]], dtype=np.int32))
            row_ids = [identifier("test", "a"), identifier("test", "b")]
            identity = write_feature_h5(
                path,
                row_ids=row_ids,
                rna=matrix,
                gene_ids=["ENSG1", "ENSG2", "ENSG3"],
                gene_names=["A", "B", "C"],
            )
            self.assertEqual(identity["shape"], [2, 3])
            with h5py.File(path, "r") as handle:
                self.assertEqual(set(handle["obs"].keys()), {
                    "atac_observed_mask", "rna_observed_mask", "rna_status", "row_id"
                })
                self.assertEqual(handle.attrs["label_visibility"], "evaluator_only")
                self.assertTrue(np.all(handle["obs/atac_observed_mask"][:] == 0))


if __name__ == "__main__":
    unittest.main()
