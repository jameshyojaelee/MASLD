from __future__ import annotations

import gzip
from hashlib import sha256
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy import sparse

from scripts.materialize_peakvi_cross_cohort_cells import LINEAGES, WindowIndex
from scripts.materialize_peakvi_gse244832_target_cells import (
    DATASET_ID,
    PeakVICellMaterializationError,
    process_fragment_file,
    select_cells,
)


class MaterializePeakVIGSE244832TargetCellsTests(unittest.TestCase):
    def _fixture(self) -> tuple[list[dict[str, str]], list[dict[str, str]], np.ndarray]:
        units: list[dict[str, str]] = []
        membership: list[dict[str, str]] = []
        states = np.full((18, 4), 3, dtype=np.uint8)
        for offset in range(72):
            donor_index, lineage_index = divmod(offset, 4)
            donor = f"D{donor_index + 1:02d}"
            lineage = LINEAGES[lineage_index]
            state = 0 if offset < 48 else (1 if offset < 52 else 3)
            cells = 50 if state == 0 else (0 if state == 1 else 4)
            states[donor_index, lineage_index] = state
            units.append(
                {
                    "donor_index": str(donor_index),
                    "donor_id": donor,
                    "lineage_index": str(lineage_index),
                    "lineage_id": lineage,
                    "outer_fold": str(donor_index % 5),
                    "cells": str(cells),
                    "evidence_state": (
                        "observed" if state == 0 else "structurally_missing" if state == 1 else "below_qc"
                    ),
                    "missing_state_code": str(state),
                    "eligible_min_50_cells": "true" if state == 0 else "false",
                }
            )
            for cell in range(cells):
                barcode = (
                    "RECURRENT"
                    if cell == 0 and lineage_index == 0
                    else f"BC-{lineage_index}-{cell:03d}"
                )
                membership.append(
                    {
                        "dataset_id": DATASET_ID,
                        "donor_id": donor,
                        "raw_barcode": barcode,
                        "lineage_id": lineage,
                        "analysis_role": "primary",
                        "outer_fold": str(donor_index % 5),
                    }
                )
        return units, membership, states

    def test_registered_mask_is_applied_before_cell_selection(self) -> None:
        units, membership, states = self._fixture()
        selected = select_cells(membership, units, states, cap=64, seed=20260825)
        self.assertEqual(len(selected), 48 * 50)
        self.assertEqual(
            {(row["donor_id"], row["lineage_id"]) for row in selected},
            {
                (row["donor_id"], row["lineage_id"])
                for offset, row in enumerate(units)
                if states.reshape(-1)[offset] == 0
            },
        )
        self.assertEqual(len({row["cell_hash"] for row in selected}), len(selected))

    def test_recurrent_barcodes_are_globally_namespaced(self) -> None:
        units, membership, states = self._fixture()
        selected = select_cells(membership, units, states, cap=64, seed=20260825)
        recurrent = [row for row in selected if row["raw_barcode"] == "RECURRENT"]
        self.assertGreater(len(recurrent), 1)
        self.assertEqual(len({row["cell_hash"] for row in recurrent}), len(recurrent))

    def test_selection_rejects_ineligible_unit_count_drift(self) -> None:
        units, membership, states = self._fixture()
        units[-1]["cells"] = "5"
        with self.assertRaises(PeakVICellMaterializationError):
            select_cells(membership, units, states, cap=64, seed=20260825)

    def test_custom_bowtie_fragment_uses_plus4_minus5_geometry(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fragment = root / "D01_fragments.tsv.gz"
            with gzip.open(fragment, "wt", encoding="utf-8", newline="") as handle:
                handle.write("chr1\t100\t120\tAA\t1\n")
            expected_hash = sha256(fragment.read_bytes()).hexdigest()
            index = WindowIndex(
                rows=(("chr1", 104, 105, "left"), ("chr1", 115, 116, "right")),
                by_contig={"chr1": ((104, 115), (105, 116), (0, 1))},
            )
            output = root / "shard"
            receipt = process_fragment_file(
                donor_id="D01",
                path=str(fragment),
                expected_size=fragment.stat().st_size,
                expected_sha256=expected_hash,
                selected_lookup={"AA": 0},
                windows={"valid": index, "test": index},
                reference_lengths={"chr1": 1000},
                output=str(output),
                n_rows=1,
            )
            matrix = sparse.load_npz(output / "valid.counts.npz").tocsr()
            self.assertEqual(matrix.shape, (1, 16000))
            self.assertEqual(matrix[0, 0], 1)
            self.assertEqual(matrix[0, 1], 1)
            self.assertEqual(int(matrix.sum()), 2)
            self.assertEqual(receipt["fragment_cut_sites"], "start_plus_4_and_end_minus_5")


if __name__ == "__main__":
    unittest.main()
