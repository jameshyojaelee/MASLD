from __future__ import annotations

from copy import deepcopy
import inspect
import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np
from scipy import sparse

from masld_bench.observed_multiome_materialization import (
    ObservedMultiomeMaterializationError,
    donor_lineage_aggregation,
    selectively_read_csr_columns,
)
from scripts import materialize_gse296875_observed_multiome_model_inputs as materializer


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_model_input_materialization_20260825.json"


class ObservedMultiomeModelInputMaterializationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_config_passes_and_firewall_mutation_fails(self) -> None:
        resolved = materializer.validate_config(ROOT, deepcopy(self.config))
        self.assertEqual(set(resolved), {"materialization_contract", "mask_plan", "input_h5"})
        mutated = deepcopy(self.config)
        mutated["firewall"]["target_atac_values_read"] = True
        with self.assertRaises(materializer.ModelInputMaterializationError):
            materializer.validate_config(ROOT, mutated)

    def test_selective_sparse_read_only_returns_requested_columns(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "matrix.h5"
            with h5py.File(path, "w") as handle:
                group = handle.create_group("csr")
                group.create_dataset("data", data=np.asarray([5, 6, 7, 8, 9]))
                group.create_dataset("indices", data=np.asarray([0, 2, 1, 3, 2]))
                group.create_dataset("indptr", data=np.asarray([0, 2, 4, 5]))
                group.create_dataset("shape", data=np.asarray([3, 4]))
            with h5py.File(path, "r") as handle:
                matrix, positions = selectively_read_csr_columns(
                    handle["csr"], [2, 0]
                )
            np.testing.assert_array_equal(positions, np.asarray([0, 1, 4]))
            np.testing.assert_array_equal(
                matrix.toarray(), np.asarray([[6, 5], [0, 0], [9, 0]])
            )

    def test_donor_lineage_aggregation_rejects_cell_census_drift(self) -> None:
        planned = [
            {"donor_hash": "d1", "lineage": "Hepatocyte", "nuclei": "2"},
            {"donor_hash": "d2", "lineage": "Stellate", "nuclei": "1"},
        ]
        aggregation, keys = donor_lineage_aggregation(
            cell_unit_keys=[
                ("d1", "Hepatocyte"),
                ("d2", "Stellate"),
                ("d1", "Hepatocyte"),
            ],
            planned_units=planned,
        )
        np.testing.assert_array_equal(aggregation.sum(axis=1), [[2], [1]])
        self.assertEqual(keys, ["d1\0Hepatocyte", "d2\0Stellate"])
        with self.assertRaises(ObservedMultiomeMaterializationError):
            donor_lineage_aggregation(
                cell_unit_keys=[("d1", "Hepatocyte"), ("d2", "Stellate")],
                planned_units=planned,
            )

    def test_written_fold_has_no_target_values_or_raw_donor_field(self) -> None:
        planned = [
            {
                "donor_hash": "a" * 64,
                "lineage": "Hepatocyte",
                "donor_fold": "0",
                "nuclei": "2",
            }
        ]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fold.h5"
            materializer._write_fold_h5(
                path=path,
                held_fold=0,
                planned_units=planned,
                row_hashes=["b" * 64],
                rna_counts=sparse.csr_matrix([[2, 3]]),
                ensembl_ids=["ENSG1", "ENSG2"],
                gene_names=["G1", "G2"],
                atac_counts=sparse.csr_matrix([[4, 5]]),
                input_rows=[
                    {
                        "input_hash": "c" * 64,
                        "chromosome": "chr2",
                        "bed_start_0based": "10",
                        "bed_end_half_open": "20",
                    },
                    {
                        "input_hash": "d" * 64,
                        "chromosome": "chr2",
                        "bed_start_0based": "30",
                        "bed_end_half_open": "40",
                    },
                ],
                target_rows=[
                    {
                        "target_hash": "e" * 64,
                        "chromosome": "chr1",
                        "bed_start_0based": "50",
                        "bed_end_half_open": "60",
                    }
                ],
            )
            with h5py.File(path, "r") as handle:
                self.assertFalse(bool(handle.attrs["target_atac_values_present"]))
                self.assertIn("unit_hash", handle["rows"])
                self.assertNotIn("donor_id", handle["rows"])
                self.assertNotIn("donor_hash", handle["rows"])
                self.assertEqual(
                    set(handle["targets_without_values"]),
                    {"target_hash", "chromosome", "bed_start_0based", "bed_end_half_open"},
                )
                self.assertFalse(
                    bool(handle["targets_without_values"].attrs["atac_values_present"])
                )

    def test_source_never_binds_full_atac_values(self) -> None:
        source = inspect.getsource(materializer)
        self.assertEqual(source.count("read_full_csr("), 1)
        self.assertIn('read_full_csr(handle["rna/counts_csr"])', source)
        self.assertNotIn('read_full_csr(handle["atac', source)
        self.assertNotIn('atac_group["data"][:]', source)
        self.assertNotIn("target_atac_values=", source)


if __name__ == "__main__":
    unittest.main()
