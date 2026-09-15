from __future__ import annotations

import gzip
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from masld_bench.artifacts import verify_frozen_tree

from scripts.materialize_gse281367_label_free_query_atac import (
    QueryATACMaterializationError,
    WindowIndex,
    process_donor,
    validate_config,
    windows_at,
    write_role_output,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse281367_label_free_query_atac_materialization_20260825.json"


class MaterializeGSE281367LabelFreeQueryATACTests(unittest.TestCase):
    def _windows(self) -> dict[str, WindowIndex]:
        return {
            "valid": WindowIndex(
                (("chr1", 100, 1100, "valid_w1"),),
                {"chr1": ((100,), (1100,), (0,))},
            ),
            "test": WindowIndex(
                (("chr2", 100, 1100, "test_w1"),),
                {"chr2": ((100,), (1100,), (0,))},
            ),
        }

    def test_real_config_binds_label_free_axis_and_no_outcome_input(self) -> None:
        import json

        config = json.loads(CONFIG.read_text(encoding="utf-8"))
        validate_config(ROOT, config)
        self.assertFalse(config["evaluator_outcome_artifact_read"])
        self.assertFalse(config["condition_values_read"])
        self.assertFalse(config["firewall"]["champion_claim_allowed"])

    def test_half_open_window_lookup(self) -> None:
        index = self._windows()["valid"]
        self.assertEqual(list(windows_at(index, "chr1", 100)), [0])
        self.assertEqual(list(windows_at(index, "chr1", 1099)), [0])
        self.assertEqual(list(windows_at(index, "chr1", 1100)), [])

    def test_each_output_contains_only_its_context_contig_values(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fragments.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("chr1\t100\t120\tcell_a\t2\n")
                handle.write("chr2\t200\t220\tcell_a\t3\n")
            result = process_donor(
                "Z01",
                path,
                path.stat().st_size,
                {"cell_a": "hepatocyte"},
                {"cell_a"},
                self._windows(),
                {"chr1": 2000, "chr2": 2000},
            )
            self.assertEqual(int(result["counts"]["valid"].sum()), 2)
            self.assertEqual(int(result["counts"]["test"].sum()), 2)
            self.assertEqual(int(result["counts"]["valid"][:, 1:].sum()), 0)
            self.assertEqual(int(result["counts"]["test"][:, 1:].sum()), 0)

    def test_primary_barcode_absence_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fragments.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("chr1\t100\t120\tcell_a\t1\n")
            with self.assertRaises(QueryATACMaterializationError):
                process_donor(
                    "Z01",
                    path,
                    path.stat().st_size,
                    {"cell_a": "hepatocyte", "cell_b": "hepatocyte"},
                    {"cell_a", "cell_b"},
                    self._windows(),
                    {"chr1": 2000, "chr2": 2000},
                )

    def test_output_freezes_explicit_axes_and_input_authorities(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "query"
            authorities = {
                "campaign_id": "fixture",
                "campaign_sha256": "a" * 64,
                "axis_artifacts_sha256": "b" * 64,
                "source_artifacts_sha256": "c" * 64,
                "source_inventory_sha256": "d" * 64,
                "reference_fai_sha256": "e" * 64,
                "builder_sha256": "f" * 64,
                "fragment_sha256_by_donor": {f"Z{index:02d}": str(index) * 64 for index in range(1, 13)},
                "condition_or_phenotype_authority_included": False,
                "evaluator_outcome_authority_included": False,
            }
            write_role_output(
                output,
                "valid",
                self._windows()["valid"],
                np.zeros((12, 4, 1, 20), dtype=np.uint32),
                np.zeros((12, 4, 1), dtype=np.uint32),
                [
                    {
                        "donor_id": f"Z{index:02d}",
                        "selected_fragment_rows": 0,
                        "primary_barcodes": 1,
                        "all_membership_barcodes": 1,
                        "gzip_crc_verified_to_eof": True,
                    }
                    for index in range(1, 13)
                ],
                authorities,
            )
            manifest = verify_frozen_tree(output)
            contract = json.loads((output / "query_contract.json").read_text(encoding="utf-8"))
            observed_authorities = json.loads((output / "input_authorities.json").read_text(encoding="utf-8"))
            axis_lines = (output / "query_axis.tsv").read_text(encoding="utf-8").splitlines()
            self.assertEqual(manifest["metadata"]["artifact_class"], "gse281367_label_free_query_atac")
            self.assertEqual(observed_authorities, authorities)
            self.assertEqual(len(axis_lines), 49)
            self.assertEqual(contract["audit_only_files"], ["processing_qc.tsv"])
            self.assertNotIn("processing_qc.tsv", contract["model_input_files"])


if __name__ == "__main__":
    unittest.main()
