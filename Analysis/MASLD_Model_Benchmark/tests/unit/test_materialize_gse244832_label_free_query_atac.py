from __future__ import annotations

import gzip
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from masld_bench.artifacts import verify_frozen_tree
from scripts.materialize_gse244832_label_free_query_atac import (
    QueryATACMaterializationError,
    WindowIndex,
    process_donor,
    tn5_positions,
    validate_config,
    write_role_output,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse244832_label_free_query_atac_materialization_20260825.json"


class MaterializeGSE244832LabelFreeQueryATACTests(unittest.TestCase):
    def _windows(self) -> dict[str, WindowIndex]:
        return {
            "valid": WindowIndex(
                (("chr1", 104, 1104, "valid_w1"),),
                {"chr1": ((104,), (1104,), (0,))},
            ),
            "test": WindowIndex(
                (("chr2", 104, 1104, "test_w1"),),
                {"chr2": ((104,), (1104,), (0,))},
            ),
        }

    def test_real_config_binds_guarded_axis_and_no_outcome_input(self) -> None:
        config = json.loads(CONFIG.read_text(encoding="utf-8"))
        validate_config(ROOT, config)
        self.assertFalse(config["evaluator_outcome_artifact_read"])
        self.assertFalse(config["rna_assay_opened"])
        self.assertFalse(config["sequence_execution_authorized"])

    def test_custom_bowtie2_tn5_geometry(self) -> None:
        self.assertEqual(tn5_positions(100, 120), (104, 115))
        with self.assertRaises(QueryATACMaterializationError):
            tn5_positions(100, 109)

    def test_each_output_contains_only_its_context_contig_values(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fragments.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("chr1\t100\t120\tcell_a\t2\n")
                handle.write("chr2\t100\t120\tcell_a\t3\n")
            result = process_donor(
                "D01",
                path,
                path.stat().st_size,
                {"cell_a": "hepatocyte"},
                {"cell_a"},
                self._windows(),
                {"chr1": 2000, "chr2": 2000},
            )
            self.assertEqual(int(result["counts"]["valid"].sum()), 2)
            self.assertEqual(int(result["counts"]["test"].sum()), 2)

    def test_output_freezes_observed_atac_only_contract(self) -> None:
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
                "fragment_sha256_by_donor": {f"D{index:02d}": str(index) * 64 for index in range(1, 19)},
                "condition_or_phenotype_authority_included": False,
                "rna_or_cross_assay_authority_included": False,
                "evaluator_outcome_authority_included": False,
                "sequence_execution_authorized": False,
            }
            write_role_output(
                output,
                "valid",
                self._windows()["valid"],
                np.zeros((18, 4, 1, 20), dtype=np.uint32),
                np.zeros((18, 4, 1), dtype=np.uint32),
                [
                    {
                        "donor_id": f"D{index:02d}",
                        "selected_fragment_rows": 0,
                        "primary_barcodes": 1,
                        "all_membership_barcodes": 1,
                        "membership_barcodes_seen": 1,
                        "gzip_crc_verified_to_eof": True,
                    }
                    for index in range(1, 19)
                ],
                authorities,
            )
            metadata = verify_frozen_tree(output)["metadata"]
            contract = json.loads((output / "query_contract.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["artifact_class"], "gse244832_label_free_query_atac")
            self.assertFalse(metadata["sequence_execution_authorized"])
            self.assertFalse(contract["rna_assay_opened"])
            self.assertNotIn("processing_qc.tsv", contract["model_input_files"])


if __name__ == "__main__":
    unittest.main()
