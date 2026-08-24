from __future__ import annotations

import csv
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

from scipy import sparse


MODULE = (
    Path(__file__).parents[2] / "scripts" / "scbasset_build_fold_inputs.py"
)
SPEC = importlib.util.spec_from_file_location("scbasset_fold_inputs_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
scbasset = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(scbasset)


class ScBassetFoldInputTest(unittest.TestCase):
    def test_production_wrapper_leaves_builder_output_absent(self) -> None:
        wrapper = (
            Path(__file__).parents[2]
            / "slurm"
            / "scbasset_build_donor0_genomic0_inputs.sbatch"
        ).read_text(encoding="utf-8")
        install_block = wrapper.split("install -d -m 0750", 1)[1].split(
            "module load", 1
        )[0]
        self.assertNotIn('"${STAGE}/inputs"', install_block)
        self.assertIn('--output "${STAGE}/inputs"', wrapper)
        self.assertIn('mv "${STAGE}" "${failed}"', wrapper)
        self.assertNotIn("rm -rf", wrapper)

    def test_held_region_identifier_preserves_canonical_ccre_join_key(self) -> None:
        self.assertEqual(
            scbasset._canonical_region_id(
                name="fixed_valid|ccre_000123|PLS",
                role="valid",
                contig="chr2",
                center=1500,
            ),
            "ccre_000123",
        )
        self.assertEqual(
            scbasset._canonical_region_id(
                name="train_macs2|peak_17",
                role="train",
                contig="chr3",
                center=500,
            ),
            "chr3:500:train_macs2|peak_17",
        )
        with self.assertRaisesRegex(
            scbasset.ScBassetInputError, "fixed cCRE identifier"
        ):
            scbasset._canonical_region_id(
                name="fixed_test|not_canonical|PLS",
                role="test",
                contig="chr1",
                center=1500,
            )

    def test_fragment_overlap_is_binary_and_excludes_genomic_test(self) -> None:
        regions = [
            {
                "region_id": "train_a",
                "block_id": "chr1",
                "role": "train",
                "contig": "chr1",
                "target_start": 100,
                "target_end": 200,
                "sequence_start": 0,
                "sequence_end": 1344,
            },
            {
                "region_id": "valid_a",
                "block_id": "chr2",
                "role": "valid",
                "contig": "chr2",
                "target_start": 300,
                "target_end": 400,
                "sequence_start": 0,
                "sequence_end": 1344,
            },
            {
                "region_id": "test_a",
                "block_id": "chr3",
                "role": "test",
                "contig": "chr3",
                "target_start": 500,
                "target_end": 600,
                "sequence_start": 0,
                "sequence_end": 1344,
            },
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fragment = root / "hepatocyte.fragments.tsv.gz"
            with gzip.open(fragment, "wt", encoding="utf-8", newline="") as handle:
                handle.write("chr1\t90\t110\tcell_a\t99\n")
                handle.write("chr1\t150\t180\tcell_a\t7\n")
                handle.write("chr1\t190\t210\tcell_b\t1\n")
                handle.write("chr2\t350\t360\tcell_b\t10\n")
                handle.write("chr3\t550\t560\tcell_a\t3\n")
            output = root / "hepatocyte.npz"
            result = scbasset._build_lineage_matrix(
                {
                    "lineage": "hepatocyte",
                    "source": str(fragment),
                    "output": str(output),
                    "cells": ["cell_a", "cell_b"],
                    "interval_index": scbasset._interval_index(regions),
                    "test_contigs": ["chr3"],
                    "n_regions": 3,
                }
            )
            matrix = sparse.load_npz(output).toarray()
            self.assertEqual(matrix.tolist(), [[1, 1], [0, 1], [0, 0]])
            self.assertEqual(result["binary_accessible_pairs"], 3)
            self.assertEqual(
                result["genomic_test_fragments_parsed_and_excluded"], 1
            )
            self.assertFalse(result["read_support_used_as_signal_weight"])

    def test_training_roster_excludes_both_nested_held_folds(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            membership = Path(directory) / "membership.tsv.gz"
            with gzip.open(membership, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=scbasset.MEMBERSHIP_FIELDS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for fold in range(5):
                    for lineage in scbasset.LINEAGES:
                        writer.writerow(
                            {
                                "well_id": "well1",
                                "raw_barcode": f"raw_{fold}_{lineage}",
                                "cell_id": f"cell_{fold}_{lineage}",
                                "donor_id": f"donor_{fold}",
                                "source_label": lineage,
                                "lineage_id": lineage,
                                "analysis_role": "primary",
                                "outer_fold": fold,
                            }
                        )
            rows, by_lineage = scbasset._training_roster(membership, (2, 3, 4))
            self.assertEqual({int(row["outer_fold"]) for row in rows}, {2, 3, 4})
            self.assertEqual(len(rows), 15)
            self.assertTrue(all(len(by_lineage[lineage]) == 3 for lineage in scbasset.LINEAGES))

    def test_intermediate_full_axis_matrices_are_removed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            lineage_matrices = output / "lineage_matrices"
            lineage_matrices.mkdir()
            for lineage in scbasset.LINEAGES:
                (lineage_matrices / f"{lineage}.npz").write_bytes(b"matrix")
                (lineage_matrices / f"{lineage}.json").write_text(
                    "{}\n", encoding="utf-8"
                )
            scbasset._remove_intermediate_lineage_matrices(output)
            self.assertFalse(lineage_matrices.exists())


if __name__ == "__main__":
    unittest.main()
