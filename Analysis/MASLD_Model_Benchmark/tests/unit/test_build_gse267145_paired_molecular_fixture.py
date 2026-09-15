from __future__ import annotations

import csv
import gzip
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.build_gse267145_paired_molecular_fixture import (
    PairedMolecularFixtureError,
    build_fixture,
)


PARTICIPANTS = ("P01", "P02", "P03", "P04", "P05")


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_inputs(root: Path, *, integer_rna: bool = False, fractional_h3: bool = False) -> tuple[Path, Path, Path]:
    matrix = root / "matrix"
    (matrix / "raw").mkdir(parents=True)
    join = root / "participant_join.tsv"
    crosswalk = root / "rna_gene_crosswalk.tsv"
    join_fields = (
        "participant_id",
        "rna_gsm",
        "h3k27ac_gsm",
        "stage",
        "fibrosis",
        "sex",
        "pairing",
        "outer_fold",
        "rna_measurement_state",
        "h3k27ac_measurement_state",
    )
    write_tsv(
        join,
        join_fields,
        [
            {
                "participant_id": participant,
                "rna_gsm": f"RNA{index}",
                "h3k27ac_gsm": f"H3{index}",
                "stage": "hidden_value",
                "fibrosis": "hidden_value",
                "sex": "hidden_value",
                "pairing": "same_sample_different_aliquot",
                "outer_fold": str(index),
                "rna_measurement_state": "observed_unit_semantics_unresolved",
                "h3k27ac_measurement_state": "observed_integer_counts",
            }
            for index, participant in enumerate(PARTICIPANTS)
        ],
    )
    write_tsv(
        crosswalk,
        ("matrix_gene_id", "stable_gene_id", "mapping_state", "allowed_project_input"),
        [
            {
                "matrix_gene_id": f"ENSG{index}",
                "stable_gene_id": f"ENSG{index}",
                "mapping_state": (
                    "stable_id_exact_v98_and_v49" if index != 2 else "ensembl98_only_retired_or_absent_v49"
                ),
                "allowed_project_input": "true" if index != 2 else "false",
            }
            for index in range(4)
        ],
    )
    with gzip.open(matrix / "raw/GSE269412_RNA.txt.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("ensembl_gene_id", *PARTICIPANTS))
        for gene in range(4):
            values = [str(gene + participant + 1) for participant in range(5)]
            if gene == 1 and not integer_rna:
                values[3] = "2.75"
            writer.writerow((f"ENSG{gene}", *values))
    with gzip.open(matrix / "raw/GSE267119_H3K27ac.txt.gz", "wt", encoding="utf-8") as handle:
        labels = [f'"{participant}_hidden_H3K27ac_cutrun"' for participant in reversed(PARTICIPANTS)]
        handle.write(" ".join(labels) + "\n")
        for feature in ("chr1:1-3", "opaque_key", "chrQ:7-9"):
            values = [str(index + 1) for index in range(5)]
            if feature == "opaque_key" and fractional_h3:
                values[2] = "1.5"
            handle.write(f'"{feature}" ' + " ".join(values) + "\n")
    return matrix, join, crosswalk


class GSE267145PairedMolecularFixtureTests(unittest.TestCase):
    def build(self, root: Path, **input_options: bool) -> tuple[Path, dict]:
        matrix, join, crosswalk = make_inputs(root, **input_options)
        output = root / "fixture"
        receipt = build_fixture(
            matrix_root=matrix,
            join_path=join,
            crosswalk_path=crosswalk,
            output=output,
            expected_participants=5,
            expected_rna_source_features=4,
            expected_rna_allowed_features=3,
            expected_h3_features=3,
            expected_fold_counts={0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
        )
        return output, receipt

    def test_exact_values_axes_pairing_masks_and_folds(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output, receipt = self.build(Path(temporary))
            rna = np.load(output / "molecular/rna_values.npy", allow_pickle=False)
            h3 = np.load(output / "molecular/h3k27ac_counts.npy", allow_pickle=False)
            self.assertEqual(rna.shape, (5, 3))
            self.assertEqual(h3.shape, (5, 3))
            self.assertEqual(rna.dtype, np.float64)
            self.assertEqual(h3.dtype, np.uint32)
            self.assertEqual(rna[3, 1], 2.75)
            np.testing.assert_array_equal(h3[:, 0], np.array([5, 4, 3, 2, 1], dtype=np.uint32))
            self.assertTrue(np.load(output / "molecular/rna_observed_mask.npy", allow_pickle=False).all())
            self.assertTrue(np.load(output / "molecular/h3k27ac_observed_mask.npy", allow_pickle=False).all())
            keys = (output / "molecular/h3k27ac_feature_axis.tsv").read_text(encoding="utf-8")
            self.assertIn("chr1:1-3", keys)
            self.assertIn("opaque_key", keys)
            self.assertFalse(receipt["h3k27ac"]["coordinate_semantics_inferred"])
            folds = (output / "folds/participant_outer_folds.tsv").read_text(encoding="utf-8")
            self.assertNotIn("hidden_value", folds)

    def test_outcomes_and_outcome_bearing_source_labels_are_not_emitted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output, receipt = self.build(Path(temporary))
            participant_text = (output / "molecular/participant_axis.tsv").read_text(encoding="utf-8")
            self.assertNotIn("hidden_value", participant_text)
            self.assertNotIn("_hidden_H3K27ac", participant_text)
            self.assertNotIn("stage", participant_text)
            self.assertNotIn("fibrosis", participant_text)
            self.assertNotIn("sex", participant_text)
            self.assertFalse(receipt["feature_firewall"]["outcome_values_accessed"])
            self.assertFalse(receipt["feature_firewall"]["outcomes_used_for_feature_construction"])

    def test_integer_only_rna_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            matrix, join, crosswalk = make_inputs(Path(temporary), integer_rna=True)
            with self.assertRaises(PairedMolecularFixtureError):
                build_fixture(
                    matrix_root=matrix,
                    join_path=join,
                    crosswalk_path=crosswalk,
                    output=Path(temporary) / "fixture",
                    expected_participants=5,
                    expected_rna_source_features=4,
                    expected_rna_allowed_features=3,
                    expected_h3_features=3,
                    expected_fold_counts={0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
                )

    def test_fractional_h3_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            matrix, join, crosswalk = make_inputs(Path(temporary), fractional_h3=True)
            with self.assertRaises(PairedMolecularFixtureError):
                build_fixture(
                    matrix_root=matrix,
                    join_path=join,
                    crosswalk_path=crosswalk,
                    output=Path(temporary) / "fixture",
                    expected_participants=5,
                    expected_rna_source_features=4,
                    expected_rna_allowed_features=3,
                    expected_h3_features=3,
                    expected_fold_counts={0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
                )

    def test_participant_pairing_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            matrix, join, crosswalk = make_inputs(root)
            text = join.read_text(encoding="utf-8").replace("P05\tRNA4", "PX5\tRNA4")
            join.write_text(text, encoding="utf-8")
            with self.assertRaises(PairedMolecularFixtureError):
                build_fixture(
                    matrix_root=matrix,
                    join_path=join,
                    crosswalk_path=crosswalk,
                    output=root / "fixture",
                    expected_participants=5,
                    expected_rna_source_features=4,
                    expected_rna_allowed_features=3,
                    expected_h3_features=3,
                    expected_fold_counts={0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
                )


if __name__ == "__main__":
    unittest.main()
