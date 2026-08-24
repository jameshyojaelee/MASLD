from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np


MODULE = (
    Path(__file__).parents[2]
    / "src"
    / "masld_bench"
    / "adapters"
    / "scbasset_sequence_profile.py"
)
SPEC = importlib.util.spec_from_file_location("scbasset_sequence_profile_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
scbasset = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(scbasset)


def cells() -> list[dict[str, str]]:
    rows = []
    for fold in (2, 3, 4):
        for lineage in scbasset.LINEAGES:
            rows.append(
                {
                    "cell_id": f"fold{fold}_{lineage}",
                    "donor_id": f"donor_{fold}",
                    "lineage": lineage,
                    "outer_fold": str(fold),
                }
            )
    return rows


def regions() -> list[dict[str, str]]:
    return [
        {
            "region_id": "ccre_a",
            "block_id": "chr2",
            "role": "valid",
            "sequence_start": "1000",
            "sequence_end": "2344",
        },
        {
            "region_id": "ccre_b",
            "block_id": "chr8",
            "role": "valid",
            "sequence_start": "5000",
            "sequence_end": "6344",
        },
        {
            "region_id": "ccre_c",
            "block_id": "chr13",
            "role": "valid",
            "sequence_start": "9000",
            "sequence_end": "10344",
        },
    ]


class ScBassetSequenceProfileTest(unittest.TestCase):
    def test_aggregates_training_cells_and_normalizes_each_lineage(self) -> None:
        roster = cells()
        matrix = np.linspace(
            0.01, 0.99, num=3 * len(roster), dtype=np.float64
        ).reshape(3, len(roster))
        observed = scbasset.aggregate_training_cell_predictions(matrix, roster)
        self.assertEqual(observed.shape, (5, 3))
        np.testing.assert_allclose(observed.sum(axis=1), 1.0, atol=1e-12)
        self.assertTrue(np.all(observed > 0))

    def test_rejects_held_atac_and_held_cell_columns(self) -> None:
        roster = cells()
        matrix = np.full((3, len(roster)), 0.5, dtype=np.float32)
        with self.assertRaisesRegex(scbasset.ScBassetProfileError, "held-donor ATAC"):
            scbasset.aggregate_training_cell_predictions(
                matrix, roster, held_atac=np.ones((2, 3))
            )
        with self.assertRaisesRegex(scbasset.ScBassetProfileError, "held-cell"):
            scbasset.aggregate_training_cell_predictions(
                matrix, roster, held_cell_ids=["query_cell"]
            )

    def test_rejects_donor_or_genomic_partition_leakage(self) -> None:
        roster = cells()
        roster[0] = {**roster[0], "outer_fold": "1"}
        matrix = np.full((3, len(roster)), 0.5, dtype=np.float32)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(scbasset.ScBassetProfileError, "training cells use"):
                scbasset.export_prediction_bundle(
                    sequence_predictions=matrix,
                    training_cells=roster,
                    regions=regions(),
                    held_donors=[
                        {
                            "donor_id": "held_valid",
                            "outer_fold": "1",
                            "evaluation_role": "valid",
                        }
                    ],
                    output=Path(directory) / "out",
                    run_id="a" * 64,
                    namespace="fixture",
                    donor_test_fold=0,
                    donor_valid_fold=1,
                    model_artifact_sha256="b" * 64,
                )
        bad_regions = regions()
        bad_regions[0] = {**bad_regions[0], "role": "train"}
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(scbasset.ScBassetProfileError, "valid or test"):
                scbasset.export_prediction_bundle(
                    sequence_predictions=np.full(
                        (3, len(cells())), 0.5, dtype=np.float32
                    ),
                    training_cells=cells(),
                    regions=bad_regions,
                    held_donors=[
                        {
                            "donor_id": "held_valid",
                            "outer_fold": "1",
                            "evaluation_role": "valid",
                        }
                    ],
                    output=Path(directory) / "out",
                    run_id="a" * 64,
                    namespace="fixture",
                    donor_test_fold=0,
                    donor_valid_fold=1,
                    model_artifact_sha256="b" * 64,
                )

    def test_bundle_repeats_sequence_only_profile_across_held_donors(self) -> None:
        roster = cells()
        matrix = np.linspace(
            0.02, 0.98, num=3 * len(roster), dtype=np.float32
        ).reshape(3, len(roster))
        held = [
            {
                "donor_id": "valid_a",
                "outer_fold": "1",
                "evaluation_role": "valid",
            },
            {
                "donor_id": "valid_b",
                "outer_fold": "1",
                "evaluation_role": "valid",
            },
        ]
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "bundle"
            bundle_path = scbasset.export_prediction_bundle(
                sequence_predictions=matrix,
                training_cells=roster,
                regions=regions(),
                held_donors=held,
                output=output,
                run_id="a" * 64,
                namespace="fixture",
                donor_test_fold=0,
                donor_valid_fold=1,
                model_artifact_sha256="b" * 64,
            )
            bundle = json.loads(bundle_path.read_text())
            self.assertTrue(bundle["metadata"]["profile_fixture_passed"])
            self.assertFalse(bundle["metadata"]["held_atac_input_exposed"])
            self.assertFalse(bundle["metadata"]["held_cell_embedding_available"])
            self.assertEqual(bundle["n_predictions"], 2 * 5 * 3)
            with (output / "predictions.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            by_donor: dict[str, dict[tuple[str, str], float]] = {}
            region_by_row = {}
            for donor in held:
                donor_hash = scbasset._join_hash(
                    "fixture", "unit", donor["donor_id"]
                )
                by_donor[donor_hash] = {}
            for donor in held:
                for lineage in scbasset.LINEAGES:
                    for region in regions():
                        row_hash = scbasset._join_hash(
                            "fixture",
                            "row",
                            f"{donor['donor_id']}\0{lineage}\0{region['region_id']}",
                        )
                        region_by_row[row_hash] = (lineage, region["region_id"])
            for row in rows:
                by_donor[row["donor_hash"]][region_by_row[row["row_hash"]]] = float(
                    row["predicted"]
                )
            values = list(by_donor.values())
            self.assertEqual(values[0], values[1])


if __name__ == "__main__":
    unittest.main()
