from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from masld_bench.artifacts import freeze_tree
import scripts.aggregate_gse267145_histology_production_predictions as aggregator


class GSE267145ProductionAggregationTests(unittest.TestCase):
    def _write_tsv(self, path: Path, fields, rows) -> None:
        with path.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)

    def _prediction_row(self, participant: str, fold: int, seed: int) -> dict[str, str]:
        shift = (seed - 1700) / 100_000.0
        stage = [0.6 - shift, 0.3, 0.1 + shift]
        fibrosis = [0.5, 0.3 + shift, 0.2 - shift]
        return {
            "participant_id": participant,
            "outer_fold": str(fold),
            "probability_NOR": format(stage[0], ".17g"),
            "probability_NAFL": format(stage[1], ".17g"),
            "probability_NASH": format(stage[2], ".17g"),
            "predicted_stage3": "NOR",
            "predicted_nash_crn_component_sum": format(1.0 + shift, ".17g"),
            "predicted_fibrosis_cumulative_expected": format(0.8 + shift, ".17g"),
            "predicted_fibrosis_regression": format(0.7 + shift, ".17g"),
            "probability_fibrosis_F0": format(fibrosis[0], ".17g"),
            "probability_fibrosis_F1": format(fibrosis[1], ".17g"),
            "probability_fibrosis_F2_3": format(fibrosis[2], ".17g"),
            "predicted_fibrosis_group3": "F0",
        }

    def _fixture(
        self,
        root: Path,
        *,
        omit: tuple[int, int] | None = None,
        metric_unit: tuple[int, int] | None = None,
    ) -> tuple[Path, Path, list[str]]:
        folds = root / "folds"
        folds.mkdir()
        participants = [f"P{index:03d}" for index in range(99)]
        fold_rows = [
            {"participant_id": participant, "outer_fold": str(index % 5)}
            for index, participant in enumerate(participants)
        ]
        self._write_tsv(
            folds / "participant_outer_folds.tsv",
            ("participant_id", "outer_fold"),
            fold_rows,
        )
        freeze_tree(folds, {"fixture": True})
        units = root / "units"
        units.mkdir()
        for outer_fold in aggregator.OUTER_FOLDS:
            fold_root = units / f"outer_{outer_fold}"
            fold_root.mkdir()
            held = [
                row["participant_id"]
                for row in fold_rows
                if int(row["outer_fold"]) == outer_fold
            ]
            for seed in aggregator.SEEDS:
                if omit == (outer_fold, seed):
                    continue
                unit = fold_root / f"seed_{seed}"
                predictions = unit / "fit/predictions"
                predictions.mkdir(parents=True)
                for model_id in aggregator.MODEL_IDS:
                    self._write_tsv(
                        predictions / f"{model_id}.tsv",
                        aggregator.PREDICTION_FIELDS,
                        [
                            self._prediction_row(participant, outer_fold, seed)
                            for participant in held
                        ],
                    )
                bad_metrics = metric_unit == (outer_fold, seed)
                freeze_tree(
                    predictions,
                    {
                        "artifact_class": "gse267145_histology_production_unit_predictions",
                        "outer_fold": outer_fold,
                        "seed": seed,
                        "prediction_files": 11,
                        "outcomes_read": False,
                        "metrics_calculated": bad_metrics,
                        "status": "passed_unscored",
                    },
                )
                (unit / "unit_receipt.json").write_text(
                    json.dumps(
                        {
                            "outer_fold": outer_fold,
                            "seed": seed,
                            "model_family_fits": 11,
                            "prediction_files": 11,
                            "outer_test_outcomes_read": False,
                            "outer_test_metrics_calculated": False,
                            "production_unit_complete": True,
                        },
                        sort_keys=True,
                    )
                    + "\n",
                    encoding="utf-8",
                )
                freeze_tree(
                    unit,
                    {
                        "artifact_class": "gse267145_histology_production_outer_seed_unit",
                        "outer_fold": outer_fold,
                        "seed": seed,
                        "outcomes_read": False,
                        "metrics_calculated": False,
                        "status": "passed_unscored",
                    },
                )
        return units, folds, participants

    def test_assembles_55_seed_files_and_11_ensembles_in_axis_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            units, folds, participants = self._fixture(root)
            output = root / "aggregate"
            receipt = aggregator.aggregate(units=units, folds=folds, output=output)
            self.assertEqual(receipt["logical_outer_seed_units"], 25)
            self.assertEqual(receipt["model_family_unit_fits"], 275)
            self.assertEqual(receipt["prediction_files"], 66)
            self.assertFalse(receipt["metrics_calculated"])
            self.assertFalse(receipt["seeds_are_biological_replicates"])
            self.assertEqual(len(list((output / "predictions").glob("*.tsv"))), 66)
            fields, rows = aggregator.read_tsv(
                output / "predictions/rna_hvg_pca_elastic_net.tsv"
            )
            self.assertEqual(fields, aggregator.PREDICTION_FIELDS)
            self.assertEqual([row["participant_id"] for row in rows], participants)
            expected = sum((seed - 1700) / 100_000.0 for seed in aggregator.SEEDS) / 5
            self.assertAlmostEqual(float(rows[0]["probability_NASH"]), 0.1 + expected)
            _, unit_rows = aggregator.read_tsv(output / "unit_manifest.tsv")
            self.assertEqual(len(unit_rows), 25)

    def test_missing_unit_fails_closed_before_prediction_publication(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            units, folds, _ = self._fixture(root, omit=(4, 1733))
            with self.assertRaises((aggregator.AggregationError, FileNotFoundError)):
                aggregator.aggregate(units=units, folds=folds, output=root / "aggregate")

    def test_metric_marked_unit_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            units, folds, _ = self._fixture(root, metric_unit=(2, 1721))
            with self.assertRaisesRegex(aggregator.AggregationError, "metadata differs"):
                aggregator.aggregate(units=units, folds=folds, output=root / "aggregate")

    def test_mutated_frozen_prediction_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            units, folds, _ = self._fixture(root)
            target = units / "outer_0/seed_1701/fit/predictions/training_stage_distribution.tsv"
            target.write_text("mutated\n", encoding="utf-8")
            with self.assertRaisesRegex(Exception, "size mismatch|checksum mismatch"):
                aggregator.aggregate(units=units, folds=folds, output=root / "aggregate")


if __name__ == "__main__":
    unittest.main()
