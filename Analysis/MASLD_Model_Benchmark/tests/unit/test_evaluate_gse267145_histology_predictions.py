from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from scripts.evaluate_gse267145_histology_predictions import (
    FIBROSIS_GROUP3,
    MODEL_IDS,
    PREDICTION_FIELDS,
    SEEDS,
    STAGE3,
    canonical_json,
    evaluate_campaign,
    write_tsv,
)


ROOT = Path(__file__).resolve().parents[2]
OUTCOMES = ROOT / "executions/model-data-061-21079623/activation/participant_endpoints.tsv"
FOLDS = ROOT / "executions/model-data-064-21079902/fixture/folds/participant_outer_folds.tsv"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class GSE267145HistologyEvaluatorTests(unittest.TestCase):
    def test_canonical_json_is_stable(self) -> None:
        self.assertEqual(canonical_json({"b": 2, "a": 1}), '{"a":1,"b":2}')

    def test_independent_evaluator_reports_paired_differences_pareto_and_limitations(self) -> None:
        endpoint_rows = read_rows(OUTCOMES)
        prediction_rows = []
        for row in endpoint_rows:
            stage = row["stage3"]
            fibrosis = int(row["fibrosis"])
            group = "F0" if fibrosis == 0 else "F1" if fibrosis == 1 else "F2_3"
            prediction_rows.append(
                {
                    "participant_id": row["participant_id"],
                    "outer_fold": row["outer_fold"],
                    **{f"probability_{label}": "1" if label == stage else "0" for label in STAGE3},
                    "predicted_stage3": stage,
                    "predicted_nash_crn_component_sum": row["nash_crn_component_sum"],
                    "predicted_fibrosis_cumulative_expected": row["fibrosis"],
                    "predicted_fibrosis_regression": row["fibrosis"],
                    **{f"probability_fibrosis_{label}": "1" if label == group else "0" for label in FIBROSIS_GROUP3},
                    "predicted_fibrosis_group3": group,
                }
            )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            predictions = root / "predictions"
            predictions.mkdir()
            for model_id in MODEL_IDS:
                write_tsv(predictions / f"{model_id}.tsv", PREDICTION_FIELDS, prediction_rows)
                for seed in SEEDS:
                    write_tsv(
                        predictions / f"{model_id}--seed-{seed}.tsv",
                        PREDICTION_FIELDS,
                        prediction_rows,
                    )
            output = root / "evaluation"
            receipt = evaluate_campaign(
                predictions_root=predictions,
                outcomes_path=OUTCOMES,
                folds_path=FOLDS,
                output=output,
                bootstrap_replicates=20,
            )
            self.assertEqual(receipt["models"], 11)
            self.assertEqual(len(read_rows(output / "paired_differences.tsv")), 99)
            self.assertEqual(len(read_rows(output / "pareto_ranking.tsv")), 11)
            disposition = json.loads((output / "disposition.json").read_text())
            self.assertEqual(
                disposition["fibrosis"]["limitation_code"],
                "sparse_f3_inner_training_instability",
            )
            self.assertFalse(disposition["recorded_sex"]["sex_fairness_claim_allowed"])
            male = [
                row for row in read_rows(output / "sex_error_audit.tsv")
                if row["recorded_sex"] == "M"
            ]
            self.assertTrue(all(row["stage3_macro_f1"].startswith("not_applicable") for row in male))

    def test_evaluator_does_not_import_fitter(self) -> None:
        source = (ROOT / "scripts/evaluate_gse267145_histology_predictions.py").read_text()
        self.assertNotIn("fit_gse267145_histology_baselines", source)


if __name__ == "__main__":
    unittest.main()
