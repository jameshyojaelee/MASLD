from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest

from scripts.evaluate_gse267145_histology_predictions import PREDICTION_FIELDS
from scripts.validate_gse267145_histology_preflight_predictions import (
    MODEL_IDS,
    validate_preflight,
)


ROOT = Path(__file__).resolve().parents[2]
OUTCOMES_ROOT = ROOT / "executions/model-data-061-21079623/activation"
FOLDS_ROOT = ROOT / "executions/model-data-064-21079902/fixture/folds"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class GSE267145PreflightValidationTests(unittest.TestCase):
    def test_preflight_validator_joins_after_prediction_validation_without_metrics(self) -> None:
        held = [row for row in read_rows(FOLDS_ROOT / "participant_outer_folds.tsv") if row["outer_fold"] == "0"]
        rows = [
            {
                "participant_id": row["participant_id"],
                "outer_fold": "0",
                "probability_NOR": "0.34",
                "probability_NAFL": "0.33",
                "probability_NASH": "0.33",
                "predicted_stage3": "NOR",
                "predicted_nash_crn_component_sum": "3",
                "predicted_fibrosis_cumulative_expected": "1",
                "predicted_fibrosis_regression": "1",
                "probability_fibrosis_F0": "0.34",
                "probability_fibrosis_F1": "0.33",
                "probability_fibrosis_F2_3": "0.33",
                "predicted_fibrosis_group3": "F0",
            }
            for row in held
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            predictions = root / "predictions"
            predictions.mkdir()
            from scripts.evaluate_gse267145_histology_predictions import write_tsv

            for model_id in MODEL_IDS:
                write_tsv(predictions / f"{model_id}.tsv", PREDICTION_FIELDS, rows)
            output = root / "validation"
            result = validate_preflight(
                predictions=predictions,
                outcomes=OUTCOMES_ROOT,
                folds=FOLDS_ROOT,
                output=output,
                outer_fold=0,
                seed=1701,
            )
            self.assertFalse(result["metrics_calculated"])
            self.assertFalse(result["f3_limitation"]["champion_gate_eligible"])
            self.assertFalse((output / "metrics.tsv").exists())


if __name__ == "__main__":
    unittest.main()
