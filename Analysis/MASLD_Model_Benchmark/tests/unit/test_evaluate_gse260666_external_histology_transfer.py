from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from scripts.evaluate_gse260666_external_histology_transfer import (
    CLASSES,
    ExternalTransferEvaluatorError,
    MANDATORY_BASELINES,
    evaluate,
)


class GSE260666ExternalEvaluatorTests(unittest.TestCase):
    def _write_tsv(self, path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def _labels(self, root: Path) -> Path:
        source = ["healthy_control"] * 6 + ["non-alcoholic_fatty_liver_disease_(nafld)"] * 6 + ["non-alcoholic_steatohepatitis_(nash)"] * 4
        mapping = {
            "healthy_control": "NOR",
            "non-alcoholic_fatty_liver_disease_(nafld)": "NAFL",
            "non-alcoholic_steatohepatitis_(nash)": "NASH",
        }
        rows = [
            {
                "row_id": f"row_{index:02d}",
                "source_label": label,
                "evaluation_stage3": mapping[label],
            }
            for index, label in enumerate(source)
        ]
        path = root / "labels.tsv"
        self._write_tsv(path, tuple(rows[0]), rows)
        return path

    def _predictions(self, root: Path, model_id: str, *, correct: bool) -> Path:
        truth = [0] * 6 + [1] * 6 + [2] * 4
        rows = []
        for index, observed in enumerate(truth):
            predicted = observed if correct else (observed + 1) % 3
            probability = [0.05, 0.05, 0.05]
            probability[predicted] = 0.90
            rows.append(
                {
                    "row_id": f"row_{index:02d}",
                    "model_id": model_id,
                    "probability_NOR": probability[0],
                    "probability_NAFL": probability[1],
                    "probability_NASH": probability[2],
                    "predicted_stage3": CLASSES[predicted],
                }
            )
        path = root / f"{model_id}.tsv"
        self._write_tsv(path, tuple(rows[0]), rows)
        return path

    def test_complete_baseline_and_candidate_evaluation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            labels = self._labels(root)
            predictions = [
                self._predictions(root, model_id, correct=index == 1)
                for index, model_id in enumerate(MANDATORY_BASELINES)
            ]
            predictions.append(self._predictions(root, "candidate_model", correct=True))
            output = root / "evaluation"
            receipt = evaluate(
                labels_path=labels,
                prediction_paths=predictions,
                output=output,
                bootstrap_replicates=100,
                bootstrap_seed=17,
            )
            self.assertEqual(receipt["strongest_baseline"], MANDATORY_BASELINES[1])
            self.assertFalse(receipt["nafl_nash_pooled"])
            self.assertFalse(receipt["champion_claim_eligible"])
            self.assertTrue(receipt["external_development_only"])
            self.assertFalse(receipt["external_outcomes_used_for_model_repair"])
            self.assertFalse(receipt["nominal_or_confirmatory_p_values_computed"])
            self.assertFalse(receipt["diagnostic_or_prognostic_claim_eligible"])
            self.assertFalse(receipt["confirmatory_claim_eligible"])
            class_rows = (output / "per_class.tsv").read_text(encoding="utf-8")
            self.assertIn("candidate_model\tNAFL\t6", class_rows)
            gain = (output / "candidate_vs_strongest_baseline.tsv").read_text(encoding="utf-8")
            self.assertIn("candidate_model", gain)
            class_ci = (output / "per_class_f1_bootstrap_intervals.tsv").read_text(
                encoding="utf-8"
            )
            self.assertIn("candidate_model\tNASH\t1.0", class_ci)
            confusion = json.loads((output / "confusion_matrices.json").read_text())
            self.assertEqual(confusion["candidate_model"]["NASH"]["NASH"], 4)

    def test_missing_mandatory_baseline_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            labels = self._labels(root)
            prediction = self._predictions(root, MANDATORY_BASELINES[0], correct=True)
            with self.assertRaises(ExternalTransferEvaluatorError):
                evaluate(
                    labels_path=labels,
                    prediction_paths=[prediction],
                    output=root / "evaluation",
                    bootstrap_replicates=10,
                )


if __name__ == "__main__":
    unittest.main()
