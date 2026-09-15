from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.commit_gse244832_atac_transport_predictions import validate_evidence_contract
from scripts.peakvi_gse244832_inference_only import (
    DONORS,
    LINEAGES,
    aggregate_probabilities,
    make_prediction_manifest,
    read_target_cells,
)


class PeakVIGSE244832InferenceOnlyTests(unittest.TestCase):
    def _states(self) -> np.ndarray:
        states = np.full((18, 4), 3, dtype=np.uint8)
        states.reshape(-1)[:48] = 0
        states.reshape(-1)[48:52] = 1
        return states

    def _cells(self, states: np.ndarray) -> list[dict[str, str]]:
        output: list[dict[str, str]] = []
        for donor_index, lineage_index in zip(*np.nonzero(states == 0), strict=True):
            output.append(
                {
                    "donor_id": DONORS[donor_index],
                    "lineage_id": LINEAGES[lineage_index],
                }
            )
        return output

    def test_masked_aggregation_uses_nan_for_ineligible_units(self) -> None:
        states = self._states()
        cells = self._cells(states)
        probabilities = np.full((48, 16_000), 0.25, dtype=np.float32)
        result = aggregate_probabilities(probabilities, cells, states)
        self.assertEqual(result.shape, (18, 4, 16_000))
        self.assertTrue(np.all(result[states == 0] == np.float32(0.25)))
        self.assertTrue(np.all(np.isnan(result[states != 0])))

    def test_ineligible_cell_cannot_enter_aggregation(self) -> None:
        states = self._states()
        cells = self._cells(states)
        cells[0] = {"donor_id": "D13", "lineage_id": "cholangiocyte"}
        with self.assertRaises(ValueError):
            aggregate_probabilities(
                np.full((48, 16_000), 0.5, dtype=np.float32), cells, states
            )

    def test_target_axis_requires_all_and_only_48_observed_units(self) -> None:
        states = self._states()
        fields = (
            "cell_index", "cell_hash", "donor_id", "lineage_id", "outer_fold",
            "selection_rank", "available_nuclei",
        )
        groups = self._cells(states)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cell_axis.tsv"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                for offset in range(3049):
                    group = groups[offset % len(groups)]
                    writer.writerow(
                        {
                            "cell_index": offset,
                            "cell_hash": f"{offset:064x}",
                            "donor_id": group["donor_id"],
                            "lineage_id": group["lineage_id"],
                            "outer_fold": "0",
                            "selection_rank": offset // len(groups),
                            "available_nuclei": "64",
                        }
                    )
            rows = read_target_cells(path, states)
            self.assertEqual(len(rows), 3049)

    def test_prediction_manifest_passes_commit_evidence_firewall(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            predictions = root / "predictions.float32.npy"
            missing = root / "missing_state.uint8.npy"
            model = root / "model_execution_manifest.json"
            np.save(predictions, np.ones((18, 4, 16_000), dtype=np.float32))
            np.save(missing, np.zeros((18, 4, 16_000), dtype=np.uint8))
            model.write_text(json.dumps({"model_id": "peakvi"}) + "\n", encoding="utf-8")
            config = {
                "seeds": [1103, 2207, 3301],
                "exchange_axis": {"artifacts_sha256": "a" * 64},
                "execution_registration": {"artifacts_sha256": "b" * 64},
                "query_authorities": {
                    "valid": {"artifacts_sha256": "c" * 64},
                    "test": {"artifacts_sha256": "d" * 64},
                },
            }
            manifest = make_prediction_manifest(
                config=config,
                config_sha256="e" * 64,
                rotation_id="valid_context_test_target",
                model_manifest=model,
                query_authority=root,
                predictions=predictions,
                missing_state=missing,
                output=root,
            )
            validate_evidence_contract(manifest)
            self.assertEqual(manifest["model_manifest"]["path"], "model_execution_manifest.json")

    def test_runner_has_no_training_optimizer(self) -> None:
        source = Path(__file__).parents[2] / "scripts/peakvi_gse244832_inference_only.py"
        text = source.read_text(encoding="utf-8")
        self.assertNotIn("torch.optim", text)
        self.assertNotIn("fit_peakvae", text)
        self.assertNotIn("fit_decoder", text)


if __name__ == "__main__":
    unittest.main()
