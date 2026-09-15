from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path
import tempfile
import unittest

import torch

from scripts.materialize_gse281364_sequence_control_training_targets import (
    aggregate_training_targets,
)
from scripts.train_predict_gse281364_task_native_sequence_controls import (
    SequenceTrainingError,
    composite_loss_rows,
    fit_standardization,
    select_epoch_one_standard_error,
    target_arrays,
)


class GSE281364SequenceControlProductionTests(unittest.TestCase):
    def test_composite_loss_has_frozen_weights_and_signed_delta(self) -> None:
        predicted_reference = torch.tensor([[0.0, 1.0]])
        predicted_alternative = torch.tensor([[1.0, 0.0]])
        target_reference = torch.zeros((1, 2))
        target_alternative = torch.zeros((1, 2))
        observed = composite_loss_rows(
            predicted_reference,
            predicted_alternative,
            target_reference,
            target_alternative,
        )
        expected = 0.25 * torch.tensor([[0.0, 0.5]])
        expected += 0.25 * torch.tensor([[0.5, 0.0]])
        expected += torch.tensor([[0.5, 0.5]])
        self.assertTrue(torch.equal(observed, expected))

    def test_one_standard_error_selects_earliest_eligible_epoch(self) -> None:
        history = [
            {
                "epoch": float(epoch),
                "block_mean_loss": loss,
                "block_loss_standard_error": error,
            }
            for epoch, loss, error in (
                (9, 0.80, 0.02),
                (10, 0.72, 0.02),
                (11, 0.68, 0.03),
                (12, 0.66, 0.05),
                (13, 0.64, 0.05),
            )
        ]
        selected = select_epoch_one_standard_error(history, 10)
        self.assertEqual(selected["leader_epoch"], 13)
        self.assertEqual(selected["selected_epoch_count"], 11)
        self.assertAlmostEqual(selected["one_standard_error_threshold"], 0.69)

    def test_standardization_uses_only_passed_training_indices(self) -> None:
        reference = torch.tensor([[0.0, 2.0], [2.0, 4.0], [999.0, 999.0]])
        alternative = torch.tensor([[2.0, 4.0], [4.0, 6.0], [-999.0, -999.0]])
        mean, sd = fit_standardization(
            reference, alternative, torch.tensor([0, 1])
        )
        self.assertTrue(torch.equal(mean, torch.tensor([2.0, 4.0])))
        self.assertTrue(torch.equal(sd, torch.tensor([math.sqrt(2.0)] * 2)))

    def test_target_array_rejects_extra_held_element(self) -> None:
        entries = [{"element_id": "train"}]
        targets = [
            {
                "element_id": element,
                "context_id": context,
                "reference_activity": "0",
                "alternative_activity": "1",
            }
            for element in ("train", "held")
            for context in ("HepG2_control", "HepG2_PAOA")
        ]
        with self.assertRaisesRegex(SequenceTrainingError, "held or unused"):
            target_arrays(entries, targets)

    def test_materializer_skips_held_counts_before_numeric_parsing(self) -> None:
        fold_counts = [208, 208, 207, 208, 202]
        manifest = []
        for fold_index, count in enumerate(fold_counts):
            manifest.extend(
                {
                    "element_id": f"element-{fold_index}-{index}",
                    "outer_fold": f"fold-{fold_index}",
                }
                for index in range(count)
            )
        fields = (
            "element_id",
            "allele",
            "context_id",
            "experimental_replicate",
            "sample_id",
            "DNA",
            "RNA",
            "assay_state",
            "missing_reason",
            "pairing",
            "biological_unit",
            "donor_id",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "outcomes.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                for entry in manifest:
                    held = entry["outer_fold"] == "fold-0"
                    for context in ("HepG2_control", "HepG2_PAOA"):
                        for allele in ("ref", "alt"):
                            for replicate in range(1, 5):
                                writer.writerow(
                                    {
                                        "element_id": entry["element_id"],
                                        "allele": allele,
                                        "context_id": context,
                                        "experimental_replicate": replicate,
                                        "sample_id": f"{context}_r{replicate}",
                                        "DNA": "withheld" if held else "10",
                                        "RNA": "withheld" if held else "20",
                                        "assay_state": "withheld" if held else "observed",
                                        "missing_reason": "withheld" if held else "not_applicable",
                                        "pairing": "withheld" if held else "same_sample_different_aliquot",
                                        "biological_unit": "withheld" if held else "experimental_replicate",
                                        "donor_id": "withheld" if held else "not_applicable",
                                    }
                                )
            rows, receipt = aggregate_training_targets(path, manifest, "fold-0")
        self.assertEqual(len(rows), 825 * 2)
        self.assertEqual(receipt["held_test_elements"], 208)
        self.assertEqual(receipt["held_test_target_rows_emitted"], 0)
        self.assertEqual(
            receipt["held_test_replicate_rows_skipped_before_count_parsing"],
            208 * 2 * 2 * 4,
        )


if __name__ == "__main__":
    unittest.main()
