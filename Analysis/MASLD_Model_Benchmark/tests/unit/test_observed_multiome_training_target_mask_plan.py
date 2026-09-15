from __future__ import annotations

from copy import deepcopy
import inspect
import json
from pathlib import Path
import unittest

from scripts.plan_gse296875_observed_multiome_training_targets import (
    TrainingTargetMaskPlanError,
    select_training_targets,
    validate_config,
)
from scripts import plan_gse296875_observed_multiome_training_targets as planner


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_training_target_mask_plan_20260825.json"


class ObservedMultiomeTrainingTargetMaskPlanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_config_passes_without_count_access(self) -> None:
        resolved = validate_config(ROOT, deepcopy(self.config))
        self.assertTrue(resolved["input_h5"].is_file())
        self.assertFalse(self.config["selection"]["selection_uses_count_values"])

    def test_synthetic_selection_is_buffered_balanced_and_deterministic(self) -> None:
        peaks = []
        assignment = {}
        for fold in range(3):
            chromosome = f"chr{fold + 1}"
            assignment[chromosome] = fold
            for index in range(12):
                start = index * 2000
                peaks.append({"peak_id": f"p{fold}_{index}", "chromosome": chromosome, "bed_start_0based": start, "bed_end_half_open": start + 100})
        inputs = [{"chromosome": "chr2", "bed_start_0based": "0", "bed_end_half_open": "100", "held_genomic_fold": "0"}]
        selected, census = select_training_targets(peaks, assignment, inputs, [], folds=3, per_source_fold=2, buffer_bp=500, namespace="test", seed=7)
        self.assertEqual(len(selected), 12)
        self.assertEqual(len(census), 6)
        self.assertTrue(all(row["source_genomic_fold"] != row["held_genomic_fold"] for row in selected))
        self.assertTrue(all(row["nearest_observed_input_distance_bp"] in {-1} or row["nearest_observed_input_distance_bp"] >= 500 for row in selected))
        replay, _ = select_training_targets(peaks, assignment, inputs, [], folds=3, per_source_fold=2, buffer_bp=500, namespace="test", seed=7)
        self.assertEqual(selected, replay)

    def test_insufficient_buffered_candidates_fail_closed(self) -> None:
        peaks = [
            {"peak_id": "p1", "chromosome": "chr1", "bed_start_0based": 0, "bed_end_half_open": 100},
            {"peak_id": "p2", "chromosome": "chr2", "bed_start_0based": 0, "bed_end_half_open": 100},
        ]
        inputs = [{"chromosome": "chr2", "bed_start_0based": "0", "bed_end_half_open": "100", "held_genomic_fold": "0"}]
        with self.assertRaises(TrainingTargetMaskPlanError):
            select_training_targets(peaks, {"chr1": 0, "chr2": 1}, inputs, [], folds=2, per_source_fold=1, buffer_bp=500, namespace="test", seed=7)

    def test_count_selected_plan_is_rejected(self) -> None:
        mutated = deepcopy(self.config)
        mutated["selection"]["selection_uses_count_values"] = True
        with self.assertRaises(TrainingTargetMaskPlanError):
            validate_config(ROOT, mutated)

    def test_source_never_opens_count_or_donor_datasets(self) -> None:
        source = inspect.getsource(planner)
        self.assertNotIn("counts_csr", source)
        self.assertNotIn('obs/donor_id', source)
        self.assertNotIn('obs/broad_label', source)


if __name__ == "__main__":
    unittest.main()
