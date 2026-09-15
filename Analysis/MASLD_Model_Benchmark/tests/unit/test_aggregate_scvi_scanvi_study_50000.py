#!/usr/bin/env python3
"""Contract checks for the three-seed scVI/scANVI bundle and ensemble."""

from __future__ import annotations

from pathlib import Path
import unittest

from scripts import aggregate_scvi_scanvi_study_50000 as aggregate


ROOT = Path(__file__).parents[2]
WRAPPER = ROOT / "slurm/run_scvi_scanvi_study_50000_seed_bundle.sbatch"
SCRIPT = ROOT / "scripts/aggregate_scvi_scanvi_study_50000.py"
SCORER = ROOT / "scripts/score_scvi_scanvi_study_50000.py"
SCORER_WRAPPER = ROOT / "slurm/score_scvi_scanvi_study_50000.sbatch"


class AggregateSCVISCANVIStudy50000Tests(unittest.TestCase):
    def test_schema_and_roster_are_frozen(self) -> None:
        self.assertEqual(
            aggregate.MODEL_IDS, ("scvi_baseline", "scanvi_baseline")
        )
        self.assertEqual(
            aggregate.SEEDS, (20260824, 20260825, 20260826)
        )
        self.assertEqual(len(aggregate.FIELDS), 10)
        self.assertEqual(aggregate.FIELDS[0], "row_id")
        self.assertEqual(aggregate.FIELDS[4], "predicted_class")

    def test_one_gpu_bundle_contains_all_seeds_and_aggregation(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        header = "\n".join(
            line for line in text.splitlines() if line.startswith("#SBATCH")
        )
        self.assertIn("--job-name=model-work-151", header)
        self.assertIn("--partition=gpu", header)
        self.assertIn("--account=nslab", header)
        self.assertIn("--qos=nslab", header)
        self.assertIn("--gres=gpu:b6k:1", header)
        self.assertIn("--time=72:00:00", header)
        self.assertNotIn("--array", header)
        self.assertEqual(text.count("for SEED in 20260824 20260825 20260826"), 1)
        self.assertIn("aggregate_scvi_scanvi_study_50000.py", text)
        self.assertFalse(
            any(line.lstrip().startswith("sbatch ") for line in text.splitlines())
        )

    def test_aggregator_never_scores_or_reads_labels(self) -> None:
        text = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("broad_label", text)
        self.assertNotIn("macro_f1", text.lower())
        self.assertNotIn("brier", text.lower())
        self.assertIn('"metrics_calculated": False', text)
        self.assertIn('"sealed_outcomes_read": False', text)

    def test_independent_scorer_is_roster_only_adapter(self) -> None:
        text = SCORER.read_text(encoding="utf-8")
        self.assertIn('MODEL_IDS = ("scvi_baseline", "scanvi_baseline")', text)
        self.assertIn("evaluator.run(", text)
        self.assertNotIn("fit(", text)
        self.assertNotIn("LogisticRegression", text)
        self.assertIn("10,000-donor-bootstrap", text)
        wrapper = SCORER_WRAPPER.read_text(encoding="utf-8")
        header = "\n".join(
            line for line in wrapper.splitlines() if line.startswith("#SBATCH")
        )
        self.assertIn("--partition=cpu", header)
        self.assertIn("--cpus-per-task=4", header)
        self.assertIn("--mem=32G", header)
        self.assertNotIn("--gres=gpu", header)
        self.assertIn("PREDICTIONS_ROOT", wrapper)
        self.assertIn("PREDICTIONS_SHA256", wrapper)


if __name__ == "__main__":
    unittest.main()
