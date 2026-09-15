from __future__ import annotations

from pathlib import Path
import unittest

import scripts.run_gse267145_histology_production_bundle as runner


ROOT = Path(__file__).resolve().parents[2]
WRAPPER = ROOT / "slurm/run_gse267145_histology_production_bundle_cpu.sbatch"


class GSE267145ProductionBundleContractTests(unittest.TestCase):
    def test_slurm_resources_and_stable_campaign_are_locked(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        for line in (
            "#SBATCH --job-name=model-training-069",
            "#SBATCH --partition=cpu",
            "#SBATCH --cpus-per-task=16",
            "#SBATCH --mem=32G",
            "#SBATCH --time=04:00:00",
            "#SBATCH --qos=nslab",
            "#SBATCH --account=nslab",
            "CAMPAIGN=${ROOT}/executions/model-training-069-production-v2",
            "DIAGNOSTIC=${ROOT}/executions/model-check-070-21082520",
            "export OMP_NUM_THREADS=3",
            "export MKL_NUM_THREADS=3",
            "export OPENBLAS_NUM_THREADS=3",
            "--fold-workers 5",
            "--blas-threads 3",
        ):
            self.assertIn(line, text)
        self.assertNotIn("#SBATCH --array", text)
        self.assertNotIn("--partition=gpu", text)
        self.assertNotIn("innovation", text)

    def test_bundle_has_exact_25_unit_roster_and_atomic_publication(self) -> None:
        self.assertEqual(runner.OUTER_FOLDS, tuple(range(5)))
        self.assertEqual(runner.SEEDS, (1701, 1709, 1721, 1723, 1733))
        self.assertEqual(len(runner.OUTER_FOLDS) * len(runner.SEEDS), 25)
        self.assertEqual(runner.MODEL_COUNT, 11)
        source = (ROOT / "scripts/run_gse267145_histology_production_bundle.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("ThreadPoolExecutor(max_workers=arguments.fold_workers)", source)
        self.assertIn("publish_directory_noreplace", source)
        self.assertIn("freeze_tree", source)
        self.assertIn("verify_frozen_tree", source)
        self.assertIn("LOCK_EX | fcntl.LOCK_NB", source)

    def test_wrapper_has_no_label_source_or_scoring_command(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        self.assertNotIn("evaluate_gse267145", text)
        self.assertNotIn("validate_gse267145_histology_preflight_predictions", text)
        self.assertNotIn("--outcomes", text)
        self.assertNotIn("--workers", text)

    def test_validated_one_unit_hashes_are_bound(self) -> None:
        self.assertEqual(
            runner.FITTER_SHA256,
            "b9832e9ddb1fe794103167db23467099288e3d9c351971716e043d6405d9dd27",
        )
        self.assertEqual(
            runner.DIAGNOSTIC_ARTIFACTS_SHA256,
            "da2e5cb0cd6da1d7f7ecb15cf9e38cd0d11a7ea5a8a0e8c989cf4f1174b6f6ef",
        )
        self.assertEqual(
            runner.PREFLIGHT_ARTIFACTS_SHA256,
            "2b95408bc92719706b728fd795ae918bb03280ff99ff8c549965318fee5d5bf3",
        )
        self.assertEqual(
            runner.FIT_VIEWS_ARTIFACTS_SHA256,
            "8254274f1b11089656f4b6d5ec446a6ba491e2f28ae8abcc088fc6e78ac75d6a",
        )


if __name__ == "__main__":
    unittest.main()
