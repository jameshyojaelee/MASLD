from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import unittest
from unittest.mock import patch

import scripts.run_gse267145_histology_production_bundle_v4 as runner


ROOT = Path(__file__).resolve().parents[2]
WRAPPER = ROOT / "slurm/run_gse267145_histology_production_bundle_v4_cpu.sbatch"


def diagnostic_receipt() -> dict:
    return {
        "schema_version": "masld-bench-gse267145-production-v4-validation-v1",
        "status": "passed_unscored",
        "revision_id": runner.REVISION_ID,
        "delegated_algorithm_revision_id": runner.DELEGATED_ALGORITHM_REVISION_ID,
        "actual_child_environment_self_test_passed": True,
        "unit_tests_passed": True,
        "v2_outer_1_seed_1701_predictions_byte_identical": True,
        "units": {
            "1/1701": {"activation_count": 0, "activated_model_ids": []},
            "1/1721": {
                "activation_count": 1,
                "activated_model_ids": ["h3_variance_pca_knn"],
            },
        },
        "primary_endpoint_code_path_changed": False,
        "secondary_endpoint_code_paths_changed": False,
        "outer_test_outcomes_read": False,
        "metrics_calculated": False,
        "scorer_called": False,
        "production_submission_authorized": False,
    }


class GSE267145ProductionBundleV4ContractTests(unittest.TestCase):
    def test_v4_sources_and_delegated_algorithm_are_pinned(self) -> None:
        for relative, expected in (
            (
                "config/evaluation/gse267145_histology_production_v4.json",
                runner.SOFTWARE_CONTRACT_SHA256,
            ),
            (
                "scripts/fit_gse267145_histology_baselines_v4.py",
                runner.V4_FITTER_SHA256,
            ),
            (
                "scripts/fit_gse267145_histology_baselines_v3.py",
                runner.v3.V3_FITTER_SHA256,
            ),
            (
                "config/evaluation/gse267145_histology_production_v3_reuse_manifest.json",
                runner.v3.REUSE_MANIFEST_SHA256,
            ),
        ):
            self.assertEqual(sha256((ROOT / relative).read_bytes()).hexdigest(), expected)

    def test_exact_reuse_and_recompute_rosters_are_unchanged(self) -> None:
        reused, recomputed = runner.v3._reuse_manifest()
        self.assertEqual(len(reused), 22)
        self.assertEqual(recomputed, {(1, 1721), (1, 1723), (1, 1733)})

    def test_v4_parity_receipt_requires_child_execution_and_exact_fallback(self) -> None:
        runner._validate_diagnostic_receipt(diagnostic_receipt())
        invalid = diagnostic_receipt()
        invalid["actual_child_environment_self_test_passed"] = False
        with self.assertRaises(runner.v3.v2.ProductionBundleError):
            runner._validate_diagnostic_receipt(invalid)

    def test_bundle_validation_exports_remain_required(self) -> None:
        with patch.dict("os.environ", {}, clear=True):
            with self.assertRaises(runner.v3.v2.ProductionBundleError):
                runner._bundle_validation_environment()

    def test_slurm_wrapper_is_new_immutable_target_without_scoring(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        for line in (
            "#SBATCH --job-name=model-training-077",
            "#SBATCH --partition=cpu",
            "#SBATCH --cpus-per-task=16",
            "#SBATCH --mem=32G",
            "#SBATCH --time=04:00:00",
            "CAMPAIGN=${ROOT}/executions/model-training-069-production-v4",
            "DIAGNOSTIC=${ROOT}/executions/model-check-075-21088989",
            "FITTER=${ROOT}/scripts/fit_gse267145_histology_baselines_v4.py",
            "RUNNER=${ROOT}/scripts/run_gse267145_histology_production_bundle_v4.py",
            "--fold-workers 5",
            "--blas-threads 3",
        ):
            self.assertIn(line, text)
        self.assertNotIn("#SBATCH --array", text)
        self.assertNotIn("evaluate_gse267145", text)
        self.assertNotIn("--outcomes", text)


if __name__ == "__main__":
    unittest.main()
