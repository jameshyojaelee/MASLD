from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import unittest
from unittest.mock import patch

import scripts.run_gse267145_histology_production_bundle_v3 as runner


ROOT = Path(__file__).resolve().parents[2]
WRAPPER = ROOT / "slurm/run_gse267145_histology_production_bundle_v3_cpu.sbatch"


def diagnostic_receipt() -> dict:
    return {
        "schema_version": "masld-bench-gse267145-production-v3-validation-v1",
        "status": "passed_unscored",
        "revision_id": runner.REVISION_ID,
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
        "other_secondary_endpoint_code_paths_changed": False,
        "outer_test_outcomes_read": False,
        "metrics_calculated": False,
        "scorer_called": False,
        "production_submission_authorized": False,
    }


def fit_receipt() -> dict:
    return {
        "schema_version": (
            "masld-bench-gse267145-histology-baseline-fit-preflight-v3"
        ),
        "revision_id": runner.REVISION_ID,
        "status": "passed_timing_predictions_unscored",
        "primary_endpoint_code_path_changed": False,
        "other_secondary_endpoint_code_paths_changed": False,
        "outer_test_outcomes_read": False,
        "outer_test_metrics_calculated": False,
        "fibrosis_group3_fallback": {
            "fallback_id": runner.FALLBACK_ID,
            "activation_count": 1,
            "activated_model_ids": ["h3_variance_pca_knn"],
            "eligible_endpoint": "fibrosis_group3",
            "eligible_endpoint_role": "secondary_only",
            "outer_test_outcomes_used": False,
            "metrics_calculated": False,
        },
    }


def selection_receipts() -> dict:
    return {
        "h3_variance_pca_knn": {
            "stage3": {},
            "nash_crn_component_sum": {},
            "fibrosis_exact_regression": {},
            "fibrosis_cumulative": {},
            "fibrosis_group3": {
                "fallback_triggered": True,
                "fallback_id": runner.FALLBACK_ID,
                "parameters": {"fallback_id": runner.FALLBACK_ID},
                "valid_candidate_count": 0,
                "candidate_count": 12,
                "candidate_failure_count": 12,
                "candidate_failure_reasons": {runner.CONVERGENCE_FAILURE: 12},
                "one_standard_error_applied": False,
                "one_standard_error_not_applicable_no_fallback_hyperparameters": True,
                "outer_test_features_used": False,
                "outer_test_outcomes_used": False,
                "recorded_sex_used": False,
                "source_stage5_used": False,
                "randomness_used": False,
            },
        }
    }


class GSE267145ProductionBundleV3ContractTests(unittest.TestCase):
    def test_frozen_axes_and_parent_sources_are_preserved(self) -> None:
        self.assertEqual(runner.v2.OUTER_FOLDS, tuple(range(5)))
        self.assertEqual(runner.v2.SEEDS, (1701, 1709, 1721, 1723, 1733))
        self.assertEqual(runner.v2.MODEL_COUNT, 11)
        for relative, expected in (
            (
                "scripts/run_gse267145_histology_production_bundle.py",
                runner.BASE_RUNNER_SHA256,
            ),
            (
                "scripts/fit_gse267145_histology_baselines_v3.py",
                runner.V3_FITTER_SHA256,
            ),
            (
                "config/evaluation/gse267145_histology_production_v3.json",
                runner.CONTRACT_SHA256,
            ),
            (
                "config/evaluation/gse267145_histology_production_v3_reuse_manifest.json",
                runner.REUSE_MANIFEST_SHA256,
            ),
        ):
            observed = sha256((ROOT / relative).read_bytes()).hexdigest()
            self.assertEqual(observed, expected)

    def test_exactly_22_v2_units_are_reused_and_three_are_recomputed(self) -> None:
        reused, recomputed = runner._reuse_manifest()
        self.assertEqual(len(reused), 22)
        self.assertEqual(recomputed, {(1, 1721), (1, 1723), (1, 1733)})
        self.assertEqual(
            set(reused).union(recomputed),
            {
                (outer, seed)
                for outer in runner.v2.OUTER_FOLDS
                for seed in runner.v2.SEEDS
            },
        )

    def test_validation_receipt_requires_control_parity_and_failed_unit_activation(self) -> None:
        runner._validate_diagnostic_receipt(diagnostic_receipt())
        invalid = diagnostic_receipt()
        invalid["units"]["1/1721"]["activation_count"] = 0
        with self.assertRaises(runner.v2.ProductionBundleError):
            runner._validate_diagnostic_receipt(invalid)

    def test_unit_receipt_accepts_only_secondary_group_fallback(self) -> None:
        runner._validate_v3_fit_receipt(fit_receipt(), selection_receipts())
        invalid = deepcopy(selection_receipts())
        invalid["h3_variance_pca_knn"]["stage3"]["fallback_triggered"] = True
        with self.assertRaises(runner.v2.ProductionBundleError):
            runner._validate_v3_fit_receipt(fit_receipt(), invalid)

    def test_bundle_validation_exports_are_required(self) -> None:
        with patch.dict("os.environ", {}, clear=True):
            with self.assertRaises(runner.v2.ProductionBundleError):
                runner._bundle_validation_environment()

    def test_slurm_wrapper_locks_v3_resources_and_has_no_scoring_path(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        for line in (
            "#SBATCH --job-name=model-training-073",
            "#SBATCH --partition=cpu",
            "#SBATCH --cpus-per-task=16",
            "#SBATCH --mem=32G",
            "#SBATCH --time=04:00:00",
            "CAMPAIGN=${ROOT}/executions/model-training-069-production-v3",
            "DIAGNOSTIC=${ROOT}/executions/model-check-072-21088262",
            "FITTER=${ROOT}/scripts/fit_gse267145_histology_baselines_v3.py",
            "RUNNER=${ROOT}/scripts/run_gse267145_histology_production_bundle_v3.py",
            "REUSE_MANIFEST=${ROOT}/config/evaluation/gse267145_histology_production_v3_reuse_manifest.json",
            "--fold-workers 5",
            "--blas-threads 3",
            ': "${BUNDLE_VALIDATION_ROOT:?submit with the frozen bundle-validation root export}"',
            ': "${BUNDLE_VALIDATION_SHA256:?submit with the frozen bundle-validation SHA-256 export}"',
            "export BUNDLE_VALIDATION_ROOT",
            "export BUNDLE_VALIDATION_SHA256",
        ):
            self.assertIn(line, text)
        self.assertNotIn("#SBATCH --array", text)
        self.assertNotIn("evaluate_gse267145", text)
        self.assertNotIn("--outcomes", text)
        self.assertNotIn("scorer", text.lower())


if __name__ == "__main__":
    unittest.main()
