#!/usr/bin/env python3
"""Focused checks for the GSE267145 paired molecular baseline campaign."""

from __future__ import annotations

import ast
import csv
from pathlib import Path
import tempfile
import tomllib
import unittest

import numpy as np

from scripts.build_gse267145_paired_baseline_views import build_views
from scripts.build_gse267145_paired_molecular_fixture import build_fixture
from scripts.evaluate_gse267145_paired_baseline_predictions import (
    BOOTSTRAP_REPLICATES,
    multinomial_deviance,
    tie_aware_retrieval,
)
from scripts.fit_predict_gse267145_paired_baselines import (
    h3_hellinger,
    library_log1p,
    normalize_profile,
    reduced_rank_predict,
    ridge_predict,
    sparse_cca_loadings,
)
from tests.unit.test_build_gse267145_paired_molecular_fixture import make_inputs


ROOT = Path(__file__).parents[2]
FITTER = ROOT / "scripts/fit_predict_gse267145_paired_baselines.py"
EVALUATOR = ROOT / "scripts/evaluate_gse267145_paired_baseline_predictions.py"
VIEW_WRAPPER = ROOT / "slurm/build_gse267145_paired_baseline_views_cpu.sbatch"
FIT_WRAPPER = ROOT / "slurm/run_gse267145_paired_baselines_cpu.sbatch"
EVALUATE_WRAPPER = ROOT / "slurm/evaluate_gse267145_paired_baselines_cpu.sbatch"
TASK = ROOT / "config/evaluation/paired_bulk_rna_h3k27ac_task.toml"
GATE = ROOT / "config/evaluation/paired_bulk_rna_h3k27ac_promotion_gate.json"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class GSE267145PairedBaselineCampaignTests(unittest.TestCase):
    def test_registered_task_is_participant_safe_and_development_only(self) -> None:
        task = tomllib.loads(TASK.read_text(encoding="utf-8"))
        import json

        gate = json.loads(GATE.read_text(encoding="utf-8"))
        self.assertEqual(task["unit_of_inference"], "participant")
        self.assertEqual(
            task["required_pairing_levels"], ["same_sample_different_aliquot"]
        )
        self.assertEqual(
            task["baseline_model_ids"],
            [
                "training_mean_h3_profile",
                "pca_ridge",
                "reduced_rank_regression",
                "pls2",
                "sparse_cca",
            ],
        )
        contract = "\n".join(task["admission_gates"]).lower()
        for blocked in ("fibrosis", "sex", "opaque", "held-participant h3k27ac"):
            self.assertIn(blocked, contract)
        self.assertFalse(gate["champion_eligible"])
        self.assertEqual(gate["claim_mode"], "development_only")

    def test_fold_views_structurally_withhold_query_h3(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            matrix, join, crosswalk = make_inputs(root)
            fixture = root / "fixture"
            build_fixture(
                matrix_root=matrix,
                join_path=join,
                crosswalk_path=crosswalk,
                output=fixture,
                expected_participants=5,
                expected_rna_source_features=4,
                expected_rna_allowed_features=3,
                expected_h3_features=3,
                expected_fold_counts={0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
            )
            output = root / "views"
            result = build_views(
                molecular=fixture / "molecular",
                folds=fixture / "folds",
                output=output,
                task_spec_sha256="a" * 64,
                promotion_gate_sha256="b" * 64,
                expected_participants=5,
                expected_rna_features=3,
                expected_h3_features=3,
                expected_fold_counts={0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
            )
            self.assertTrue(result["held_participant_h3_is_evaluator_only"])
            for fold in range(5):
                fold_root = output / f"model/outer_{fold}"
                self.assertTrue((fold_root / "training_h3k27ac_counts.npy").is_file())
                self.assertFalse(any("query_h3" in path.name for path in fold_root.iterdir()))
                training = rows(fold_root / "training_participants.tsv")
                query = rows(fold_root / "query_participants.tsv")
                self.assertEqual(len(training), 4)
                self.assertEqual(len(query), 1)
                self.assertFalse(
                    {row["participant_id"] for row in training}
                    & {row["participant_id"] for row in query}
                )
            self.assertTrue((output / "evaluator/observed_h3k27ac_counts.npy").is_file())
            self.assertFalse(any("rna" in path.name for path in (output / "evaluator").iterdir()))

    def test_assay_native_transforms_and_predictions_are_valid(self) -> None:
        rna = np.asarray([[1.0, 3.0, 0.0], [2.0, 1.0, 4.0]])
        transformed = library_log1p(rna)
        self.assertEqual(transformed.shape, rna.shape)
        self.assertTrue(np.all(np.isfinite(transformed)))
        h3, totals = h3_hellinger(np.asarray([[1, 3, 0], [0, 2, 2]], dtype=np.uint32))
        np.testing.assert_allclose(np.sum(h3**2, axis=1), np.ones(2))
        np.testing.assert_array_equal(totals, np.asarray([4.0, 4.0]))
        profile = normalize_profile(np.asarray([[1.0, -1.0, 0.0]]))
        self.assertTrue(np.all(profile > 0.0))
        np.testing.assert_allclose(profile.sum(axis=1), np.ones(1))

    def test_ridge_and_reduced_rank_predict_without_observed_query_h3(self) -> None:
        rng = np.random.default_rng(17)
        train_x = rng.normal(size=(12, 4))
        query_x = rng.normal(size=(3, 4))
        train_y = rng.normal(size=(12, 7))
        ridge = ridge_predict(train_x, query_x, train_y, alpha=1.0)
        reduced = reduced_rank_predict(train_x, query_x, train_y)
        self.assertEqual(ridge.shape, (3, 7))
        self.assertEqual(reduced.shape, (3, 7))
        self.assertTrue(np.all(np.isfinite(ridge)))
        self.assertTrue(np.all(np.isfinite(reduced)))

    def test_sparse_cca_loadings_are_sparse_and_finite(self) -> None:
        rng = np.random.default_rng(19)
        covariance = rng.normal(size=(20, 30))
        x, y, iterations = sparse_cca_loadings(
            covariance,
            components=3,
            x_nonzero=5,
            y_nonzero=7,
        )
        self.assertEqual(x.shape, (20, 3))
        self.assertEqual(y.shape, (30, 3))
        self.assertEqual(len(iterations), 3)
        self.assertTrue(np.all(np.count_nonzero(x, axis=0) <= 5))
        self.assertTrue(np.all(np.count_nonzero(y, axis=0) <= 7))
        self.assertTrue(np.all(np.isfinite(x)))
        self.assertTrue(np.all(np.isfinite(y)))

    def test_deviance_and_retrieval_are_evaluator_only_metrics(self) -> None:
        observed = np.asarray([2, 3, 5], dtype=np.uint32)
        perfect = observed / observed.sum()
        self.assertAlmostEqual(multinomial_deviance(observed, perfect), 0.0, places=10)
        self.assertGreater(
            multinomial_deviance(observed, np.asarray([0.8, 0.1, 0.1])), 0.0
        )
        reciprocal, top1 = tie_aware_retrieval(np.eye(3))
        np.testing.assert_array_equal(reciprocal, np.ones(3))
        np.testing.assert_array_equal(top1, np.ones(3))
        tied_reciprocal, tied_top1 = tie_aware_retrieval(np.ones((3, 3)))
        np.testing.assert_allclose(tied_reciprocal, np.full(3, 0.5))
        np.testing.assert_allclose(tied_top1, np.full(3, 1.0 / 3.0))

    def test_prediction_and_metric_code_are_independent(self) -> None:
        fitter_text = FITTER.read_text(encoding="utf-8")
        evaluator_text = EVALUATOR.read_text(encoding="utf-8")
        ast.parse(fitter_text)
        ast.parse(evaluator_text)
        self.assertNotIn("evaluate_gse267145_paired_baseline_predictions", fitter_text)
        self.assertNotIn("fit_predict_gse267145_paired_baselines", evaluator_text)
        self.assertNotIn("observed_h3k27ac_counts.npy", fitter_text)
        self.assertIn('"held_participant_h3_read": False', fitter_text)
        self.assertIn('"metrics_calculated": False', fitter_text)
        self.assertIn('"models_fit_inside_evaluator": False', evaluator_text)
        self.assertIn('"champion_promotion_performed": False', evaluator_text)
        self.assertEqual(BOOTSTRAP_REPLICATES, 10_000)

    def test_wrappers_are_cpu_nslab_bundled_and_not_self_submitting(self) -> None:
        expectations = (
            (VIEW_WRAPPER, "model-data-087"),
            (FIT_WRAPPER, "model-cpu-fit-510"),
            (EVALUATE_WRAPPER, "model-cpu-score-511"),
        )
        for path, job_name in expectations:
            text = path.read_text(encoding="utf-8")
            header = "\n".join(
                line for line in text.splitlines() if line.startswith("#SBATCH")
            )
            self.assertIn(f"--job-name={job_name}", header)
            self.assertIn("--partition=cpu", header)
            self.assertIn("--account=nslab", header)
            self.assertIn("--qos=nslab", header)
            self.assertNotIn("--array", header)
            self.assertNotIn("innovation", text.lower())
            self.assertFalse(
                any(line.lstrip().startswith("sbatch ") for line in text.splitlines())
            )
        self.assertIn("for OUTER_FOLD in 0 1 2 3 4", FIT_WRAPPER.read_text())


if __name__ == "__main__":
    unittest.main()
