from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

from masld_bench.evaluators.rna_atac_cobolt_5seed_development import (
    CAMPAIGN_ID,
    EXPECTED_RUNS,
    FOLDS,
    MODEL_IDS,
    RNAATACCoboltFiveSeedEvaluationError,
    SEEDS,
    comparison_summary,
    discover_verified_run_grid,
    ensemble_mean,
    select_strongest_task_native_baseline,
)
from masld_bench.hashing import sha256_file


class RNAATACCoboltFiveSeedEvaluatorTests(unittest.TestCase):
    @staticmethod
    def _fixture(root: Path) -> tuple[Path, Path, str, list[dict[str, object]]]:
        candidate = root / "candidate"
        execution = root / "executions"
        candidate.mkdir()
        execution.mkdir()
        (candidate / "ARTIFACTS.json").write_text("{}\n", encoding="utf-8")
        runs: list[dict[str, object]] = []
        for model_id in MODEL_IDS:
            for seed in SEEDS:
                for fold in FOLDS:
                    run_id = sha256(f"{model_id}:{seed}:{fold}".encode()).hexdigest()
                    runs.append(
                        {
                            "task_id": "rna_conditioned_atac",
                            "model_id": model_id,
                            "seed": seed,
                            "fold": fold,
                            "run_id": run_id,
                            "action": ["prepare", "fit", "predict"],
                        }
                    )
                    bundle = (
                        execution
                        / "runs"
                        / run_id
                        / "attempt-001"
                        / "adapter_actions"
                        / "003-predict"
                        / "prediction_bundle.json"
                    )
                    bundle.parent.mkdir(parents=True)
                    bundle.write_text("{}\n", encoding="utf-8")
        (candidate / "plan.json").write_text(
            json.dumps(
                {"campaign": {"campaign_id": CAMPAIGN_ID}, "runs": runs}
            )
            + "\n",
            encoding="utf-8",
        )
        return candidate, execution, sha256_file(candidate / "ARTIFACTS.json"), runs

    def test_exact_grid_verifies_every_attempt_before_returning_paths(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate, execution, candidate_hash, _ = self._fixture(Path(temporary))
            verified: list[Path] = []

            def verifier(path: Path, *, require_succeeded: bool) -> dict[str, object]:
                self.assertTrue(require_succeeded)
                verified.append(path)
                return {"status": "succeeded"}

            grid = discover_verified_run_grid(
                candidate=candidate,
                execution_root=execution,
                candidate_sha256=candidate_hash,
                verifier=verifier,
            )
            self.assertEqual(len(verified), EXPECTED_RUNS)
            self.assertEqual(set(grid), set(MODEL_IDS))
            for model_id in MODEL_IDS:
                self.assertEqual(set(grid[model_id]), set(SEEDS))
                for seed in SEEDS:
                    self.assertEqual(set(grid[model_id][seed]), set(FOLDS))

    def test_missing_or_duplicate_grid_cell_fails_before_verification(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate, execution, candidate_hash, runs = self._fixture(Path(temporary))
            calls = 0

            def verifier(path: Path, *, require_succeeded: bool) -> dict[str, object]:
                nonlocal calls
                calls += 1
                return {"status": "succeeded"}

            runs.pop()
            (candidate / "plan.json").write_text(
                json.dumps(
                    {"campaign": {"campaign_id": CAMPAIGN_ID}, "runs": runs}
                )
                + "\n",
                encoding="utf-8",
            )
            with self.assertRaises(RNAATACCoboltFiveSeedEvaluationError):
                discover_verified_run_grid(
                    candidate=candidate,
                    execution_root=execution,
                    candidate_sha256=candidate_hash,
                    verifier=verifier,
                )
            self.assertEqual(calls, 0)

        with tempfile.TemporaryDirectory() as temporary:
            candidate, execution, candidate_hash, runs = self._fixture(Path(temporary))
            runs[-1] = dict(runs[0])
            (candidate / "plan.json").write_text(
                json.dumps(
                    {"campaign": {"campaign_id": CAMPAIGN_ID}, "runs": runs}
                )
                + "\n",
                encoding="utf-8",
            )
            with self.assertRaises(RNAATACCoboltFiveSeedEvaluationError):
                discover_verified_run_grid(
                    candidate=candidate,
                    execution_root=execution,
                    candidate_sha256=candidate_hash,
                    verifier=lambda *_args, **_kwargs: {"status": "succeeded"},
                )

    def test_five_seed_ensemble_is_an_arithmetic_mean_with_no_seed_inference(self) -> None:
        values = ensemble_mean(
            [
                [1.0, 2.0],
                [2.0, 3.0],
                [3.0, 4.0],
                [4.0, 5.0],
                [5.0, 6.0],
            ]
        )
        self.assertEqual(values, [3.0, 4.0])
        with self.assertRaises(RNAATACCoboltFiveSeedEvaluationError):
            ensemble_mean([[1.0], [2.0]])
        with self.assertRaises(RNAATACCoboltFiveSeedEvaluationError):
            ensemble_mean([[1.0], [2.0], [3.0], [4.0], [0.0]])

    def test_every_model_is_compared_to_one_frozen_task_native_baseline(self) -> None:
        ensembles = {
            model_id: {
                "total_multinomial_deviance": 100.0 + index,
                "lineage_mean_deviance": {
                    lineage: 10.0 + index for lineage in (
                        "cholangiocyte",
                        "fibroblast",
                        "hepatocyte",
                        "macrophage",
                        "t_cell",
                    )
                },
            }
            for index, model_id in enumerate(MODEL_IDS)
        }
        ensembles["mean_track"]["total_multinomial_deviance"] = 80.0
        ensembles["mean_track"]["lineage_mean_deviance"] = {
            lineage: 8.0 for lineage in (
                "cholangiocyte",
                "fibroblast",
                "hepatocyte",
                "macrophage",
                "t_cell",
            )
        }
        per_seed = {
            model_id: {
                seed: {
                    "total_multinomial_deviance": float(
                        ensembles[model_id]["total_multinomial_deviance"]
                    )
                }
                for seed in SEEDS
            }
            for model_id in MODEL_IDS
        }
        strongest = select_strongest_task_native_baseline(ensembles)
        self.assertEqual(strongest, "mean_track")
        summary = comparison_summary(
            per_seed=per_seed,
            ensembles=ensembles,
            strongest_baseline=strongest,
        )
        self.assertEqual(
            set(
                summary[
                    "ensemble_relative_deviance_reduction_vs_strongest_task_native"
                ]
            ),
            set(MODEL_IDS),
        )
        self.assertEqual(
            summary["positive_gain_seed_count"]["mean_track"], 0
        )


if __name__ == "__main__":
    unittest.main()
