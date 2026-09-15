from __future__ import annotations

import inspect
from pathlib import Path
import unittest

from scripts import evaluate_scbasset_five_seed_rectangle as evaluator


ROOT = Path(__file__).resolve().parents[2]


class ScBassetFiveSeedRectangleEvaluationTests(unittest.TestCase):
    def test_contract_locks_exact_native_rectangle_and_claim_limits(self) -> None:
        contract = evaluator.load_evaluator_contract(
            ROOT
            / "config/evaluation/scbasset_five_seed_rectangular_evaluator_20260825.json",
            ROOT,
        )
        self.assertEqual(contract["rectangle"]["expected_fits"], 25)
        self.assertEqual(contract["rectangle"]["fixed_seeds"], list(evaluator.SEEDS))
        self.assertFalse(contract["claims"]["rna_conditioned_claim_allowed"])
        self.assertFalse(contract["claims"]["cross_task_family_ranking_allowed"])

    def test_all_prediction_sources_precede_outcome_authority_resolution(self) -> None:
        source = inspect.getsource(evaluator.evaluate)
        self.assertLess(source.index("load_locked_sources"), source.index("_verify_outcome_authorities"))

    def test_summary_ensembles_seeds_but_does_not_treat_them_as_replicates(self) -> None:
        units = []
        secondary = []
        for fold in evaluator.FOLDS:
            for lineage in evaluator.LINEAGES:
                for donor in ("d0", "d1"):
                    secondary.extend(
                        [
                            {
                                "outer_fold": fold,
                                "model_id": "training_lineage_mean",
                                "seed": "training_only",
                                "donor_hash": f"f{fold}-{donor}",
                                "lineage_id": lineage,
                                "peak_auprc": "0.5",
                                "profile_spearman": "0.4",
                            },
                            {
                                "outer_fold": fold,
                                "model_id": "training_global_mean",
                                "seed": "training_only",
                                "donor_hash": f"f{fold}-{donor}",
                                "lineage_id": lineage,
                                "peak_auprc": "0.4",
                                "profile_spearman": "0.3",
                            },
                            {
                                "outer_fold": fold,
                                "model_id": "scbasset",
                                "seed": "ensemble",
                                "donor_hash": f"f{fold}-{donor}",
                                "lineage_id": lineage,
                                "peak_auprc": "0.6",
                                "profile_spearman": "0.5",
                            },
                        ]
                    )
                    for block in ("b0", "b1"):
                        common = {
                            "outer_fold": fold,
                            "donor_hash": f"f{fold}-{donor}",
                            "block_hash": f"f{fold}-{block}",
                            "lineage_id": lineage,
                            "observed_insertions": 10,
                            "regions": 2,
                        }
                        units.extend(
                            [
                                {
                                    **common,
                                    "model_id": "training_lineage_mean",
                                    "seed": "training_only",
                                    "deviance_per_insertion": "1.0",
                                },
                                {
                                    **common,
                                    "model_id": "training_global_mean",
                                    "seed": "training_only",
                                    "deviance_per_insertion": "1.2",
                                },
                                {
                                    **common,
                                    "model_id": "scbasset",
                                    "seed": "ensemble",
                                    "deviance_per_insertion": "0.8",
                                },
                            ]
                        )
                        units.extend(
                            {
                                **common,
                                "model_id": "scbasset",
                                "seed": str(seed),
                                "deviance_per_insertion": "0.8",
                            }
                            for seed in evaluator.SEEDS
                        )
        result = evaluator.summarize(
            units, secondary, n_resamples=100, bootstrap_seed=17
        )
        self.assertEqual(
            result["strongest_training_only_baseline"], "training_lineage_mean"
        )
        self.assertAlmostEqual(
            result["five_seed_ensemble_relative_deviance_reduction"], 0.2
        )
        self.assertEqual(result["positive_gain_seeds"], 5)
        self.assertTrue(result["seed_stability_requirement_met"])
        self.assertEqual(result["paired_two_way_bootstrap"]["biological_unit"], "donor")

    def test_partial_paired_unit_universe_is_rejected(self) -> None:
        candidate = [
            {
                "outer_fold": 0,
                "donor_hash": "d0",
                "block_hash": "b0",
                "lineage_id": "hepatocyte",
                "deviance_per_insertion": "0.8",
            }
        ]
        baseline = [
            {
                "outer_fold": 0,
                "donor_hash": "d1",
                "block_hash": "b0",
                "lineage_id": "hepatocyte",
                "deviance_per_insertion": "1.0",
            }
        ]
        with self.assertRaises(evaluator.ScBassetFiveSeedEvaluationError):
            evaluator._paired_skill_rows(candidate, baseline)


if __name__ == "__main__":
    unittest.main()
