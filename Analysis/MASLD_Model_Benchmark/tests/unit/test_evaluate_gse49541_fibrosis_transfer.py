from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.evaluate_gse49541_fibrosis_transfer import (
    CLASSES,
    ExternalTransferEvaluatorError,
    average_precision,
    brier,
    evaluate_gate,
    load_outcomes,
    load_predictions,
    load_seed_scores,
    macro_f1,
    permutation_null,
    seed_direction_consistency,
    write_tsv,
)


GATE = {
    "gate_id": "gse49541_fibrosis_external_development_v1",
    "primary_endpoint": "participant_auprc_advanced_f3_f4",
    "external_development_only": True,
    "champion_eligible": False,
    "development_advance_conditions": {
        "brier_score_degradation_maximum": 0.01,
        "gain_direction_required_bootstrap_lower_bound_above_zero": True,
        "minimum_absolute_auprc_gain_over_strongest_baseline": 0.03,
        "minimum_seed_direction_consistency": "4_of_5",
    },
}


def write_labels(path: Path, mild: int = 40, advanced: int = 32) -> list[str]:
    rows = []
    for index in range(mild):
        rows.append({"row_id": f"p{index:03d}", "fibrosis_stage_group": "mild_f0_f1"})
    for index in range(advanced):
        rows.append(
            {"row_id": f"p{mild + index:03d}", "fibrosis_stage_group": "advanced_f3_f4"}
        )
    write_tsv(path, ("row_id", "fibrosis_stage_group"), rows)
    return [row["row_id"] for row in rows]


class NullReferenceTests(unittest.TestCase):
    def test_random_scorer_null_mean_exceeds_prevalence(self) -> None:
        rng = np.random.default_rng(5)
        observed = np.asarray([0] * 40 + [1] * 32, dtype=np.int64)
        scores = rng.uniform(size=72)
        null = permutation_null(observed, scores, replicates=2000, seed=99)
        prevalence = float(np.mean(observed))
        self.assertAlmostEqual(prevalence, 32 / 72, places=6)
        # The whole point: a random scorer does not average the prevalence.
        self.assertGreater(null["auprc_null_mean"], prevalence)
        self.assertGreater(null["auprc_null_p95"], null["auprc_null_mean"])
        self.assertAlmostEqual(null["auroc_null_mean"], 0.5, delta=0.02)

    def test_constant_scorer_has_a_degenerate_null(self) -> None:
        observed = np.asarray([0] * 40 + [1] * 32, dtype=np.int64)
        scores = np.full(72, 0.0444)
        null = permutation_null(observed, scores, replicates=500, seed=7)
        self.assertEqual(null["distinct_score_values"], 1)
        # Every permutation gives the same value, so the null has no spread.
        self.assertAlmostEqual(null["auprc_null_mean"], null["auprc_null_p95"], places=9)
        self.assertAlmostEqual(null["auprc_null_mean"], float(np.mean(observed)), places=9)

    def test_low_cardinality_scorer_null_is_tighter_than_continuous(self) -> None:
        rng = np.random.default_rng(11)
        observed = np.asarray([0] * 40 + [1] * 32, dtype=np.int64)
        continuous = rng.uniform(size=72)
        binary = (continuous > 0.5).astype(np.float64)
        wide = permutation_null(observed, continuous, replicates=2000, seed=3)
        tight = permutation_null(observed, binary, replicates=2000, seed=3)
        self.assertGreater(wide["auprc_null_p95"], tight["auprc_null_p95"])


class MetricTests(unittest.TestCase):
    def test_constant_prevalence_scorer_scores_prevalence_on_auprc(self) -> None:
        observed = np.asarray([0] * 40 + [1] * 32, dtype=np.int64)
        scores = np.full(72, 0.25)
        self.assertAlmostEqual(average_precision(observed, scores), 32 / 72, places=9)

    def test_low_constant_scorer_calls_everything_mild(self) -> None:
        observed = np.asarray([0] * 40 + [1] * 32, dtype=np.int64)
        scores = np.full(72, 0.044)
        value = macro_f1(observed, scores)
        # Only the mild class earns any F1; the advanced class earns zero.
        self.assertGreater(value, 0.0)
        self.assertLess(value, 0.5)

    def test_brier_is_the_binary_form(self) -> None:
        observed = np.asarray([0, 1], dtype=np.int64)
        self.assertAlmostEqual(brier(observed, np.asarray([0.0, 1.0])), 0.0)
        self.assertAlmostEqual(brier(observed, np.asarray([1.0, 0.0])), 1.0)
        self.assertAlmostEqual(brier(observed, np.asarray([0.5, 0.5])), 0.25)


class LoadingTests(unittest.TestCase):
    def test_label_count_mismatch_raises(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "labels.tsv"
            write_labels(path, mild=39, advanced=32)
            with self.assertRaises(ExternalTransferEvaluatorError):
                load_outcomes(path)

    def test_group_census_mismatch_raises(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "labels.tsv"
            write_labels(path, mild=41, advanced=31)
            with self.assertRaises(ExternalTransferEvaluatorError):
                load_outcomes(path)

    def test_predictions_must_sum_to_one(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "predictions.tsv"
            write_tsv(
                path,
                (
                    "row_id",
                    "model_id",
                    "predicted_fibrosis_group",
                    "probability_mild_f0_f1",
                    "probability_advanced_f3_f4",
                ),
                [
                    {
                        "row_id": "p000",
                        "model_id": "m",
                        "predicted_fibrosis_group": "mild_f0_f1",
                        "probability_mild_f0_f1": "0.3",
                        "probability_advanced_f3_f4": "0.3",
                    }
                ],
            )
            with self.assertRaises(ExternalTransferEvaluatorError):
                load_predictions(path, ["p000"])

    def test_seed_scores_must_be_rectangular(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "seed_scores.tsv"
            write_tsv(
                path,
                ("row_id", "model_seed", "probability_advanced_f3_f4"),
                [
                    {"row_id": "p000", "model_seed": "1", "probability_advanced_f3_f4": "0.1"},
                    {"row_id": "p001", "model_seed": "1", "probability_advanced_f3_f4": "0.2"},
                    {"row_id": "p000", "model_seed": "2", "probability_advanced_f3_f4": "0.3"},
                ],
            )
            with self.assertRaises(ExternalTransferEvaluatorError):
                load_seed_scores(path, ["p000", "p001"])


class SeedConsistencyTests(unittest.TestCase):
    def test_identical_seeds_are_reported_as_one_deterministic_fit(self) -> None:
        observed = np.asarray([0] * 40 + [1] * 32, dtype=np.int64)
        rng = np.random.default_rng(2)
        one = rng.uniform(size=72)
        seeds = np.tile(one, (5, 1))
        summary = seed_direction_consistency(observed, seeds, 0.0)
        self.assertEqual(summary["distinct_seed_prediction_vectors"], 1)
        self.assertIs(summary["deterministic_single_fit"], True)
        self.assertIs(summary["consistency_is_evidence_of_stability"], False)
        self.assertEqual(summary["seeds_above_comparator"], 5)

    def test_genuinely_varying_seeds_count_as_replicates(self) -> None:
        observed = np.asarray([0] * 40 + [1] * 32, dtype=np.int64)
        rng = np.random.default_rng(4)
        seeds = rng.uniform(size=(5, 72))
        summary = seed_direction_consistency(observed, seeds, 0.0)
        self.assertEqual(summary["distinct_seed_prediction_vectors"], 5)
        self.assertIs(summary["consistency_is_evidence_of_stability"], True)


class GateTests(unittest.TestCase):
    def _seed_summary(self, above: int, distinct: int) -> dict:
        return {
            "seeds": 5,
            "distinct_seed_prediction_vectors": distinct,
            "seed_auprc": [0.6] * 5,
            "seeds_above_comparator": above,
            "deterministic_single_fit": distinct == 1,
            "consistency_is_evidence_of_stability": distinct > 1,
            "interpretation": "test",
        }

    def test_all_four_conditions_met_passes(self) -> None:
        verdict = evaluate_gate(
            gate=GATE,
            candidate_id="c",
            comparator_id="b",
            candidate_metrics={
                "participant_auprc_advanced_f3_f4": 0.60,
                "participant_brier_fibrosis_group": 0.19,
            },
            comparator_metrics={
                "participant_auprc_advanced_f3_f4": 0.50,
                "participant_brier_fibrosis_group": 0.20,
            },
            paired_gain={"ci_low": 0.02, "ci_high": 0.18},
            seed_summary=self._seed_summary(5, 5),
        )
        self.assertEqual(verdict["verdict"], "PASS")
        self.assertEqual(verdict["conditions_met"], 4)

    def test_gain_below_threshold_fails(self) -> None:
        verdict = evaluate_gate(
            gate=GATE,
            candidate_id="c",
            comparator_id="b",
            candidate_metrics={
                "participant_auprc_advanced_f3_f4": 0.52,
                "participant_brier_fibrosis_group": 0.19,
            },
            comparator_metrics={
                "participant_auprc_advanced_f3_f4": 0.50,
                "participant_brier_fibrosis_group": 0.20,
            },
            paired_gain={"ci_low": 0.005, "ci_high": 0.04},
            seed_summary=self._seed_summary(5, 5),
        )
        self.assertEqual(verdict["verdict"], "FAIL")

    def test_deterministic_seeds_do_not_satisfy_the_consistency_condition(self) -> None:
        verdict = evaluate_gate(
            gate=GATE,
            candidate_id="c",
            comparator_id="b",
            candidate_metrics={
                "participant_auprc_advanced_f3_f4": 0.60,
                "participant_brier_fibrosis_group": 0.19,
            },
            comparator_metrics={
                "participant_auprc_advanced_f3_f4": 0.50,
                "participant_brier_fibrosis_group": 0.20,
            },
            paired_gain={"ci_low": 0.02, "ci_high": 0.18},
            seed_summary=self._seed_summary(5, 1),
        )
        seed_condition = next(
            row
            for row in verdict["conditions"]
            if row["condition_id"] == "minimum_seed_direction_consistency"
        )
        self.assertIs(seed_condition["met"], False)
        self.assertEqual(verdict["verdict"], "FAIL")

    def test_brier_degradation_beyond_the_cap_fails(self) -> None:
        verdict = evaluate_gate(
            gate=GATE,
            candidate_id="c",
            comparator_id="b",
            candidate_metrics={
                "participant_auprc_advanced_f3_f4": 0.60,
                "participant_brier_fibrosis_group": 0.25,
            },
            comparator_metrics={
                "participant_auprc_advanced_f3_f4": 0.50,
                "participant_brier_fibrosis_group": 0.20,
            },
            paired_gain={"ci_low": 0.02, "ci_high": 0.18},
            seed_summary=self._seed_summary(5, 5),
        )
        self.assertEqual(verdict["verdict"], "FAIL")

    def test_bootstrap_lower_bound_at_zero_fails(self) -> None:
        verdict = evaluate_gate(
            gate=GATE,
            candidate_id="c",
            comparator_id="b",
            candidate_metrics={
                "participant_auprc_advanced_f3_f4": 0.60,
                "participant_brier_fibrosis_group": 0.19,
            },
            comparator_metrics={
                "participant_auprc_advanced_f3_f4": 0.50,
                "participant_brier_fibrosis_group": 0.20,
            },
            paired_gain={"ci_low": -0.01, "ci_high": 0.22},
            seed_summary=self._seed_summary(5, 5),
        )
        self.assertEqual(verdict["verdict"], "FAIL")


if __name__ == "__main__":
    unittest.main()


class UnfittableBaselineTests(unittest.TestCase):
    """A missing baseline must be declared, never merely absent."""

    def _run(self, *, present, declared, candidate="gene_median_linear_svm"):
        from scripts.evaluate_gse49541_fibrosis_transfer import evaluate

        rng = np.random.default_rng(13)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            row_ids = write_labels(root / "labels.tsv")
            gate_path = root / "gate.json"
            gate_path.write_text(json.dumps(GATE), encoding="utf-8")
            roots = []
            for model_id in present:
                bundle = root / model_id
                bundle.mkdir()
                scores = rng.uniform(size=72)
                write_tsv(
                    bundle / "predictions.tsv",
                    (
                        "row_id",
                        "model_id",
                        "predicted_fibrosis_group",
                        "probability_mild_f0_f1",
                        "probability_advanced_f3_f4",
                    ),
                    [
                        {
                            "row_id": row_id,
                            "model_id": model_id,
                            "predicted_fibrosis_group": (
                                "advanced_f3_f4" if value > 0.5 else "mild_f0_f1"
                            ),
                            "probability_mild_f0_f1": format(1.0 - value, ".17g"),
                            "probability_advanced_f3_f4": format(value, ".17g"),
                        }
                        for row_id, value in zip(row_ids, scores)
                    ],
                )
                write_tsv(
                    bundle / "seed_scores.tsv",
                    ("row_id", "model_seed", "probability_advanced_f3_f4"),
                    [
                        {
                            "row_id": row_id,
                            "model_seed": seed,
                            "probability_advanced_f3_f4": format(value, ".17g"),
                        }
                        for seed in (1, 2)
                        for row_id, value in zip(row_ids, scores)
                    ],
                )
                roots.append(bundle)
            return evaluate(
                labels_path=root / "labels.tsv",
                prediction_roots=roots,
                gate_path=gate_path,
                candidate_model_id=candidate,
                output=root / "evaluation",
                arm_id="test",
                gate_eligible_arm=True,
                bootstrap_replicates=200,
                permutation_replicates=200,
                unfittable_baselines=declared,
            )

    def test_undeclared_missing_baseline_raises(self) -> None:
        present = [
            "gene_median_elastic_net",
            "gene_median_linear_svm",
            "per_array_rank_elastic_net",
            "training_prevalence",
        ]
        with self.assertRaises(ExternalTransferEvaluatorError):
            self._run(present=present, declared=None)

    def test_declared_missing_baseline_is_recorded(self) -> None:
        present = [
            "gene_median_elastic_net",
            "gene_median_linear_svm",
            "per_array_rank_elastic_net",
            "training_prevalence",
        ]
        receipt = self._run(
            present=present, declared=["gene_median_pca_elastic_net"]
        )
        self.assertEqual(
            receipt["unfittable_baselines_declared"], ["gene_median_pca_elastic_net"]
        )
        self.assertNotIn("gene_median_pca_elastic_net", receipt["scored_baselines"])

    def test_declaring_a_baseline_that_is_actually_present_raises(self) -> None:
        present = [
            "gene_median_elastic_net",
            "gene_median_linear_svm",
            "gene_median_pca_elastic_net",
            "per_array_rank_elastic_net",
            "training_prevalence",
        ]
        with self.assertRaises(ExternalTransferEvaluatorError):
            self._run(present=present, declared=["gene_median_pca_elastic_net"])

    def test_candidate_cannot_be_the_unfittable_one(self) -> None:
        present = [
            "gene_median_elastic_net",
            "gene_median_linear_svm",
            "per_array_rank_elastic_net",
            "training_prevalence",
        ]
        with self.assertRaises(ExternalTransferEvaluatorError):
            self._run(
                present=present,
                declared=["gene_median_pca_elastic_net"],
                candidate="gene_median_pca_elastic_net",
            )
