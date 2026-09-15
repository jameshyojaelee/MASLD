"""Contract and estimator tests for the GSE267145 paired development contrasts."""

from __future__ import annotations

import json
from pathlib import Path
import unittest

import numpy as np

from masld_bench.hashing import canonical_sha256, sha256_file

import scripts.analyze_gse267145_histology_paired_contrasts as analyzer


ROOT = Path(__file__).resolve().parents[2]
PRESPEC_PATH = ROOT / "config/evaluation/gse267145_histology_paired_contrast_prespec_v1.json"
SEAL_SBATCH = ROOT / "slurm/seal_gse267145_histology_paired_contrast_prespec_cpu.sbatch"
RUN_SBATCH = ROOT / "slurm/run_gse267145_histology_paired_contrasts_cpu.sbatch"
FROZEN_METRICS = (
    ROOT / "executions/model-scoring-078-21092130/scores/standardized_metrics.tsv"
)


class PrespecBindingTest(unittest.TestCase):
    def setUp(self) -> None:
        self.prespec = json.loads(PRESPEC_PATH.read_text(encoding="utf-8"))

    def test_prespec_digest_is_bound_identically_in_both_wrappers(self) -> None:
        digest = sha256_file(PRESPEC_PATH)
        for wrapper in (SEAL_SBATCH, RUN_SBATCH):
            text = wrapper.read_text(encoding="utf-8")
            self.assertIn(f"PRESPEC_SHA256={digest}", text)

    def test_prespec_declares_no_confirmatory_inference(self) -> None:
        estimator = self.prespec["estimator"]
        boundary = self.prespec["claim_boundary"]
        self.assertFalse(estimator["p_values_calculated"])
        self.assertFalse(estimator["bh_adjustment_calculated"])
        self.assertFalse(boundary["model_ranking_allowed"])
        self.assertFalse(boundary["champion_claim_allowed"])
        self.assertFalse(boundary["confirmatory_inference_allowed"])
        self.assertEqual(
            estimator["multiplicity_family_declared"],
            "not_applicable_no_confirmatory_lock",
        )

    def test_declared_contrast_totals_are_internally_consistent(self) -> None:
        totals = self.prespec["declared_contrast_totals"]
        families = self.prespec["contrast_families"]
        declared = (
            families["modality_matched_algorithm"]["contrast_count"]
            + families["modality_mean_across_pairs"]["contrast_count"]
            + families["fusion_versus_parents"]["contrast_count"]
            + families["versus_training_stage_distribution"]["declared_contrast_count"]
        )
        self.assertEqual(declared, totals["declared"])
        self.assertEqual(totals["declared"], 189)
        self.assertEqual(totals["computable"], 169)
        self.assertEqual(totals["not_applicable"], 20)
        self.assertEqual(
            totals["computable"] + totals["not_applicable"], totals["declared"]
        )

    def test_both_failure_readings_are_declared_before_the_run(self) -> None:
        cells = self.prespec["declared_outcome_cells"]
        for required in ("M3", "M4", "F2", "F3", "B2"):
            self.assertIn(required, cells)
        self.assertIn("within noise", cells["M3"])
        self.assertIn("selection artifact", cells["M4"])
        self.assertIn("ordinal expected-value head", cells["F3"])

    def test_rejected_engine_is_named_with_both_reasons(self) -> None:
        estimator = self.prespec["estimator"]
        self.assertEqual(
            estimator["rejected_engine"],
            "masld_bench.evaluators.stats.paired_cluster_bootstrap",
        )
        self.assertEqual(len(estimator["rejected_engine_reasons"]), 2)


class NotApplicableCellTest(unittest.TestCase):
    def test_baseline_fibrosis_spearman_is_absent_not_zero(self) -> None:
        """The frozen scoring run must record the two cells as not estimable."""

        absent = []
        for line in FROZEN_METRICS.read_text(encoding="utf-8").splitlines()[1:]:
            fields = line.split("\t")
            if fields[1] != "training_stage_distribution" or fields[8] != "observed":
                continue
            if fields[5] in {
                "fibrosis_cumulative_spearman",
                "fibrosis_regression_spearman",
            }:
                absent.append((fields[5], fields[11], fields[14]))
        self.assertEqual(len(absent), 2)
        for _metric_id, estimate, valid_replicates in absent:
            self.assertEqual(estimate, "not_estimable")
            self.assertEqual(valid_replicates, "0")

    def test_not_applicable_count_is_ten_models_by_two_metrics(self) -> None:
        prespec = json.loads(PRESPEC_PATH.read_text(encoding="utf-8"))
        family = prespec["contrast_families"]["versus_training_stage_distribution"]
        self.assertEqual(len(family["arm_a_roster"]), 10)
        self.assertEqual(len(family["not_applicable_cells"]["metrics"]), 2)
        self.assertEqual(family["not_applicable_contrast_count"], 10 * 2)


class BootstrapReuseTest(unittest.TestCase):
    def test_regenerated_indices_match_the_frozen_scoring_digest(self) -> None:
        indices = np.random.default_rng(analyzer.BOOTSTRAP_SEED).integers(
            0,
            analyzer.EXPECTED_PARTICIPANTS,
            size=(analyzer.BOOTSTRAP_REPLICATES, analyzer.EXPECTED_PARTICIPANTS),
            endpoint=False,
        )
        self.assertEqual(
            canonical_sha256(indices.tolist()), analyzer.BOOTSTRAP_INDICES_SHA256
        )

    def test_analyzer_does_not_use_the_cluster_bootstrap_helper(self) -> None:
        source = (
            ROOT / "scripts/analyze_gse267145_histology_paired_contrasts.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("from masld_bench.evaluators.stats import", source)
        self.assertNotIn("paired_cluster_bootstrap(", source)


class SignedDifferenceTest(unittest.TestCase):
    """A positive signed difference must always favor arm A, in both directions."""

    def _contrast(self, metric_id: str, a_point: float, b_point: float):
        rng = np.random.default_rng(11)
        a_draws = a_point + rng.normal(0.0, 1e-6, size=256)
        b_draws = b_point + rng.normal(0.0, 1e-6, size=256)
        row, difference = analyzer._contrast(
            contrast_id=f"T--{metric_id}",
            family_id="T",
            metric_id=metric_id,
            arm_a="arm_a",
            arm_b="arm_b",
            arm_a_point=a_point,
            arm_b_point=b_point,
            arm_a_draws=a_draws,
            arm_b_draws=b_draws,
        )
        return row, difference

    def test_higher_is_better_metric_favors_the_larger_arm(self) -> None:
        row, _ = self._contrast("stage3_macro_f1", 0.70, 0.60)
        self.assertEqual(row["direction"], "maximize")
        self.assertGreater(float(row["signed_difference"]), 0.0)
        self.assertGreater(float(row["probability_improvement"]), 0.99)

    def test_lower_is_better_metric_favors_the_smaller_arm(self) -> None:
        row, _ = self._contrast("fibrosis_cumulative_ordinal_mae", 0.60, 0.70)
        self.assertEqual(row["direction"], "minimize")
        self.assertGreater(float(row["signed_difference"]), 0.0)
        self.assertGreater(float(row["probability_improvement"]), 0.99)

    def test_lower_is_better_metric_penalizes_the_larger_arm(self) -> None:
        row, _ = self._contrast("stage3_multiclass_brier", 0.70, 0.60)
        self.assertLess(float(row["signed_difference"]), 0.0)
        self.assertLess(float(row["probability_improvement"]), 0.01)

    def test_non_estimable_reference_yields_a_declared_not_applicable_cell(self) -> None:
        row, difference = analyzer._contrast(
            contrast_id="T--absent",
            family_id="T",
            metric_id="fibrosis_cumulative_spearman",
            arm_a="rna_hvg_pca_knn",
            arm_b=analyzer.BASELINE_MODEL_ID,
            arm_a_point=0.4,
            arm_b_point=float("nan"),
            arm_a_draws=np.full(8, 0.4),
            arm_b_draws=np.full(8, np.nan),
        )
        self.assertIsNone(difference)
        self.assertEqual(row["applicability_state"], "not_applicable")
        self.assertEqual(
            row["applicability_reason"], "not_applicable_baseline_not_estimable"
        )
        self.assertEqual(row["valid_bootstrap_replicates"], 0)

    def test_every_row_withholds_p_and_q(self) -> None:
        row, _ = self._contrast("nash_crn_component_sum_spearman", 0.70, 0.60)
        self.assertEqual(
            row["p_value"], "not_calculated_development_no_confirmatory_lock"
        )
        self.assertEqual(
            row["bh_adjusted_q_value"],
            "not_calculated_development_no_confirmatory_lock",
        )
        self.assertEqual(row["confirmatory_inference_allowed"], "false")


class MetricDirectionTest(unittest.TestCase):
    def test_directions_match_the_frozen_scoring_table(self) -> None:
        frozen: dict[str, str] = {}
        for line in FROZEN_METRICS.read_text(encoding="utf-8").splitlines()[1:]:
            fields = line.split("\t")
            frozen[fields[5]] = fields[10]
        self.assertEqual(set(frozen), set(analyzer.METRIC_DIRECTION))
        for metric_id, direction in analyzer.METRIC_DIRECTION.items():
            expected = "maximize" if direction > 0 else "minimize"
            self.assertEqual(frozen[metric_id], expected, metric_id)

    def test_matched_modality_pairs_hold_the_algorithm_fixed(self) -> None:
        for h3_id, rna_id, algorithm in analyzer.MATCHED_MODALITY_PAIRS:
            self.assertTrue(h3_id.startswith("h3_variance_pca_"))
            self.assertTrue(rna_id.startswith("rna_hvg_pca_"))
            self.assertTrue(h3_id.endswith(algorithm))
            self.assertTrue(rna_id.endswith(algorithm))


if __name__ == "__main__":
    unittest.main()
