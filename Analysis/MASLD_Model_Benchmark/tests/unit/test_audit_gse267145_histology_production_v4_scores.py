from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import tempfile
import unittest

import numpy as np

import scripts.audit_gse267145_histology_production_v4_scores as auditor
import scripts.score_gse267145_histology_production_v4 as scorer
from tests.unit.test_score_gse267145_histology_production_v4 import (
    perfect_rows,
    write_complete_prediction_roster,
)


ROOT = Path(__file__).resolve().parents[2]
OUTCOMES_ROOT = ROOT / "executions/model-data-061-21079623/activation"
FOLDS_ROOT = ROOT / "executions/model-data-064-21079902/fixture/folds"


class GSE267145ProductionV4ScoreAuditorTests(unittest.TestCase):
    def test_reference_evaluator_rederives_all_nine_point_estimates(self) -> None:
        participants, outer_folds, outcomes = auditor._load_outcomes(
            OUTCOMES_ROOT, FOLDS_ROOT
        )
        self.assertEqual(len(participants), 99)
        self.assertEqual(len(outer_folds), 99)
        points = auditor._reference_points(
            outcomes, auditor._prediction_arrays(perfect_rows())
        )
        self.assertEqual(set(points), {item[0] for item in auditor.METRIC_SPECS})
        self.assertEqual(points["stage3_macro_f1"], "1")
        self.assertEqual(points["stage3_multiclass_brier"], "0")

    def test_metric_audit_enforces_endpoint_applicability_and_secondary_lock(self) -> None:
        participants, outer_folds, outcomes = auditor._load_outcomes(
            OUTCOMES_ROOT, FOLDS_ROOT
        )
        prediction_rows = perfect_rows()
        predictions = scorer._load_predictions(prediction_rows)
        indices = np.random.default_rng(scorer.BOOTSTRAP_SEED).integers(
            0, 99, size=(25, 99), endpoint=False
        )
        points = scorer._point_metrics(outcomes, predictions)
        distributions = scorer._bootstrap_metrics(outcomes, predictions, indices)
        observed = scorer._metric_rows(
            model_id="training_stage_distribution",
            endpoint_id="fibrosis_group3",
            set_id="training_stage_distribution--endpoint-fibrosis_group3",
            points=points,
            distributions=distributions,
        )
        formatted_points = auditor._reference_points(
            outcomes, auditor._prediction_arrays(prediction_rows)
        )
        computed, not_applicable = auditor._validate_metric_rows(
            rows=observed,
            model_id="training_stage_distribution",
            endpoint_id="fibrosis_group3",
            set_id="training_stage_distribution--endpoint-fibrosis_group3",
            points=formatted_points,
        )
        self.assertEqual((computed, not_applicable), (1, 8))
        changed = deepcopy(observed)
        changed[2]["p_value"] = "0.01"
        with self.assertRaises(auditor.HistologyProductionScoreAuditError):
            auditor._validate_metric_rows(
                rows=changed,
                model_id="training_stage_distribution",
                endpoint_id="fibrosis_group3",
                set_id="training_stage_distribution--endpoint-fibrosis_group3",
                points=formatted_points,
            )

    def test_independent_auditor_rederives_exact_five_seed_ensemble(self) -> None:
        participants, outer_folds, _ = auditor._load_outcomes(
            OUTCOMES_ROOT, FOLDS_ROOT
        )
        with tempfile.TemporaryDirectory() as temporary:
            prediction_root = Path(temporary) / "predictions"
            write_complete_prediction_roster(prediction_root, perfect_rows())
            ensembles, lineage = auditor._verify_five_seed_ensembles(
                prediction_root=prediction_root,
                participants=participants,
                outer_folds=outer_folds,
            )
            self.assertEqual(set(ensembles), set(auditor.reference.MODEL_IDS))
            self.assertEqual(len(lineage), 11)
            self.assertTrue(
                all(len(value["seed_prediction_sha256"]) == 5 for value in lineage.values())
            )

    def test_constant_native_spearman_is_locked_not_estimable(self) -> None:
        participants, outer_folds, outcomes = auditor._load_outcomes(
            OUTCOMES_ROOT, FOLDS_ROOT
        )
        prediction_rows = perfect_rows()
        for row in prediction_rows:
            row["predicted_fibrosis_cumulative_expected"] = "0"
        predictions = scorer._load_predictions(prediction_rows)
        indices = np.random.default_rng(scorer.BOOTSTRAP_SEED).integers(
            0, 99, size=(25, 99), endpoint=False
        )
        points = scorer._point_metrics(outcomes, predictions)
        distributions = scorer._bootstrap_metrics(outcomes, predictions, indices)
        observed = scorer._metric_rows(
            model_id="training_stage_distribution",
            endpoint_id="fibrosis_cumulative",
            set_id="training_stage_distribution--endpoint-fibrosis_cumulative",
            points=points,
            distributions=distributions,
        )
        formatted_points = auditor._reference_points(
            outcomes, auditor._prediction_arrays(prediction_rows)
        )
        spearman = next(
            row
            for row in observed
            if row["metric_id"] == "fibrosis_cumulative_spearman"
        )
        self.assertEqual(spearman["estimate"], "not_estimable")
        self.assertEqual(spearman["ci95_low"], "not_estimable")
        self.assertEqual(spearman["ci95_high"], "not_estimable")
        self.assertEqual(spearman["valid_bootstrap_replicates"], 0)
        auditor._validate_metric_rows(
            rows=observed,
            model_id="training_stage_distribution",
            endpoint_id="fibrosis_cumulative",
            set_id="training_stage_distribution--endpoint-fibrosis_cumulative",
            points=formatted_points,
        )
        changed = deepcopy(observed)
        next(
            row
            for row in changed
            if row["metric_id"] == "fibrosis_cumulative_spearman"
        )["estimate"] = "0"
        with self.assertRaises(auditor.HistologyProductionScoreAuditError):
            auditor._validate_metric_rows(
                rows=changed,
                model_id="training_stage_distribution",
                endpoint_id="fibrosis_cumulative",
                set_id="training_stage_distribution--endpoint-fibrosis_cumulative",
                points=formatted_points,
            )

    def test_auditor_does_not_import_scorer_or_fit_source(self) -> None:
        source = (
            ROOT / "scripts/audit_gse267145_histology_production_v4_scores.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("import scripts.score_gse267145_histology_production_v4", source)
        self.assertNotIn("fit_gse267145_histology_baselines", source)
        self.assertIn("bootstrap_interval_values_rederived\": False", source)


if __name__ == "__main__":
    unittest.main()
