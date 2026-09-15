from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from masld_bench.artifacts import verify_frozen_tree
from masld_bench.contracts import PredictionBundle
from masld_bench.hashing import sha256_file

import scripts.score_gse267145_histology_production_v4 as scorer


ROOT = Path(__file__).resolve().parents[2]
OUTCOMES = ROOT / "executions/model-data-061-21079623/activation/participant_endpoints.tsv"
FOLDS = ROOT / "executions/model-data-064-21079902/fixture/folds/participant_outer_folds.tsv"
CONTRACT = ROOT / "config/evaluation/gse267145_histology_production_v4_scoring.json"


def perfect_rows() -> list[dict[str, str]]:
    _, endpoints = scorer._read_tsv(OUTCOMES)
    rows: list[dict[str, str]] = []
    for endpoint in endpoints:
        fibrosis = int(endpoint["fibrosis"])
        fibrosis_group = "F0" if fibrosis == 0 else "F1" if fibrosis == 1 else "F2_3"
        rows.append(
            {
                "participant_id": endpoint["participant_id"],
                "outer_fold": endpoint["outer_fold"],
                **{
                    f"probability_{label}": "1" if label == endpoint["stage3"] else "0"
                    for label in scorer.reference.STAGE3
                },
                "predicted_stage3": endpoint["stage3"],
                "predicted_nash_crn_component_sum": endpoint[
                    "nash_crn_component_sum"
                ],
                "predicted_fibrosis_cumulative_expected": endpoint["fibrosis"],
                "predicted_fibrosis_regression": endpoint["fibrosis"],
                **{
                    f"probability_fibrosis_{label}": "1"
                    if label == fibrosis_group
                    else "0"
                    for label in scorer.reference.FIBROSIS_GROUP3
                },
                "predicted_fibrosis_group3": fibrosis_group,
            }
        )
    return rows


def write_complete_prediction_roster(root: Path, rows: list[dict[str, str]]) -> None:
    root.mkdir()
    for filename in sorted(scorer._production_prediction_files()):
        scorer._write_tsv(root / filename, scorer.reference.PREDICTION_FIELDS, rows)


class GSE267145ProductionV4ScorerTests(unittest.TestCase):
    def test_contract_accepts_only_root_bound_production_hash(self) -> None:
        contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
        production_sha256 = "a" * 64
        contract["status"] = "bound_root_verified_production_artifacts"
        contract["production_campaign"]["artifacts_sha256"] = production_sha256
        contract["independent_source_validation"] = {
            "path": "executions/model-check-079-test",
            "artifacts_sha256": "c" * 64,
            "status": "passed_unscored_source_and_unit_validation",
            "production_scores_calculated": False,
        }
        scorer._validate_contract(
            contract, production_artifacts_sha256=production_sha256
        )
        changed = deepcopy(contract)
        changed["claim_boundary"]["champion_claim_allowed"] = True
        with self.assertRaises(scorer.HistologyProductionScoringError):
            scorer._validate_contract(
                changed, production_artifacts_sha256=production_sha256
            )

    def test_freezes_55_endpoint_bundles_and_11_non_bundle_dispositions(self) -> None:
        rows = perfect_rows()
        _, fold_rows = scorer._read_tsv(FOLDS)
        participants = [row["participant_id"] for row in fold_rows]
        outer_folds = [row["outer_fold"] for row in fold_rows]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            write_complete_prediction_roster(source, rows)
            output = root / "bundles"
            records, artifacts_sha256 = scorer._freeze_prediction_bundles(
                source_root=source,
                output=output,
                production_sha256="b" * 64,
                participants=participants,
                outer_folds=outer_folds,
            )
            self.assertEqual(len(records), 55)
            self.assertTrue(all("--seed-" not in set_id for set_id in records))
            self.assertEqual(sha256_file(output / "ARTIFACTS.json"), artifacts_sha256)
            self.assertFalse(verify_frozen_tree(output)["metadata"]["outcomes_read"])
            first = output / "training_stage_distribution--endpoint-stage3"
            bundle = PredictionBundle.load_json(first / "prediction_bundle.json")
            bundle.validate_artifacts(first)
            self.assertEqual(bundle.biological_unit, "participant")
            self.assertEqual(bundle.n_predictions, 99)
            self.assertEqual(bundle.metadata["endpoint_id"], "stage3")
            self.assertEqual(bundle.metadata["prediction_aggregation"], scorer.SEED_AGGREGATION)
            self.assertEqual(
                tuple(bundle.metadata["model_seed_roster"]), scorer.reference.SEEDS
            )
            self.assertFalse(bundle.metadata["champion_claim_allowed"])
            scoring_predictions = scorer._load_frozen_endpoint_predictions(
                prediction_root=output,
                prediction_records=records,
                production_sha256="b" * 64,
                participants=participants,
                outer_folds=outer_folds,
            )
            self.assertEqual(set(scoring_predictions), set(scorer.reference.MODEL_IDS))
            self.assertEqual(
                scorer._point_metrics(
                    scorer._load_outcomes(OUTCOMES.parent, participants, outer_folds),
                    scoring_predictions["training_stage_distribution"],
                )["stage3_macro_f1"],
                1.0,
            )

            dispositions = root / "dispositions"
            disposition_rows, _ = scorer._freeze_source_stage5_dispositions(
                output=dispositions,
                production_sha256="b" * 64,
            )
            self.assertEqual(len(disposition_rows), 11)
            self.assertFalse(any(dispositions.rglob("prediction_bundle.json")))
            document = json.loads(
                (
                    dispositions
                    / "training_stage_distribution/endpoint_disposition.json"
                ).read_text(encoding="utf-8")
            )
            self.assertFalse(document["prediction_bundle_created"])
            self.assertEqual(document["disposition"], scorer.SOURCE_STAGE5_DISPOSITION)

    def test_five_seed_ensemble_must_be_exact_before_outcome_read(self) -> None:
        rows = perfect_rows()
        _, fold_rows = scorer._read_tsv(FOLDS)
        participants = [row["participant_id"] for row in fold_rows]
        outer_folds = [row["outer_fold"] for row in fold_rows]
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source"
            write_complete_prediction_roster(source, rows)
            changed = deepcopy(rows)
            changed[0]["predicted_nash_crn_component_sum"] = "0.5"
            bad = source / "training_stage_distribution.tsv"
            bad.unlink()
            scorer._write_tsv(bad, scorer.reference.PREDICTION_FIELDS, changed)
            with self.assertRaises(scorer.HistologyProductionScoringError):
                scorer._verify_five_seed_ensembles(
                    source_root=source,
                    participants=participants,
                    outer_folds=outer_folds,
                )

    def test_rectangular_metric_mapping_is_495_99_396(self) -> None:
        rows = perfect_rows()
        _, fold_rows = scorer._read_tsv(FOLDS)
        participants = [row["participant_id"] for row in fold_rows]
        outer_folds = [row["outer_fold"] for row in fold_rows]
        outcomes = scorer._load_outcomes(OUTCOMES.parent, participants, outer_folds)
        predictions = scorer._load_predictions(rows)
        indices = np.random.default_rng(scorer.BOOTSTRAP_SEED).integers(
            0, 99, size=(25, 99), endpoint=False
        )
        points = scorer._point_metrics(outcomes, predictions)
        distributions = scorer._bootstrap_metrics(outcomes, predictions, indices)
        rows_per_model = []
        for endpoint_id in scorer.PREDICTED_ENDPOINTS:
            metric_rows = scorer._metric_rows(
                model_id="training_stage_distribution",
                endpoint_id=endpoint_id,
                set_id=f"training_stage_distribution--endpoint-{endpoint_id}",
                points=points,
                distributions=distributions,
            )
            self.assertEqual(len(metric_rows), 9)
            rows_per_model.extend(metric_rows)
        self.assertEqual(
            sum(row["applicability_state"] == "observed" for row in rows_per_model),
            9,
        )
        self.assertEqual(
            sum(
                row["applicability_state"] == "not_applicable"
                for row in rows_per_model
            ),
            36,
        )
        secondary = [row for row in rows_per_model if row["endpoint_role"] == "secondary"]
        self.assertTrue(
            all(row["multiplicity_family"] == scorer.MULTIPLICITY_FAMILY for row in secondary)
        )
        self.assertTrue(
            all(row["p_value"] == scorer.NO_CONFIRMATORY_STATE for row in secondary)
        )
        self.assertTrue(
            all(
                row["bh_adjusted_q_value"] == scorer.NO_CONFIRMATORY_STATE
                for row in secondary
            )
        )
        self.assertEqual(55 * 9, 495)
        self.assertEqual(11 * 9, 99)
        self.assertEqual(11 * 36, 396)

    def test_scorer_has_no_fit_import_or_ranking_surface(self) -> None:
        source = (
            ROOT / "scripts/score_gse267145_histology_production_v4.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("fit_gse267145_histology_baselines", source)
        self.assertNotIn("pareto", source.lower())
        self.assertNotIn("strongest_by_metric", source)
        post_outcome_source = source.split(
            '_verify_file(outcomes / "ARTIFACTS.json"', maxsplit=1
        )[1]
        self.assertNotIn('production_state["prediction_root"]', post_outcome_source)
        self.assertIn(
            '"metric_scoring_prediction_source": "55_frozen_endpoint_prediction_bundles"',
            post_outcome_source,
        )


if __name__ == "__main__":
    unittest.main()
