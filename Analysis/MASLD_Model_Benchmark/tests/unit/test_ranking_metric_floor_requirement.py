"""Lock the ranking-metric floor requirement and the tied-baseline annotations.

Both are forward-looking. The requirement binds metrics that have not been
scored yet; the annotations are record-only and alter no reported number. These
tests exist to keep that boundary from eroding.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
REQUIREMENT = ROOT / "config/campaigns/ranking_metric_floor_requirement_20260825.json"
ANNOTATIONS = ROOT / "config/campaigns/tied_baseline_floor_annotations_20260825.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text())


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


class FloorRequirementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = _load(REQUIREMENT)

    def test_it_is_prospective_and_alters_no_frozen_metric(self) -> None:
        self.assertIs(self.spec["applies_prospectively_only"], True)
        self.assertIs(self.spec["no_frozen_metric_is_altered"], True)
        self.assertIs(self.spec["supersedes_nothing"], True)
        self.assertIs(self.spec["metrics_calculated"], False)
        self.assertIs(self.spec["outcomes_read"], False)

    def test_both_references_are_required(self) -> None:
        requirement = self.spec["requirement"]
        self.assertIs(requirement["must_report_both_references"], True)
        self.assertIs(requirement["neither_reference_substitutes_for_the_other"], True)
        self.assertIn("reference_1_continuous_random_score", requirement)
        self.assertIn("reference_2_label_permutation", requirement)

    def test_prevalence_as_baseline_is_prohibited(self) -> None:
        self.assertIs(self.spec["requirement"]["prevalence_as_baseline_prohibited"], True)

    def test_n_and_positive_count_must_reach_the_receipt(self) -> None:
        recorded = self.spec["requirement"]["receipt_must_record"]
        for field in ("n", "positive_count", "negative_count"):
            self.assertIn(field, recorded)
        self.assertIn("distinct_score_values_emitted_by_the_scorer", recorded)

    def test_donor_peak_auprc_is_covered_before_it_is_scored(self) -> None:
        self.assertIn("donor_peak_auprc", self.spec["requirement"]["applies_to_metrics"])
        roster = self.spec["authorities"]["roster_activation_cited_not_modified"]
        self.assertEqual(roster["task_status_at_citation"], "candidate")
        self.assertEqual(roster["profile_fixture_passed_models_at_citation"], [])
        self.assertIs(roster["modified"], False)

    def test_gain_metrics_are_preferred_over_absolute_ones(self) -> None:
        preferred = self.spec["preferred_metric_form"]
        self.assertIn("gain metric", preferred["statement"])
        self.assertIs(preferred["absolute_values_require_both_references"], True)

    def test_prior_art_is_cited(self) -> None:
        prior = self.spec["prior_art_in_this_repository"]
        self.assertEqual(prior["flag"], "full_positive_set_background_AUPRC_prohibited")
        self.assertIs(prior["value"], True)
        self.assertEqual(len(prior["artifacts_sha256"]), 64)


class TiedBaselineAnnotationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = _load(ANNOTATIONS)

    def test_annotations_are_record_only(self) -> None:
        self.assertIs(self.spec["record_only"], True)
        self.assertIs(self.spec["consumed_by_runtime"], False)
        self.assertIs(self.spec["reported_metrics_unchanged"], True)
        self.assertIs(self.spec["supersedes_nothing"], True)
        self.assertIs(self.spec["metrics_calculated"], False)

    def test_every_annotated_baseline_is_a_comparator_not_a_floor(self) -> None:
        for entry in self.spec["annotated_baselines"]:
            self.assertIs(entry["is_a_legitimate_baseline"], True)
            self.assertIs(entry["is_a_valid_chance_level_floor"], False)
            self.assertIn("constant", entry["score_vector_form"])

    def test_the_exposure_direction_is_stated(self) -> None:
        direction = self.spec["defect_annotated"]["direction_matters"]
        self.assertIn("REFERENCE side", direction)
        self.assertIn("understates the floor", direction)
        self.assertIn("does not manufacture a p-value", direction)

    def test_no_interpretation_changes_today(self) -> None:
        self.assertIn(
            "protects the next result",
            self.spec["defect_annotated"]["no_interpretation_changes_today"],
        )
        for entry in self.spec["annotated_baselines"]:
            if "candidates_on_this_task_are_at_ceiling" in entry:
                self.assertIs(entry["candidates_on_this_task_are_at_ceiling"], True)

    def test_gse274114_class_counts_match_the_label_authority(self) -> None:
        entries = {
            entry["task_id"]: entry
            for entry in self.spec["annotated_baselines"]
            if "task_id" in entry
        }
        hiseq = entries["gse274114_hiseq_healthy_vs_hbv"]
        self.assertEqual(hiseq["class_counts"], {"CTRL": 9, "ENEG": 11})
        self.assertEqual(hiseq["n"], 20)
        novaseq = entries["gse274114_novaseq_mash_vs_mash_hbv"]
        self.assertEqual(novaseq["class_counts"], {"NASH": 10, "ENEG_NASH": 9})
        self.assertEqual(novaseq["n"], 19)

    def test_gse267145_has_no_auprc_to_annotate(self) -> None:
        entry = next(
            e for e in self.spec["annotated_baselines"]
            if e["model_id"] == "training_stage_distribution"
        )
        self.assertIn("no reference of any kind", entry["note"])
        self.assertIn("baseline_comparison_performed = false", entry["note"])


class CitedArtifactTests(unittest.TestCase):
    def test_every_cited_artifact_still_verifies(self) -> None:
        cited: list[tuple[str, str]] = []
        for path in (REQUIREMENT, ANNOTATIONS):
            def walk(node: object) -> None:
                if isinstance(node, dict):
                    if "path" in node and "artifacts_sha256" in node:
                        cited.append((str(node["path"]), str(node["artifacts_sha256"])))
                    for value in node.values():
                        walk(value)
                elif isinstance(node, list):
                    for value in node:
                        walk(value)

            walk(_load(path))
        self.assertGreaterEqual(len(cited), 6)
        for relative, expected in cited:
            manifest = ROOT / relative / "ARTIFACTS.json"
            self.assertTrue(manifest.is_file(), f"missing cited artifact: {relative}")
            self.assertEqual(
                _digest(manifest), expected, f"cited artifact changed: {relative}"
            )


if __name__ == "__main__":
    unittest.main()


SIDEDNESS = ROOT / "config/campaigns/phenotype_null_sidedness_annotation_20260825.json"


class NullSidednessAnnotationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = _load(SIDEDNESS)

    def test_annotation_is_record_only(self) -> None:
        self.assertIs(self.spec["record_only"], True)
        self.assertIs(self.spec["consumed_by_runtime"], False)
        self.assertIs(self.spec["reported_metrics_unchanged"], True)
        self.assertIs(self.spec["metrics_calculated"], False)

    def test_the_two_sides_are_named(self) -> None:
        defect = self.spec["defect"]
        self.assertEqual(
            defect["continuous_random_score_reference_side"], "absolute_two_sided"
        )
        self.assertEqual(
            defect["label_permutation_null_side"], "signed_one_sided"
        )
        self.assertNotEqual(
            defect["continuous_random_score_reference_null_percentile_95"],
            defect["label_permutation_null_null_percentile_95"],
        )

    def test_no_conclusion_changes_and_the_reason_is_checkable(self) -> None:
        impact = self.spec["impact_on_conclusions"]
        self.assertIs(impact["any_conclusion_changes"], False)
        self.assertIs(impact["any_q_value_changes"], False)
        self.assertIn("0.230", impact["reason"])

    def test_the_frozen_artifact_it_annotates_still_verifies(self) -> None:
        defect = self.spec["defect"]
        manifest = ROOT / defect["artifact"] / "ARTIFACTS.json"
        self.assertEqual(_digest(manifest), defect["artifacts_sha256"])

    def test_source_hash_drift_is_explained_not_hidden(self) -> None:
        note = self.spec["correction"]["source_hash_drift_note"]
        self.assertIn("no longer matches", note)
        self.assertIn("records what ran", note)
