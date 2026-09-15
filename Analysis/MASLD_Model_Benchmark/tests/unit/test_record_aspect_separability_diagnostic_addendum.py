"""An addendum recorded late must say it was recorded late.

The whole value of this record is that it does not pretend to have come first.
These tests pin the properties that make it honest rather than decorative: it
states the ordering it actually had, it cannot alter a frozen verdict, it
registers no gate, and its declared decision role stays asymmetric so it can
never rescue a failed criterion.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "record_aspect_separability_diagnostic_addendum.py"
SPEC = importlib.util.spec_from_file_location("addendum_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
addendum_module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(addendum_module)


def _fake_result(decision: str = "STOP") -> dict:
    reliability = {
        aspect: {"spearman_brown_full_length_reliability": value}
        for aspect, value in (
            ("steatosis", 0.553),
            ("ballooning", 0.479),
            ("lobular_inflammation", 0.518),
        )
    }
    return {
        "decision": decision,
        "criteria": {
            "c1_labels_are_not_redundant": {
                "pairs": [
                    {"pair": ["steatosis", "ballooning"], "spearman": 0.68},
                    {"pair": ["steatosis", "lobular_inflammation"], "spearman": 0.7284},
                    {
                        "pair": ["ballooning", "lobular_inflammation"],
                        "spearman": 0.8181,
                    },
                ]
            }
        },
        "diagnostics_not_criteria": {
            "split_half_reliability": {"per_aspect": reliability},
            "contrast_axes_through_the_same_pipeline": {
                "pairs": [
                    {
                        "contrast_axis": "fibrosis",
                        "aspect": "steatosis",
                        "spearman_of_association_vectors": 0.5665,
                    }
                ]
            },
        },
    }


def _build(decision: str = "STOP"):
    return addendum_module.build(
        prespec_digest="a" * 64,
        result_digest="b" * 64,
        result=_fake_result(decision),
        contrast_reliability={
            "fibrosis": {"spearman_brown_full_length_reliability": 0.51},
            "lobular_necrosis": {"spearman_brown_full_length_reliability": 0.44},
        },
        analysis_job_id=21181334,
        analysis_end="2026-08-27T15:12:41",
    )


class OrderingHonestyTests(unittest.TestCase):
    def setUp(self) -> None:
        self.addendum, self.observations = _build()

    def test_it_does_not_claim_to_have_come_first(self) -> None:
        ordering = self.addendum["honesty_about_ordering"]
        self.assertIs(ordering["frozen_before_the_analysis"], False)
        self.assertIs(ordering["recorded_after_the_analysis"], True)
        self.assertIs(ordering["in_force_for_this_decision"], False)
        self.assertEqual(ordering["analysis_job_id"], 21181334)

    def test_it_explains_why_it_could_not_have_changed_the_decision(self) -> None:
        ordering = self.addendum["honesty_about_ordering"]
        self.assertIs(ordering["could_it_have_changed_this_decision"], False)
        self.assertIn("STOP", ordering["why_it_could_not"])

    def test_it_alters_no_frozen_verdict_and_registers_no_gate(self) -> None:
        self.assertIs(self.addendum["frozen_verdict_altered"], False)
        self.assertIs(self.addendum["promotion_gate_registered"], False)
        self.assertIs(self.addendum["record_is_additive_overlay"], True)

    def test_it_names_the_artifacts_it_covers_by_digest(self) -> None:
        covered = self.addendum["frozen_artifacts_this_record_covers"]
        self.assertEqual(len(covered["prespecification_artifacts_sha256"]), 64)
        self.assertEqual(len(covered["result_artifacts_sha256"]), 64)


class DecisionRoleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.addendum, _ = _build()

    def test_the_role_is_asymmetric(self) -> None:
        role = self.addendum["declared_decision_role"]
        self.assertIs(role["may_downgrade_a_go_to_indeterminate"], True)
        self.assertIs(role["may_never_upgrade_a_stop_to_a_go"], True)

    def test_the_threshold_is_marked_a_judgment_call(self) -> None:
        self.assertIs(
            self.addendum["declared_decision_role"]["threshold_is_a_judgment_call"],
            True,
        )
        self.assertIs(
            self.addendum["reliability_gate"]["threshold_is_a_judgment_call"], True
        )

    def test_the_proposed_form_is_recorded_verbatim_not_silently_dropped(self) -> None:
        gate = self.addendum["reliability_gate"]
        self.assertIn("geometric mean", gate["proposed_form_recorded_verbatim"])
        self.assertIn("attenuation ceiling", gate["why_the_proposed_form_was_not_adopted"])

    def test_the_counterexample_is_measured_not_asserted(self) -> None:
        example = self.addendum["reliability_gate"]["measured_counterexample"]
        expected = (0.51 * 0.553) ** 0.5
        self.assertAlmostEqual(
            example["geometric_mean_ceiling_for_fibrosis_vs_steatosis"],
            expected,
            places=12,
        )
        self.assertGreater(abs(example["observed_fibrosis_vs_steatosis"]), expected)
        self.assertIs(example["observed_exceeds_the_proposed_ceiling"], True)

    def test_severity_residualization_is_reported_never_gated(self) -> None:
        role = self.addendum["severity_residualized_sensitivity_role"]
        self.assertIs(role["reported"], True)
        self.assertIs(role["gated_on"], False)


class PostHocObservationTests(unittest.TestCase):
    def setUp(self) -> None:
        _, self.observations = _build()

    def test_it_is_marked_post_hoc_and_authorises_nothing(self) -> None:
        self.assertIs(self.observations["explicitly_post_hoc"], True)
        self.assertIs(self.observations["derived_by_looking_at_the_answer"], True)
        self.assertIs(self.observations["authorises_nothing"], True)
        self.assertIs(self.observations["no_work_was_built_toward_it"], True)

    def test_it_records_the_pairwise_values_it_is_reading(self) -> None:
        values = self.observations["observation"]["c1_pairwise_spearman"]
        self.assertEqual(values["ballooning|lobular_inflammation"], 0.8181)
        self.assertEqual(len(values), 3)

    def test_it_defers_the_decision_to_open_it(self) -> None:
        self.assertIn(
            "campaign lead", self.observations["observation"]["who_decides_whether_to_open_it"]
        )
        self.assertIn(
            "own prespecification", self.observations["observation"]["why_it_is_not_a_plan"]
        )


if __name__ == "__main__":
    unittest.main()
