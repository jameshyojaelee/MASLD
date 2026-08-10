from __future__ import annotations

import sys
import unittest
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_ROOT))

from clinical_response_common import (  # noqa: E402
    ASSAY_IDS, CENSUS_METRICS, classify_portfolio, classify_response,
)


class ClinicalResponseIntakeTests(unittest.TestCase):
    @staticmethod
    def answers(value: str = "yes") -> dict[str, str]:
        return {f"CF{i:02d}": value for i in range(1, 17)}

    @staticmethod
    def assays(value: str = "yes") -> dict[str, str]:
        return {assay: value for assay in ASSAY_IDS}

    @staticmethod
    def census() -> dict[str, int]:
        values = {metric: 20 for metric in CENSUS_METRICS}
        values["unique_participants"] = 30
        values["histologic_improvers"] = 8
        values["histologic_non_improvers"] = 12
        return values

    def test_complete_response_passes_census_and_nested_gate(self) -> None:
        result = classify_response(self.answers(), self.assays(), self.census(), True, True)
        self.assertEqual(result["feasibility_verdict"], "pass_to_blinded_census")
        self.assertTrue(result["full_goal_nested_candidate"])
        self.assertFalse(result["molecular_outcome_access_authorized"])

    def test_missing_orthogonal_assay_preserves_human_gate_only(self) -> None:
        assays = self.assays()
        assays["secreted_proteomics"] = "no"
        answers = self.answers()
        answers["CF10"] = "no"
        result = classify_response(answers, assays, self.census(), True, True)
        self.assertEqual(result["feasibility_verdict"], "pass_to_blinded_census")
        self.assertFalse(result["full_goal_nested_candidate"])

    def test_one_explicitly_recoverable_gate_is_conditional(self) -> None:
        answers = self.answers()
        answers["CF08"] = "unresolved"
        result = classify_response(
            answers, self.assays(), self.census(), True, True,
            missing_gate_recoverable=True, missing_gate_id="CF08",
        )
        self.assertEqual(result["feasibility_verdict"], "conditional_missing_one_gate")
        self.assertFalse(result["advance_to_overlap_power_audit"])

    def test_missing_rna_is_supportive_only_when_orthogonal_data_exist(self) -> None:
        answers = self.answers()
        answers["CF05"] = "no"
        assays = self.assays()
        assays["bulk_rna"] = "no"
        result = classify_response(answers, assays, self.census(), True, True)
        self.assertEqual(result["feasibility_verdict"], "supportive_only")
        self.assertFalse(result["pass_to_blinded_census"])

    def test_absent_response_direction_fails(self) -> None:
        census = self.census()
        census["histologic_non_improvers"] = 0
        result = classify_response(self.answers(), self.assays(), census, True, True)
        self.assertEqual(result["feasibility_verdict"], "fail_inaccessible_or_unpaired")

    def test_two_independent_passes_are_still_outcome_locked(self) -> None:
        cohort_a = {
            "candidate_cohort_role": "candidate_cohort_A", "route_id": "a",
            "independence_group": "group_a", "institution_uid": "inst_a",
            "trial_uid": "trial_a", "pass_to_blinded_census": True,
            "full_goal_nested_candidate": False,
        }
        cohort_b = {
            "candidate_cohort_role": "candidate_cohort_B", "route_id": "b",
            "independence_group": "group_b", "institution_uid": "inst_b",
            "trial_uid": "trial_b", "pass_to_blinded_census": True,
            "full_goal_nested_candidate": True,
        }
        independence = {
            "institutionally_independent": True, "trial_independent": True,
            "participant_independent": True, "overlap_resolved": True,
            "cohort_b_locked_before_cohort_a_outcomes": True,
        }
        result = classify_portfolio(cohort_a, cohort_b, independence)
        self.assertTrue(result["two_cohort_source_gate_pass"])
        self.assertTrue(result["orthogonal_nested_source_available"])
        self.assertFalse(result["molecular_outcome_access_authorized"])

    def test_unresolved_overlap_fails_portfolio(self) -> None:
        cohort_a = {
            "candidate_cohort_role": "candidate_cohort_A", "route_id": "a",
            "independence_group": "ga", "institution_uid": "ia", "trial_uid": "ta",
            "pass_to_blinded_census": True, "full_goal_nested_candidate": True,
        }
        cohort_b = {
            "candidate_cohort_role": "candidate_cohort_B", "route_id": "b",
            "independence_group": "gb", "institution_uid": "ib", "trial_uid": "tb",
            "pass_to_blinded_census": True, "full_goal_nested_candidate": True,
        }
        independence = {
            "institutionally_independent": True, "trial_independent": True,
            "participant_independent": True, "overlap_resolved": False,
            "cohort_b_locked_before_cohort_a_outcomes": True,
        }
        self.assertFalse(classify_portfolio(cohort_a, cohort_b, independence)["two_cohort_source_gate_pass"])


if __name__ == "__main__":
    unittest.main()

