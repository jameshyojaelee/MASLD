from __future__ import annotations

import copy
import unittest

from scripts.build_gse83452_baseline_selection import (
    BaselineSelectionError,
    select_baseline,
)


def _record(*, index: int, participant: str, timepoint: str, topology: str) -> dict:
    return {
        "cohort_family_id": "antwerp_inserm_shared",
        "series": "GSE83452",
        "participant_id": participant,
        "sample_accession": f"GSM{2203254 + index}",
        "timepoint": timepoint,
        "repeat_topology": topology,
        "paired_baseline_accession": f"GSM{2203254 + index}",
        "outer_group_key": f"antwerp_inserm_shared::{participant}",
        "participant_has_repeat_records": topology != "followup_only_unpaired",
        "prediction_eligible": True,
        "nash_status": None,
        "nash_status_observed": False,
        "age": None,
        "age_observed": False,
        "sex": None,
        "sex_observed": False,
        "intervention": None,
        "intervention_observed": False,
    }


def make_masks(*, baseline: int = 152, paired: int = 60, unpaired: int = 19) -> dict:
    """Mirror the deposited topology: 231 records over 171 participants.

    152 baseline anchors, 60 one-year follow-ups that reuse a baseline
    participant, and 19 follow-up-only participants with no baseline array.
    """

    mask_records = []
    index = 0
    baseline_participants = [f"gse83452::GSM{2203254 + i}" for i in range(baseline)]
    for participant in baseline_participants:
        mask_records.append(
            _record(
                index=index,
                participant=participant,
                timepoint="baseline",
                topology="baseline_only_or_pair_anchor",
            )
        )
        index += 1
    for offset in range(paired):
        mask_records.append(
            _record(
                index=index,
                participant=baseline_participants[offset],
                timepoint="follow-up",
                topology="paired_baseline_one_year_followup",
            )
        )
        index += 1
    for offset in range(unpaired):
        mask_records.append(
            _record(
                index=index,
                participant=f"gse83452::followup{offset:03d}",
                timepoint="follow-up",
                topology="followup_only_unpaired",
            )
        )
        index += 1
    return {
        "status": "pass_label_blind_participant_masks",
        "series": "GSE83452",
        "cohort_family_id": "antwerp_inserm_shared",
        "labels_read": False,
        "outcome_columns_read": False,
        "withheld_fields_imputed": False,
        "family_firewall_reason": "one cohort family",
        "mask_records": mask_records,
    }


class BaselineSelectionTests(unittest.TestCase):
    def test_selection_returns_one_array_per_participant(self) -> None:
        selection = select_baseline(make_masks())
        self.assertEqual(selection["records"], 152)
        self.assertEqual(selection["participants"], 152)
        self.assertEqual(selection["distinct_outer_groups"], 152)
        self.assertEqual(selection["followup_records_excluded"], 79)
        self.assertEqual(selection["arrays_per_participant_after_selection"], 1)

    def test_selection_reads_only_timepoint_and_topology(self) -> None:
        selection = select_baseline(make_masks())
        self.assertEqual(
            selection["selection_fields_read"], ["timepoint", "repeat_topology"]
        )
        self.assertTrue(selection["selection_is_label_blind"])
        for record in selection["mask_records"]:
            for field in ("nash_status", "age", "sex", "intervention"):
                self.assertIsNone(record[field])
                self.assertFalse(record[f"{field}_observed"])

    def test_a_filled_outcome_is_rejected(self) -> None:
        masks = make_masks()
        masks["mask_records"][0]["nash_status"] = "nash"
        with self.assertRaises(BaselineSelectionError):
            select_baseline(masks)

    def test_a_covariate_declared_observed_is_rejected(self) -> None:
        masks = make_masks()
        masks["mask_records"][3]["age_observed"] = True
        with self.assertRaises(BaselineSelectionError):
            select_baseline(masks)

    def test_two_baseline_arrays_on_one_participant_are_rejected(self) -> None:
        """Otherwise the downstream participant bootstrap is an array bootstrap."""

        masks = make_masks()
        # Point a second baseline array at participant 0.  Participant 1 keeps
        # its identity through its own paired follow-up row, so 231 records over
        # 171 participants still holds and the guard under test is the only one
        # that can fire.
        masks["mask_records"][1]["participant_id"] = masks["mask_records"][0][
            "participant_id"
        ]
        masks["mask_records"][1]["outer_group_key"] = masks["mask_records"][0][
            "outer_group_key"
        ]
        self.assertEqual(
            len({record["participant_id"] for record in masks["mask_records"]}), 171
        )
        with self.assertRaises(BaselineSelectionError):
            select_baseline(masks)

    def test_an_outer_group_that_is_not_family_scoped_is_rejected(self) -> None:
        masks = make_masks()
        masks["mask_records"][0]["outer_group_key"] = "gse83452::alone"
        with self.assertRaises(BaselineSelectionError):
            select_baseline(masks)

    def test_a_mask_that_is_not_label_blind_is_rejected(self) -> None:
        masks = make_masks()
        masks["labels_read"] = True
        with self.assertRaises(BaselineSelectionError):
            select_baseline(masks)

    def test_a_short_baseline_census_is_rejected(self) -> None:
        masks = make_masks(baseline=151, unpaired=20)
        with self.assertRaises(BaselineSelectionError):
            select_baseline(masks)

    def test_selection_does_not_mutate_its_input(self) -> None:
        masks = make_masks()
        before = copy.deepcopy(masks)
        select_baseline(masks)
        self.assertEqual(masks, before)


if __name__ == "__main__":
    unittest.main()
