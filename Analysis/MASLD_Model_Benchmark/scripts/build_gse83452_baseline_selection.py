#!/usr/bin/env python3
"""Select the GSE83452 baseline arrays from the frozen label-blind mask.

The primary transfer endpoint is baseline NASH status.  GSE83452 also deposits
79 one-year follow-up arrays taken after a bariatric-surgery or diet
intervention; those belong to a separate frozen paired/intervention task and are
not this lane's deliverable.

Selection reads ``timepoint`` and ``repeat_topology`` and nothing else, so it is
label-blind by construction: the mask this consumes has already replaced
``nash_status``, ``age``, ``sex`` and ``intervention`` with explicit nulls, and
this module re-checks that rather than trusting it.

Cohort-family grouping survives the selection.  GSE83452 and GSE106737 are one
project-exposed family sharing 78 fingerprint-confirmed arrays and 41
participant aliases, so the outer group stays ``cohort_family_id::participant_id``
and every alias of one participant sits in one group.  On the baseline subset
that grouping is checked to be one array per group, which is what makes an
unstratified participant bootstrap downstream a participant bootstrap rather
than an array bootstrap.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Any

SERIES = "GSE83452"
COHORT_FAMILY_ID = "antwerp_inserm_shared"
ACCESSION_ALIASES = ("GSE106737", "GSE83452")
BASELINE_TIMEPOINT = "baseline"
EXPECTED_ALL_RECORDS = 231
EXPECTED_ALL_PARTICIPANTS = 171
EXPECTED_BASELINE_RECORDS = 152
STRUCTURALLY_WITHHELD = ("nash_status", "age", "sex", "intervention")
FORBIDDEN_VALUE_TOKENS = ("nash", "no_nash", "undefined", "BS", "Diet")


class BaselineSelectionError(RuntimeError):
    """Raised when baseline selection would leak a label or break the grouping."""


def select_baseline(masks: dict[str, Any]) -> dict[str, Any]:
    if (
        masks.get("status") != "pass_label_blind_participant_masks"
        or masks.get("labels_read") is not False
        or masks.get("outcome_columns_read") is not False
        or masks.get("withheld_fields_imputed") is not False
        or masks.get("series") != SERIES
        or masks.get("cohort_family_id") != COHORT_FAMILY_ID
    ):
        raise BaselineSelectionError("participant mask is not the label-blind GSE83452 mask")
    records = masks.get("mask_records")
    if not isinstance(records, list) or len(records) != EXPECTED_ALL_RECORDS:
        raise BaselineSelectionError("participant mask record census differs")
    if len({record["participant_id"] for record in records}) != EXPECTED_ALL_PARTICIPANTS:
        raise BaselineSelectionError("participant mask participant census differs")

    for record in records:
        for field in STRUCTURALLY_WITHHELD:
            if record.get(field, "sentinel") is not None:
                raise BaselineSelectionError(f"{field} is filled in the mask")
            if record.get(f"{field}_observed", "sentinel") is not False:
                raise BaselineSelectionError(f"{field} is declared observed in the mask")
        values = {value for value in record.values() if isinstance(value, str)}
        if set(FORBIDDEN_VALUE_TOKENS) & values:
            raise BaselineSelectionError("participant mask leaked an outcome value")

    topology = Counter(record["repeat_topology"] for record in records)
    baseline = sorted(
        (record for record in records if record["timepoint"] == BASELINE_TIMEPOINT),
        key=lambda record: record["participant_id"],
    )
    if len(baseline) != EXPECTED_BASELINE_RECORDS:
        raise BaselineSelectionError(
            f"baseline selection returned {len(baseline)} records, not "
            f"{EXPECTED_BASELINE_RECORDS}"
        )
    participants = [record["participant_id"] for record in baseline]
    if len(set(participants)) != len(participants):
        raise BaselineSelectionError(
            "a participant contributes more than one baseline array; the "
            "participant bootstrap downstream would be an array bootstrap"
        )
    groups = [record["outer_group_key"] for record in baseline]
    if any(
        group != f"{COHORT_FAMILY_ID}::{participant}"
        for group, participant in zip(groups, participants, strict=True)
    ):
        raise BaselineSelectionError("baseline outer group key is not family-scoped")
    if len(set(groups)) != len(groups):
        raise BaselineSelectionError("two baseline arrays share one outer group")
    excluded = [
        record for record in records if record["timepoint"] != BASELINE_TIMEPOINT
    ]
    if {record["timepoint"] for record in excluded} != {"follow-up"}:
        raise BaselineSelectionError("an excluded record is not a follow-up visit")

    selected = []
    for record in baseline:
        projected = dict(record)
        projected["baseline_selected"] = True
        projected["prediction_eligible"] = True
        selected.append(projected)

    return {
        "schema_version": "masld-bench-gse83452-baseline-selection-v1",
        "status": "pass_label_blind_participant_masks",
        "series": SERIES,
        "cohort_family_id": COHORT_FAMILY_ID,
        "accession_aliases": list(ACCESSION_ALIASES),
        "independent_evaluation_between_accession_aliases": False,
        "family_firewall_reason": masks["family_firewall_reason"],
        "selection_rule": "timepoint == 'baseline'",
        "selection_fields_read": ["timepoint", "repeat_topology"],
        "selection_is_label_blind": True,
        "records_before_selection": len(records),
        "participants_before_selection": EXPECTED_ALL_PARTICIPANTS,
        "repeat_topology_before_selection": {
            key: int(topology[key]) for key in sorted(topology)
        },
        "records": len(selected),
        "participants": len(set(participants)),
        "arrays_per_participant_after_selection": 1,
        "followup_records_excluded": len(excluded),
        "followup_exclusion_reason": (
            "the one-year post-intervention visits belong to the separate frozen "
            "gse83452_paired_intervention_stress task and are not this endpoint"
        ),
        "outer_group_key_rule": "cohort_family_id::participant_id",
        "distinct_outer_groups": len(set(groups)),
        "withheld_from_preprocessing": list(STRUCTURALLY_WITHHELD),
        "withheld_fields_imputed": False,
        "outcome_columns_read": False,
        "labels_read": False,
        "mask_records": selected,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--masks", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--output-columns", required=True, type=Path)
    arguments = parser.parse_args()
    masks = json.loads(arguments.masks.read_text(encoding="utf-8"))
    selection = select_baseline(masks)
    arguments.output.write_text(
        json.dumps(selection, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    arguments.output_columns.write_text(
        "\n".join(record["sample_accession"] for record in selection["mask_records"])
        + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {key: selection[key] for key in selection if key != "mask_records"},
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
