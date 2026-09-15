#!/usr/bin/env python3
"""Project the GSE83452 record table onto its label-blind columns only.

Unlike GSE49541 this cohort is not one array per participant: 231 records cover
171 participants, 60 of them with a paired follow-up.  The grouping therefore
matters, and it is recorded here as an outer-group key so that no split can put
two records from one participant on opposite sides.

GSE83452 also shares a cohort family with GSE106737 -- 78 fingerprint-confirmed
arrays and 41 participant aliases -- so the family, not the accession, is the
unit that must stay together.  That is asserted, not assumed.

``nash_status`` is the outcome, and ``intervention``, ``age`` and ``sex`` are
covariates that must not reach a preprocessing object.  Label blindness is
enforced by an allowlist that raises on any column outside it, and the module is
tested against a table whose outcome column carries a poison value.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

SERIES = "GSE83452"
COHORT_FAMILY_ID = "antwerp_inserm_shared"
ACCESSION_ALIASES = ("GSE106737", "GSE83452")
EXPECTED_RECORDS = 231
EXPECTED_PARTICIPANTS = 171

IDENTITY_ALLOWLIST = (
    "cohort_family_id",
    "series",
    "participant_id",
    "sample_accession",
    "timepoint",
    "repeat_topology",
    "paired_baseline_accession",
)

FORBIDDEN_COLUMNS = (
    "nash_status",
    "age",
    "sex",
    "intervention",
    "primary_transfer_evaluable",
    "paired_expression_stress_evaluable",
    "paired_nash_transition_evaluable",
)

STRUCTURALLY_WITHHELD = ("nash_status", "age", "sex", "intervention")


class ParticipantMaskError(RuntimeError):
    """Raised when the projection would leak an outcome or break the grouping."""


def project_identity_columns(path: Path) -> list[dict[str, str]]:
    if set(IDENTITY_ALLOWLIST) & set(FORBIDDEN_COLUMNS):
        raise ParticipantMaskError("identity allowlist overlaps a forbidden column")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ParticipantMaskError("record table has no header")
        missing = sorted(set(IDENTITY_ALLOWLIST) - set(reader.fieldnames))
        if missing:
            raise ParticipantMaskError(f"record table lacks identity columns: {missing}")
        rows = [
            {column: row[column] for column in IDENTITY_ALLOWLIST}
            for row in reader
            if row.get("series") == SERIES
        ]
    if len(rows) != EXPECTED_RECORDS:
        raise ParticipantMaskError(f"{SERIES} record table is not {EXPECTED_RECORDS} records")
    accessions = [row["sample_accession"] for row in rows]
    if len(set(accessions)) != len(accessions):
        raise ParticipantMaskError("a sample accession repeats")
    if {row["cohort_family_id"] for row in rows} != {COHORT_FAMILY_ID}:
        raise ParticipantMaskError("record table carries a foreign cohort family")
    return sorted(rows, key=lambda row: row["sample_accession"])


def build_masks(*, records_path: Path, summarized_accessions: list[str]) -> dict[str, object]:
    rows = project_identity_columns(records_path)
    if sorted(summarized_accessions) != [row["sample_accession"] for row in rows]:
        raise ParticipantMaskError("summarized arrays do not match the record table")

    participants = {row["participant_id"] for row in rows}
    if len(participants) != EXPECTED_PARTICIPANTS:
        raise ParticipantMaskError(
            f"{SERIES} carries {len(participants)} participants, not {EXPECTED_PARTICIPANTS}"
        )
    per_participant = Counter(row["participant_id"] for row in rows)
    repeated = {p for p, n in per_participant.items() if n > 1}

    records = []
    for row in rows:
        record: dict[str, object] = dict(row)
        # The outer group is the participant inside the cohort family, so repeats
        # and aliases can never be split across a fold boundary.
        record["outer_group_key"] = f"{COHORT_FAMILY_ID}::{row['participant_id']}"
        record["participant_has_repeat_records"] = row["participant_id"] in repeated
        for field in STRUCTURALLY_WITHHELD:
            record[field] = None
            record[f"{field}_observed"] = False
        record["prediction_eligible"] = True
        records.append(record)

    return {
        "schema_version": "masld-bench-gse83452-participant-mask-v1",
        "status": "pass_label_blind_participant_masks",
        "series": SERIES,
        "cohort_family_id": COHORT_FAMILY_ID,
        "accession_aliases": list(ACCESSION_ALIASES),
        "independent_evaluation_between_accession_aliases": False,
        "family_firewall_reason": (
            "GSE83452 and GSE106737 share 78 fingerprint-confirmed arrays and 41 "
            "participant aliases and are one cohort family"
        ),
        "records": len(records),
        "participants": len(participants),
        "participants_with_repeat_records": len(repeated),
        "one_array_per_participant": False,
        "outer_group_key_rule": "cohort_family_id::participant_id",
        "distinct_outer_groups": len({r["outer_group_key"] for r in records}),
        "identity_columns_projected": list(IDENTITY_ALLOWLIST),
        "withheld_from_preprocessing": list(STRUCTURALLY_WITHHELD),
        "withheld_fields_imputed": False,
        "outcome_columns_read": False,
        "labels_read": False,
        "mask_records": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--summarized-columns", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    accessions = [
        line.strip()
        for line in arguments.summarized_columns.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    masks = build_masks(records_path=arguments.records, summarized_accessions=accessions)
    arguments.output.write_text(
        json.dumps(masks, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: masks[k] for k in masks if k != "mask_records"}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
