#!/usr/bin/env python3
"""Project the GSE49541 participant table onto its label-blind columns only.

The frozen TaskSpec keeps exact fibrosis stage, age, and sex as structurally
missing and forbids imputing them from the group label, publication summaries,
array intensity, or another cohort.  The deposited fibrosis group itself is an
outcome and must not reach preprocessing.

Label blindness is enforced here rather than promised: the reader carries an
allowlist of identity columns, and any attempt to project an outcome-bearing or
undeclared column raises instead of returning data.  Structurally missing fields
are emitted as an explicit mask with a null value, never as a filled cell.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

EXPECTED_PARTICIPANTS = 72
SERIES = "GSE49541"
COHORT_FAMILY_ID = "gse31803_gse49541_fibrosis_array"

# The only columns preprocessing may see.  Everything else in the participant
# table is an outcome, a derived eligibility flag, or a structurally missing
# field, and none of those may enter a preprocessing object.
IDENTITY_ALLOWLIST = (
    "cohort_family_id",
    "series",
    "participant_id",
    "sample_accession",
    "timepoint",
    "repeat_topology",
)

# Named so that a future edit that tries to widen the allowlist fails loudly.
FORBIDDEN_COLUMNS = (
    "fibrosis_stage_group",
    "exact_fibrosis_stage",
    "nash_status",
    "age",
    "sex",
    "intervention",
    "primary_transfer_evaluable",
    "paired_expression_stress_evaluable",
    "paired_nash_transition_evaluable",
)

STRUCTURALLY_MISSING = ("exact_fibrosis_stage", "age", "sex")


class ParticipantMaskError(RuntimeError):
    """Raised when the participant projection would leak an outcome."""


def project_identity_columns(path: Path) -> list[dict[str, str]]:
    if set(IDENTITY_ALLOWLIST) & set(FORBIDDEN_COLUMNS):
        raise ParticipantMaskError("identity allowlist overlaps a forbidden column")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ParticipantMaskError("participant table has no header")
        missing = sorted(set(IDENTITY_ALLOWLIST) - set(reader.fieldnames))
        if missing:
            raise ParticipantMaskError(f"participant table lacks identity columns: {missing}")
        rows = [
            {column: row[column] for column in IDENTITY_ALLOWLIST}
            for row in reader
            if row.get("series") == SERIES
        ]
    if len(rows) != EXPECTED_PARTICIPANTS:
        raise ParticipantMaskError(
            f"{SERIES} participant table is not {EXPECTED_PARTICIPANTS} records"
        )
    accessions = [row["sample_accession"] for row in rows]
    participants = [row["participant_id"] for row in rows]
    if len(set(accessions)) != len(accessions):
        raise ParticipantMaskError("a sample accession repeats")
    if len(set(participants)) != len(participants):
        raise ParticipantMaskError("a participant identifier repeats")
    if {row["cohort_family_id"] for row in rows} != {COHORT_FAMILY_ID}:
        raise ParticipantMaskError("participant table carries a foreign cohort family")
    return sorted(rows, key=lambda row: row["sample_accession"])


def build_masks(*, participants_path: Path, summarized_accessions: list[str]) -> dict[str, object]:
    rows = project_identity_columns(participants_path)
    if sorted(summarized_accessions) != [row["sample_accession"] for row in rows]:
        raise ParticipantMaskError("summarized arrays do not match the participant table")

    records = []
    for row in rows:
        record: dict[str, object] = dict(row)
        for field in STRUCTURALLY_MISSING:
            record[field] = None
            record[f"{field}_observed"] = False
        record["prediction_eligible"] = True
        records.append(record)

    return {
        "schema_version": "masld-bench-gse49541-participant-mask-v1",
        "status": "pass_label_blind_participant_masks",
        "series": SERIES,
        "cohort_family_id": COHORT_FAMILY_ID,
        "participants": len(records),
        "arrays": len(summarized_accessions),
        "one_array_per_participant": True,
        "identity_columns_projected": list(IDENTITY_ALLOWLIST),
        "structurally_missing_fields": list(STRUCTURALLY_MISSING),
        "structurally_missing_imputed": False,
        "outcome_columns_read": False,
        "labels_read": False,
        "records": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--participants", type=Path, required=True)
    parser.add_argument("--summarized-columns", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    accessions = [
        line.strip()
        for line in arguments.summarized_columns.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    masks = build_masks(
        participants_path=arguments.participants, summarized_accessions=accessions
    )
    arguments.output.write_text(
        json.dumps(masks, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {key: masks[key] for key in masks if key != "records"}, indent=2, sort_keys=True
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
