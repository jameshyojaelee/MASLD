#!/usr/bin/env python3
"""Validate a newly created blank Plan 46 clinical response scaffold."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from clinical_response_common import CENSUS_METRICS, read_tsv, sha256_file


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    root = Path(parser.parse_args().root).resolve()
    marker = json.loads((root / "DRAFT_EDITABLE.json").read_text(encoding="utf-8"))
    if marker.get("status") != "editable_metadata_only_clinical_response":
        raise RuntimeError("Invalid clinical response scaffold")
    for field in (
        "participant_rows_present", "clinical_outcomes_present",
        "molecular_outcomes_present", "permission_assumed",
        "analysis_authorized", "paper_promotion_authorized",
    ):
        if marker.get(field) is not False:
            raise RuntimeError(f"Clinical response firewall crossed: {field}")
    for relative, expected in marker["initial_template_sha256"].items():
        path = root / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Initial clinical response scaffold drift: {relative}")
    questions = read_tsv(root / "clinical_cohort_feasibility_questions.tsv")
    if [row["question_id"] for row in questions] != [f"CF{i:02d}" for i in range(1, 17)]:
        raise RuntimeError("Clinical response question universe drift")
    if any(row["answer"] or row["evidence_path"] or row["owner"] or row["signed_utc"] for row in questions):
        raise RuntimeError("Clinical response question fields were prefilled")
    census = read_tsv(root / "blinded_cohort_census.tsv")
    if [row["metric"] for row in census] != list(CENSUS_METRICS):
        raise RuntimeError("Clinical response census universe drift")
    if any(row["value"] or row["evidence_path"] for row in census):
        raise RuntimeError("Clinical response census was prefilled")
    identity = read_tsv(root / "clinical_response_identity.tsv")
    if len(identity) != 1 or identity[0]["response_status"] != "editable_draft":
        raise RuntimeError("Clinical response identity drift")
    if read_tsv(root / "evidence_file_register.tsv"):
        raise RuntimeError("Clinical response evidence register was prefilled")
    print(
        "CLINICAL_RESPONSE_SCAFFOLD_VALIDATION_PASS "
        f"route={marker['route_id']} questions=16 census_metrics=7 outcomes=false"
    )


if __name__ == "__main__":
    main()

