#!/usr/bin/env python3
"""Freeze an outcome-blind sensitivity membership for pathological Scrublet thresholds."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path

from masld_bench.artifacts import verify_frozen_tree


EXPECTED_QC_ARTIFACTS = (
    "3499eed58c6da15e644f5114da33b2874e3252a86fad6d6937960424e781ecbc"
)
EXPECTED_DOUBLEt_RATE = 0.10
PATHOLOGICAL_MULTIPLIER = 3.0
PATHOLOGICAL_FLOOR = 0.25


class GSE256398QCSensitivityError(RuntimeError):
    """Raised when the frozen QC output file or sensitivity policy differs."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def pathological_auto_call(predicted: int, total: int) -> bool:
    if predicted < 0 or total < 1 or predicted > total:
        raise GSE256398QCSensitivityError("invalid doublet census")
    fraction = predicted / total
    return (
        fraction > PATHOLOGICAL_FLOOR
        and fraction > EXPECTED_DOUBLEt_RATE * PATHOLOGICAL_MULTIPLIER
    )


def capped_retained_row_ids(
    rows: list[dict[str, str]], *, expected_rate: float = EXPECTED_DOUBLEt_RATE
) -> set[str]:
    basic = [row for row in rows if row["basic_qc_retained"] == "true"]
    if not basic:
        raise GSE256398QCSensitivityError("donor has no basic-QC nuclei")
    if not 0.0 < expected_rate < 1.0:
        raise GSE256398QCSensitivityError("expected doublet rate differs")
    excluded = math.ceil(expected_rate * len(basic))
    ranked = sorted(
        basic,
        key=lambda row: (-float(row["scrublet_score"]), row["row_id"]),
    )
    return {row["row_id"] for row in ranked[excluded:]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qc-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise GSE256398QCSensitivityError("sensitivity output exists")
    verify_frozen_tree(args.qc_root)
    if sha256_file(args.qc_root / "ARTIFACTS.json") != EXPECTED_QC_ARTIFACTS:
        raise GSE256398QCSensitivityError("QC ARTIFACTS SHA-256 differs")
    receipt = json.loads((args.qc_root / "receipt.json").read_text(encoding="utf-8"))
    if (
        receipt["status"] != "pass"
        or receipt["donors"] != 26
        or receipt["nuclei_input"] != 197_942
        or receipt["final_retained"] != 165_372
        or receipt["seed"] != 0
        or receipt["disease_age_sex_or_histology_read"]
    ):
        raise GSE256398QCSensitivityError("frozen QC receipt differs")
    with (args.qc_root / "nuclei_qc.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 197_942 or len({row["row_id"] for row in rows}) != len(rows):
        raise GSE256398QCSensitivityError("nucleus row census differs")
    by_donor: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_donor.setdefault(row["source_sample_id"], []).append(row)
    if len(by_donor) != 26:
        raise GSE256398QCSensitivityError("donor census differs")

    args.output.mkdir(parents=True, exist_ok=False)
    membership_path = args.output / "expected_rate_capped_membership.tsv"
    audit_path = args.output / "donor_policy_audit.tsv"
    retained_total = 0
    pathological_donors: list[str] = []
    with membership_path.open("x", newline="", encoding="utf-8") as membership_handle, audit_path.open(
        "x", newline="", encoding="utf-8"
    ) as audit_handle:
        membership_fields = ("row_id", "gsm", "source_sample_id", "barcode")
        membership_writer = csv.DictWriter(
            membership_handle, fieldnames=membership_fields, delimiter="\t", lineterminator="\n"
        )
        membership_writer.writeheader()
        audit_fields = (
            "source_sample_id",
            "nuclei_input",
            "basic_qc_retained",
            "automatic_doublets",
            "automatic_doublet_fraction",
            "automatic_threshold",
            "pathological_auto_call",
            "sensitivity_policy",
            "sensitivity_retained",
        )
        audit_writer = csv.DictWriter(
            audit_handle, fieldnames=audit_fields, delimiter="\t", lineterminator="\n"
        )
        audit_writer.writeheader()
        donor_receipts = {
            row["source_sample_id"]: row for row in receipt["donor_receipts"]
        }
        for donor in sorted(by_donor, key=lambda value: int(value.removeprefix("S"))):
            donor_rows = by_donor[donor]
            donor_receipt = donor_receipts[donor]
            predicted = sum(row["scrublet_predicted_doublet"] == "true" for row in donor_rows)
            pathological = pathological_auto_call(predicted, len(donor_rows))
            if pathological:
                pathological_donors.append(donor)
                retained = capped_retained_row_ids(donor_rows)
                policy = "top_10pct_scrublet_score_excluded_after_basic_qc"
            else:
                retained = {
                    row["row_id"]
                    for row in donor_rows
                    if row["final_retained"] == "true"
                }
                policy = "source_exact_automatic_scrublet_then_basic_qc"
            for row in donor_rows:
                if row["row_id"] in retained:
                    membership_writer.writerow({key: row[key] for key in membership_fields})
            retained_total += len(retained)
            audit_writer.writerow(
                {
                    "source_sample_id": donor,
                    "nuclei_input": len(donor_rows),
                    "basic_qc_retained": sum(
                        row["basic_qc_retained"] == "true" for row in donor_rows
                    ),
                    "automatic_doublets": predicted,
                    "automatic_doublet_fraction": format(predicted / len(donor_rows), ".17g"),
                    "automatic_threshold": format(
                        float(donor_receipt["scrublet_threshold"]), ".17g"
                    ),
                    "pathological_auto_call": str(pathological).lower(),
                    "sensitivity_policy": policy,
                    "sensitivity_retained": len(retained),
                }
            )
    if pathological_donors != ["S35"] or retained_total != 172_997:
        raise GSE256398QCSensitivityError(
            f"sensitivity census differs: {pathological_donors}, {retained_total}"
        )
    output_receipt = {
        "schema_version": "masld-bench-gse256398-qc-sensitivity-v1",
        "status": "pass_with_prespecified_qc_sensitivity",
        "dataset_id": "gse256398",
        "donors": 26,
        "nuclei_input": 197_942,
        "source_exact_retained": 165_372,
        "sensitivity_retained": retained_total,
        "pathological_auto_call_donors": pathological_donors,
        "pathological_rule": "predicted_fraction_gt_0.25_and_gt_3x_expected_rate",
        "expected_doublet_rate": EXPECTED_DOUBLEt_RATE,
        "flagged_donor_policy": "exclude_highest_10pct_scrublet_scores_among_basic_qc_nuclei",
        "primary_membership": "source_exact_remains_primary_until_development_sensitivity_is_compared",
        "required_sensitivity": "repeat_every_GSE256398_result_with_expected_rate_capped_and_leave_S35_out_memberships",
        "selection_outcomes_read": False,
        "disease_age_sex_or_histology_read": False,
        "sealed_outcomes_read": False,
        "model_training_activated": False,
        "qc_artifacts_sha256": EXPECTED_QC_ARTIFACTS,
    }
    (args.output / "receipt.json").write_text(
        json.dumps(output_receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(output_receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
