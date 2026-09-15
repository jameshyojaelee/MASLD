#!/usr/bin/env python3
"""Freeze source-native GSE256398 donor metadata without stage invention."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path


EXPECTED_SOURCE_ARTIFACTS = (
    "f54ac799b22cb02912b82432929255e6400d871cc0fd99dbf9b4756f8f977a20"
)


class GSE256398MetadataError(RuntimeError):
    """Raised when donor metadata differ from the source-native requirements."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def fibrosis_metadata(source_title: str, disease_group: str) -> dict[str, object]:
    if disease_group == "masld_f0":
        return {
            "fibrosis_source_label": "F0",
            "fibrosis_state": "observed",
            "fibrosis_min": 0,
            "fibrosis_max": 0,
            "fibrosis_numeric_use": "exact",
        }
    if disease_group == "mash_fibrosis":
        if "F2-3" in source_title:
            label, minimum, maximum, use = "F2-3", 2, 3, "interval"
        elif "F3" in source_title:
            label, minimum, maximum, use = "F3", 3, 3, "exact"
        else:
            raise GSE256398MetadataError("MASH fibrosis source label differs")
        return {
            "fibrosis_source_label": label,
            "fibrosis_state": "observed",
            "fibrosis_min": minimum,
            "fibrosis_max": maximum,
            "fibrosis_numeric_use": use,
        }
    if disease_group in {"mash_cirrhosis", "alcohol_associated_cirrhosis"}:
        return {
            "fibrosis_source_label": "cirrhosis",
            "fibrosis_state": "observed",
            "fibrosis_min": "",
            "fibrosis_max": "",
            "fibrosis_numeric_use": "not_numeric_without_source_stage",
        }
    return {
        "fibrosis_source_label": "",
        "fibrosis_state": "structurally_missing",
        "fibrosis_min": "",
        "fibrosis_max": "",
        "fibrosis_numeric_use": "not_available",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise GSE256398MetadataError("metadata output exists")
    if sha256_file(args.source / "ARTIFACTS.json") != EXPECTED_SOURCE_ARTIFACTS:
        raise GSE256398MetadataError("source ARTIFACTS SHA-256 differs")
    with (args.source / "audit/donors.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(source_rows) != 26:
        raise GSE256398MetadataError("donor census differs")
    rows: list[dict[str, object]] = []
    for source in source_rows:
        if source["age_state"] != "observed" or source["sex_state"] != "observed":
            raise GSE256398MetadataError("age or recorded sex completeness differs")
        fibrosis = fibrosis_metadata(source["source_title"], source["disease_group"])
        rows.append(
            {
                "person_key": source["person_key"],
                "gsm": source["gsm"],
                "source_sample_id": source["source_sample_id"],
                "disease_group": source["disease_group"],
                "development_role": source["development_role"],
                "age_years": source["age_years"],
                "age_state": source["age_state"],
                "recorded_sex": source["recorded_sex"],
                "sex_state": source["sex_state"],
                **fibrosis,
                "nas_value": "",
                "nas_state": "structurally_missing",
                "other_histology_state": "structurally_missing",
            }
        )
    if len({row["person_key"] for row in rows}) != 26:
        raise GSE256398MetadataError("person key is not unique")
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / "donor_metadata.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    receipt = {
        "schema_version": "masld-bench-gse256398-donor-metadata-v1",
        "status": "pass",
        "dataset_id": "gse256398",
        "donors": 26,
        "masld_relevant_donors": sum(
            row["development_role"] == "masld_relevant" for row in rows
        ),
        "alcohol_ood_donors": sum(
            row["development_role"] == "etiology_ood" for row in rows
        ),
        "age_observed_donors": sum(row["age_state"] == "observed" for row in rows),
        "recorded_sex_observed_donors": sum(
            row["sex_state"] == "observed" for row in rows
        ),
        "recorded_sex_counts": dict(
            sorted(Counter(str(row["recorded_sex"]) for row in rows).items())
        ),
        "fibrosis_source_label_counts": dict(
            sorted(
                Counter(str(row["fibrosis_source_label"]) for row in rows).items()
            )
        ),
        "nas_observed_donors": 0,
        "healthy_as_fibrosis_zero_forbidden": True,
        "cirrhosis_as_f4_without_source_forbidden": True,
        "alcohol_as_masld_negative_forbidden": True,
        "unit_of_replication": "donor",
        "labels_are_development_visible": True,
        "sealed_or_final_outcomes_read": False,
        "model_training_activated": False,
        "source_artifacts_sha256": EXPECTED_SOURCE_ARTIFACTS,
    }
    (args.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
