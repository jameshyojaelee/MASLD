#!/usr/bin/env python3
"""Diagnose deposited visit-token differences for frozen array aliases."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path


class CrossAccessionTimepointError(RuntimeError):
    """Raised when the label-blind timepoint diagnostic differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def token_disposition(left: str, right: str) -> str:
    if left == right:
        return "exact_same_source_token"
    left_punctuation_free = left.replace("-", "").replace("_", "")
    right_punctuation_free = right.replace("-", "").replace("_", "")
    if left_punctuation_free == right_punctuation_free:
        return "encoding_semantics_punctuation_only"
    return "semantic_mismatch_review_required"


def read_prior_matches(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "gse106737_sample",
            "gse83452_best_sample",
            "pearson_r",
            "reciprocal_best",
            "duplicate_threshold_pass",
            "n_fixed_probes",
            "probe_selection",
        }
        if not required.issubset(reader.fieldnames or ()):
            raise CrossAccessionTimepointError("prior fingerprint fields differ")
        rows = [row for row in reader if row["duplicate_threshold_pass"] == "TRUE"]
    if len(rows) != 78 or any(row["reciprocal_best"] != "TRUE" for row in rows):
        raise CrossAccessionTimepointError("prior reciprocal fingerprint count differs")
    return rows


def read_gse106737(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"dataset_id", "sample_id", "biological_unit_id", "timepoint"}
        if not required.issubset(reader.fieldnames or ()):
            raise CrossAccessionTimepointError("GSE106737 manifest fields differ")
        rows = [row for row in reader if row["dataset_id"] == "GSE106737"]
    result = {row["sample_id"]: row for row in rows}
    if len(rows) != 111 or len(result) != 111:
        raise CrossAccessionTimepointError("GSE106737 sample census differs")
    return result


def read_gse83452(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_accession", "participant_id", "timepoint"}
        if not required.issubset(reader.fieldnames or ()):
            raise CrossAccessionTimepointError("GSE83452 manifest fields differ")
        rows = [dict(row) for row in reader]
    result = {row["sample_accession"]: row for row in rows}
    if len(rows) != 231 or len(result) != 231:
        raise CrossAccessionTimepointError("GSE83452 sample census differs")
    return result


def diagnose(
    fingerprint: Path, gse106737_manifest: Path, gse83452_records: Path
) -> tuple[list[dict[str, object]], dict[str, object]]:
    matches = read_prior_matches(fingerprint)
    left = read_gse106737(gse106737_manifest)
    right = read_gse83452(gse83452_records)
    rows: list[dict[str, object]] = []
    for match in matches:
        left_sample = match["gse106737_sample"]
        right_sample = match["gse83452_best_sample"]
        if left_sample not in left or right_sample not in right:
            raise CrossAccessionTimepointError("fingerprint sample alias is absent")
        left_record = left[left_sample]
        right_record = right[right_sample]
        disposition = token_disposition(
            left_record["timepoint"], right_record["timepoint"]
        )
        rows.append(
            {
                "cohort_family_id": "antwerp_inserm_shared",
                "gse106737_sample_accession": left_sample,
                "gse83452_sample_accession": right_sample,
                "gse106737_participant_alias": (
                    f"gse106737::{left_record['biological_unit_id']}"
                ),
                "gse83452_participant_alias": right_record["participant_id"],
                "gse106737_source_timepoint": left_record["timepoint"],
                "gse83452_source_timepoint": right_record["timepoint"],
                "timepoint_token_disposition": disposition,
                "pearson_r": match["pearson_r"],
                "reciprocal_best": True,
                "duplicate_threshold_pass": True,
                "n_fixed_probes": int(match["n_fixed_probes"]),
                "probe_selection": match["probe_selection"],
                "bad_match_indicated_by_timepoint": (
                    disposition == "semantic_mismatch_review_required"
                ),
            }
        )
    participant_pairs = {
        (row["gse106737_participant_alias"], row["gse83452_participant_alias"])
        for row in rows
    }
    states = Counter(str(row["timepoint_token_disposition"]) for row in rows)
    pairs = Counter(
        (
            str(row["gse106737_source_timepoint"]),
            str(row["gse83452_source_timepoint"]),
        )
        for row in rows
    )
    if len(rows) != 78 or len(participant_pairs) != 41:
        raise CrossAccessionTimepointError("diagnostic alias census differs")
    if states.get("semantic_mismatch_review_required", 0) != 0:
        status = "blocked_semantic_timepoint_mismatch"
        conclusion = "bad_match_or_visit_mismatch_requires_review"
    else:
        status = "pass_deposited_visit_tokens_are_exact_or_punctuation_only"
        conclusion = "encoding_semantics_not_bad_match"
    summary = {
        "schema_version": "masld-bench-cross-accession-timepoint-diagnostic-v1",
        "status": status,
        "cohort_family_id": "antwerp_inserm_shared",
        "array_aliases_audited": len(rows),
        "participant_aliases_audited": len(participant_pairs),
        "source_token_pairs": {
            f"{left_value}|{right_value}": count
            for (left_value, right_value), count in sorted(pairs.items())
        },
        "disposition_counts": dict(sorted(states.items())),
        "conclusion": conclusion,
        "source_timepoints_retained_verbatim": True,
        "timepoints_coerced_or_overwritten": False,
        "labels_read": False,
        "expression_values_read": False,
        "normalization_run": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
    }
    return rows, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fingerprint", type=Path, required=True)
    parser.add_argument("--fingerprint-sha256", required=True)
    parser.add_argument("--gse106737-manifest", type=Path, required=True)
    parser.add_argument("--gse106737-manifest-sha256", required=True)
    parser.add_argument("--gse83452-records", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if sha256_file(arguments.fingerprint) != arguments.fingerprint_sha256:
        raise CrossAccessionTimepointError("prior fingerprint SHA-256 differs")
    if sha256_file(arguments.gse106737_manifest) != arguments.gse106737_manifest_sha256:
        raise CrossAccessionTimepointError("GSE106737 manifest SHA-256 differs")
    arguments.output.mkdir(parents=True, exist_ok=False)
    rows, summary = diagnose(
        arguments.fingerprint,
        arguments.gse106737_manifest,
        arguments.gse83452_records,
    )
    with (arguments.output / "cross_accession_timepoint_mismatches.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    (arguments.output / "cross_accession_timepoint_diagnostic.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
