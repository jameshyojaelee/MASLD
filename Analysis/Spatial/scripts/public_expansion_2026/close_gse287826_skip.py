#!/usr/bin/env python3
"""Close the GSE287826 branch without reading expression outcomes."""

from __future__ import annotations

import argparse
import csv
import hashlib
from datetime import datetime, timezone
from pathlib import Path


EMPTY_SCHEMAS = {
    "sample_manifest.tsv": [
        "dataset_id", "sample_id", "biological_unit_id", "technical_unit_type",
        "condition", "included", "exclusion_reason",
    ],
    "gene_mapping_audit.tsv": [
        "dataset_id", "source_feature", "gene_symbol", "mapping_status",
    ],
    "design_audit.tsv": [
        "dataset_id", "contrast", "n_biological_units", "design_rank",
        "estimable", "failure_reason",
    ],
    "program_testability.tsv": [
        "dataset_id", "program_uid", "n_mapped_genes", "retained_l1_weight",
        "testability_status", "reason",
    ],
    "per_sample_program_scores.tsv": [
        "dataset_id", "biological_unit_id", "program_uid", "score",
        "score_transform",
    ],
    "program_effects.tsv": [
        "dataset_id", "program_uid", "contrast", "effect", "se", "ci_low",
        "ci_high", "p", "q", "status",
    ],
    "sensitivity.tsv": [
        "dataset_id", "program_uid", "sensitivity_id", "effect", "direction_agree",
    ],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_gate(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one source-gate row, found {len(rows)}")
    return rows[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-gate", required=True, type=Path)
    parser.add_argument("--download-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    gate = read_gate(args.source_gate)
    if gate.get("status") != "skipped_no_donor_key":
        raise ValueError(
            "Refusing to close as skipped: source gate is " + gate.get("status", "missing")
        )
    if gate.get("inference_authorized", "").lower() == "true":
        raise ValueError("Refusing skip closure: source gate authorizes inference")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filename, columns in EMPTY_SCHEMAS.items():
        with (args.output_dir / filename).open("w", newline="") as handle:
            csv.writer(handle, delimiter="\t", lineterminator="\n").writerow(columns)

    closed_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    gate_columns = [
        "dataset_id", "status", "source_gate_sha256", "inference_authorized",
        "outcomes_read", "reason", "closed_utc",
    ]
    gate_row = [
        "GSE287826", "skipped_no_donor_key", sha256(args.source_gate), "false",
        "false", gate.get("reason", "no_explicit_public_donor_key"), closed_at,
    ]
    with (args.output_dir / "gate_status.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(gate_columns)
        writer.writerow(gate_row)

    inputs = [args.source_gate, args.download_manifest, Path(__file__).resolve()]
    with (args.output_dir / "execution_manifest.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["role", "path", "bytes", "sha256", "recorded_utc"])
        for role, path in zip(("source_gate", "download_manifest", "producer"), inputs):
            writer.writerow([role, str(path), path.stat().st_size, sha256(path), closed_at])


if __name__ == "__main__":
    main()
