#!/usr/bin/env python3
"""Build the shared gene axis and label-blind target roster for the NAS transfer.

Both sides of this lane are RNA-seq with a stable-gene-id feature axis, so the
shared axis is a plain intersection of two frozen feature axes.  Neither the
source outcome nor the target outcome is opened here.

The target roster is projected from the GSE267145 molecular participant axis,
which carries identity and assay-observation state only.  The projection asserts
that no histology column is present rather than trusting that it is not, because
the NASH-CRN component sum is the evaluator's alone.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Sequence

# Any of these reaching the model side means the evaluator-only outcome leaked.
FORBIDDEN_TARGET_COLUMNS = (
    "nash_crn_component_sum",
    "steatosis",
    "ballooning",
    "lobular_inflammation",
    "lobular_necrosis",
    "fibrosis",
    "stage3",
    "stage5",
    "nas_score",
    "recorded_sex",
)
TARGET_IDENTITY_ALLOWLIST = (
    "participant_index",
    "participant_id",
    "rna_source_sample_accession",
    "h3k27ac_source_sample_accession",
    "pairing",
    "rna_observation_state",
    "h3k27ac_observation_state",
)


class NasActivationError(RuntimeError):
    """Raised when the activation would leak an outcome or lose a participant."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise NasActivationError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def read_feature_axis(path: Path) -> list[str]:
    fields, rows = read_tsv(path)
    if "stable_gene_id" not in fields:
        raise NasActivationError(f"feature axis lacks stable_gene_id: {path}")
    ids = [row["stable_gene_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise NasActivationError(f"feature axis repeats a gene: {path}")
    if any("." in value for value in ids):
        raise NasActivationError(f"feature axis carries a version suffix: {path}")
    return ids


def project_target_roster(path: Path) -> list[dict[str, str]]:
    fields, rows = read_tsv(path)
    leaked = sorted(set(FORBIDDEN_TARGET_COLUMNS) & set(fields))
    if leaked:
        raise NasActivationError(f"target participant axis carries outcomes: {leaked}")
    missing = sorted({"participant_id"} - set(fields))
    if missing:
        raise NasActivationError(f"target participant axis lacks {missing}")
    projected = [
        {key: row[key] for key in TARGET_IDENTITY_ALLOWLIST if key in row}
        for row in rows
    ]
    ids = [row["participant_id"] for row in projected]
    if len(ids) != len(set(ids)):
        raise NasActivationError("a target participant identifier repeats")
    return sorted(projected, key=lambda row: row["participant_id"])


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[dict]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-feature-axis", required=True, type=Path)
    parser.add_argument("--target-feature-axis", required=True, type=Path)
    parser.add_argument("--target-participant-axis", required=True, type=Path)
    parser.add_argument("--expected-source-genes", required=True, type=int)
    parser.add_argument("--expected-target-genes", required=True, type=int)
    parser.add_argument("--expected-target-participants", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()

    source_ids = read_feature_axis(arguments.source_feature_axis)
    target_ids = read_feature_axis(arguments.target_feature_axis)
    if (
        len(source_ids) != arguments.expected_source_genes
        or len(target_ids) != arguments.expected_target_genes
    ):
        raise NasActivationError("declared gene census differs from disk")
    shared = sorted(set(source_ids) & set(target_ids))
    if not shared:
        raise NasActivationError("source and target gene axes do not intersect")
    roster = project_target_roster(arguments.target_participant_axis)
    if len(roster) != arguments.expected_target_participants:
        raise NasActivationError("declared target participant census differs from disk")

    arguments.output.mkdir(parents=True)
    write_tsv(
        arguments.output / "common_stable_gene_axis.tsv",
        ("feature_index", "stable_gene_id"),
        [
            {"feature_index": index, "stable_gene_id": value}
            for index, value in enumerate(shared)
        ],
    )
    write_tsv(
        arguments.output / "target_participant_roster.tsv",
        tuple(roster[0]),
        roster,
    )
    receipt = {
        "schema_version": "masld-bench-gse267145-nas-transfer-activation-v1",
        "status": "pass_label_blind_shared_axis_and_target_roster",
        "source_series": "GSE135251",
        "target_series": "GSE267145",
        "source_genes": len(source_ids),
        "target_genes": len(target_ids),
        "shared_genes": len(shared),
        "target_participants": len(roster),
        "axis_order": "ascending_stable_gene_id",
        "target_identity_columns_projected": [
            key for key in TARGET_IDENTITY_ALLOWLIST if key in roster[0]
        ],
        "target_outcome_columns_read": False,
        "source_outcome_columns_read": False,
        "expression_values_read": False,
        "shared_axis_sha256": hashlib.sha256(
            "\n".join(shared).encode("utf-8")
        ).hexdigest(),
        "source_feature_axis_sha256": sha256_file(arguments.source_feature_axis),
        "target_feature_axis_sha256": sha256_file(arguments.target_feature_axis),
    }
    with (arguments.output / "activation_receipt.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
