#!/usr/bin/env python3
"""Independently validate candidate-only ATAC v3 integration products."""

from __future__ import annotations

import argparse
import csv
import hashlib
from collections import Counter
from pathlib import Path


RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
ROOT = Path(__file__).resolve().parents[4]
EXPECTED = (
    ROOT / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID
).resolve()
V3_FIELDS = (
    "atac_v3_release_id",
    "atac_v3_promoter_da_states",
    "atac_v3_promoter_da_source_dependent",
    "atac_v3_program_uids",
    "atac_v3_program_promoter_coverage_states",
    "atac_v3_program_score_states",
    "atac_v3_program_contrast_states",
    "atac_v3_program_cross_cohort_states",
    "atac_v3_testability_reasons",
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, default=EXPECTED)
    parser.add_argument("--canonical-before-sha256", required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or ()), list(reader)


def require_unique(table: list[dict[str, str]], fields: tuple[str, ...], label: str) -> None:
    keys = [tuple(row[field] for field in fields) for row in table]
    if len(keys) != len(set(keys)):
        raise RuntimeError(f"duplicate {label} keys")


def main() -> None:
    args = arguments()
    supplied = args.candidate_root
    root = supplied.resolve()
    if root != EXPECTED or supplied.is_symlink():
        raise RuntimeError(f"unsafe candidate root: {root}")
    out = root / "integration/INTEGRATION_READY"
    if out.exists():
        raise RuntimeError(f"refusing to overwrite integration seal: {out}")

    manifest = rows(root / "integration/integration_manifest.tsv")
    for row in manifest:
        if row["release_id"] != RELEASE_ID:
            raise RuntimeError("integration manifest release mismatch")
        if row["role"] == "integration_artifact":
            artifact = root / "integration" / row["artifact"]
        elif row["role"] == "sealed_non_genetic_input_manifest":
            artifact = root / row["artifact"]
        elif row["role"] == "post_gate_integration_producer":
            artifact = ROOT / row["artifact"]
        else:
            raise RuntimeError(f"unknown integration manifest role: {row['role']}")
        if not artifact.is_file() or sha256(artifact) != row["sha256"]:
            raise RuntimeError(f"integration hash mismatch: {row['artifact']}")

    integration = root / "integration"
    original = (integration / "fig4f_source_matrix_original.tsv").read_bytes()
    extended = (integration / "fig4f_source_matrix_with_atac_v3.tsv").read_bytes()
    if extended[: len(original)] != original:
        raise RuntimeError("protected Figure 4F source prefix changed")
    original_rows = rows(integration / "fig4f_source_matrix_original.tsv")
    extended_rows = rows(integration / "fig4f_source_matrix_with_atac_v3.tsv")
    if len(extended_rows) - len(original_rows) != 234:
        raise RuntimeError("Figure 4F source does not append exactly 234 ATAC v3 rows")
    appended = extended_rows[len(original_rows):]
    require_unique(appended, ("dataset_id", "program_uid"), "Figure 4F ATAC")
    if Counter(row["dataset_id"] for row in appended) != {
        "GSE244832": 117, "GSE281367": 117,
    }:
        raise RuntimeError("Figure 4F ATAC cohort family size mismatch")

    states = rows(integration / "fig4f_atac_v3_state_contract.tsv")
    if len(states) != 234:
        raise RuntimeError("Figure 4F state contract must contain 234 rows")
    require_unique(states, ("cohort", "program_uid"), "Figure 4F state")
    if Counter(row["cohort"] for row in states) != {
        "GSE244832": 117, "GSE281367": 117,
    }:
        raise RuntimeError("Figure 4F state cohort family size mismatch")
    allowed = {
        "promoter_measurement_state": {"observed", "unobserved"},
        "program_score_state": {"testable", "untestable"},
        "contrast_state": {"testable", "untestable"},
    }
    for field, values in allowed.items():
        if not {row[field] for row in states}.issubset(values):
            raise RuntimeError(f"invalid {field}")

    promoter = rows(integration / "gene_catalog_promoter_da_adapter.tsv")
    programs = rows(integration / "gene_catalog_program_atac_adapter.tsv")
    require_unique(promoter, ("gene_symbol", "lineage", "promoter_da_state"), "promoter adapter")
    require_unique(programs, ("gene_symbol", "program_uid", "cohort"), "program adapter")
    if not {row["promoter_da_state"] for row in promoter}.issubset(
        {"supported", "discordant", "source_dependent", "indeterminate", "untestable"}
    ):
        raise RuntimeError("invalid promoter DA adapter state")

    before = root / "compatibility/default_before/l8_atac_columns.csv"
    after = root / "compatibility/default_after/l8_atac_columns.csv"
    if before.read_bytes() != after.read_bytes():
        raise RuntimeError("default integration output changed after optional adapter patch")
    base_fields, base_rows = csv_rows(after)
    v3_fields, v3_rows = csv_rows(root / "compatibility/with_v3/l8_atac_columns.csv")
    if any(field in base_fields for field in V3_FIELDS):
        raise RuntimeError("default integration unexpectedly contains ATAC v3 fields")
    if v3_fields != base_fields + list(V3_FIELDS):
        raise RuntimeError("candidate integration field order mismatch")
    if len(base_rows) != len(v3_rows) or len(v3_rows) != 27187:
        raise RuntimeError("candidate catalog row count mismatch")
    for base, v3 in zip(base_rows, v3_rows):
        if any(base[field] != v3[field] for field in base_fields):
            raise RuntimeError("candidate adapter changed a protected default L8 value")
        if v3["atac_v3_release_id"] not in {"", RELEASE_ID}:
            raise RuntimeError("candidate catalog release field mismatch")
    if not any(row["atac_v3_release_id"] == RELEASE_ID for row in v3_rows):
        raise RuntimeError("candidate catalog contains no ATAC v3 mappings")

    canonical = ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
    canonical_after = sha256(canonical)
    if canonical_after != args.canonical_before_sha256:
        raise RuntimeError("canonical Gene Catalog changed during candidate integration")

    producer = ROOT / "Analysis/ATAC/Integration/scripts/35_atac_integration.py"
    fields = ("release_id", "gate", "status", "artifact", "sha256")
    seal_paths = (
        integration / "integration_manifest.tsv",
        integration / "fig4f_atac_v3_state_contract.tsv",
        root / "compatibility/with_v3/l8_atac_columns.csv",
    )
    seal_rows = [
        {
            "release_id": RELEASE_ID,
            "gate": "INTEGRATION_READY",
            "status": "READY",
            "artifact": str(path.relative_to(root)),
            "sha256": sha256(path),
        }
        for path in seal_paths
    ]
    seal_rows.append({
        "release_id": RELEASE_ID,
        "gate": "INTEGRATION_READY",
        "status": "READY",
        "artifact": str(producer.relative_to(ROOT)),
        "sha256": sha256(producer),
    })
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(seal_rows)
    print(
        f"INTEGRATION_READY: {len(states)} program states, {len(promoter)} promoter rows, "
        f"{len(programs)} program-gene rows"
    )


if __name__ == "__main__":
    main()
