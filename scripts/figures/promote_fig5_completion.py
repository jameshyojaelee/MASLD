#!/usr/bin/env python3
"""Promote the validated Figure 5 completion candidate with recoverable archival."""

from __future__ import annotations

import csv
import hashlib
import os
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CANDIDATE = ROOT / "figures/candidates/fig5-completion-2026-08-24-v6"
CANONICAL = ROOT / "figures/main/fig5_molecular_context"
ARCHIVE = CANONICAL / "_legacy/2026-08-24_pre_fig5_completion"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv_atomic(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    temporary = path.with_name(f".{path.name}.pending.{os.getpid()}")
    if temporary.exists():
        raise RuntimeError(f"pending manifest exists: {temporary}")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def copy_atomic(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.pending.{os.getpid()}")
    if temporary.exists():
        raise RuntimeError(f"pending file exists: {temporary}")
    shutil.copy2(source, temporary)
    os.replace(temporary, destination)


def main() -> None:
    if not (CANDIDATE / "VALIDATED").is_file():
        raise RuntimeError("Figure 5 candidate lacks VALIDATED gate")
    if ARCHIVE.exists():
        raise RuntimeError(f"refusing to overwrite promotion archive: {ARCHIVE}")
    candidate_manifest = {row["callout"]: row for row in read_tsv(CANDIDATE / "CANDIDATE_MAIN_PANELS.tsv")}
    canonical_manifest = read_tsv(CANONICAL / "CANONICAL_MAIN_PANELS.tsv")
    if set(candidate_manifest) != {f"5{letter}" for letter in "ABCDEF"}:
        raise RuntimeError("candidate manifest does not contain Figure 5A-F")
    for row in canonical_manifest:
        path = CANONICAL / row["filename"]
        if not path.is_file() or sha256(path) != row["sha256"] or path.stat().st_size != int(row["bytes"]):
            raise RuntimeError(f"canonical panel drift before promotion: {row['callout']}")
    for row in candidate_manifest.values():
        path = CANDIDATE / row["filename"]
        if not path.is_file() or sha256(path) != row["sha256"] or path.stat().st_size != int(row["bytes"]):
            raise RuntimeError(f"candidate panel drift before promotion: {row['callout']}")

    ARCHIVE.mkdir(parents=True)
    archive_files = [
        CANONICAL / "CANONICAL_MAIN_PANELS.tsv",
        CANONICAL / "CANDIDATE_PANEL_INDEX.tsv",
        CANONICAL / "panels/fig5a_input_firewall.pdf",
        CANONICAL / "panels/fig5c_mrna_protein_composite.pdf",
        CANONICAL / "panels/fig5d_snatac_accessibility.pdf",
        CANONICAL / "panels/data/fig5a_input_firewall.tsv",
        CANONICAL / "data/composite_mrna_protein_corrected.csv",
        CANONICAL / "data/composite_protein_histology_partial_corrected.csv",
        CANONICAL / "data/fig5d_snatac_accessibility_source.tsv",
    ]
    for source in archive_files:
        if source.is_file():
            relative = source.relative_to(CANONICAL)
            destination = ARCHIVE / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)

    promoted = {
        "5A": ("figure5/panels/fig5a_input_firewall.pdf", "panels/fig5a_input_firewall.pdf"),
        "5C": ("figure5/panels/fig5c_mrna_protein_composite.pdf", "panels/fig5c_mrna_protein_composite.pdf"),
        "5D": ("figure5/panels/fig5d_snatac_accessibility.pdf", "panels/fig5d_snatac_accessibility.pdf"),
    }
    for _, (source_relative, destination_relative) in promoted.items():
        copy_atomic(CANDIDATE / source_relative, CANONICAL / destination_relative)
    copy_atomic(
        CANDIDATE / "figure5/panels/data/fig5a_input_firewall.tsv",
        CANONICAL / "panels/data/fig5a_input_firewall.tsv",
    )
    copy_atomic(
        CANDIDATE / "figure5/data/composite_mrna_protein_corrected.csv",
        CANONICAL / "data/composite_mrna_protein_corrected.csv",
    )
    copy_atomic(
        CANDIDATE / "figure5/data/composite_protein_histology_partial_corrected.csv",
        CANONICAL / "data/composite_protein_histology_partial_corrected.csv",
    )
    copy_atomic(
        CANDIDATE / "figure5/panels/data/fig5d_variant_consistent_accessibility.tsv",
        CANONICAL / "panels/data/fig5d_variant_consistent_accessibility.tsv",
    )
    stale_source = CANONICAL / "data/fig5d_snatac_accessibility_source.tsv"
    if stale_source.exists():
        stale_source.unlink()
    copy_atomic(
        CANDIDATE / "figure5/fig5_molecular_context.pdf",
        CANONICAL / "fig5_molecular_context.pdf",
    )
    copy_atomic(CANDIDATE / "ASSEMBLED_FIGURE.tsv", CANONICAL / "ASSEMBLED_FIGURE.tsv")

    source_updates = {
        "5A": "scripts/figures/fig5a_input_firewall.R",
        "5C": "scripts/figures/composite_mrna_protein.R",
        "5D": "scripts/figures/fig5d_variant_consistent_accessibility.R",
    }
    input_updates = {
        "5D": "ATAC v3 FULL_READY; genetic_lineage_context_primary_pairs.tsv; promoted 2026-08-17 COLOC",
    }
    status_updates = {
        "5A": "accepted_current_source_wording",
        "5C": "accepted_fixed_selection_conditioned_cornerstone",
        "5D": "accepted_variant_consistent_promoted_coloc_context",
    }
    output_manifest = []
    for row in canonical_manifest:
        callout = row["callout"]
        if callout in promoted:
            path = CANONICAL / row["filename"]
            row["source"] = source_updates[callout]
            row["primary_inputs"] = input_updates.get(callout, row["primary_inputs"])
            row["sha256"] = sha256(path)
            row["bytes"] = str(path.stat().st_size)
            row["status"] = status_updates[callout]
        output_manifest.append(row)
    write_tsv_atomic(
        CANONICAL / "CANONICAL_MAIN_PANELS.tsv",
        output_manifest,
        ["callout", "filename", "source", "primary_inputs", "sha256", "bytes", "status"],
    )

    candidate_index = read_tsv(CANONICAL / "CANDIDATE_PANEL_INDEX.tsv")
    state_updates = {
        "5A": ("promoted_2026-08-24_current_wording", "accepted_current_source_wording"),
        "5C": ("promoted_2026-08-24_fixed_selection_conditioned", "accepted_cornerstone"),
        "5D": ("promoted_2026-08-24_variant_consistent", "accepted_promoted_coloc_variant_consistency"),
    }
    for row in candidate_index:
        if row["callout"] in state_updates:
            row["state"], row["gate"] = state_updates[row["callout"]]
            row["generator_or_source"] = source_updates[row["callout"]]
    write_tsv_atomic(
        CANONICAL / "CANDIDATE_PANEL_INDEX.tsv",
        candidate_index,
        ["callout", "candidate_relative_path", "generator_or_source", "state", "gate"],
    )

    assembled = CANONICAL / "fig5_molecular_context.pdf"
    promotion = [{
        "promotion_id": "fig5-completion-2026-08-24-v6",
        "candidate_root": CANDIDATE.relative_to(ROOT).as_posix(),
        "archive_root": ARCHIVE.relative_to(ROOT).as_posix(),
        "assembled_pdf": assembled.relative_to(ROOT).as_posix(),
        "assembled_sha256": sha256(assembled),
        "assembled_bytes": assembled.stat().st_size,
        "status": "promoted_pending_postcopy_validation",
    }]
    write_tsv_atomic(
        CANONICAL / "PROMOTION_2026-08-24.tsv",
        promotion,
        list(promotion[0]),
    )
    print("PROMOTED_PENDING_POSTCOPY_VALIDATION")


if __name__ == "__main__":
    main()
