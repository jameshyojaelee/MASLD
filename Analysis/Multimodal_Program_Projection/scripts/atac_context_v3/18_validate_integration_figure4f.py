#!/usr/bin/env python3
"""Independently validate the post-gate candidate Figure 4F ATAC panel."""

from __future__ import annotations

import csv
import hashlib
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
ROOT = Path(__file__).resolve().parents[4]
CANDIDATE = (
    ROOT / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID
).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    integration_ready = read(CANDIDATE / "integration/INTEGRATION_READY")
    if not integration_ready or any(row["status"] != "READY" for row in integration_ready):
        raise RuntimeError("INTEGRATION_READY is required")
    out = CANDIDATE / "integration/FIG4F_READY"
    if out.exists():
        raise RuntimeError(f"refusing to overwrite Figure 4F seal: {out}")

    states = read(CANDIDATE / "integration/fig4f_atac_v3_state_contract.tsv")
    source = read(CANDIDATE / "integration/figures/fig4f_atac_v3_contract_source.tsv")
    if len(states) != 234 or len(source) != 30:
        raise RuntimeError("Figure 4F source family size mismatch")
    if Counter(row["cohort"] for row in states) != {"GSE244832": 117, "GSE281367": 117}:
        raise RuntimeError("Figure 4F program cohort family mismatch")

    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in states:
        grouped[(row["cohort"], row["lineage"])].append(row)
    expected = {}
    for key, group in grouped.items():
        measured = sum(int(row["n_promoter_measured_genes"]) for row in group)
        members = sum(int(row["n_program_genes"]) for row in group)
        expected[key + ("Promoter-gene coverage",)] = (measured, members, measured / members)
        testable = sum(row["program_score_state"] == "testable" for row in group)
        expected[key + ("Program score testable",)] = (testable, len(group), testable / len(group))
        contrast = sum(row["contrast_state"] == "testable" for row in group)
        expected[key + ("MASH contrast testable",)] = (contrast, len(group), contrast / len(group))
    if len(expected) != 30:
        raise RuntimeError("Figure 4F expected metric grid is incomplete")
    for row in source:
        key = (row["cohort"], row["lineage"], row["metric"])
        if key not in expected:
            raise RuntimeError(f"unexpected Figure 4F metric key: {key}")
        numerator, denominator, proportion = expected.pop(key)
        if int(float(row["numerator"])) != numerator or int(float(row["denominator"])) != denominator:
            raise RuntimeError(f"Figure 4F count mismatch: {key}")
        if abs(float(row["proportion"]) - proportion) > 1e-12:
            raise RuntimeError(f"Figure 4F proportion mismatch: {key}")
    if expected:
        raise RuntimeError("Figure 4F source is missing metric rows")

    figure_dir = CANDIDATE / "integration/figures"
    pdf = figure_dir / "fig4f_atac_v3_contract.pdf"
    info = subprocess.run(
        ["pdfinfo", str(pdf)], check=True, capture_output=True, text=True
    ).stdout
    if "Pages:           1" not in info:
        raise RuntimeError("Figure 4F PDF must contain one page")
    if shutil.which("pdffonts"):
        fonts = subprocess.run(
            ["pdffonts", str(pdf)], check=True, capture_output=True, text=True
        ).stdout.lower()
        if "type 3" in fonts or "type3" in fonts:
            raise RuntimeError("Figure 4F PDF contains a Type 3 font")
    elif b"/Subtype /Type3" in pdf.read_bytes():
        raise RuntimeError("Figure 4F PDF contains a Type 3 font")

    manifest = read(figure_dir / "figure_manifest.tsv")
    if len(manifest) != 1:
        raise RuntimeError("Figure 4F manifest must contain one row")
    row = manifest[0]
    source_path = figure_dir / row["source_table"]
    if row["pdf_sha256"] != sha256(pdf) or row["source_sha256"] != sha256(source_path):
        raise RuntimeError("Figure 4F manifest hash mismatch")

    artifacts = (
        figure_dir / "figure_manifest.tsv",
        source_path,
        pdf,
        figure_dir / "sessionInfo.txt",
    )
    with out.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("release_id", "gate", "status", "artifact", "sha256"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for artifact in artifacts:
            writer.writerow({
                "release_id": RELEASE_ID,
                "gate": "FIG4F_READY",
                "status": "READY",
                "artifact": str(artifact.relative_to(CANDIDATE)),
                "sha256": sha256(artifact),
            })
    print("FIG4F_READY: 234 programs, 30 explicit observability/testability cells")


if __name__ == "__main__":
    main()
