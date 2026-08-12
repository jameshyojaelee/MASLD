#!/usr/bin/env python3
"""Independent validation and readiness sealing for promoted-COLOC context."""

from __future__ import annotations

import csv
import hashlib
import math
import os
import subprocess
import tempfile
from collections import Counter
from pathlib import Path


RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
ROOT = Path(__file__).resolve().parents[4]
CANDIDATE = (ROOT / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID).resolve()


def read(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def classify(mapped: float, shared: float, left: float, right: float) -> str:
    if mapped < 0.95:
        return "untestable"
    if shared >= 0.5:
        return "replicated_accessible"
    if max(left, right) >= 0.5:
        return "source_dependent"
    if max(shared, left, right) > 0:
        return "partial"
    return "indeterminate"


def atomic(path: Path, columns: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    if path.exists():
        raise RuntimeError(f"refusing to overwrite genetics validator output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False, encoding="utf-8", newline="") as handle:
        temp = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp, path)


def main() -> None:
    if not (CANDIDATE / "NON_GENETIC_READY").is_file():
        raise RuntimeError("NON_GENETIC_READY is required before genetics sealing")
    if not (CANDIDATE / "genetics/PROMOTED_UPSTREAM_GATE.tsv").is_file():
        raise RuntimeError("promoted COLOC gate is absent")
    out = CANDIDATE / "genetics/validation"
    if out.exists() or (CANDIDATE / "GENETIC_READY").exists() or (CANDIDATE / "FULL_READY").exists():
        raise RuntimeError("refusing to overwrite genetics readiness outputs")
    all_rows = read(CANDIDATE / "genetics/context/genetic_lineage_context_all_pairs.tsv")
    primary = read(CANDIDATE / "genetics/context/genetic_lineage_context_primary_pairs.tsv")
    if not all_rows or not primary:
        raise RuntimeError("genetic lineage context is empty")
    keys = Counter((row["gwas_name"], row["ensembl"], row["signal_pair_index"]) for row in all_rows)
    if any(value != 5 for value in keys.values()):
        raise RuntimeError("each signal pair must have exactly five lineage rows")
    primary_keys = Counter((row["gwas_name"], row["ensembl"]) for row in primary)
    if any(value != 5 for value in primary_keys.values()):
        raise RuntimeError("each gene-study pair must have one primary pair across five lineages")
    plan = read(CANDIDATE / "genetics/replay_plan.tsv")
    expected_pairs = {(row["gwas_name"], row["ensembl"]): float(row["promoted_pp_h4_susie"]) for row in plan}
    if set(primary_keys) != set(expected_pairs):
        raise RuntimeError("replay did not return every promoted supported gene-study pair")
    replay_max: dict[tuple[str, str], float] = {}
    for row in all_rows:
        key = (row["gwas_name"], row["ensembl"])
        replay_max[key] = max(replay_max.get(key, float("-inf")), float(row["pp_h4"]))
    for key, expected in expected_pairs.items():
        if not math.isclose(replay_max[key], expected, rel_tol=1e-7, abs_tol=1e-9):
            raise RuntimeError(f"replayed maximum PP.H4 does not reproduce promoted aggregate: {key}")
    for row in all_rows:
        mapped = float(row["mapped_posterior_mass"])
        lost = float(row["lost_posterior_mass"])
        left = float(row["gse244832_accessible_mass"])
        right = float(row["gse281367_accessible_mass"])
        shared = float(row["shared_accessible_mass"])
        any_mass = float(row["any_accessible_mass"])
        pp4 = float(row["pp_h4"])
        if not math.isclose(mapped + lost, 1.0, abs_tol=1e-8):
            raise RuntimeError("mapped and lost posterior mass do not sum to one")
        if shared > min(left, right) + 1e-10 or any_mass + 1e-10 < max(left, right):
            raise RuntimeError("lineage posterior mass set relation failed")
        if not math.isclose(float(row["joint_any_accessible_mass"]), pp4 * any_mass, abs_tol=1e-10):
            raise RuntimeError("joint accessible mass mismatch")
        if row["evidence_state"] != classify(mapped, shared, left, right):
            raise RuntimeError("genetic evidence state mismatch")
    audit = read(CANDIDATE / "genetics/context/variant_liftover_audit.tsv")
    posterior_sums: Counter[tuple[str, str, str]] = Counter()
    for row in audit:
        posterior_sums[(row["gwas_name"], row["ensembl"], row["signal_pair_index"])] += float(row["snp_pp_h4"])
    if any(not math.isclose(value, 1.0, abs_tol=1e-6) for value in posterior_sums.values()):
        raise RuntimeError("pre-liftover signal-pair posterior does not sum to one")
    figure_manifest = read(CANDIDATE / "genetics/figures/figure_manifest.tsv")
    if len(figure_manifest) != 1:
        raise RuntimeError("genetic Figure 4D candidate manifest must contain one panel")
    figure = figure_manifest[0]
    pdf = CANDIDATE / figure["pdf"]
    source = CANDIDATE / figure["source_table"]
    if sha256(pdf) != figure["pdf_sha256"] or sha256(source) != figure["source_sha256"]:
        raise RuntimeError("genetic Figure 4D checksum mismatch")
    pdfinfo = subprocess.run(["pdfinfo", str(pdf)], text=True, capture_output=True, check=True).stdout
    if "Pages:           1" not in pdfinfo:
        raise RuntimeError("genetic Figure 4D candidate is not a one-page PDF")
    validation = [{
        "release_id": RELEASE_ID,
        "check": "variant_consistent_lineage_context",
        "status": "PASS",
        "n_signal_pairs": len(keys),
        "n_primary_gene_study_pairs": len(primary_keys),
        "n_variant_rows": len(audit),
    }]
    atomic(
        out / "validation_results.tsv",
        ("release_id", "check", "status", "n_signal_pairs", "n_primary_gene_study_pairs", "n_variant_rows"),
        validation,
    )
    artifacts = (
        CANDIDATE / "genetics/context/variant_liftover_audit.tsv",
        CANDIDATE / "genetics/context/genetic_lineage_context_all_pairs.tsv",
        CANDIDATE / "genetics/context/genetic_lineage_context_primary_pairs.tsv",
        CANDIDATE / "genetics/figures/figure_manifest.tsv",
        out / "validation_results.tsv",
    )
    seal = [{
        "release_id": RELEASE_ID, "gate": "GENETIC_READY", "status": "READY",
        "artifact": path.relative_to(CANDIDATE).as_posix(), "sha256": sha256(path),
    } for path in artifacts]
    atomic(CANDIDATE / "GENETIC_READY", ("release_id", "gate", "status", "artifact", "sha256"), seal)
    full = [{
        "release_id": RELEASE_ID,
        "gate": "FULL_READY",
        "status": "READY",
        "non_genetic_ready_sha256": sha256(CANDIDATE / "NON_GENETIC_READY"),
        "genetic_ready_sha256": sha256(CANDIDATE / "GENETIC_READY"),
    }]
    atomic(
        CANDIDATE / "FULL_READY",
        ("release_id", "gate", "status", "non_genetic_ready_sha256", "genetic_ready_sha256"),
        full,
    )
    print("GENETIC_READY\nFULL_READY")


if __name__ == "__main__":
    main()
