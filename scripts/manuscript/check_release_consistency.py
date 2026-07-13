#!/usr/bin/env python3
"""Fail when the active manuscript diverges from the frozen release outputs."""

from __future__ import annotations

import json
import os
import re
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
RELEASE_ID = os.environ.get("MANUSCRIPT_RELEASE_ID", "2026-07-10-r1")
OUT = ROOT / "RNA-seq/results/manuscript_release" / RELEASE_ID
DOC_RELEASE = ROOT / "docs/manuscript/release"

ACTIVE_DOCS = [
    ROOT / "docs/paper_outline.md",
    ROOT / "docs/manuscript/working/abstract.md",
    ROOT / "docs/manuscript/working/01_intro.md",
    ROOT / "docs/manuscript/working/fig2.md",
    ROOT / "docs/manuscript/working/fig3.md",
    ROOT / "docs/manuscript/working/fig4.md",
    ROOT / "docs/manuscript/working/fig5_discussion.md",
    ROOT / "docs/manuscript/working/METHODS.md",
]


class Audit:
    def __init__(self) -> None:
        self.errors: list[str] = []
        self.notes: list[str] = []

    def require(self, condition: bool, message: str) -> None:
        if not condition:
            self.errors.append(message)

    def note(self, message: str) -> None:
        self.notes.append(message)


def one_row(df: pd.DataFrame, **filters: object) -> pd.Series:
    hit = df.copy()
    for key, value in filters.items():
        hit = hit[hit[key].eq(value)]
    if len(hit) != 1:
        raise AssertionError(f"Expected one row for {filters}, observed {len(hit)}")
    return hit.iloc[0]


def has_unsupported_assertion(text: str, phrase: str) -> list[str]:
    """Return phrase hits that are not explicitly negative, retired, or historical."""
    flagged = []
    negative_markers = (
        "retired", "historical", "do not", "does not", "did not", "not support",
        "cannot", "prohibited", "remove", "without implying", "rather than",
        "excluded", "failed", "no longer", "not a", "avoid",
    )
    for line in text.splitlines():
        if phrase.lower() not in line.lower():
            continue
        if not any(marker in line.lower() for marker in negative_markers):
            flagged.append(line.strip())
    return flagged


def main() -> None:
    audit = Audit()
    required = [
        OUT / "orthogonality_audit.tsv",
        OUT / "evidence_class_table.tsv",
        OUT / "evidence_class_validation_summary.tsv",
        OUT / "acceptance_gates.tsv",
        OUT / "strict_drug_verdict.json",
        OUT / "gsmap_trait_replication.tsv",
        OUT / "credible_set_coding_summary.tsv",
        DOC_RELEASE / "analysis_release_manifest.tsv",
        DOC_RELEASE / "claim_ledger.tsv",
    ]
    for path in required:
        audit.require(path.exists(), f"Missing required release artifact: {path}")
    if audit.errors:
        raise SystemExit("\n".join(audit.errors))

    orth = pd.read_csv(OUT / "orthogonality_audit.tsv", sep="\t")
    primary = one_row(
        orth, trait_scope="all", genetic_definition="susie", pp4_threshold=0.5
    )
    union = one_row(
        orth, trait_scope="all", genetic_definition="susie_or_abf", pp4_threshold=0.5
    )
    direct = one_row(
        orth, trait_scope="direct_disease",
        genetic_definition="susie_or_abf", pp4_threshold=0.5,
    )
    expected = {
        "primary": (473, 447, 34),
        "union": (1030, 971, 74),
        "direct": (172, 160, 11),
    }
    for name, row in [("primary", primary), ("union", union), ("direct", direct)]:
        observed = (
            int(row["n_genetic_total"]),
            int(row["n_genetic_joint"]),
            int(row["n_overlap"]),
        )
        audit.require(observed == expected[name], f"{name} counts drifted: {observed}")
        audit.require(
            float(row["pct_genetic_not_deg"]) >= 85,
            f"{name} no longer passes the non-overlap gate",
        )

    gates = pd.read_csv(OUT / "acceptance_gates.tsv", sep="\t")
    gate_map = dict(zip(gates["gate"], gates["passed"].astype(str).str.upper().eq("TRUE")))
    audit.require(gate_map.get("orthogonality_core") is True, "Orthogonality gate must pass")
    for gate in [
        "class_validation_article", "drug_calibration_article",
        "gsmap_direct_disease_article", "coding_architecture_coverage",
    ]:
        audit.require(gate_map.get(gate) is False, f"Frozen release expects failed gate: {gate}")

    with open(OUT / "strict_drug_verdict.json", encoding="utf-8") as handle:
        drug = json.load(handle)
    audit.require(drug["headline_pass"] is False, "Drug headline unexpectedly passed")
    audit.require(drug["n_approved_in_ranking_universe"] == 1, "Approved coverage drifted")
    audit.require(drug["n_matched_positive"] == 22, "Strict matched-positive count drifted")

    replication = pd.read_csv(OUT / "gsmap_trait_replication.tsv", sep="\t")
    audit.require(
        int(replication["replicated_both_cohorts"].astype(bool).sum()) == 0,
        "A gsMap trait now replicates; rerun claim review before changing the manuscript",
    )

    coding = pd.read_csv(OUT / "credible_set_coding_summary.tsv", sep="\t")
    coding_all = one_row(coding, trait_scope="all")
    audit.require(int(coding_all["n_coloc_genes"]) == 75, "Credible-set linked-gene count drifted")
    audit.require(
        float(coding_all["pct_primary_genes_linked"]) < 80,
        "Coding coverage now passes; rerun claim review before promoting the result",
    )

    combined = "\n".join(path.read_text(encoding="utf-8") for path in ACTIVE_DOCS)
    abstract = (ROOT / "docs/manuscript/working/abstract.md").read_text(encoding="utf-8")
    for token in ["473", "447", "34", "92.4%"]:
        audit.require(token in abstract, f"Abstract is missing frozen token {token}")
    audit.require("1,031" not in abstract, "Abstract contains retired blank-row union count 1,031")

    unsupported_phrases = [
        "recovers both approved",
        "full drug-development gradient",
        "tissue-real only",
        "convergence is required",
        "therapeutic windows",
        "multicellular cascade",
        "prospectively nominated",
    ]
    for phrase in unsupported_phrases:
        hits = has_unsupported_assertion(combined, phrase)
        audit.require(not hits, f"Unsupported positive assertion for '{phrase}': {hits[:3]}")

    figure_paths = [
        ROOT / "figures/main/fig1_atlas_overview/fig1b_complementary_maps.pdf",
        ROOT / "figures/main/fig1_atlas_overview/fig1c_genetic_trait_scope.pdf",
        ROOT / "figures/main/fig4_validation/fig4a_evidence_class_positive_rates.pdf",
        ROOT / "figures/main/fig4_validation/fig4b_evidence_class_adjusted_or.pdf",
    ]
    for path in figure_paths:
        audit.require(path.exists() and path.stat().st_size > 1000, f"Missing/empty figure: {path}")

    ledger = pd.read_csv(DOC_RELEASE / "claim_ledger.tsv", sep="\t")
    audit.require(ledger["analysis_release_id"].eq(RELEASE_ID).all(), "Claim ledger release ID drift")
    audit.note(
        f"Orthogonality: {int(primary.n_overlap)}/{int(primary.n_genetic_joint)} overlap; "
        f"{primary.pct_genetic_not_deg:.1f}% non-DE"
    )
    audit.note("Article escalation gates passed: 0/4")

    report = [
        f"# Manuscript Release Consistency Report: {RELEASE_ID}",
        "",
        f"**Status:** {'PASS' if not audit.errors else 'FAIL'}",
        "",
        "## Checks",
        "",
        *[f"- {note}" for note in audit.notes],
        f"- Active documents checked: {len(ACTIVE_DOCS)}",
        f"- Required release artifacts checked: {len(required)}",
        f"- Figure outputs checked: {len(figure_paths)}",
    ]
    if audit.errors:
        report.extend(["", "## Errors", "", *[f"- {e}" for e in audit.errors]])
    report_path = DOC_RELEASE / "consistency_report.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))
    if audit.errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
