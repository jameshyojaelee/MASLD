#!/usr/bin/env python3
"""Validate the timestamped complete Figure 5 candidate without promotion."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import re
import subprocess
from collections import Counter
from pathlib import Path

from PyPDF2 import PdfReader


ROOT = Path(__file__).resolve().parents[2]
CANONICAL = ROOT / "figures/main/fig5_molecular_context"
ATAC = ROOT / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def pdf_text(path: Path) -> str:
    return subprocess.run(
        ["pdftotext", str(path), "-"],
        text=True,
        capture_output=True,
        check=True,
    ).stdout


def check(condition: bool, name: str, detail: str, report: list[dict[str, str]]) -> None:
    report.append({"check": name, "status": "PASS" if condition else "FAIL", "detail": detail})
    if not condition:
        raise RuntimeError(f"{name}: {detail}")


def main() -> None:
    args = arguments()
    candidate = args.candidate_root.resolve()
    try:
        candidate.relative_to((ROOT / "figures/candidates").resolve())
    except ValueError as error:
        raise RuntimeError("candidate root must be under figures/candidates") from error
    figure = candidate / "figure5"
    panels = figure / "panels"
    report: list[dict[str, str]] = []

    expected = {
        "5A": "fig5a_input_firewall.pdf",
        "5B": "fig5b_protein_triage.pdf",
        "5C": "fig5c_mrna_protein_composite.pdf",
        "5D": "fig5d_snatac_accessibility.pdf",
        "5E": "fig5e_multimodal_program_summary.pdf",
        "5F": "fig5f_spatial_program_calibration.pdf",
    }
    check(
        all((panels / filename).is_file() for filename in expected.values()),
        "six_panel_roster",
        "all six Figure 5 panel PDFs exist",
        report,
    )
    for callout, filename in expected.items():
        path = panels / filename
        reader = PdfReader(str(path))
        check(
            len(reader.pages) == 1 and path.stat().st_size > 1000,
            f"{callout}_one_page_nonempty",
            f"{filename}: {path.stat().st_size} bytes",
            report,
        )
        check(
            b"/Subtype /Type3" not in path.read_bytes(),
            f"{callout}_no_type3",
            f"{filename} has no Type 3 font",
            report,
        )

    text_a = pdf_text(panels / expected["5A"])
    check(
        "prespecified display set" in text_a.lower() and "22 frozen programs" not in text_a.lower(),
        "5A_current_wording",
        "5A uses prespecified display-set wording",
        report,
    )
    firewall = read_tsv(panels / "data/fig5a_input_firewall.tsv")
    check(
        len(firewall) == 3
        and {row["branch"] for row in firewall}
        == {"prioritized_gene_context", "prespecified_program_projection", "fixed_protein_display"},
        "5A_input_firewall",
        "three fixed input branches are explicit",
        report,
    )

    canonical_manifest = {row["callout"]: row for row in read_tsv(CANONICAL / "CANONICAL_MAIN_PANELS.tsv")}
    # FIG5B_BULK_DEG: 5B is rendered from that bulk fit instead of copied, so its
    # printed numbers are rederived from the per-gene rows rather than hash-matched.
    rendered_5b = bool(os.environ.get("FIG5B_BULK_DEG"))
    if rendered_5b:
        genes = read_tsv(panels / "data/fig5b_protein_triage_genes.tsv")
        summary_5b = read_tsv(panels / "data/fig5b_protein_triage_summary.tsv")[0]
        n_conf = sum(row["class"] == "confirmed" for row in genes)
        n_disc = sum(row["class"] == "discordant" for row in genes)
        n_sig = n_conf + n_disc
        text_b = pdf_text(panels / expected["5B"])
        check(
            summary_5b["bulk_deg_path"] == str(Path(os.environ["FIG5B_BULK_DEG"]).resolve())
            and int(summary_5b["n_prioritized_measured"]) == len(genes)
            and int(summary_5b["n_protein_significant"]) == n_sig
            and int(summary_5b["n_significant_concordant"]) == n_conf
            and int(summary_5b["n_significant_discordant"]) == n_disc
            and f"(n = {n_sig:,})" in text_b
            and f"{n_conf} ({100 * n_conf / n_sig:.1f}%)" in text_b
            and f"{n_disc} ({100 * n_disc / n_sig:.1f}%)" in text_b
            and f"ρ = {float(summary_5b['rho_prioritized']):.2f}" in text_b,
            "5B_rendered_corrected_bulk",
            f"5B rendered from FIG5B_BULK_DEG: {n_sig} protein-significant, "
            f"{n_conf} concordant, {n_disc} discordant rederived and printed",
            report,
        )
    for callout in ("5E", "5F") if rendered_5b else ("5B", "5E", "5F"):
        observed = sha256(panels / expected[callout])
        check(
            observed == canonical_manifest[callout]["sha256"],
            f"{callout}_canonical_identity",
            f"accepted {callout} is byte-identical to canonical",
            report,
        )

    protein_rows = read_csv(figure / "data/composite_mrna_protein_corrected.csv")
    histology_rows = read_csv(figure / "data/composite_protein_histology_partial_corrected.csv")
    fixed_rows = read_tsv(
        ROOT / "Analysis/Multimodal_Program_Projection/results/proteomics/fixed_panel4c_rows.tsv"
    )
    fixed_genes = {row["gene"] for row in fixed_rows}
    check(
        len(protein_rows) == 25
        and len({row["gene"] for row in protein_rows}) == 25
        and {row["gene"] for row in protein_rows} == fixed_genes,
        "5C_fixed_rows",
        "5C contains the frozen 25-protein display without reselection",
        report,
    )
    check(
        len(histology_rows) == 25
        and {row["gene"] for row in histology_rows} == fixed_genes
        and {"Steatosis", "Ballooning", "Inflammation", "Fibrosis", "NAS"}.issubset(histology_rows[0]),
        "5C_histology_family",
        "5C contains all 25 proteins by five histology features",
        report,
    )
    text_c = pdf_text(panels / expected["5C"])
    check(
        all(term in text_c for term in ("Protein z-score", "Partial Spearman rho", "mRNA vs protein")),
        "5C_artwork_content",
        "5C retains abundance, histology, and RNA-protein views",
        report,
    )

    full_ready = read_tsv(ATAC / "FULL_READY")
    genetic_ready = ATAC / "GENETIC_READY"
    check(
        len(full_ready) == 1
        and full_ready[0]["status"] == "READY"
        and full_ready[0]["genetic_ready_sha256"] == sha256(genetic_ready),
        "5D_atac_full_ready",
        "ATAC v3 promoted-COLOC replay remains FULL_READY",
        report,
    )
    genetic_manifest = {row["artifact"]: row["sha256"] for row in read_tsv(genetic_ready)}
    context_relative = "genetics/context/genetic_lineage_context_primary_pairs.tsv"
    context_path = ATAC / context_relative
    check(
        genetic_manifest.get(context_relative) == sha256(context_path),
        "5D_sealed_input",
        "5D primary-pair input matches the sealed ATAC artifact",
        report,
    )

    # FIG5D_PRIMARY_PAIRS: the sealed context restricted to COLOC-eligible signal
    # pairs by atac_context_v3/24_restrict_primary_pairs_to_eligible.py.
    expected_rows, expected_pairs, expected_genes = 4080, 816, 462
    restricted = os.environ.get("FIG5D_PRIMARY_PAIRS")
    if restricted:
        context_path = Path(restricted)
        restriction = {row["role"]: row["sha256"]
                       for row in read_tsv(context_path.parent / "restriction_manifest.tsv")}
        check(
            restriction.get("restricted_primary_pairs") == sha256(context_path)
            and restriction.get("sealed_primary_pairs") == genetic_manifest.get(context_relative)
            and restriction.get("sealed_all_pairs")
            == genetic_manifest.get("genetics/context/genetic_lineage_context_all_pairs.tsv"),
            "5D_restricted_input",
            "5D primary pairs are the sealed ATAC context restricted to COLOC-eligible pairs",
            report,
        )
        census = {row["metric"]: row["r2"]
                  for row in read_tsv(context_path.parent / "restriction_summary.tsv")}
        expected_rows = int(census["lineage_rows"])
        expected_pairs = int(census["gene_study_pairs"])
        expected_genes = int(census["genes"])
    primary = read_tsv(context_path)
    source_d = read_tsv(panels / "data/fig5d_variant_consistent_accessibility.tsv")
    observed = Counter(
        (row["trait_class"], row["lineage"], row["evidence_state"]) for row in primary
    )
    denominators = Counter(
        (row["trait_class"], row["gwas_name"], row["ensembl"])
        for row in primary
        if row["lineage"] == "hepatocyte"
    )
    expected_denominators = Counter(key[0] for key in denominators)
    check(
        len(primary) == expected_rows
        and len({(row["gwas_name"], row["ensembl"]) for row in primary}) == expected_pairs
        and len({row["ensembl"] for row in primary}) == expected_genes,
        "5D_primary_cardinality",
        f"{expected_rows:,} lineage rows represent {expected_pairs:,} gene-study pairs "
        f"and {expected_genes:,} genes",
        report,
    )
    check(
        len(source_d) == 50,
        "5D_complete_state_grid",
        "2 trait classes x 5 lineages x 5 evidence states",
        report,
    )
    for row in source_d:
        key = (row["trait_class"], row["lineage"], row["evidence_state"])
        check(
            int(row["n_pairs"]) == observed[key]
            and int(row["denominator"]) == expected_denominators[row["trait_class"]]
            and math.isclose(float(row["proportion"]), int(row["n_pairs"]) / int(row["denominator"]), abs_tol=1e-12),
            f"5D_cell_{'_'.join(key)}",
            f"{row['n_pairs']}/{row['denominator']} rederived",
            report,
        )
    text_d = pdf_text(panels / expected["5D"])
    check(
        all(term in text_d for term in ("MASLD/liver fat", "Liver enzyme", "Both sources", "Partial mass")),
        "5D_artwork_scope",
        "5D separates trait classes and retains unresolved context",
        report,
    )

    assembled = figure / "fig5_molecular_context.pdf"
    assembled_reader = PdfReader(str(assembled))
    page = assembled_reader.pages[0]
    width = float(page.mediabox.width)
    height = float(page.mediabox.height)
    assembled_text = pdf_text(assembled)
    check(
        len(assembled_reader.pages) == 1
        and math.isclose(width, 7.2 * 72, abs_tol=1.0)
        and math.isclose(height, 9.2 * 72, abs_tol=1.0),
        "assembled_page_contract",
        f"one-page vector assembly, {width:.1f} x {height:.1f} pt",
        report,
    )
    check(
        all(re.search(rf"(^|\n){label}($|\n)", assembled_text) for label in "ABCDEF"),
        "assembled_callouts",
        "assembled page contains callouts A-F",
        report,
    )
    check(
        b"/Subtype /Type3" not in assembled.read_bytes(),
        "assembled_no_type3",
        "assembled PDF has no Type 3 font",
        report,
    )
    subprocess.run(
        ["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", str(assembled)],
        check=True,
    )
    check(True, "assembled_ghostscript", "Ghostscript parsed the complete PDF", report)

    manifest_rows = []
    sources = {
        "5A": "scripts/figures/fig5a_input_firewall.R",
        "5B": "scripts/figures/fig5b_protein_triage.R",
        "5C": "scripts/figures/composite_mrna_protein.R",
        "5D": "scripts/figures/fig5d_variant_consistent_accessibility.R",
        "5E": "scripts/figures/fig5e_multimodal_program_summary.R",
        "5F": "scripts/figures/fig5f_spatial_program_calibration.py",
    }
    statuses = {
        "5A": "validated_current_source_wording",
        "5B": "validated_rendered_from_corrected_bulk" if os.environ.get("FIG5B_BULK_DEG")
              else "accepted_canonical_identity",
        "5C": "validated_fixed_selection_conditioned_cornerstone",
        "5D": "validated_variant_consistent_promoted_coloc_context",
        "5E": "accepted_canonical_identity",
        "5F": "accepted_canonical_identity",
    }
    for callout, filename in expected.items():
        path = panels / filename
        manifest_rows.append({
            "callout": callout,
            "filename": f"figure5/panels/{filename}",
            "source": sources[callout],
            "sha256": sha256(path),
            "bytes": str(path.stat().st_size),
            "status": statuses[callout],
        })
    validation_dir = candidate / "validation"
    validation_dir.mkdir(parents=True, exist_ok=False)
    with (candidate / "CANDIDATE_MAIN_PANELS.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)
    with (candidate / "ASSEMBLED_FIGURE.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ("filename", "sha256", "bytes", "page_width_pt", "page_height_pt", "status")
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({
            "filename": "figure5/fig5_molecular_context.pdf",
            "sha256": sha256(assembled),
            "bytes": assembled.stat().st_size,
            "page_width_pt": f"{width:.3f}",
            "page_height_pt": f"{height:.3f}",
            "status": "validated_candidate_not_promoted",
        })
    with (validation_dir / "validation_report.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("check", "status", "detail"), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(report)
    print(f"PASS: {len(report)} Figure 5 completion checks")


if __name__ == "__main__":
    main()
