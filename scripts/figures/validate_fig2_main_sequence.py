#!/usr/bin/env python3
"""Validate the indexed working-main Figure 2 sequence or an A-I candidate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PANELS = ROOT / "figures/main/fig2_genetics/panels"
SIZE_SPEC = ROOT / "figures/layout_specs/figure2_panel_sizes.tsv"
LEGEND = ROOT / "docs/manuscript/05_figure_legends.md"
REINDEX_MANIFEST = PANELS / "Figure2_reindex_manifest.tsv"
PROHIBITED = ("ABF-fallback", "exploratory", "causal gene", "ancestry-specific")
TOLERANCE_IN = 0.015


def pdf_info(path: Path) -> tuple[int, float, float]:
    text = subprocess.check_output(["pdfinfo", str(path)], text=True)
    pages = int(re.search(r"^Pages:\s+(\d+)", text, re.M).group(1))
    size = re.search(r"^Page size:\s+([\d.]+) x ([\d.]+) pts", text, re.M)
    return pages, float(size.group(1)) / 72, float(size.group(2)) / 72


def pdf_text(path: Path, errors: list[str]) -> str | None:
    """Extract PDF text, recording a FAIL instead of aborting the whole run.

    A partial candidate render (source sidecar present, PDF absent) otherwise
    kills validation with a CalledProcessError traceback before the collected
    failures are printed.
    """
    if not path.exists():
        errors.append(f"missing panel PDF for text check: {path.name}")
        return None
    try:
        return subprocess.check_output(["pdftotext", str(path), "-"], text=True)
    except subprocess.CalledProcessError as exc:
        errors.append(f"pdftotext failed on {path.name}: {exc}")
        return None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-dir", type=Path, default=PANELS)
    parser.add_argument(
        "--layout-only",
        action="store_true",
        help="validate the A-I assembled layout candidate; omit the auxiliary J panel and promotion manifest",
    )
    parser.add_argument(
        "--expected-json",
        type=Path,
        help="r2_summary.json from GWAS/finemapping/src/perf/11_assemble_coloc_r2.py; the "
        "multi-signal, single-signal-only and union counts come from its fig2 block. "
        "Without it the adopted 2026-08-17 release counts (462 / 551 / 1,013) are checked.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    panel_dir = args.panel_dir.resolve()
    errors: list[str] = []
    # Gene counts of the COLOC release the panels were rendered from. The class,
    # ancestry and roster breakdowns below are literal only for the adopted release;
    # for another release they must sum to these counts.
    release_counts = args.expected_json is not None
    multi, single, union = 462, 551, 1013
    if release_counts:
        fig2 = json.loads(args.expected_json.read_text())["fig2"]
        multi, single = fig2["multi_signal"], fig2["single_signal_only"]
        union = fig2["multi_or_single_signal_coloc_gene_union"]
        if multi + single != union:
            errors.append(f"expected-json counts do not add up: {multi} + {single} != {union}")
    with SIZE_SPEC.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected_all = {
        Path(row["pdf"]).name: (float(row["width_in"]), float(row["height_in"]))
        for row in rows
    }
    # Select the A-I sequence by identity, not by position: slicing [:9] silently
    # validates a different nine panels if the size index is reordered.
    expected = (
        {name: size for name, size in expected_all.items() if not name.startswith("Fig2J_")}
        if args.layout_only
        else expected_all
    )
    indexed = list(expected)
    expected_n = 9 if args.layout_only else 10
    if len(indexed) != expected_n:
        errors.append(f"size index selected {len(indexed)} panels, expected {expected_n}")

    for name, (want_w, want_h) in expected.items():
        path = panel_dir / name
        if not path.exists():
            errors.append(f"missing indexed panel {name}")
            continue
        pages, width, height = pdf_info(path)
        if pages != 1:
            errors.append(f"{name}: {pages} pages")
        if abs(width - want_w) > TOLERANCE_IN or abs(height - want_h) > TOLERANCE_IN:
            errors.append(
                f"{name}: {width:.3f}x{height:.3f} in, expected {want_w:.2f}x{want_h:.2f}"
            )
        artwork = subprocess.check_output(["pdftotext", str(path), "-"], text=True)
        for term in PROHIBITED:
            if term.lower() in artwork.lower():
                errors.append(f"{name}: prohibited artwork label {term!r}")
        if b"/FontFile" not in path.read_bytes():
            errors.append(f"{name}: no embedded font object detected")
        if b"Helvetica" not in path.read_bytes():
            errors.append(f"{name}: embedded Helvetica name not detected")
        if subprocess.run(
            ["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", str(path)],
            check=False,
        ).returncode:
            errors.append(f"{name}: Ghostscript parse failed")

    required_sources = [
        "fig2A_gwas_cascade_source.tsv",
        "Fig2B_PIP_vs_SuSiE-coloc_source.csv",
        "Fig2C_coding_noncoding_source.csv",
        "Fig2D_pip_architecture_by_trait_directness_source.tsv",
        "Fig2E_ancestry_unique_coloc_GWS_source.tsv",
        "Fig2F_crossancestry_coloc_source.tsv",
        "Fig2G_RORA_NAFLD_source.csv",
        "Fig2H_FABP1_ALT_source.csv",
    ]
    if not args.layout_only:
        required_sources.append("Fig2J_phenotype_provenance_source.csv")
    for name in required_sources:
        if not (panel_dir / name).exists():
            errors.append(f"missing source sidecar {name}")

    fig2a_source = panel_dir / "fig2A_gwas_cascade_source.tsv"
    if fig2a_source.exists():
        with fig2a_source.open(newline="") as handle:
            fig2a_rows = list(csv.DictReader(handle, delimiter="\t"))
        study_rows = [row for row in fig2a_rows if row["record_type"] == "study_registry"]
        if len(study_rows) != 35 or len({row["study"] for row in study_rows}) != 35:
            errors.append("Fig2A source does not contain 35 unique GWAS strata")
        summaries = {
            row["metric"]: int(row["value"])
            for row in fig2a_rows if row["record_type"] == "global_summary"
        }
        expected_summaries = {
            "gwas_strata": 35,
            "globally_unique_fine_mapped_loci": 265,
            "multi_or_single_signal_coloc_gene_union": union,
        }
        if summaries != expected_summaries:
            errors.append(f"Fig2A global summaries disagree: {summaries}")
        band_totals: dict[str, int] = {}
        for model in ("multi_signal", "single_signal_only"):
            band_totals[model] = sum(
                int(row["value"]) for row in fig2a_rows
                if row["record_type"] == "posterior_band" and row["signal_model"] == model
            )
        if band_totals != {"multi_signal": multi, "single_signal_only": single}:
            errors.append(f"Fig2A signal-model totals disagree: {band_totals}")
        routed_totals: dict[str, int] = {}
        for model in ("multi_signal", "single_signal_only"):
            routed_totals[model] = sum(
                int(row["value"]) for row in fig2a_rows
                if row["record_type"] == "ancestry_trait_posterior_band"
                and row["signal_model"] == model
            )
        if routed_totals != band_totals:
            errors.append(
                f"Fig2A routed PP.H4 bands {routed_totals} disagree with endpoints {band_totals}"
            )
        fig2a_text = pdf_text(panel_dir / "fig2A_gwas_cascade.pdf", errors)
        if fig2a_text is not None and f"{union:,}" not in fig2a_text:
            errors.append(f"Fig2A artwork does not display the {union:,}-gene union")

    fig2c_source = panel_dir / "Fig2C_coding_noncoding_source.csv"
    if fig2c_source.exists():
        with fig2c_source.open(newline="") as handle:
            fig2c_rows = list(csv.DictReader(handle))
        class_counts = {row["consequence"]: int(row["n_genes"]) for row in fig2c_rows}
        expected_classes = {
            "coding": 36,
            "UTR": 73,
            "promoter": 107,
            "intron": 455,
            "intergenic": 342,
        }
        if release_counts:
            if set(class_counts) != set(expected_classes) or sum(class_counts.values()) != union:
                errors.append(f"Fig2C classes do not partition the {union:,}-gene union: {class_counts}")
        elif class_counts != expected_classes:
            errors.append(f"Fig2C promoted-union classes disagree: {class_counts}")
        fig2c_text = pdf_text(panel_dir / "Fig2C_coding_noncoding_split.pdf", errors)
        if fig2c_text is not None:
            if str(union) not in fig2c_text.replace(",", ""):
                errors.append(f"Fig2C artwork does not display the {union:,}-gene union")
            if union != 1030 and "1030" in fig2c_text.replace(",", ""):
                errors.append("Fig2C artwork still displays the superseded 1,030-gene union")

    fig2e_source = panel_dir / "Fig2E_ancestry_unique_coloc_GWS_source.tsv"
    if fig2e_source.exists():
        with fig2e_source.open(newline="") as handle:
            fig2e_rows = list(csv.DictReader(handle, delimiter="\t"))
        ancestry_counts = {row["category"]: int(row["n_genes"]) for row in fig2e_rows}
        expected_ancestry = {
            "EUR_only": 352,
            "shared": 70,
            "non_EUR_sub_threshold": 26,
            "non_EUR_GWS_unique": 14,
        }
        if release_counts:
            if set(ancestry_counts) != set(expected_ancestry) or sum(ancestry_counts.values()) != multi:
                errors.append(f"Fig2E ancestry bins do not partition the {multi} multi-signal genes: {ancestry_counts}")
        elif ancestry_counts != expected_ancestry:
            errors.append(f"Fig2E promoted ancestry partition disagrees: {ancestry_counts}")

    fig2j_source = panel_dir / "Fig2J_phenotype_provenance_source.csv"
    if fig2j_source.exists():
        with fig2j_source.open(newline="") as handle:
            fig2j_rows = list(csv.DictReader(handle))
        phenotype_totals: dict[str, int] = {}
        for row in fig2j_rows:
            key = row["phenotype_stratum"]
            phenotype_totals[key] = phenotype_totals.get(key, 0) + int(row["n_studies"])
        expected_phenotypes = {
            "alt_ast_or_ggt": 20,
            "direct_masld_mash_diagnosis": 12,
            "mri_pdff_or_histologic_steatosis": 3,
        }
        if phenotype_totals != expected_phenotypes:
            errors.append(f"Fig2J phenotype totals disagree: {phenotype_totals}")

    fig2f_source = panel_dir / "Fig2F_crossancestry_coloc_source.tsv"
    fig2f_pdf = panel_dir / "Fig2F_crossancestry_coloc.pdf"
    if fig2f_source.exists() and fig2f_pdf.exists():
        with fig2f_source.open(newline="") as handle:
            fig2f_rows = list(csv.DictReader(handle, delimiter="\t"))
        expected_genes = {
            "GGT1", "PANX1", "EPHA2", "C2orf16", "ACTG1", "SHROOM3",
            "MLIP", "MSL2", "EFHD1", "GOT2", "GLDC", "HSCB", "CHEK2",
        }
        observed_genes = {row["gene"] for row in fig2f_rows}
        if release_counts:
            # Another release may change the >=4-ancestry roster; keep one cell per
            # gene and ancestry (five ancestries) and report the roster.
            expected_genes = observed_genes if observed_genes else {"<empty roster>"}
        if observed_genes != expected_genes or len(fig2f_rows) != 5 * len(expected_genes):
            errors.append(
                "Fig2F does not contain the complete promoted >=4-ancestry roster: "
                f"{sorted(observed_genes)} ({len(fig2f_rows)} cells)"
            )
        support_by_gene = {
            gene: sum(
                row["evidence_state"] in {"multi_signal", "single_signal_only"}
                for row in fig2f_rows if row["gene"] == gene
            )
            for gene in observed_genes
        }
        if set(support_by_gene.values()) != {4}:
            errors.append(f"Fig2F ancestry-support counts disagree: {support_by_gene}")
        observed_states = {row["evidence_state"] for row in fig2f_rows}
        fig2f_text = pdf_text(fig2f_pdf, errors)
        state_labels = {
            "multi_signal": "multi-signal",
            "single_signal_only": "single-signal only",
            "evaluated_no_support": "evaluated ≤ 0.5",
            "not_evaluable": "Not evaluable",
        }
        if fig2f_text is not None:
            for state, label in state_labels.items():
                present = label.lower() in fig2f_text.lower()
                if state not in observed_states and present:
                    errors.append(f"Fig2F legend includes unused evidence state {state}")
                # The reverse case matters just as much: a state the panel draws
                # but the legend omits leaves the reader an unexplained mark.
                if state in observed_states and not present:
                    errors.append(f"Fig2F draws evidence state {state} with no legend entry")

    check_promotion = not args.layout_only and panel_dir == PANELS.resolve()
    reindex_manifest = panel_dir / REINDEX_MANIFEST.name
    if check_promotion and not reindex_manifest.exists():
        errors.append("missing Figure2_reindex_manifest.tsv")
    elif check_promotion:
        with reindex_manifest.open(newline="") as handle:
            manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
        if [row["panel"] for row in manifest_rows] != list("ABCDEFGHIJ"):
            errors.append("reindex manifest does not contain ordered panels A-J")
        manifest_names = [row["canonical_pdf"] for row in manifest_rows]
        if manifest_names != indexed:
            errors.append(
                f"manifest filenames disagree with size index: {manifest_names} versus {indexed}"
            )
        for row in manifest_rows:
            path = panel_dir / row["canonical_pdf"]
            if not path.exists():
                errors.append(f"manifest target is missing: {row['canonical_pdf']}")
            elif sha256(path) != row["sha256"]:
                errors.append(f"{row['canonical_pdf']}: hash disagrees with reindex manifest")

    stale = (
        "Fig2D_ancestry_unique_coloc_GWS.pdf",
        "Fig2E_crossancestry_coloc.pdf",
        "Fig2F_RORA_NAFLD.pdf",
        "Fig2G_FABP1_ALT.pdf",
        "Fig2H_locus_legend.pdf",
        "fig2I_eqtl_observability.pdf",
        "Fig2J_pip_architecture_by_trait_directness.pdf",
        "fig2K_phenotype_provenance.pdf",
        "Fig2J_eqtl_observability.pdf",
        "Fig2J_eqtl_observability_source.csv",
        "Fig2K_phenotype_provenance.pdf",
        "Fig2K_phenotype_provenance_source.csv",
        "fig2A_gwas_alluvial.pdf",
    )
    for name in stale:
        if check_promotion and (panel_dir / name).exists():
            errors.append(f"stale pre-reindex canonical file remains: {name}")

    legend_text = LEGEND.read_text()
    for term in PROHIBITED:
        if term.lower() in legend_text.lower():
            errors.append(f"legend contains prohibited label {term!r}")

    if errors:
        print("FAIL")
        print("\n".join(f"- {error}" for error in errors))
        return 1
    print("PASS")
    print(f"indexed_panels\t{len(indexed)}")
    print("one_page_exact_dimensions\tTRUE")
    print("embedded_fonts\tTRUE")
    print("prohibited_labels_absent\tTRUE")
    print("source_sidecars_present\tTRUE")
    print(f"panel_dir\t{panel_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
