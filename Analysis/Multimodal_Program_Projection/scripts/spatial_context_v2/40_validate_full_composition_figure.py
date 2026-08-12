#!/usr/bin/env python3
"""Validate the candidate-only supplementary full-composition panel."""

from __future__ import annotations

import argparse
import importlib.util
import re
import subprocess
import sys
from pathlib import Path


def import_renderer(script_dir: Path):
    path = script_dir / "39_render_full_composition_figure.py"
    spec = importlib.util.spec_from_file_location("full_composition_figure_renderer", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import renderer: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Checks:
    def __init__(self) -> None:
        self.rows = []

    def require(self, condition: bool, check_id: str, detail: str) -> None:
        self.rows.append({"check_id": check_id, "status": "PASS" if condition else "FAIL", "detail": detail})

    @property
    def passed(self) -> bool:
        return bool(self.rows) and all(row["status"] == "PASS" for row in self.rows)


def command_output(command: list[str]) -> str:
    return subprocess.run(command, check=True, capture_output=True, text=True).stdout


def pdf_pages(path: Path) -> int:
    for line in command_output(["pdfinfo", str(path)]).splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise RuntimeError(f"pdfinfo did not report pages: {path}")


def validate(project_root: Path, input_root: Path | None, output_root: Path | None):
    import pandas as pd

    renderer = import_renderer(Path(__file__).resolve().parent)
    root = project_root.resolve()
    input_path = (input_root or renderer.default_input(root)).resolve()
    output_path = (output_root or renderer.default_output(root)).resolve()
    if not output_path.is_dir():
        raise RuntimeError(f"figure candidate is missing: {output_path}")
    if (output_path / "READY").exists() or (output_path / "validation_report.tsv").exists():
        raise RuntimeError("refusing to overwrite validated figure artifacts")
    checks = Checks()

    renderer.verify_input(input_path)
    manifest = pd.read_csv(output_path / "panel_manifest.tsv", sep="\t", dtype=str)
    checks.require(len(manifest) == 1, "one_panel_manifest_row", f"rows={len(manifest)}")
    row = manifest.iloc[0]
    panel = output_path / renderer.PANEL_NAME
    source_path = output_path / renderer.SOURCE_NAME
    checks.require(row["figure_release_id"] == renderer.FIGURE_RELEASE_ID, "figure_release", renderer.FIGURE_RELEASE_ID)
    checks.require(row["canonical_written"] == "FALSE", "canonical_write_forbidden", "FALSE")
    checks.require(
        row["input_release_id"] == renderer.INPUT_RELEASE_ID
        and row["input_ready_sha256"] == renderer.INPUT_READY_SHA256,
        "input_release_link",
        renderer.INPUT_RELEASE_ID,
    )
    checks.require(panel.is_file() and row["panel_sha256"] == renderer.sha256_file(panel), "panel_checksum", "SHA256 rederived")
    checks.require(source_path.is_file() and row["source_sha256"] == renderer.sha256_file(source_path), "source_checksum", "SHA256 rederived")

    observed = pd.read_csv(source_path, sep="\t")
    expected = renderer.build_source(input_path)
    expected_columns = list(observed.columns)
    expected = expected[expected_columns]
    try:
        pd.testing.assert_frame_equal(observed, expected, check_dtype=False, rtol=1e-12, atol=1e-14)
        exact = True
        detail = "four plotted rows rederived from validated sensitivity"
    except AssertionError as exc:
        exact = False
        detail = str(exc).splitlines()[0]
    checks.require(exact, "source_table_exact", detail)
    checks.require(
        len(observed) == 4
        and not observed.duplicated(["dataset", "program_id"]).any(),
        "complete_family",
        "2 programs x 2 datasets",
    )
    checks.require(
        set(observed.loc[observed["dataset"] == "Vu_et_al_2025", "evidence_state"]) == {"source_dependent"},
        "vu_source_dependence",
        "both Vu rows remain source-dependent",
    )
    checks.require(
        observed["claim_boundary"].eq(
            "attenuation_is_composition_linked_context_not_mediation_causality_or_cell_intrinsic_proof"
        ).all(),
        "claim_boundary",
        "no mediation, causality, or cell-intrinsic upgrade",
    )
    checks.require(observed["matched_set_hash_identical"].astype(str).str.lower().eq("true").all(), "matched_sets", "4/4 identical")
    checks.require(pdf_pages(panel) == 1, "one_page_pdf", "Pages=1")
    checks.require(
        re.search(rb"/Subtype\s*/Type3\b", panel.read_bytes()) is None,
        "no_type3_fonts",
        "raw PDF object scan",
    )
    text = command_output(["pdftotext", str(panel), "-"])
    checks.require(
        "All 16 cell types + QC" in text and "Hepatocyte + QC adjustment" in text,
        "adjustment_labels",
        "both prespecified adjustments present",
    )
    checks.require(
        "source-dependent" in text and "not evidence of mediation or causality" in text,
        "paper_claim_boundaries_visible",
        "source dependence and noncausal wording rendered",
    )
    checks.require(
        "Stromal ECM (IGFBP7)" in text and "Ductular injury (BICC1)" in text,
        "program_labels",
        "both frozen program identities visible",
    )

    report = output_path / "validation_report.tsv"
    pd.DataFrame(checks.rows).to_csv(report, sep="\t", index=False)
    if not checks.passed:
        failures = [row["check_id"] for row in checks.rows if row["status"] == "FAIL"]
        raise RuntimeError("figure validation failed: " + ", ".join(failures))
    pd.DataFrame(
        [{
            "figure_release_id": renderer.FIGURE_RELEASE_ID,
            "status": "pass_full_composition_supplementary_figure",
            "panel_sha256": renderer.sha256_file(panel),
            "source_sha256": renderer.sha256_file(source_path),
            "validation_report_sha256": renderer.sha256_file(report),
            "input_ready_sha256": renderer.INPUT_READY_SHA256,
            "canonical_written": "FALSE",
        }]
    ).to_csv(output_path / "READY", sep="\t", index=False)
    print(f"full-composition figure validated: {len(checks.rows)} checks passed", flush=True)
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        validate(args.project_root, args.input_root, args.output_root)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
